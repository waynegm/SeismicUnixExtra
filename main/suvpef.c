/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPEF: $Revision: 1.47 $ ; $Date: 2014/12/04 23:10:17 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"

/****************************** self documentation ******************************/
char *sdoc[] = {
"# SUVPEF ",
"Wiener (least squares) predictive error filtering with spatially varying lag ",
" ",
"## Usage ",
"   suvpef < stdin > stdout  [optional parameters] ",
" ",
"### Required Parameters ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| dt=       | time sampling interval (sec)                    | trcheader     |",
" ",
"### Optional Parameters ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| key=      | header word specifying lag locations            | cdp           |",
"| xlag=     | array of lag locations as given by key          |               |",
"| lag=      | array of prediction filter lags in seconds      | dt            |",
"| len=      | operator length in seconds                      | tr.ns*dt/20   |",
"| pnoise=   | relative additive noise level                   | 0.001         |",
"| mincorr=  | start of autocorrelation window in seconds      | tmin          |",
"| maxcorr=  | end of autocorrelation window in seconds        | tmax          |",
"| verbose=  | 0 - no advisory messages, 1 - for messages      | 0             |",
"                                                                               ",
"Trace header fields accessed: ns, dt                                           ",
"Trace header fields modified: none                                             ",
"                                                                               ",
"## Notes                                                                       ",
"For fixed lag omit the key= and xlag= parameters and give a single lag= value. ",
"For spatially variable lag specify key=, xlag= and lag= arrays.                ",
"Linear interpolation and constant extrapolation used to compute lag from the ",
"arrays.                                                                        ",
"                                                                               ",
"This is a simplified version of supef                                          ",                                                                               
"                                                                               ",
NULL};
/*********************** end self doc *******************************************/

segy intrace, outtrace;

int
main(int argc, char **argv)
{
    float dt;
    int ns;
    cwp_Bool seismic;

    cwp_String key;
    int nlag;
    float* xlag=NULL;
    float* lag=NULL;
    int ilag;
    float len;
    int ilen, lenbytes;
    float pnoise;
    float mincorr;
    float maxcorr;
    int imincorr, imaxcorr, ncorr, lcorr;
    int verbose;

    cwp_String keyType;
    int keyIndex;
    Value keyVal;
    float val;
    
    float* wiener=NULL;
    float* spiker=NULL;
    float* crosscorr=NULL;
    float* autocorr=NULL;
    
/* Initialize */
    initargs(argc, argv);
    requestdoc(1);

/* Get info from first trace */ 
    if (!gettr(&intrace)) err("can't get first trace");
    ns = intrace.ns;
    if (!getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
    if (!dt) err("dt field is zero and not getparred");

/* Get parameters */
    if (!getparint("verbose", &verbose)) verbose=0;
    if (!getparstring("key", &key)) key="cdp";
    nlag = countparval("xlag");
    if (nlag>0)
        if (countparname("lag")!=nlag) err("a lag value must be given for each %s in xlag=", key);
    else
        nlag = 1;
    xlag = ealloc1float(nlag);
    if (!getparfloat("xlag",xlag)) xlag[0] = intrace.cdp;
    lag = ealloc1float(nlag);
    if (!getparfloat("lag",lag)) lag[0] = dt;
    
    if (getparfloat("mincorr", &mincorr)) {
        if (mincorr < 0.0) {
            warn("mincorr=%g too small resetting to start of trace", mincorr);
            imincorr=0;
        } else
            imincorr = NINT(mincorr/dt);
    } else
        imincorr = 0;
    
    if (getparfloat("maxcorr", &maxcorr)) {
        imaxcorr = NINT(maxcorr/dt);
        if (imaxcorr > ns) {
            warn("maxcorr=%g too big resetting to end of trace", maxcorr);
            imaxcorr = ns;
        } else if (imaxcorr<imincorr) 
            err("maxcorr: %g is less than mincorr: %g", maxcorr, mincorr);
    } else
        imaxcorr = ns;
    
    if (!getparfloat("pnoise", &pnoise) pnoise = 0.001;
    if (!getparfloat("len", &len) len = (float) ns*dt/20
    ilen = NINT(len/dt);

/* Get key type and index */
    keyType = hdtype(key);
    keyIndex = getindex(key);
    
    ncorr = imaxcorr - imincorr + 1;

/* Compute byte sizes in wiener/spiker and autocorr */
    lenbytes = FSIZE*ilen;

/* Allocate memory */
    wiener	 = ealloc1float(ilen);
    spiker	 = ealloc1float(ilen);
    autocorr = ealloc1float(ncorr);

/* Main processing loop */
    do {
        seismic = ISSEISMIC(intrace.trid);
        if (!seismic) {
            if (verbose)
                warn("ignoring input trace=%d with non-seismic trcid=%d", tr.tracl, tr.trid);
            continue;
        }
        gethval(&intrace, keyIndex, &keyVal);
        val = vtof(keyType, keyVal);
        
        float filtlag = linterp(nlag, xlag, lag, val);
        ilag = NINT(filtlag/dt);
        lcorr = ilag + ilen + 1;
        crosscorr = autocorr + ilag;
        
/* zero out filter vectors */
        memset((void *) wiener, 0, lenbytes);
        memset((void *) spiker, 0, lenbytes);
        memset((void *) autocorr, 0, ncorr*FSIZE);

/* Form autocorrelation vector */
        xcor(ncorr, imincorr, intrace.data, ncorr, imincorr, intrace.data, lcorr, 0, autocorr);

/* Leave trace alone if autocorr[0] vanishes */
        if (autocorr[0] == 0.0) {
            puttr(&intrace);
            continue;
        }

/* Whiten */
        autocorr[0] *= 1.0 + pnoise;

/* Get inverse filter by Wiener-Levinson */
        stoepf(ilen, autocorr, crosscorr, wiener, spiker);

/* Convolve pefilter with trace - don't do zero multiplies */
        for (int i = 0; i < nt; ++i) {
            register int j;
            register int n = MIN(i, ilag); 
            register float sum = intrace.data[i];

            for (j = ilag; j <= n; ++j)
                sum -= wiener[j-ilag] * intrace.data[i-j];

            outtrace.data[i] = sum;
        }

/* Output filtered trace */
        memcpy( (void *) &outtrace, (const void *) &intrace, HDRBYTES);
        puttr(&outtrace);
    } while (gettr(&intrace));

    free1float(wiener);
    free1float(spiker);
    free1float(autocorr);
    free1float(xlag);
    free1float(lag);

    return(CWP_Exit());
}
