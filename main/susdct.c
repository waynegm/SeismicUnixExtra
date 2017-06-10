/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"# SUSDCT ",
"Time-frequency decomposition using the sliding discrete cosine transform ",
" ",
"## Usage ",
"   susdct < stdin > stdout ",
" ",
"### Required Parameters ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| dt=       | time sampling interval (sec)                    | trcheader     |",
"                                                                               ",
"### Optional Parameters                                                        ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| nwin=     | number (odd) of samples in SDFT window          | 31            |",
"| window=   | none - no window applied to each data segment   | none          |",
"|           | hann - Hann window                              |               |",
"|           | hamming - Hamming window                        |               |",
"|           | blackman - Blackman window                      |               |",
"| verbose=  | 0 - no advisory messages, 1 - for messages      | 0             |",
" ",
"## Notes ",
"This process calculates a time-frequency decomposition of seismic data using ",
"the sliding discrete cosine transform. ",
" ",
"## Examples ",
"   suvibro | susdct | suximage ",
" ",
NULL};

/* Author: Wayne Mogg, Apr 2017
 *
 * Trace header fields accessed: ns,dt, trid, ntr
 * Trace header fields modified: tracl, tracr, d1, f2, d2, trid, ntr
 */
/**************** end self doc ***********************************/


segy tr;

int
main(int argc, char **argv)
{
    float dt;
    int nwin;
    int verbose;
    cwp_String window;
    sux_Window iwind = None;

    int nt;
    float df;
    int tracr=0;
    
    int nf;
    int i,j;
    cwp_Bool seismic;
    hSDCT dctHandle;
    float** specbuff;
	
/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

/* Get info from first trace */ 
    if (!gettr(&tr))  err("can't get first trace");
    nt = tr.ns;
    if (!getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
    if (!dt) err("dt field is zero and not getparred");
    
/* Get parameters */
    if (!getparint("verbose", &verbose)) verbose=0;
    if (!getparint("nwin", &nwin)) nwin = 31;
    if (nwin%2==0) {
        nwin++;
        if (verbose)
            warn("adjusting nwin to be odd, was %d now %d",nwin-1, nwin);
    }
    if (!getparstring("window", &window)) window = "none";
    if      (STREQ(window, "hann")) iwind = Hann;
    else if (STREQ(window, "hamming")) iwind = Hamming;
    else if (STREQ(window, "blackman")) iwind = Blackman;
    else if (!STREQ(window, "none")) 
        err("unknown window=\"%s\", see self-doc", window);
    
    
/* Set up DCT parameters and workspaces */
    df = 1.0/(2.0*nwin*dt);
    nf = nwin;
    dctHandle = SDCT_init( nwin, nt );
    specbuff = ealloc2float(nt, nf );
    
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (!seismic) {
            if (verbose)
                warn("ignoring input trace=%d with non-seismic trcid=%d", tr.tracl, tr.trid);
            continue;
        }
        SDCT( dctHandle, iwind, tr.data, specbuff );
        tracr = 0;
        for ( i=0; i<nf; i++ ) {
            for (j=0; j<nt; j++)
                tr.data[j] = specbuff[i][j];
            tr.trid = AMPLITUDE;
            tr.d1 = dt;
            tr.tracr = ++tracr;
            tr.gx = tr.tracr;
            tr.f2 = 0.0;
            tr.d2 = df;
            puttr(&tr);
        }
    } while (gettr(&tr));
                
    SDCT_free(dctHandle);
    free2float( specbuff );

    return (CWP_Exit());
}
