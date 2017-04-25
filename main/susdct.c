/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                   ",
" SUSDCT -  Outputs a time-frequency representation of seismic data ",
"           via the sliding discrete cosine transform               ",
"                                                                   ",
" susdct <stdin >sdout                                              ",
"                                                                   ",
" Required parameters:                                              ",
" dt=(from header)  time sampling interval (sec)                    ",
"                                                                   ",
" Optional parameters:                                              ",
" nwin=31           number (odd) of samples in SDFT window          ",
" verbose=0         no advisory messages                            ",
"         1         for advisory messages                           ",
"                                                                   ",
" This process provides a time-frequency representation of seismic  ",
" data using the sliding discrete cosine transform.                 ",
"                                                                   ",
" Examples:                                                         ",
"   suvibro | susdct | suximage                                     ",
"                                                                   ",
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
    
/* Set up DCT parameters and workspaces */
    df = 1.0/(2.0*nwin*dt);
    nf = nwin;
    dctHandle = initSDCT( nwin, nt );
    specbuff = ealloc2float(nt, nf );
    
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (!seismic) {
            if (verbose)
                warn("ignoring input trace=%d with non-seismic trcid=%d", tr.tracl, tr.trid);
            continue;
        }
        SDCT( dctHandle, tr.data, specbuff );
        windowSDCT( dctHandle, Blackman, specbuff );
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
                
    destroySDCT(dctHandle);
    free2float( specbuff );

    return (CWP_Exit());
}
