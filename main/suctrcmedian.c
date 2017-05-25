/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                   ",
" SUTRCMEDIAN - rolling median filter over a panel of seismic traces",
"                                                                   ",
" sutrcmedian <stdin >sdout                                         ",
" Optional parameters:                                              ",
" ntr=5         number (odd) of traces to compute median            ",
" mode=0        output filtered trace                               ",
"     =1        output estimated noise                              ",
" verbose=0     no advisory messages                                ",
"         1     for advisory messages                               ",
"                                                                   ",
"                                                                   ",
NULL};

/* Author: Wayne Mogg, May 2017
 *
 * Trace header fields accessed: ns, trid
 */
/**************** end self doc ***********************************/

segy tr;

int
main(int argc, char **argv)
{
    int ntr;
    int mode;
    int verbose;

    int is, tcount, imed;
    int nsamples;
    cwp_Bool seismic;
    float* databuf;
    hCTB ctbHandle;
	
// Initialize
	initargs(argc, argv);
	requestdoc(1);

/* Get info from first trace */ 
    nsamples = 0;
    while (gettr(&tr)) {
        seismic = ISSEISMIC(tr.trid);
        if (seismic) {
            nsamples = tr.ns;
            break;
        } else
            warn("skipping non-seismic trace with trid=%d", tr.trid);
    }
    if (nsamples==0) 
        err("zero length traces not allowed.");
    
// Get parameters
    if (!getparint("ntr", &ntr)) ntr = 5;
    if (!ntr%2) {
        ntr++;
        if (verbose)
            warn("adjusting ntr to be odd, was %d now %d",ntr-1, ntr);
    }
    if (!getparint("mode", &mode)) mode = 0;
    if (!getparint("verbose", &verbose)) verbose=0;
    
// Set up trace buffer and work space
    databuf = ealloc1float( ntr );
    ctbHandle = CTB_init( ntr, nsamples );
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (seismic) {
            if (CTB_push(ctbHandle, &tr)) {
                tcount = CTB_traces(ctbHandle);
                imed = tcount/2;
                for (is=0; is<nsamples; is++) {
                    int icur = CTB_getSlice( ctbHandle, is, databuf );
                    float curval = databuf[icur];
                    qkfind( imed, tcount, databuf );
                    tr.data[is] = (mode==1)? curval - databuf[imed]: databuf[imed];
                }
                CTB_copyCurrentHdr( ctbHandle, &tr );
                puttr(&tr);
            }
        } else 
            if (verbose) warn("skipping non-seismic trace with trid=%d", tr.trid);
    } while (gettr(&tr));

/* Handle last traces in buffer */
    while(CTB_push(ctbHandle, 0)) {
        tcount = CTB_traces(ctbHandle);
        imed = tcount/2;
        float curval;
        for (is=0; is<nsamples; is++) {
            int icur = CTB_getSlice( ctbHandle, is, databuf );
            curval = databuf[icur];
            qkfind( imed, tcount, databuf);
            tr.data[is] = (mode==1)? curval - databuf[imed]: databuf[imed];
        }
        CTB_copyCurrentHdr( ctbHandle, &tr );
        puttr(&tr);
    };

    free1(databuf);
    CTB_free( ctbHandle );

    return EXIT_SUCCESS;
}
