/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"# SUTRCMEDIAN ",
"Rolling median filter over a panel of seismic traces by an ordered trace buffer ",
" ",
"## Usage ",
"   sutrcmedian < stdin > stdout ",
"### Optional Parameters ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| ntr=      | number (odd) of traces in filter panel          | 5             |",
"| mode=     | =0 output filtered trace, =1 output noise       | 0             |",
"| verbose=  | =0 no advisory messages, =1 for messages        | 0             |",
" ",
"## Notes ",
"This is primarily a demonstation and test platform for the ordered trace buffer implementation. ",
" ",
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
    hOTB otbHandle;
	
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
    otbHandle = OTB_init( ntr, nsamples );
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (seismic) {
            if (OTB_push( otbHandle, &tr )) {
                tcount = OTB_traces(otbHandle);
                imed = tcount/2;
                for (is=0; is<nsamples; is++) {
                    int icur = OTB_getSlice( otbHandle, is, databuf );
                    float curval = databuf[icur];
                    qkfind( imed, tcount, databuf );
                    tr.data[is] = (mode==1)? curval - databuf[imed]: databuf[imed];
                }
                OTB_copyCurrentHdr( otbHandle, &tr );
                puttr(&tr);
            }
        } else 
            if (verbose) warn("skipping non-seismic trace with trid=%d", tr.trid);
    } while (gettr(&tr));

/* Handle last traces in buffer */
    while (OTB_push( otbHandle, 0 )) {
        tcount = OTB_traces(otbHandle);
        imed = tcount/2;
        float curval;
        for (is=0; is<nsamples; is++) {
            int icur = OTB_getSlice( otbHandle, is, databuf );
            curval = databuf[icur];
            qkfind( imed, tcount, databuf);
            tr.data[is] = (mode==1)? curval - databuf[imed]: databuf[imed];
        }
        OTB_copyCurrentHdr( otbHandle, &tr );
        puttr(&tr);
    };

    free1(databuf);
    OTB_free( otbHandle );

    return EXIT_SUCCESS;
}
