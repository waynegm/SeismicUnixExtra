/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"# SULPASMOOTH                                                                  ",
"Rolling LPA filter over a panel of seismic traces                              ",
"                                                                               ",
"## Usage                                                                       ",
"   sulpasmooth < stdin > sdout [optional parameters]                           ",
"                                                                               ",
"### Optional Parameters                                                        ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| mode=     | =0 output filtered trace, =1 output noise       | 0             |",
"| verbose=  | =0 no advisory messages, =1 for messages        | 0             |",
"                                                                               ",
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
    int nsize=5;
    int ntr = nsize;
    int mode;
    int verbose;

    int is;
    int nsamples;
    cwp_Bool seismic;
    float** db;
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
    if (!getparint("mode", &mode)) mode = 0;
    if (!getparint("verbose", &verbose)) verbose=0;
    
// Set up trace buffer and work space
    db = ealloc2float( nsize, ntr );
    otbHandle = OTB_init( ntr, nsamples );
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (seismic) {
            if (OTB_push( otbHandle, &tr )) {
                for (is=0; is<nsamples; is++) {
                    OTB_getSlab(otbHandle, is, nsize, db);
                    float outval = (-13*(db[0][0]+db[4][0]+db[0][4]+db[4][4])+
                                    2*(db[1][0]+db[3][0]+db[0][1]+db[4][1]+db[0][3]+db[4][3]+db[1][4]+db[3][4])+
                                    7*(db[2][0]+db[0][2]+db[4][2]+db[2][4])+
                                   17*(db[1][1]+db[3][1]+db[1][3]+db[3][3])+
                                   22*(db[2][1]+db[2][3]+db[1][2]+db[3][2])+
                                   27*db[2][2])/175;
                    tr.data[is] = (mode==1)? db[nsize/2][nsize/2] - outval: outval;
                }
                OTB_copyCurrentHdr( otbHandle, &tr );
                puttr(&tr);
            }
        } else 
            if (verbose) warn("skipping non-seismic trace with trid=%d", tr.trid);
    } while (gettr(&tr));

/* Handle last traces in buffer */
    while (OTB_push( otbHandle, 0 )) {
        for (is=0; is<nsamples; is++) {
            OTB_getSlab(otbHandle, is, nsize, db);
            float outval = (-13*(db[0][0]+db[4][0]+db[0][4]+db[4][4])+
            2*(db[1][0]+db[3][0]+db[0][1]+db[4][1]+db[0][3]+db[4][3]+db[1][4]+db[3][4])+
            7*(db[2][0]+db[0][2]+db[4][2]+db[2][4])+
            17*(db[1][1]+db[3][1]+db[1][3]+db[3][3])+
            22*(db[2][1]+db[2][3]+db[1][2]+db[3][2])+
            27*db[2][2])/175;
            tr.data[is] = (mode==1)? db[nsize/2][nsize/2] - outval: outval;
        }
        OTB_copyCurrentHdr( otbHandle, &tr );
        puttr(&tr);
    };

    free2float(db);
    OTB_free( otbHandle );

    return EXIT_SUCCESS;
}
