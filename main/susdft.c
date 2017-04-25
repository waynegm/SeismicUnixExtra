/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                   ",
" SUSDFT -  Outputs a time-frequency representation of seismic data ",
"           via the sliding discrete fourier transform              ",
"                                                                   ",
" susdft <stdin >sdout                                              ",
"                                                                   ",
" Required parameters:                                              ",
" dt=(from header)  time sampling interval (sec)                    ",
"                                                                   ",
" Optional parameters:                                              ",
" nwin=31           number (odd) of samples in SDFT window          ",
" mode=amp          output spectral amplitude                       ",
"     =phase        output spectral phase                           ",
"     =real         output real part                                ",
"     =imag         output imaginary part                           ",
" verbose=0         no advisory messages                            ",
"         1         for advisory messages                           ",
"                                                                   ",
" This process provides a time-frequency representation of seismic  ",
" data using the sliding discrete fourier transform. A hamming      ",
" window is used                                                    ",
"                                                                   ",
" Examples:                                                         ",
"   suvibro | susdft | suximage                                     ",
"                                                                   ",
NULL};

/* Author: Wayne Mogg, Apr 2017
 *
 * Trace header fields accessed: ns,dt, trid, ntr
 * Trace header fields modified: tracl, tracr, d1, f2, d2, trid, ntr
 */
/**************** end self doc ***********************************/


#define REAL    1
#define IMAG    2
#define AMP     3
#define ARG     4

segy tr;

int
main(int argc, char **argv)
{
    float dt;
    int nwin;
    cwp_String mode;
    int imode=AMP;
    int verbose;

    int nt;
    float df;
    int tracr=0;
    float re, im;
    
    int nf;
    int i,j;
    cwp_Bool seismic;
    complex** specbuff;
    complex* cfactors;
	
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
    if (!getparstring("mode", &mode))	mode = "amp";
    
    if      (STREQ(mode, "phase")) imode = ARG;
    else if (STREQ(mode, "real"))  imode = REAL;
    else if (STREQ(mode, "imag"))  imode = IMAG;
    else if (!STREQ(mode, "amp"))
        err("unknown mode=\"%s\", see self-doc", mode);
    
/* Set up DFT parameters and workspaces */
    df = 1.0/(nwin*dt);
    nf = nwin/2+1;
    cfactors = ealloc1complex(nf);
    initSDFT(nwin, cfactors);
    specbuff = ealloc2complex(nt, nf );
    
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (!seismic) {
            if (verbose)
                warn("ignoring input trace=%d with non-seismic trcid=%d", tr.tracl, tr.trid);
            continue;
        }
        SDFT( nt, nwin, tr.data, cfactors, specbuff );
        
        tracr = 0;
        for ( i=0; i<nf; i++ ) {
            switch (imode) {
                case REAL:
                    for (j=0; j<nt; j++) 
                        tr.data[j] = specbuff[i][j].r;
                    tr.trid = REALPART;
                    break;
                case IMAG:
                    for (j=0; j<nt; j++)
                        tr.data[j] = specbuff[i][j].i;
                    tr.trid = IMAGPART;
                    break;
                case AMP:
                    for (j=0; j<nt; j++) {
                        re = specbuff[i][j].r;
                        im = specbuff[i][j].i;
                        tr.data[j] = (float) sqrt (re * re + im * im);
                    }
                    tr.trid = AMPLITUDE;
                    break;
                case ARG:
                    for (j=0; j<nt; j++) {
                        re = specbuff[i][j].r;
                        im = specbuff[i][j].i;
                        if (re*re+im*im)
                            tr.data[j] = atan2(im,re);
                        else
                            tr.data[j] = 0.0;
                    }
                    tr.trid = PHASE;
                    break;
            }
            tr.d1 = dt;
            tr.tracr = ++tracr;
            tr.gx = tr.tracr;
            tr.f2 = 0.0;
            tr.d2 = df;
            puttr(&tr);
        }
    } while (gettr(&tr));
                
    free1complex( cfactors );
    free2complex( specbuff );

    return (CWP_Exit());
}
