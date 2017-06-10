/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"# SUSDFT                                                                       ",
"Time-frequency decomposition by the sliding discrete fourier transform         ",
"                                                                               ",
"## Usage                                                                       ",
"   susdft < stdin > stdout                                                     ",
"                                                                               ",
"### Required Parameters                                                        ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| dt=       | time sampling interval (sec)                    | trcheader     |",
"                                                                               ",
"### Optional Parameters                                                        ",
"| Parameter | Description                                     | Default       |",
"|:---------:| ----------------------------------------------- |:-------------:|",
"| nwin=     | number (odd) of samples in SDFT window          | 31            |",
"| mode=     | complex - output complex spectra                | complex       |",
"|           | amp - output spectral amplitude                 |               |",
"|           | phase - output spectral phase                   |               |",
"|           | real - output real part                         |               |",
"|           | imag - output imaginary part                    |               |",
"| window=   | none - no window applied to each data segment   | none          |",
"|           | hann - Hann window                              |               |",
"|           | hamming - Hamming window                        |               |",
"|           | blackman - Blackman window                      |               |",
"| verbose=  | 0 - no advisory messages, 1 - for messages      | 0             |",
"                                                                               ",
"## Notes                                                                       ",
"This process calculates a time-frequency decomposition of seismic data using ",
"the sliding discrete fourier transform.                                        ",
"                                                                               ",
"## Examples: ",
"   suvibro | susdft nwin=51 mode=amp window=hann | suximage ",
" ",
" ![susdft example](images/susdft_2.png) ",                                                                               " ",
NULL};

/* Author: Wayne Mogg, Apr 2017
 *
 * Trace header fields accessed: ns,dt, trid, ntr
 * Trace header fields modified: tracl, tracr, d1, f2, d2, trid, ntr
 */
/**************** end self doc ***********************************/


#define CPLX    1
#define REAL    2
#define IMAG    3
#define AMP     4
#define ARG     5

segy tr;

int
main(int argc, char **argv)
{
    float dt;
    int nwin;
    cwp_String mode;
    int imode=CPLX;
    int verbose;
    cwp_String window;
    sux_Window iwind = None;
    
    int nt;
    float df;
    int tracr=0;
    float re, im;
    
    int nf;
    int i,j;
    cwp_Bool seismic;
    hSDFT dftHandle;
    complex** specbuff;
	
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
    if (!getparstring("mode", &mode))	mode = "complex";
    
    if      (STREQ(mode, "phase")) imode = ARG;
    else if (STREQ(mode, "real"))  imode = REAL;
    else if (STREQ(mode, "imag"))  imode = IMAG;
    else if (STREQ(mode, "amp"))  imode = AMP;
    else if (!STREQ(mode, "complex"))
        err("unknown mode=\"%s\", see self-doc", mode);

    if (!getparstring("window", &window)) window = "none";
    if      (STREQ(window, "hann")) iwind = Hann;
    else if (STREQ(window, "hamming")) iwind = Hamming;
    else if (STREQ(window, "blackman")) iwind = Blackman;
    else if (!STREQ(window, "none")) 
        err("unknown window=\"%s\", see self-doc", window);
    
/* Set up DFT parameters and workspaces */
    df = 1.0/(nwin*dt);
    nf = nwin/2+1;
    dftHandle = SDFT_init(nwin, nt);
    specbuff = ealloc2complex(nt, nf );
    
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (!seismic) {
            if (verbose)
                warn("ignoring input trace=%d with non-seismic trcid=%d", tr.tracl, tr.trid);
            continue;
        }
        SDFT( dftHandle, iwind, tr.data, specbuff );
        
        tracr = 0;
        for ( i=0; i<nf; i++ ) {
            tr.ns = nt;
            switch (imode) {
                case CPLX:
                    for (j=0; j<nt; j++) {
                        tr.data[2*j] = specbuff[i][j].r;
                        tr.data[2*j+1] = specbuff[i][j].i;
                    }
                    tr.trid = FUNPACKNYQ;
                    tr.ns = 2 * nt;
                    break;
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
                
    SDFT_free(dftHandle);
    free2complex( specbuff );

    return (CWP_Exit());
}
