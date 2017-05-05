/* Copyright (c) Wayne Mogg, 2017. */
/* All rights reserved.            */

/*********************** self documentation **********************/
/*************************************************************************
SDCT - sliding discrete cosine transform (Type II)

SDCT_init       initialise a SDCT transformer handle
SDCT            calculate the sliding DCT
ISDCT           calculate the inverse sliding DCT
SDCT_free       release a SDCT transformer handle
SDCT_window     apply a window to the SDCT transform output

************************************************************************** 
Author: Wayne Mogg
**************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"
#include "par.h"
#include "sux.h"

struct _SDCT {
    int ns;
    int nwin;
    float* cosfact;
    float* cosfact2;
};

hSDCT SDCT_init( int nwin, int nsamples ) {
    float fact;
    
    hSDCT handle = emalloc(sizeof(struct _SDCT));
    handle->ns = nsamples;
    handle->nwin = nwin;
    handle->cosfact = ealloc1float(nwin);
    handle->cosfact2 = ealloc1float(nwin);

    for (int i=0; i<nwin; i++) {
        fact = PI * (float)i/(float)nwin;
        handle->cosfact[i] = cos(fact/2.0);
        handle->cosfact2[i] = 2.0 * cos(fact);
    }

    return handle;
}

void SDCT_free( hSDCT handle ) {
    free1float(handle->cosfact);
    free1float(handle->cosfact2);
    free(handle);
    handle = 0;
}

void SDCT( hSDCT h, sux_Window window, float* data, float** result ){
    int i, its, ifr, hw, neg1;
    float val, fact, cosfact;
    float frp1, fr, fmr, fmrm1;

    int ns = h->ns;
    int nwin = h->nwin;
    hw = nwin/2;

/* Calculate DCT directly for first 2 positions */    
    for (ifr=0; ifr<nwin; ifr++) {
        fact = PI * (float)ifr/(float)nwin;
        for (its=0; its<=1; its++) {
            val = 0.0;
            for (i=-hw; i<=hw; i++) {
                cosfact = cos(fact * (i+hw+0.5));
                val += cosfact * ((i+its<0)? data[0] : data[i+its]);
            }
            result[ifr][its] = val;
        }
    }

/* Calaculate rest of DCT using sliding algorithm */
    for (its=2; its<ns; its++) {
        fmrm1 = (its-hw-2<0)? data[0] : data[its-hw-2];
        fmr = (its-hw-1<0)? data[0] : data[its-hw-1];
        frp1 = (its+hw>ns-1)? data[ns-1] : data[its+hw];
        fr = (its+hw-1>ns-1)? data[ns-1] : data[its+hw-1];
        neg1 = 1;
        for( ifr=0; ifr<nwin; ifr++) {
            val = fmrm1 - fmr + neg1 * (frp1 - fr); 
            result[ifr][its] = result[ifr][its-1] * h->cosfact2[ifr] - result[ifr][its-2] + h->cosfact[ifr] * val;
            neg1 *= -1;
        }
    }
    
/* Apply the window in the transform domain */
    SDCT_window( h, window, result );
    for (its=0; its<ns; its++)
        result[0][its] = result[0][its] / sqrt(2.0);
}


void SDCT_window( hSDCT h, sux_Window window, float** data ) {
    if ( window==None ) return;
    float a0, a1, a2;
    float cm2, cp2, cm4, cp4;
    int its, ifr;
    int ns = h->ns;
    int nwin = h->nwin;
    float* work = ealloc1float(nwin);
    
    switch(window) {
        case Hann:
            a0 = 0.5;
            a1 = -0.25;
            a2 = 0.0;
            break;
        case Hamming:
            a0 = 0.54;
            a1 = -0.23;
            a2 = 0.0;
            break;
        case Blackman:
            a0 = 0.42;
            a1 = -0.25;
            a2 = 0.04;
            break;
        default:
            err("unrecognised window function: %d", window);
    }
    
    for ( its=0; its<ns; its++ ) {
        memset((void*)work, 0, nwin*FSIZE);
        for ( ifr=0; ifr<nwin; ifr++ ) {
            cm2 = (ifr-2<0)? data[ifr+2][its] : data[ifr-2][its];
            cm4 = (ifr-4<0)? data[ifr+4][its] : data[ifr-4][its];
            cp2 = (ifr+2>=nwin)? data[nwin-2][its] : data[ifr+2][its];
            cp4 = (ifr+4>=nwin)? data[nwin-4][its] : data[ifr+4][its];
            work[ifr] = a0 * data[ifr][its] + a1 * (cm2 + cp2) + a2 * (cm4 + cp4);
        }
        for ( ifr=0; ifr<nwin; ifr++ )
            data[ifr][its] = work[ifr];
        
    }
    free1float(work);
}

void ISDCT( hSDCT h, float** specdata, float* result ) {
    int its, ifr, hw, neg1;
    float val;
    
    hw = h->nwin/2;
    
    for ( its=0; its<h->ns; its++ ) {
        val = 0.0;
        neg1 = -1;
        for ( ifr=1; ifr<=hw; ifr++ ) {
            val += (float)neg1 * specdata[2*ifr][its];
            neg1 *= -1;
        }
        result[its] = (specdata[0][its] * sqrt(2.0) + 2.0 * val) / (float)h->nwin;
    }
}

