/* Copyright (c) Wayne Mogg, 2017. */
/* All rights reserved.            */

/*********************** self documentation **********************/
/*************************************************************************
SDFT - sliding discrete fourier transform

SDFT_init       initialise a SDFT transformer handle
SDFT            calculate the sliding DFT
ISDCT           calculate the inverse sliding DFT
SDFT_free       release a SDFT transformer handle
SDFT_window     apply a window to the SDFT transform output

************************************************************************** 
Author: Wayne Mogg
**************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"
#include "par.h"
#include "sux.h"

struct _SDFT {
    int ns;
    int nwin;
    complex* cfact;
    complex* cfactinv;
};

hSDFT SDFT_init( int nwin, int nsamples ) {
    float fact;
    
    hSDFT h = emalloc(sizeof(struct _SDFT));
    h->ns = nsamples;
    h->nwin = nwin;
    int hw = nwin/2;
    int nf = hw + 1;
    h->cfact = ealloc1complex(nf);
    h->cfactinv = ealloc1complex(nwin);
    
    for (int i=0; i<nf; i++) {
        fact = 2.0 * PI * (float)i/(float)nwin;
        h->cfact[i] = cmplx(cos(fact),sin(fact));
    }
    
    for (int i=0; i<nwin; i++) {
        fact = 2.0 * PI * (float)i * (float)hw / (float)nwin;
        h->cfactinv[i] = cmplx(cos(fact),sin(fact));
    }
    
    return h;
}

void SDFT_free( hSDFT h ) {
    free1complex( h->cfact );
    free( h );
    h = 0;
}

void SDFT( hSDFT h, sux_Window window, float* data, complex** result ) {
    int i, its, ifr;
    complex cval;
    float oldv, newv, fact, jfact;
    
    int ns = h->ns;
    int nwin = h->nwin;
    int hw = nwin/2;
    int nf = hw + 1;
    
/* Calculate DFT directly for the first position */    
    for (ifr=0; ifr<nf; ifr++) {
        fact = -2.0 * PI * (float)ifr/(float)nwin;
        cval = cmplx(0.0,0.0);
        for (i=-hw; i<=hw; i++) {
            jfact = fact * (float)(i+hw);
            if (i<0) 
                cval = cadd(cval, crmul(cmplx(cos(jfact),sin(jfact)), data[0]));
            else
                cval = cadd(cval, crmul(cmplx(cos(jfact),sin(jfact)), data[i]));
        }
        result[ifr][0] = cval;
    }
    
/* Calculate rest of DFT using sliding algorithm */    
    for (its=1; its<ns; its++) {
        oldv = (its-hw-1<0)? data[0] : data[its-hw-1];
        newv = (its+hw>ns-1)? data[ns-1] : data[its+hw];
        cval = cmplx(newv-oldv,0.0);
        for( ifr=0; ifr<nf; ifr++)
            result[ifr][its] = cmul(cadd(result[ifr][its-1],cval),h->cfact[ifr]);
    }
    
/* Apply the window in the transform domain */
    SDFT_window( h, window, result );
}

void SDFT_window( hSDFT h, sux_Window window, complex** data ) {
    if (window==None) return;
    float a0, a1, a2;
    complex cm1, cm2, cp1, cp2;
    int its, ifr;
    
    int ns = h->ns;
    int nwin = h->nwin;
    int nf = nwin/2+1;
    complex* work = ealloc1complex(nf);
    
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
        memset((void*)work, 0, nf*CSIZE);
        for ( ifr=0; ifr<nf; ifr++ ) {
            cm1 = (ifr-1<0)?  conjg(data[1-ifr][its]): data[ifr-1][its];
            cm2 = (ifr-2<0)?  conjg(data[2-ifr][its]): data[ifr-2][its];
            cp1 = (ifr+1>=nf)? conjg(data[nwin-ifr-1][its]): data[ifr+1][its];
            cp2 = (ifr+2>=nf)? conjg(data[nwin-ifr-2][its]): data[ifr+2][its];
            work[ifr] = cadd(cadd(crmul(data[ifr][its], a0), crmul(cadd(cp1,cm1),a1)), crmul(cadd(cp2,cm2),a2));
        }
        for (ifr=0; ifr<nf; ifr++)
            data[ifr][its] = work[ifr];
    }
    free1complex(work);
}
    
void ISDFT( hSDFT h, complex** specdata, float* result ) {
    int its, ifr;
    complex val, F;
    int nwin = h->nwin;
    int nf = nwin/2+1;
    int ns = h->ns;

    for (its=0; its<ns; its++) {
        val=cmplx(0.0,0.0);
        for (ifr=0; ifr<nwin; ifr++) {
            F = (ifr>nf-1)? conjg(specdata[nwin-ifr][its]) : (specdata[ifr][its]);
            val = cadd(val, cmul(F, h->cfactinv[ifr]));
        }
        result[its] = val.r/(float)nwin;
    }
}

