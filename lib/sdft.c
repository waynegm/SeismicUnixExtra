/* Copyright (c) Wayne Mogg, 2017. */
/* All rights reserved.            */

/*********************** self documentation **********************/
/*************************************************************************
SDFT - sliding discrete fourier transform

initSDFT        calculate the twiddle factors for the DFT
SDFT            calculate the sliding DFT
ISDFT           calculate the inverse sliding DFT

************************************************************************** 
Author: Wayne Mogg
**************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"
#include "sux.h"

void initSDFT( int nwin, complex* cfactors ) {
    float fact;
    for (int i=0; i<nwin/2+1; i++) {
        fact = 2.0 * PI * (float)i/(float)nwin;
        cfactors[i] = cmplx(cos(fact),sin(fact));
    }
}

void SDFT( int ns, int nwin, float* data, complex* cfactors, complex** result ){
    int i, its, ifr, hw;
    complex cval, fm1, f;
    float val, fact, jfact;
    int nf = nwin/2+1;
    
    hw = nwin/2;
    
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
    
    for (its=1; its<ns; its++) {
        if (its<=hw)
            val = data[its+hw] - data[0];
        else if (its>=ns-hw)
            val = data[ns-1] - data[its-hw-1];
        else
            val = data[its+hw] -data[its-hw-1];
        cval = cmplx(val,0.0);
        for( ifr=0; ifr<nf; ifr++)
            result[ifr][its] = cmul(cadd(result[ifr][its-1],cval),cfactors[ifr]);
    }
    // Apply Hamming window
    for (its=0; its<ns; its++) {
        fm1 = cmplx(0.0,0.0);
        for (ifr=0; ifr<nf-1;ifr++) {
            f = crmul(cadd(result[ifr+1][its], fm1), -0.23);
            f = cadd(f, crmul(result[ifr][its], 0.54));
            fm1 = result[ifr][its];
            result[ifr][its] = f;
        }
        result[nf-1][its] = cadd( crmul(fm1, -0.23), crmul(result[nf-1][its],0.54) );
    }
}

void ISDFT( int ns, int nwin, complex** specdata, float* result ) {
}
