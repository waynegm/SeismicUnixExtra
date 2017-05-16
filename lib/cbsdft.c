/* Copyright (c) Wayne Mogg, 2017. */
/* All rights reserved.            */

/*********************** self documentation **********************/
/*************************************************************************
CBSDFT - circular buffer for sliding discrete fourier transform of multi-trace panel

CBSDFT_init     initialise a SDFT transformer handle
CBSDFT_traces   number of CBSDFT_traces
CBSDFT_samples  number of time samples per trace
CBSDFT_size     number of samples in SDFT window
CBSDFT_nfreq    number of frequencies in SDFT
CBSDFT_push     add the next trace to the buffer
CDSDFT_getframe get the spectral data for all traces at the specified time, frequency index
CDSDFT_current  index of central trace in multi-trace panel
CBSDFT_free     release a SDFT transformer handle

************************************************************************** 
Author: Wayne Mogg
**************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"
#include "par.h"
#include "sux.h"

struct _CBSDFT {
    int ns;
    int nwin;
    int ntr;
    suxWindow window;
    int ifirst;
    int ilast;
    int icur;
    hSDFT sdftH;
    complex** resbuf;
};

hCBSDFT CBSDFT_init( int ntraces, int nsamples, int nwin, suxWindow window ) {
    
    hCBSDFT h = emalloc(sizeof(struct _CBSDFT));
    h->ntr = ntraces;
    h->ns = nsamples;
    h->nwin = nwin;
    h->window = window;
    int hw = nwin/2;
    int nf = hw + 1;
    h->sdftH = SDFT_init( nwin, nsamples );
    h->resbuf = ealloc2complex( nsamples, nf );
    h->ifirst = 0;
    h->ilast = 0;
    h->icur = 0;
    return h;
}

void CBSDFT_free( hCBSDFT h ) {
    SDFT_free( h->sdftH );
    free( h );
    h = 0;
}

size_t CBSDFT_traces( hCBSDFT h ) {
    return h ? h->ntr : 0;
}

size_t CBSDFT_samples(  hCBSDFT h ) {
    return h ? h->ns : 0;
}

size_t CBSDFT_size(  hCBSDFT h ) {
    return h ? h->nwin : 0;
}

size_t CBSDFT_nfreq(  hCBSDFT h ) {
    return h ? h->nwin/2+1 : 0;
}

void CBSDFT_push( hCBSDFT h, float* data ) {
    int nf = h->nwin/2 + 1;
    memset((void*) h->resbuf, 0, nf * h->ns * CSIZE );
    
    SDFT( h->sdftH, h->window, data, SPECBUF );
}

void CBSDFT_getframe( hCBSDFT h, int isample, int ifreq, complex* data ) {
    
}

void CBSDFT_setresult( hCBSDFT h, int isample, int nfreq, complex data ) {
    h->resbuf[nfreq][isample] = data;
}

void CBSDFT_getresult( hCBSDFT h, float* data ) {
    ISDFT( h->sdftH, h->resbuf, data );
}
