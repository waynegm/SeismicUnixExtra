/* Copyright (c) Wayne Mogg, 2017. */
/* All rights reserved.            */

/*********************** self documentation **********************/
/*************************************************************************
CBSDFT - cyclic buffer for sliding discrete fourier transform of multi-trace panel

CBSDFT_init     initialise a cyclic SDFT buffer handle
CBSDFT_traces   return number of traces in the buffer
CBSDFT_samples  return number of time samples per trace
CBSDFT_size     return number of samples in SDFT window
CBSDFT_nfreq    return number of frequencies in SDFT
CBSDFT_push     add a seg y trace to the buffer
CDSDFT_getslice get the spectral data for all traces at the specified time, frequency index
CBSDFT_free     release a SDFT transformer handle

************************************************************************** 
Function Prototypes:
hCBSDFT CBSDFT_init(int ntraces, in nsamples, int nwin, suxWindow window);
void CBSDFT_free(hCBSDFT h);
int CBSDFT_traces(hCBSDFT h);
int CBSDFT_samples(hCBSDFT h);
int CBSDFT_size(hCBSDFT h);
int CBSDFT_nfreq(hCBSDFT h);
int CBSDFT_push(hCBSDFT h, const segy* const tr);
int CBSDFT_getSlice(hCBSDFT h, int isample, int ifreq, complex* const data);
void CBSDFT_setResult(hCBSDFT h, int isample, int ifreq, complex val);
void CBSDFT_getResult(hCBSDFT h, segy* const tr);

************************************************************************** 
Author: Wayne Mogg
**************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"
#include "par.h"
#include "sux.h"

typedef struct {
    unsigned char hdr[HDRBYTES];
} _HDR;

struct _CBSDFT {
    int ns;
    int nwin;
    int ntr;
    sux_Window window;
    int intr;
    int outtr;
    int trcount;
    hSDFT sdftH;
    complex** resbuf;
    _HDR* hdrs;
    complex*** specdata;
};

hCBSDFT CBSDFT_init( int ntraces, int nsamples, int nwin, sux_Window window ) {
    
    hCBSDFT h = emalloc(sizeof(struct _CBSDFT));
    h->ntr = ntraces;
    h->ns = nsamples;
    h->nwin = nwin;
    h->window = window;
    int hw = nwin/2;
    int nf = hw + 1;
    h->sdftH = SDFT_init( nwin, nsamples );
    h->resbuf = ealloc2complex( nsamples, nf );
    h->specdata = ealloc3complex( nsamples, nf, nwin );
    h->hdrs = ealloc1(ntraces, sizeof(_HDR));
    h->intr = -1;
    h->outtr = 0;
    h->trcount = 0;
    return h;
}

void CBSDFT_free( hCBSDFT h ) {
    if (h) {
        if (h->resbuf) free2complex(h->resbuf);
        if (h->specdata) free3complex(h->specdata);
        if (h->hdrs) free1(h->hdrs);
        SDFT_free(h->sdftH);
        free(h);
        h = 0;
    } else
        err("bad pointer in CBSDFT_free.");
}

int CBSDFT_traces( hCBSDFT h ) {
    return h ? h->trcount : 0;
}

int CBSDFT_samples( hCBSDFT h ) {
    return h ? h->ns : 0;
}

int CBSDFT_size( hCBSDFT h ) {
    return h ? h->nwin : 0;
}

int CBSDFT_nfreq(  hCBSDFT h ) {
    return h ? h->nwin/2+1 : 0;
}

int CBSDFT_push( hCBSDFT h, const segy* const tr ) {
    if (h) {
        if (tr) {
            h->intr = (h->intr + 1)%h->ntr;
            SDFT(h->sdftH, h->window, (float*) tr->data, h->specdata[h->intr]);
            memcpy( (void*)&(h->hdrs[h->intr]), (void*) tr, HDRBYTES );
            h->outtr = (h->trcount <= h->ntr/2)? h->outtr : (h->outtr + 1)%h->ntr;
            h->trcount = (h->trcount < h->ntr)? h->trcount+1 : h->ntr;
        } else {
            h->trcount = (h->trcount>0)? h->trcount-1 : 0;
            h->outtr = (h->outtr + 1)%h->ntr;
        }
    } else
        err("bad pointer in CBSDFT_push.");
    return h ? h->trcount > h->ntr/2 : 0;
}

int CBSDFT_getSlice( hCBSDFT h, int isample, int ifreq, complex* const data ) {
    if (h && data) {
        if (h->trcount >= h->ntr/2) {
            int spos = h->intr - h->trcount + 1;
            spos = (spos<0)? spos+h->ntr : spos;
            for (int itrc=0; itrc<h->trcount; itrc++ )
                data[itrc] = h->specdata[(spos+itrc)%h->ntr][ifreq][isample];
            spos = (spos > h->outtr)? spos-h->ntr : spos;
            return (h->trcount < h->ntr)? h->outtr - spos : h->ntr/2;
        } else {
            warn("trace buffer too empty in CBSDFT_getSlice.");
        }
    } else
        err("bad pointer in CBSDFT_getSlice.");
    return 0;
}

void CBSDFT_setResult( hCBSDFT h, int isample, int ifreq, complex data ) {
    if (h)
        h->resbuf[ifreq][isample] = data;
    else
        err("bad pointer in CBSDFT_setResult");
}

void CBSDFT_getResult( hCBSDFT h, segy* const tr ) {
    if (h && tr) {
        ISDFT( h->sdftH, h->resbuf, tr->data );
        memcpy( (void*)tr, (void*)&(h->hdrs[h->outtr]), HDRBYTES );
    } else
        err("bad pointer in CBSDFT_getResult");
    
}
