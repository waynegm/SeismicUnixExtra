/* Copyright (c) Wayne Mogg, 2017. */
/* All rights reserved.            */

/*********************** self documentation **********************/
/*************************************************************************
OTB - ordered buffer of SEG Y trace data - 

OTB_init             initialise and return a trace buffer handle
OTB_traces           return number of traces in the buffer
OTB_push             add a seg y trace to the buffer
OTB_copyCurrentHdr    get the header of the current centre trace 
OTB_getslice         get data at a particular sample for all traces
OTB_getData          return a pointer to the buffer data
OTB_free             release a trace buffer handle

************************************************************************** 
Function Prototypes:
hOTB OTB_init(int ntraces, int nsamples);
void OTB_free(hOTB h);
int OTB_traces(hOTB h);
int OTB_push(hOTB h, const seg* const tr);
void OTB_copyCurrentHdr(hOTB h, segy* const tr);
int OTB_getSlice(hOTB h, int isample, float* const data);
const float** const OTB_getData(hOTB h);

************************************************************************** 
OTB_init:
Input:
ntraces     number of traces in buffer - should be odd
nsamples    number of samples in each trace

Returned: trace buffer handle

************************************************************************** 
OTB_free:
Input:
h           trace buffer handle created by OTB_init

************************************************************************** 
OTB_traces:
Input:
h           trace buffer handle created by OTB_init

Returned:   number of traces currrently in the trace buffer

************************************************************************** 
OTB_push:
Input:
h           trace buffer handle created by OTB_init
tr          pointer to a seg Y trace to add to the buffer

Returned:   1 if sufficient number of traces in buffer for processing (ie >ntraces/2),
            0 otherwise

************************************************************************** 
OTB_copyCurrentHdr:
Input:
h           trace buffer handle created by OTB_init

Output:
tr          header from the central trace in the buffer is copied to the 
            seg Y trace
            
************************************************************************** 
OTB_getSlice:
Input:
h           trace buffer handle created by OTB_init
isample     trace sample for te data slice required

Output:
data        float array with the buffer data at the given sample

Returned:   position in data of the central trace

************************************************************************** 
OTB_getData:
Input:
h           trace buffer handle created by OTB_init

Returned:   pointer to 2D array holding the trace buffer data

************************************************************************** 
Notes:
Encapsulates the logic of a rolling window of traces over a panel of data.
For a usage example have a look at sutrcmedian.

The order of traces added to the buffer is preserved in the output of the
OTB_getSlice and OTB_getData. This convenience comes with the overhead of 
copying all traces in the buffer every time a new trace is added.

************************************************************************** 
Author: Wayne Mogg
**************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"
#include "par.h"
#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

typedef struct {
    unsigned char hdr[HDRBYTES];
} _HDR;

struct _OTB {
    int ntr;
    int ftr;
    int ltr;
    int ns;
    _HDR* hdrs;
    float** data;
};

hOTB OTB_init( int ntraces, int nsamples ) {
    
    hOTB h = emalloc(sizeof(struct _OTB));
    h->ntr = ntraces;
    h->ltr = ntraces-1;
    h->ftr= ntraces;
    h->ns = nsamples;
    h->data = ealloc2float( nsamples, ntraces );
    h->hdrs = ealloc1(ntraces, sizeof(_HDR));
    return h;
}

void OTB_free( hOTB h ) {
    if ( h ) {
        if (h->data) free2float( h->data );
        if (h->hdrs) free1( h->hdrs );
        free( h );
        h = 0;
    } else
        err("bad pointer in OTB_free.");
}

int OTB_traces( hOTB h ) {
    return h ? h->ltr - h->ftr + 1: 0;
}

int OTB_push( hOTB h, const segy* const tr ) {
    if (h) {
        int ntr = h->ntr;
        int ns = h->ns;
        float** data = h->data;
        h->ftr = (h->ftr > 0)? h->ftr-1 : 0;
        memcpy( (void*)&(data[0][0]), (void*)&(data[1][0]), (ntr-1)*ns*FSIZE );
        memcpy( (void*) &(h->hdrs[0]), (void*) &(h->hdrs[1]), (ntr-1)*HDRBYTES );
        if (tr) {
            memcpy( (void*)&(data[ntr-1][0]), (void*) tr->data, ns*FSIZE );
            memcpy( (void*)&(h->hdrs[ntr-1]), (void*) tr, HDRBYTES );
        } else
            h->ltr--;
//            h->ltr = (h->ltr>ntr/2)? h->ltr-1 : h->ltr;
    } else
        err("bad pointer in OTB_push.");
    return h->ltr >= h->ntr/2 && h->ftr <= h->ntr/2;
}

void OTB_copyCurrentHdr( hOTB h, segy* const tr ) {
    if ( h && tr ) {
        memcpy( (void*)tr, (void*)&(h->hdrs[h->ntr/2]), HDRBYTES );
    } else
        err("bad pointer in OTB_copyCurrentHdr");
}

int OTB_getSlice( hOTB h, int isample, float* const data ) {
    if (h && data) {
        if (h->ltr>=h->ntr/2) {
            int nt = OTB_traces(h);
            for (int itrc=0; itrc<nt; itrc++ )
                data[itrc] = h->data[itrc+h->ftr][isample];
            return h->ntr/2-h->ftr;
        } else {
            warn("trace buffer too empty in OTB_getSlice.");
            return 0;
        }
    } else
        err("bad pointer in OTB_getSlice.");
}

const float** const OTB_getData( hOTB h ) {
    if (h)
        return (const float** const) h->data;
    else
        err("bad pointer in OTB_getData.");
}

