/* Copyright (c) Wayne Mogg, 2017. */
/* All rights reserved.            */

/*********************** self documentation **********************/
/*************************************************************************
CTB - cyclic buffer of SEG Y trace data - 

CTB_init             initialise and return a trace buffer handle
CTB_traces           return number of traces in the buffer
CTB_push             add a seg y trace to the buffer
CTB_copyCurrentHdr    get the header of the current centre trace 
CTB_getslice         get data at a particular sample for all traces
CTB_getData          return a pointer to the buffer data
CTB_free             release a trace buffer handle

************************************************************************** 
Function Prototypes:
hCTB CTB_init(int ntraces, int nsamples);
void CTB_free(hCTB h);
int CTB_traces(hCTB h);
int CTB_push(hCTB h, const seg* const tr);
void CTB_copyCurrentHdr(hCTB h, segy* const tr);
int CTB_getSlice(hCTB h, int isample, float* const data);
const float** const CTB_getData(hOTB h);

************************************************************************** 
CTB_init:
Input:
ntraces     number of traces in buffer - should be odd
nsamples    number of samples in each trace

Returned: trace buffer handle

************************************************************************** 
CTB_free:
Input:
h           trace buffer handle created by CTB_init

************************************************************************** 
CTB_traces:
Input:
h           trace buffer handle created by CTB_init

Returned:   number of traces currrently in the trace buffer

************************************************************************** 
CTB_push:
Input:
h           trace buffer handle created by CTB_init
tr          pointer to a seg Y trace to add to the buffer

Returned:   1 if sufficient number of traces in buffer for processing (ie >ntraces/2),
            0 otherwise

************************************************************************** 
CTB_copyCurrentHdr:
Input:
h           trace buffer handle created by CTB_init

Output:
tr          header from the central trace in the buffer is copied to the 
            seg Y trace
            
************************************************************************** 
CTB_getSlice:
Input:
h           trace buffer handle created by CTB_init
isample     trace sample for te data slice required

Output:
data        float array with the buffer data at the given sample

Returned:   position in data of the central trace

************************************************************************** 
CTB_getData:
Input:
h           trace buffer handle created by CTB_init

Returned:   pointer to 2D array holding the trace buffer data

************************************************************************** 
Notes:
Encapsulates the logic of a rolling window of traces over a panel of data.
For a usage example have a look at sutrcmedian.

This data structure uses a circular buffer to hold the trace data. This has
the advantage that adding a new trace to the buffer does not cause any copying.
The disadvantage is that the traces are not stored in the buffer in the order
they were added which can make some operations more difficult.

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

struct _CTB {
    int ntr;
    int intr;
    int outtr;
    int trcount;
    int ns;
    _HDR* hdrs;
    float** data;
};

hCTB CTB_init( int ntraces, int nsamples ) {
    
    hCTB h = emalloc(sizeof(struct _CTB));
    h->ntr = ntraces;
    h->trcount = 0;
    h->intr= -1;
    h->outtr = 0;
    h->ns = nsamples;
    h->data = ealloc2float( nsamples, ntraces );
    h->hdrs = ealloc1(ntraces, sizeof(_HDR));
    return h;
}

void CTB_free( hCTB h ) {
    if ( h ) {
        if (h->data) free2float( h->data );
        if (h->hdrs) free1( h->hdrs );
        free( h );
        h = 0;
    } else
        err("bad pointer in CTB_free.");
}

int CTB_traces( hCTB h ) {
     return h ? h->trcount : 0;
}

int CTB_push( hCTB h, const segy* const tr ) {
    if (h) {
        int ntr = h->ntr;
        int ns = h->ns;
        float** data = h->data;
        if (tr) {
            h->intr = (h->intr + 1)%ntr;
            memcpy( (void*)&(data[h->intr][0]), (void*) tr->data, ns*FSIZE );
            memcpy( (void*)&(h->hdrs[h->intr]), (void*) tr, HDRBYTES );
            h->outtr = (h->trcount <= ntr/2)? h->outtr : (h->outtr + 1)%ntr;
            h->trcount = (h->trcount < ntr)? h->trcount+1 : ntr;
        } else {
            h->trcount = (h->trcount>0)? h->trcount-1 : 0;
            h->outtr = (h->outtr + 1)%ntr;
        }
    } else
        err("bad pointer in CTB_push.");
    return h ? h->trcount > h->ntr/2 : 0;
}

void CTB_copyCurrentHdr( hCTB h, segy* const tr ) {
    if ( h && tr ) {
        memcpy( (void*)tr, (void*)&(h->hdrs[h->outtr]), HDRBYTES );
    } else
        err("bad pointer in CTB_copyCurrentHdr");
}

int CTB_getSlice( hCTB h, int isample, float* const data ) {
    if (h && data) {
        if (h->trcount >= h->ntr/2) {
            int spos = h->intr - h->trcount + 1;
            spos = (spos<0)? spos+h->ntr : spos;
//            if (isample==0) {
//                warn("CTB_getSlice - intr: %d outtr: %d trcount: %d spos: %d icur: %d", h->intr, h->outtr, h->trcount, spos, (h->trcount<h->ntr)? h->outtr : h->ntr/2);
//            }
            for (int itrc=0; itrc<h->trcount; itrc++ )
                data[itrc] = h->data[(spos+itrc)%h->ntr][isample];
            spos = (spos > h->outtr)? spos-h->ntr : spos;
            return (h->trcount < h->ntr)? h->outtr - spos : h->ntr/2;
        } else {
            warn("trace buffer too empty in CTB_getSlice.");

        }
    } else
        err("bad pointer in CTB_getSlice.");
    return 0;
}

const float** const CTB_getData( hCTB h ) {
    if (h)
        return (const float** const) h->data;
    else
        err("bad pointer in OTB_getData.");
    return 0;
}

