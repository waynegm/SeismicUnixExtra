/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.           */


#ifndef SUX_H
#define SUX_H

#include "su.h"
#include "segy.h"
#include "header.h"

typedef enum WinType { None, Hann, Hamming, Blackman } sux_Window;

/* Sliding Discrete Fourier Transform */
typedef struct _SDFT *hSDFT;
hSDFT SDFT_init( int nwin, int nsamples);
void SDFT( hSDFT h, sux_Window window, float* data, complex** result );
void ISDFT( hSDFT h, complex** specdata, float* result );
void SDFT_window( hSDFT h, sux_Window window, complex** specdata );
void SDFT_free( hSDFT h );

/* Sliding Discrete Cosine Transform */
typedef struct _SDCT *hSDCT;
hSDCT SDCT_init( int nwin, int nsamples );
void SDCT( hSDCT handle, sux_Window window, float* data, float** result );
void ISDCT( hSDCT handle, float** specdata, float* result );
void SDCT_window( hSDCT handle, sux_Window window, float** specdata );
void SDCT_free( hSDCT handle );

/* Ordered Trace buffer for seg Y trace data */
typedef struct _OTB *hOTB;
hOTB OTB_init( int ntraces, int nsamples );
int OTB_traces( hOTB h );
int OTB_push( hOTB h, const segy* const tr );
void OTB_copyCurrentHdr( hOTB h, segy* const tr );
int OTB_getSlice( hOTB h, int isample, float* const data );
const float** const OTB_getData( hOTB h );
void OTB_free( hOTB h );

/* Cyclic Trace buffer for seg Y trace data */
typedef struct _CTB *hCTB;
hCTB CTB_init( int ntraces, int nsamples );
int CTB_traces( hCTB h );
int CTB_push( hCTB h, const segy* const tr );
void CTB_copyCurrentHdr( hCTB h, segy* const tr );
int CTB_getSlice( hCTB h, int isample, float* const data );
const float** const CTB_getData( hCTB h );
void CTB_free( hCTB h );

/* Circular buffer for multi-trace sliding discrete fourier transform */
/*
typedef struct _CBSDFT *hCBSDFT;
hCBSDFT CBSDFT_init( int ntraces. int nsamples, int nwin, suxWindow window );
size_t CBSDFT_traces( hCBSDFT h );
size_t CBSDFT_samples( hCBSDFT h );
size_t CBSDFT_size( hCBSDFT h );
size_t CBSDFT_nfreq( hCBSDFT h );
void   CBSDFT_push( hCBSDFT h, float* data );
void   CBSDFT_getframe(  hCBSDFT h, int isample, int ifreq, complex* data );
void   CBSDFT_setresult( hCBSDFT h, int isample, int ifreq, complex data );
void   CBSDFT_getresult( hCBSDFT h, float* data );
void CBSDFT_free( hCBSDFT );
*/
#endif /* end of SUX_H */

