/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.           */


#ifndef SUX_H
#define SUX_H

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

/* Circular buffer for real (float) arrays */
typedef struct _CBFLOAT *hCBFLOAT;
hCBFLOAT CBFLOAT_init( int ntraces, int nsamples );
size_t CBFLOAT_size(hCBFLOAT h);
size_t CBFLOAT_samples(hCBFLOAT h);
void CBFLOAT_push( hCBFLOAT h, float* data );
float* CBFLOAT_get( hCBFLOAT, int i ); 
void CBFLOAT_free( hCBFLOAT handle ); 

#endif /* end of SUX_H */

