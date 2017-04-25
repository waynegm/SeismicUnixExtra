/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.           */


#ifndef SUX_H
#define SUX_H

typedef enum WinType { None, Hann, Hamming, Blackman } sux_Window;

/* Sliding Discrete Fourier Transform */
void initSDFT( int nwin, complex* cfactors );
void SDFT( int ns, int nwin, float* data, complex* cfactors, complex** result );
void ISDFT( int ns, int nwin, complex** specdata, float* result );

/* Sliding Discrete Cosine Transform */
typedef struct _SDCT *hSDCT;
hSDCT initSDCT( int nwin, int nsamples );
void SDCT( hSDCT handle, sux_Window window, float* data, float** result );
void ISDCT( hSDCT handle, float** specdata, float* result );
void windowSDCT( hSDCT handle, sux_Window window, float** specdata );
void destroySDCT( hSDCT handle );

#endif /* end of SUX_H */

