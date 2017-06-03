/* Copyright (c) Wayne Mogg, 2017.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sux.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                   ",
" SUSDFT_DENOISE - time-frequency denoise over a panel of seismic   ",
"                  traces                                           ",
"                                                                   ",
" susdft_denoise <stdin >sdout                                      ",
"                                                                   ",
" Required parameters:                                              ",
" dt=(from header)  time sampling interval (sec)                    ",
"                                                                   ",
" Optional parameters:                                              ",
" ntr=9             number (odd) of traces to compute median        ",
" nwin=31           number (odd) of samples in sdft window          ",
" window=none       window applied to each data segment             ",
"       =hann       Hann window                                     ",
"       =hamming    Hamming window                                  ",
"       =blackman   Blackman window                                 ",
" reject=10         percentage of traces to reject                  ",
" type=swmean       replace rejected trace with mean of accepted    ",
"     =swmedian     replace rejected trace with median of accepted  ",
"     =median       output median of accepted traces                ",
"     =mean         output mean of accepted traces                  ",
" mode=0            output filtered trace                           ",
"     =1            output estimated noise                          ",
" verbose=0         no advisory messages                            ",
"         1         for advisory messages                           ",
"                                                                   ",
" This process implements the following:                            ",
" 1. Time-Frequency transform of the trace using the sliding DFT    ",
" 2. For each time and frequency sample acrosss the panel of traces ",
"    a. Sort by amplitude and discard the top reject% of values     ",
"    b. For type=mean or median output the mean or median           ",
"       respectively of the kept traces.                            ",
"    c. For type=swmean or swmedian output the mean or median       ",
"       respectively of the kept traces only for rejected traces    ",
"       otherwise the original value is passed unchanged            ",
" 3. Inverse sliding DFT and output                                 ",
"                                                                   ",
NULL};

/* Author: Wayne Mogg, May 2017
 *
 * Trace header fields accessed: ns, trid
 */
/**************** end self doc ***********************************/

segy tr;
typedef enum ProcType { SwMean, SwMedian, Median, Mean } proc_Type; 

int main(int argc, char **argv)
{
    float dt;
    int nwin;
    int ntr;
    cwp_String window;
    sux_Window iwind = None;
    float reject;
    int nkeep;
    cwp_String type;
    proc_Type itype = SwMean;
    int mode;
    int verbose;

    int is, tcount, imed, ifreq,nfreq;
    int nsamples;
    cwp_Bool seismic;
    complex* specbuf;
    float*  ampbuf;
    int*    idxbuf;
    hCBSDFT cbsdftH;
	
// Initialize
	initargs(argc, argv);
	requestdoc(1);

/* Get info from first trace */ 
    nsamples = 0;
    while (gettr(&tr)) {
        seismic = ISSEISMIC(tr.trid);
        if (seismic) {
            nsamples = tr.ns;
            break;
        } else
            warn("skipping non-seismic trace with trid=%d", tr.trid);
    }
    if (nsamples==0) 
        err("zero length traces not allowed.");
    
// Get parameters
    if (!getparint("verbose", &verbose)) verbose=0;
    if (!getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
    if (!dt) err("trace dt field is zero and not getparred");
    if (!getparint("nwin", &nwin)) nwin = 31;
    if (nwin%2==0) {
        nwin++;
        if (verbose)
            warn("adjusting nwin to be odd, was %d now %d",nwin-1, nwin);
    }
    if (!getparint("ntr", &ntr)) ntr = 9;
    if (!ntr%2) {
        ntr++;
        if (verbose)
            warn("adjusting ntr to be odd, was %d now %d",ntr-1, ntr);
    }
    if (!getparfloat("reject", &reject)) reject=10.0;
    if (reject<=0 || reject>100) {
        warn("reject out of range 0-100, reset to 10");
        reject = 10;
    }
    if (!getparstring("type", &type)) type = "swmean";
    if      (STREQ(type, "swmedian")) itype = SwMedian;
    else if (STREQ(type, "median")) itype = Median;
    else if (STREQ(type, "mean")) itype = Mean;
    else if (!STREQ(type, "swmean")) 
        err("unknown type\"%s\", see self-doc", type);
    
    if (!getparint("mode", &mode)) mode = 0;

    if (!getparstring("window", &window)) window = "none";
    if      (STREQ(window, "hann")) iwind = Hann;
    else if (STREQ(window, "hamming")) iwind = Hamming;
    else if (STREQ(window, "blackman")) iwind = Blackman;
    else if (!STREQ(window, "none")) 
        err("unknown window=\"%s\", see self-doc", window);
    
// Set up cyclic SDFT buffer and work space
    specbuf = ealloc1complex( ntr );
    ampbuf = ealloc1float( ntr );
    idxbuf = ealloc1int( ntr );
    cbsdftH = CBSDFT_init( ntr, nsamples, nwin, iwind );
    nfreq = CBSDFT_nfreq(cbsdftH);
/* Main processing loop */
    do {
        seismic = ISSEISMIC(tr.trid);
        if (seismic) {
            if (CBSDFT_push(cbsdftH, &tr)) {
                tcount = CBSDFT_traces(cbsdftH);
                nkeep = NINT((float) tcount * (100-reject)/100);
                if (nkeep > tcount)
                    nkeep = tcount;
                imed = nkeep/2;
                float inv_nkeep = 1.0/(float)nkeep;
                complex outval = cmplx(0.0,0.0);
                for (is=0; is<nsamples; is++) {
                    for (ifreq=0; ifreq<nfreq; ifreq++) {
                        int icur = CBSDFT_getSlice( cbsdftH, is, ifreq, specbuf );
                        for (int i=0; i<tcount; i++) {
                            idxbuf[i] = i;
                            ampbuf[i] = rcabs(specbuf[i]);
                        }
                        qkifind( nkeep, tcount, ampbuf, idxbuf );
                        switch(itype) {
                            case (Mean):
                                outval = cmplx(0.0,0.0);
                                for (int i=0; i<nkeep; i++)
                                    outval = cadd(outval, specbuf[idxbuf[i]]);
                                outval = crmul( outval, inv_nkeep);
                                break;
                            case (Median):
                                outval = specbuf[idxbuf[imed]];
                                break;
                            case (SwMedian):
                                outval = specbuf[idxbuf[imed]];
                                for (int i=0; i<nkeep; i++) {
                                    if (idxbuf[i] == icur) {
                                        outval = specbuf[icur];
                                        break;
                                    }
                                }
                                break;
                            default:
                                outval = cmplx(0.0,0.0);
                                for (int i=0; i<nkeep; i++)
                                    outval = cadd(outval, specbuf[idxbuf[i]]);
                                outval = crmul( outval, inv_nkeep);
                                for (int i=0; i<nkeep; i++) {
                                    if (idxbuf[i] == icur) {
                                        outval = specbuf[icur];
                                        break;
                                    }
                                }
                        };
                        outval = (mode==1)? csub(specbuf[icur],outval): outval;
                        CBSDFT_setResult( cbsdftH, is, ifreq, outval );
                    }
                }
                CBSDFT_getResult( cbsdftH, &tr );
                puttr(&tr);
            }
        } else 
            if (verbose) warn("skipping non-seismic trace with trid=%d", tr.trid);
    } while (gettr(&tr));

/* Handle last traces in buffer */
    while(CBSDFT_push(cbsdftH, 0)) {
        tcount = CBSDFT_traces(cbsdftH);
        nkeep = NINT((float) tcount * (100-reject)/100);
        if (nkeep > tcount)
            nkeep = tcount;
        imed = nkeep/2;
        float inv_nkeep = 1.0/(float)nkeep;
        complex outval = cmplx(0.0,0.0);
        for (is=0; is<nsamples; is++) {
            for (ifreq=0; ifreq<nfreq; ifreq++) {
                int icur = CBSDFT_getSlice( cbsdftH, is, ifreq, specbuf );
                for (int i=0; i<tcount; i++) {
                    idxbuf[i] = i;
                    ampbuf[i] = rcabs(specbuf[i]);
                }
                qkifind( nkeep, tcount, ampbuf, idxbuf );
                switch(itype) {
                    case (Mean):
                        outval = cmplx(0.0,0.0);
                        for (int i=0; i<nkeep; i++)
                            outval = cadd(outval, specbuf[idxbuf[i]]);
                        outval = crmul( outval, inv_nkeep);
                        break;
                    case (Median):
                        outval = specbuf[idxbuf[imed]];
                        break;
                    case (SwMedian):
                        outval = specbuf[idxbuf[imed]];
                        for (int i=0; i<nkeep; i++) {
                            if (idxbuf[i] == icur) {
                                outval = specbuf[icur];
                                break;
                            }
                        }
                        break;
                    default:
                        outval = cmplx(0.0,0.0);
                        for (int i=0; i<nkeep; i++)
                            outval = cadd(outval, specbuf[idxbuf[i]]);
                        outval = crmul( outval, inv_nkeep);
                        for (int i=0; i<nkeep; i++) {
                            if (idxbuf[i] == icur) {
                                outval = specbuf[icur];
                                break;
                            }
                        }
                }
                outval = (mode==1)? csub(specbuf[icur],outval): outval;
                CBSDFT_setResult( cbsdftH, is, ifreq, outval );
            }
        }
        CBSDFT_getResult( cbsdftH, &tr );
        puttr(&tr);
    };

    free1(specbuf);
    free1float(ampbuf);
    free1int(idxbuf);
    CBSDFT_free( cbsdftH );

    return EXIT_SUCCESS;
}
