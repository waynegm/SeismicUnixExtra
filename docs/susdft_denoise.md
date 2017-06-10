# SUSDFT_DENOISE 
Time-frequency denoise over a panel of seismic traces using the sliding discrete fourier transform 
 
## Usage 
   susdft_denoise < stdin > stdout 
 
### Required Parameters 
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| dt=       | time sampling interval (sec)                    | trcheader     |
                                                                               
### Optional Parameters                                                        
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| ntr=      | number (odd) of traces in analysis panel        | 9             |
| nwin=     | number (odd) of samples in SDFT window          | 31            |
| window=   | none - no window applied to each data segment   | none          |
|           | hann - Hann window                              |               |
|           | hamming - Hamming window                        |               |
|           | blackman - Blackman window                      |               |
| reject=   | percentage of high trace amplitudes to reject   | 10            |
| type=     | swmean - only replace rejected trace with mean of accepted  | swmean |
|           | swmedian - only replace rejected trace with median of accepted |     |
|           | median - output median of accepted traces       |                    |
|           | mean - output mean of accepted traces           |                    |
| mode=     | 0 - output filtered, 1 - output noise           | 0             |
| verbose=  | 0 - no advisory messages, 1 - for messages      | 0             |
 
## Notes 
This process implements the following algorithm: 
 
1. Time-Frequency transform of the trace using the sliding DFT 
2. For each time and frequency sample acrosss the panel of traces 
    - Sort by amplitude and discard the top reject% of values 
    - For type=mean or median output the mean or median respectively of the kept traces. 
    - For type=swmean or swmedian output the mean or median respectively of the kept traces       only if the current trace value was rejected, otherwise the original value is passed unchanged. 
3. Inverse sliding DFT and output 
 
