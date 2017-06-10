# SUSDCT 
Time-frequency decomposition using the sliding discrete cosine transform 
 
## Usage 
   susdct < stdin > stdout 
 
### Required Parameters 
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| dt=       | time sampling interval (sec)                    | trcheader     |
                                                                               
### Optional Parameters                                                        
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| nwin=     | number (odd) of samples in SDFT window          | 31            |
| window=   | none - no window applied to each data segment   | none          |
|           | hann - Hann window                              |               |
|           | hamming - Hamming window                        |               |
|           | blackman - Blackman window                      |               |
| verbose=  | 0 - no advisory messages, 1 - for messages      | 0             |
 
## Notes 
This process calculates a time-frequency decomposition of seismic data using 
the sliding discrete cosine transform. 
 
## Examples 
   suvibro | susdct | suximage 
 
