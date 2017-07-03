# SUVPEF 
Wiener predictive error filtering with spatially varying lag 
 
## Usage 
   suvpef < stdin > stdout  [optional parameters] 
 
### Required Parameters 
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| dt=       | time sampling interval (sec)                    | trcheader     |
 
### Optional Parameters 
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| key=      | header word specifying lag locations            | cdp           |
| xlag=     | array of lag locations as given by key          |               |
| lag=      | array of prediction filter lags in seconds      | dt            |
| len=      | operator length in seconds                      | tr.ns*dt/20   |
| pnoise=   | relative additive noise level                   | 0.001         |
| mincorr=  | start of autocorrelation window in seconds      | tmin          |
| maxcorr=  | end of autocorrelation window in seconds        | tmax          |
| verbose=  | 0 - no advisory messages, 1 - for messages      | 0             |
                                                                               
Trace header fields accessed: ns, dt                                           
Trace header fields modified: none                                             
                                                                               
## Notes                                                                       
- For fixed lag omit the key= and xlag= parameters and give a single lag= value. 
- For spatially variable lag specify key=, xlag= and lag= arrays.                
- Linear interpolation and constant extrapolation used to compute lag from the 
arrays.                                                                        
                                                                               
This is a simplified version of supef                                          
                                                                               
