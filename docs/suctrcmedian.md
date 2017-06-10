# SUCTRCMEDIAN 
Rolling median filter over a panel of seismic traces by a cyclic trace buffer 
 
## Usage 
   suctrcmedian < stdin > stdout 
 
### Optional Parameters 
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| ntr=      | number (odd) of traces in filter panel          | 5             |
| mode=     | =0 output filtered trace, =1 output noise       | 0             |
| verbose=  | =0 no advisory messages, =1 for messages        | 0             |
 
## Notes 
This is primarily a demonstation and test platform for the cyclic trace buffer implementation. 
 
