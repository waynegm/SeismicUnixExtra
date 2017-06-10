# SUTRCMEDIAN 
Rolling median filter over a panel of seismic traces by an ordered trace buffer 
 
## Usage 
   sutrcmedian < stdin > stdout 
### Optional Parameters 
| Parameter | Description                                     | Default       |
|:---------:| ----------------------------------------------- |:-------------:|
| ntr=      | number (odd) of traces in filter panel          | 5             |
| mode=     | =0 output filtered trace, =1 output noise       | 0             |
| verbose=  | =0 no advisory messages, =1 for messages        | 0             |
 
## Notes 
This is primarily a demonstation and test platform for the ordered trace buffer implementation. 
 
