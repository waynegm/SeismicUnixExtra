# susdct
Time-frequency of seismic data via the sliding discrete cosine transform 

## Usage
    susdct <stdin >sdout

Required parameters:

| Parameter        | Values | Description          | Default |
|:----------------:| ------ | -------------------- |:-------:|
| dt= || time sampling interval (sec) | from header |

Optional parameters:

| Parameter        | Values | Description | Default |
|:----------------:| ------ | ----------- |:-------:|
| nwin=            |        | number (odd) of samples in SDFT window | 31 |
| window=          | none<br>hann<br>hamming<br>blackman | window applied to each data segment | none |
| verbose= | 0<br>1 | no advisory messages<br>for advisory messages | 0 |

## Examples
    suvibro | susdct | suximage
