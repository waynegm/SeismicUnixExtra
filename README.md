# SeismicUnixExtra
My addons to the [Seismic Unix](http://www.cwp.mines.edu/cwpcodes/) seismic processing package

## Functionality

### Transforms

| Program | Description                                |
| ------- | ------------------------------------------ |
| susdct  | Sliding discrete cosine transform          |
| suisdct | Inverse sliding discrete cosine transform  |
| susdft  | Sliding discrete fourier transform         |
| suisfdt | Inverse sliding discrete fourier transform | 

## Installation
  * Install and build Seismic Unix if you haven't already done so in a folder with write access for the user building this software.
  * Make sure the CWPROOT environment variable is set
  * Download and extract the GitHub repository
  * Change to the directory where the SeismicUnixExtra files were extracted
  * Type "make"

This should compile and install all programs into your seismic unix installation.
