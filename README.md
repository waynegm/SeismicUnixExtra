# SeismicUnixExtra
My addons to the [Seismic Unix](http://www.cwp.mines.edu/cwpcodes/) seismic processing package. 

Access the [documentation](https://waynegm.github.io/SeismicUnixExtra/) for more information.

## Installation
  * Install and build Seismic Unix if you haven't already done so in a folder with write access for the user building this software.
  * Make sure the CWPROOT environment variable is set
  * Download and extract the GitHub repository
  * Change to the directory where the SeismicUnixExtra files were extracted
  * Type "make"

This should compile and install all programs into your seismic unix installation.

## Functionality

| Program | Description                                |
| ------- | -------------------------------------------|
| [susdft](docs/susdft.md) | Time-frequency decomposition by the sliding discrete fourier transform |
| [suisdft](docs/suisdft.md) | Inverse sliding discrete fourier tranform. |
| [susdct](docs/susdct.md) | Time-frequency decomposition using the sliding discrete cosine transform |
| [suisdct](docs/suisdct.md) | Inverse sliding discrete cosine tranform |
| [sutrcmedian](docs/sutrcmedian.md) | Rolling median filter over a panel of seismic traces by an ordered trace buffer |
| [suctrcmedian](docs/suctrcmedian.md) | Rolling median filter over a panel of seismic traces by a cyclic trace buffer |
| [susdft_denoise](docs/susdft_denoise.md) | Time-frequency denoise over a panel of seismic traces using the sliding discrete fourier transform |
| [sulpasmooth](docs/sulpasmooth.md) | Rolling LPA filter over a panel of seismic traces |
