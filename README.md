:hotel: [Return to Home Page](https://github.com/geophydog/geophydog.github.io/blob/master/README.md)

***

# NOTICE!!! THIS PROGRAM COMES WITH NO WARRANTY, TO THE EXTEND PERMITED BY LAW!

***

## ABC: AmBient noise and Coda
- To extract empirical green's funcions(EGFs) from seismic ambient noise or coda with seismic interferomatry method.
- Actually, it's to compute cross-correlation time functios of one station-pair.

***

### :one: raw data, here 86400 seconds long (one day) and sampling rate of 5 Hz
![raw](https://github.com/geophydog/ABC/blob/master/images/raw-sac.jpg)

***

### :two: cut data from 10k senconds to 40k seconds
![cut](https://github.com/geophydog/ABC/blob/master/images/cut.jpg)

***

### :three: band-pass filtering from 50 seconds to 50 seconds
![bp](https://github.com/geophydog/ABC/blob/master/images/bp.jpg)

***

### :four: normalization in time domain with run-absolute-mean method
![norm](https://github.com/geophydog/ABC/blob/master/images/norm.png)

***

### :five: spectral whitening with run-absolute-mean method
#### results in time domain  
  - ![whi](https://github.com/geophydog/ABC/blob/master/images/whi.jpg)
#### results in frequency domain
  - ![whi_fft](https://github.com/geophydog/ABC/blob/master/images/whi_fft.png)
  
***

### :six: compute cross-correlation time function
![cor](https://github.com/geophydog/ABC/blob/master/images/cor.jpg)

***

## Denpendencies
- Linux or Mac OS platform;
- FFTW3 (Fast Foureier Transform in the West) ;

***

## Usage

abc_egf file.lst

- file.lst: sac1 sac2 year mon day hour min sec start0 cut_npts cor_name lag_time  

| parameter | meaning  |
| --------- | -------- |
|  sac1     | SAC file1|
|  sac2     | SAC file2|
|  year     | year of seismic event happening|
|  mon      | month of seismic event happening|
|  day      | day of seismic event happening|
|  hour     | hour of seismic event happening|
|  min      | minute of seismic event happening|
|  sec:     | second of seismic event happening|
|  start0   | time of after seismic event happening|
|  cut_npts | points number of cutting|
|  cor_name | name of cross-correlation|
|  lag_time | lag time of cross-correlation|

***

## Contribution
- Author: Xuping Feng
- Email: geophydogvon@gmail.com

***

## Reference
```
Bensen, G. D., et al. "Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements."  -
Geophysical Journal International 169.3 (2007): 1239-1260.
```
