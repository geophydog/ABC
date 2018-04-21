:hotel: [Return to Home Page](https://github.com/geophydog/geophydog.github.io/blob/master/README.md)

***

## ABC: AmBient noise and Coda Cross-correlation
- To extract empirical green's funcions from seismic ambient noise with seismic interferomatry.

***

### :one: raw data
![raw](https://github.com/geophydog/ABC/blob/master/images/raw-sac.jpg)

***

### :two: cut data
![cut](https://github.com/geophydog/ABC/blob/master/images/cut.jpg)

***

### :three: band-pass filtering
![bp](https://github.com/geophydog/ABC/blob/master/images/bp.jpg)

***

### :four: normalization in time domain
![norm](https://github.com/geophydog/ABC/blob/master/images/norm.png)

***

### :five: spectral whitening
#### results in time domain  
  - ![whi](https://github.com/geophydog/ABC/blob/master/images/whi.jpg)
#### results in frequency domain
  - ![whi_fft](https://github.com/geophydog/ABC/blob/master/images/whi_fft.png)
  
***

### :six: cross-correlation
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
