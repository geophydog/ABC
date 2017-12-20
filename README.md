## ABC: AmBient noise Cross-correlation
- To extract empirical green's funcions from seismic ambient noise with seismic interferomatry.

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
