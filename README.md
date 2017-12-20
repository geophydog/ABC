## ABC: Ambient Noise Cross-correlation
- To extract empirical green's funcions from seismic ambient noise with seismic interferomatry.

***

## Denpendencies
- Linux or Mac OS platform;
- FFTW3 (Fast Foureier Transform in the West) ;

***

## Usage

abc_egf file.lst

- file.lst: sac1 sac2 year mon day hour min sec start0 cut_npts cor_name lag_time
- [x] sac1: SAC file1;
- [x] sac2: SAC file2;
- [x] year: year of seismic event happening;
- [x] mon: month of seismic event happening;
- [x] day: day of seismic event happening;
- [x] hour: hour of seismic event happening;
- [x] min: minute of seismic event happening;
- [x] sec: second of seismic event happening;
- [x] start0: time of after seismic event happening;
- [x] cut_npts: points number of cutting;
- [x] cor_name: name of cross-correlation;
- [x] lag_time: lag time of cross-correlation.
***

## Contribution
- Author: Xuping Feng
- Email: geophydogvon@gmail.com
