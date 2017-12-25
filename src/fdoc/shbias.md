# SHBias

Calculate the (cross-)power spectrum expectation of a windowed function from its global spectrum.

# Usage

call SHBias (`shh`, `lwin`, `incspectra`, `ldata`, `outcspectra`, `save_cg`, `exitstatus`)

# Parameters

`shh` : input, real\*8, dimension (`lwin`+1)
:   The power spectrum of the localizing window.

`lwin` : input, integer
:   The spherical harmonic bandwidth of the localizing window.

`incspectra` : input, real\*8, dimension (`ldata`+1)
:   The global unwindowed (cross-)power spectrum.

`ldata` : input, integer
:   The maximum degree of the global unwindowed power spectrum.

`outcspectra` : output, real\*8, dimension (`ldata`+`lwin`+1)
:   The expectation of the localized (cross-)power spectrum.

`save_cg` : optional, input, integer, default = 0
:   If set equal to 1, the Clebsch-Gordon coefficients will be precomputed and saved for future use (if `lwin` or `ldata` change, this will be recomputed). To deallocate the saved memory, set this parameter equal to -1. If set equal to 0 (default), the Clebsch-Gordon coefficients will be recomputed for each call.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHBias` will calculate the (cross-)power spectrum expectation of a function multiplied by a localizing window. This is given by equation 35 of Wieczorek and Simons (2005) and equation 2.11 of Wieczorek and Simons (2007),

`<SFG> = Sum_{j=0}^L Shh Sum_{i=|l-j|}^{|l+j|} Sfg (C_{j0i0}^{l0})^2`

where `<SFG>` is the expectation of the localized (cross-)power spectrum, `Shh` is the power spectrum of the window bandlimited to degree `L`, `Sfg` is the global unwindowed (cross-)power spectrum, and `C` is a Clebsch-Gordan coefficient. The Clebsch-Gordan coefficients are calculated using a simple relationship to the Wigner 3-j symbols. The maximum calculated degree of the windowed power spectrum expectation corresponds to the smaller of (`ldata+lwin`) and `size(outcspectra)-1`. It is implicitly assumed that the power spectrum of `inspectrum` is zero beyond degree `ldata.`. If this is not the case, the ouput power spectrum should be considered valid only for the degrees up to and including `ldata` - `lwin`.

If this routine is to be called several times using the same values of `lwin` and `ldata`, then the Clebsch-Gordon coefficients can be precomputed and saved by setting the optional parameter `save_cg` equal to 1. To deallocate the saved memory, which is a matrix of size (`lwin+ldata,lwin,2*lwin+ldata+1`), set `save_cg` equal to -1.

# References

Wieczorek, M. A. and F. J. Simons, Localized spectral analysis on the sphere, 
Geophys. J. Int., 162, 655-675, doi:10.1111/j.1365-246X.2005.02687.x, 2005.

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on a sphere, J. Fourier Anal. Appl., 13, 665-692, doi:10.1007/s00041-006-6904-1, 2007.

# See also

[shpowerspectrum](shpowerspectrum.html),[shcrosspowerspectrum](shcrosspowerspectrum.html), [wigner3j](wigner3j.html), [shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html), [shbiasadmitcorr](shbiasadmitcorr.html)
