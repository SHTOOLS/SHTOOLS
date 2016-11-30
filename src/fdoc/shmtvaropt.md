# SHMTVarOpt

Calculate the minimum variance and corresponding optimal weights of a localized multitaper spectral estimate.

# Usage

call SHMTVarOpt (`l`, `tapers`, `taper_order`, `lwin`, `kmax`, `sff`, `var_opt`, `var_unit`, `weight_opt`, `unweighted_covar`, `nocross`, `exitstatus`)

# Parameters

`l` : input, integer
:   The angular degree to determine the minimum variance and optimal weights.

`tapers` : input, real\*8, dimension (`lwin`+1, `kmax`)
:   A matrix of localization functions obtained from `SHReturnTapers` or `SHReturnTapersM`.

`taper_order` : input, integer, dimension (`kmax`)
:   The angular order of the windowing coefficients in TAPERS. If this matrix was created using `SHReturnTapersM`, then this array must be composed of zeros.

`lwin` : input, integer
:   The spherical harmonic bandwidth of the localizing windows.

`kmax` : input, integer
:   The maximum number of tapers to be used when calculating the minimum variance and optimal weights.

`sff` : input, real\*8, dimension (`l`+`lwin`+1)
:   The global unwindowed power spectrum of the function to be localized.

`var_opt` : output, real\*8, dimension (`kmax`)
:   The minimum variance of the multitaper spectral estimate for degree `l` using 1 through `kmax` tapers.

`var_unit` : output, real\*8, dimension (`kmax`)
:   The variance of the multitaper spectral estimate using equal weights for degree `l` using 1 through `kmax` tapers.

`weight_opt` : optional, output, real\*8, dimension (`kmax`, `kmax`)
:   The optimal weights (in columns) that minimize the multitaper spectral estimate's variance using 1 through `kmax` tapers.

`unweighted_covar` : optional, output, real\*8, dimension (`kmax`, `kmax`)
:   The unweighted covariance matrix of the `kmax` tapers (i.e., Fij in Wieczorek and Simons 2007).

`nocross` : optional, input, integer, default = 0
:   If 1, only the diagonal terms of the covariance matrix Fij will be computed. If 0, all terms will be computed.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHMTVarOpt` will determine the minimum variance that can be achieved by a weighted multitaper spectral analysis, as is described by Wieczorek and Simons (2007). The minimum variance is output as a function of the number of tapers utilized, from 1 to a maximum of `kmax`, and the corresponding variance using equal weights is output for comparison. The windowing functions are assumed to be solutions to the spherical-cap concentration problem, as determined by a call to `SHReturnTapers` or `SHReturnTapersM`. The minimum variance and weights are dependent upon the form of the global unwindowed power spectrum, `Sff`.

If the optional argument `weight_opt` is specified, then the optimal weights will be returned as a function of the number of tapers employed, from 1 to `kmax`. If `unweighted_covar` is specified, then the unweighted covariance matrix of the `kmax` tapers (i.e., Fij) will be output. If the optional argument `nocross` is set to 1, then only the diagnonal terms of `Fij` will be computed.

# References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

# See also

[shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html), [shmultitaperse](shmultitaperse.html), [shmultitapercse](shmultitapercse.html); [shlocalizedadmitcorr](shlocalizedadmitcorr.html), [shbiasadmitcorr](shbiasadmitcorr.html), [shbiask](shbiask.html), [shmtdebias](shmtdebias.html)
