# SHMTCouplingMatrix

This routine returns the multitaper coupling matrix for a given set of power spectra of arbitrary localization windows. This matrix relates the global power spectrum to the expectation of the localized multitaper spectrum.
# Usage

`Mmt` = SHMTCouplingMatrix (`lmax`, `tapers_power`, [`lwin`, `k`, `taper_wt`])

# Returns

`Mmt` : float, dimension (`lmax`+`lwin`+1, `lmax`+1)
:   The full multitaper coupling matrix that relates the global power spectrum to the expectation of the localized multitaper spectrum.

# Parameters

`lmax` : integer
:   The spherical harmonic bandwidth of the global power spectrum.

`tapers_power` : float, dimension (`lwinin`+1, `kin`)
:   An array of power spectra of the k windowing functions, arranged in columns.

`lwin` : optional, integer, default = `lwinin`
:   The spherical harmonic bandwidth of the windowing functions in the array `tapers`.

`k` : optional, integer, default = `kin`
:   The number of tapers utilized in the multitaper spectral analysis.

`taper_wt` : optional, float, dimension (`kin`)
:   The weights used in calculating the multitaper spectral estimates. Optimal values of the weights (for a known global power spectrum) can be obtained from the routine `SHMTVarOpt`.

# Description

`SHMTCouplingMatrix` returns the multitaper coupling matrix that relates the global power spectrum (assumed to be stationary) to the expectation of the localized multitaper spectrum. This is given by eqs 4.5 and 4.6 in Wieczorek and Simons (2007):

`< S_{Phi Phi}^(mt) > = M^(mt) S_{ff}`

where `S_{Phi Phi}` is a vector containing the `lmax+lwin+1` localized multitaper power spectral estiamtes, `S_{ff}` is a vector of the global power spectrum up to degree `lmax`, and `< ... >` is the expectation operator. The coupling matrix is given explicitly by

`M_{ij} = Sum_{l=0}^L Sum_{k=1}^K a_k S_{hh}^{k}(l) [ C_{l0j0}^{i0} ]^2`

where `a_k` are the taper weights, `S_{hh}` is the power of the window, and `C` is a Clebsch-Gordon coefficient.

Note that this routine returns the "full" coupling matrix of dimension (`lmax` + `lwin` + 1, `lmax` + 1). When multiplied by a global input power spectrum with bandwidth `lmax`, it returns the output power spectrum with a bandwidth of `lmax` + `lwin`. In doing so, it is implicitly assumed that input power spectrum is exactly zero for all degrees greater than lmax. If this is not the case, the ouput power spectrum should be considered valid only for the degrees up to and including `lmax` - `lwin`.

# References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, 665-692, doi:10.1007/s00041-006-6904-1, 2007.

# See also

[shmultitaperse](pyshmultitaperse.html), [shmultitapercse](pyshmultitapercse.html), [shreturntapers](pyshreturntapers.html), [shmtvaropt](pyshmtvaropt.html), [shmtdebias](pyshmtdebias.html)
