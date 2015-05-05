# SHMTCouplingMatrix

Returns the spherical harmonics coupling matrix for a set of tapers

#USAGE

`Mmt` = SHMTCouplingMatrix(`lmax`,`tapers`,[`lwin`,`k`,`taper_wt`])

#Returns

`Mmt` : float, dimension (lmax+1,lmax+lwin+1)
:   The spherical harmonics coupling matrix that relates input power (dim 1) with output power (dim 2)

#Parameters

`lmax`: integer
:   The input power bandwidth (matrix size)

`tapers` : float, dimension (`lwinin`+1, `kin`)
:   An array of the k windowing functions, arranged in columns, obtained from a call to `SHReturnTapers`. 

`lwin` : optional, integer, default = `lwinin`
:   The spherical harmonic bandwidth of the windowing functions in the array `tapers`.

`k` : optional, integer, default = `kin`
:   The number of tapers utilized in the multitaper spectral analysis.

`taper_wt` : optional, float, dimension (`kin`)
:   The weights used in calculating the multitaper spectral estimates. Optimal values of the weights (for a known global power spectrum) can be obtained from the routine `SHMTVarOpt`.


# Description

I<SHMTCouplingMatrix> returns the spherical harmonics coupling matrix, as
described in Wieczorek and Simons (2007). It maps input power into expected
output power after windowing, under the assumption that window and model
coefficients are uncorrelated.

# References

Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed., Cambridge Univ. Press, Cambridge, UK, 1992.

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, 665-692, doi:10.1007/s00041-006-6904-1, 2007.

# See also

[shmultitaperse](pyshmultitaperse.html), [shmultitapercse](pyshmultitapercse.html), [shreturntapers](pyshreturntapers.html), [shmtvaropt](pyshmtvaropt.html), [shmtdebias](pyshmtdebias.html)
