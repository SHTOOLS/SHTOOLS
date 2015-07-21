# SHCrossPowerDensityLC 

Compute the cross-power spectral density of two complex functions for a single spherical harmonic degree.

# Usage

`cpsd` = pyshtools.SHCrossPowerDensityLC (`cilm1`, `cilm2`, `l`)

# Returns

`cpsd` : complex
:   The cross-power spectral density of the two complex functions for spherical harmonic degree `l`.

# Parameters

`cilm1` : complex, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The first complex function expressed in complex spherical harmonics.

`cilm2` : complex, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The second complex function expressed in complex spherical harmonics.
	
`l` : integer
:   The spherical harmonic degree. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

# Description

`SHCrossPowerDensityLC` will calculate the cross-power spectral density of two complex functions expressed in complex spherical harmonics for a single spherical harmonic degree `l`. This is calculated as:

`cpsd = Sum_{i=0}^1 Sum_{m=0}^l cilm1[i, l, m] * conjg[cilm2[i, l, m]] / (2l + 1)`.

# See also 

[shpowerlc](pyshpowerlc.html), [shpowerdensitylc](pyshpowerdensitylc.html), [shcrosspowerlc](pyshcrosspowerlc.html), [shpowerspectrumc](pyshpowerspectrumc.html), [shpowerspectrumdensityc](pyshpowerspectrumdensityc.html), [shcrosspowerspectrumc](pyshcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](pyshcrosspowerspectrumdensityc.html)
