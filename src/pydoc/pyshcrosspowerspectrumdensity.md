# SHCrossPowerSpectrumDensity

Compute the cross-power spectral density of two real functions.

# Usage

`cspectrum` = pyshtools.SHCrossPowerSpectrumDensity (`cilm1`, `cilm2`, [`lmax`])

# Returns

`cspectrum` : float, dimension (`lmax`+1)
:   The cross-power spectral density of the two functions.

# Parameters

`cilm1` : float, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The first function expressed in real spherical harmonics.

`cilm2` : float, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The second function expressed in real spherical harmonics.
	
`lmax` : optional, integer, default = min(`lmaxin1`, `lmaxin2`)
:   The maximum spherical harmonic degree to calculate the cross-power spectral density. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

# Description

`SHCrossPowerSpectrumDensity` will calculate the cross-power spectral density of two functions expressed in real spherical harmonics. For a given spherical harmonic degree `l`, this is explicitly calculated as:

`cspectrum(l) = Sum_{i=0}^1 Sum_{m=0}^l cilm1[i, l, m] * cilm2[i, l, m] / (2l + 1)`.

# See also

[shpowerl](pyshpowerl.html), [shpowerdensityl](pyshpowerdensityl.html), [shcrosspowerl](pyshcrosspowerl.html), [shcrosspowerdensityl](pyshcrosspowerdensityl.html), [shpowerspectrum](pyshpowerspectrum.html), [shpowerspectrumdensity](pyshpowerspectrumdensity.html), [shcrosspowerspectrum](pyshcrosspowerspectrum.html), [shadmitcorr](pyshadmitcorr.html)
