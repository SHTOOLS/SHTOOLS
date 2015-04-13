# SHCrossPowerSpectrum 

Compute the cross-power spectrum of two real functions.

# Usage

`cspectrum` = pyshtools.SHCrossPowerSpectrum (`cilm1`, `cilm2`, [`lmax`])

# Returns

`cspectrum` : float, dimension (`lmax`+1)
:   The cross-power spectrum of the two functions.

# Parameters

`cilm1` : float, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The first function expressed in real spherical harmonics.

`cilm2` : float, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The second function expressed in real spherical harmonics.
	
`lmax` : optional, integer, default = min(`lmaxin1`, `lmaxin2`)
:   The maximum spherical harmonic degree to calculate the cross-power spectrum. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

# Description

`SHCrossPowerSpectrum` will calculate the cross-power spectrum of two real functions expressed in real spherical harmonics. For a given degree spherical harmonic degree `l`, this is explicitly calculated as:

`cspectrum(l) = Sum_{i=0}^1 Sum_{m=0}^l cilm1[i, l, m] * cilm2[i, l, m]`.

# See also

[shpowerl](pyshpowerl.html), [shpowerdensityl](pyshpowerdensityl.html), [shcrosspower](pyshcrosspower.html), [shcrosspowerdensityl](pyshcrosspowerdensityl.html), [shpowerspectrum](pyshpowerspectrum.html), [shpowerspectrumdensity](pyshpowerspectrumdensity.html), [shcrosspowerspectrumdensity](pyshcrosspowerspectrumdensity.html), [shadmitcorr](pyshadmitcorr.html)
