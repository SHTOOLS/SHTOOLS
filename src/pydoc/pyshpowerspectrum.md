# SHPowerSpectrum 

Compute the power spectrum of a real function.

# Usage

`pspectrum` = pyshtools.SHPowerSpectrum (`cilm`, [`lmax`])

# Returns

`pspectrum` : float, dimension (`lmax`+1)
:   The power spectrum of the function.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The function expressed in real spherical harmonics.
	
`lmax` : integer, optional, default = `lmaxin`
:   The maximum spherical harmonic degree used in calculating the power spectrum. This must be less than or equal to `lmaxin`.

# Description

`SHPowerSpectrum` will calculate the power spectrum of a function expressed in real spherical harmonics. For a given spherical harmonic degree `l`, this is calculated explicitly as:

`pspectrum(l) = Sum_{i=0}^1 Sum_{m=0}^l cilm[i, l, m]**2`.

# See also

[shpowerl](pyshpowerl.html), [shpowerdensityl](pyshpowerdensityl.html), [shcrosspowerl](pyshcrosspowerl.html), [shcrosspowerdensityl](pyshcrosspowerdensityl.html), [shpowerspectrumdensity](pyshpowerspectrumdensity.html), [shcrosspowerspectrum](pyshcrosspowerspectrum.html), [shcrosspowerspectrumdensity](pyshcrosspowerspectrumdensity.html), [shadmitcorr](pyshadmitcorr.html)
