# SHPowerSpectrumDensity

Compute the power spectral density of a real function.

# Usage

call SHPowerSpectrumDensity (`cilm`, `lmax`, `pspectrum`)

# Parameters

`cilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real function expressed in real spherical harmonics.
	
`lmax` : input, integer
:   The maximum spherical harmonic degree used in calculating the power spectrum. This must be less than or equal to `lmaxin`.

`pspectrum` : output, real\*8, dimension (`lmax`+1)
:   The power spectral density of the function.

# Description

`SHPowerSpectrumDensity` will calculate the power spectral density of a function expressed in real spherical harmonics. For a given spherical harmonic degree `l`, this is explicitly calculated as:

`pspectrum(l) = Sum_{i=1}^2 Sum_{m=0}^l cilm(i, l+1, m+1)**2 / (2l + 1)`.

# See also

[shpowerl](shpowerl.html), [shpowerdensityl](shpowerdensityl.html), [shcrosspowerl](shcrosspowerl.html), [shcrosspowerdensityl](shcrosspowerdensityl.html), [shpowerspectrum](shpowerspectrum.html), [shcrosspowerspectrum](shcrosspowerspectrum.html), [shcrosspowerspectrumdensity](shcrosspowerspectrumdensity.html), [shadmitcorr](shadmitcorr.html)
