# SHPowerSpectrumDensityC

Compute the power spectral density of a complex function.

# Usage

call SHPowerSpectrumDensityC (`cilm`, `lmax`, `pspectrum`)

# Parameters

`cilm` : input, complex\*16, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex function expressed in complex spherical harmonics.
	
`lmax` : input, integer
:   The maximum spherical harmonic degree of the power spectral density. This must be less than or equal to `lmaxin`.

`pspectrum` : output, real\*8, dimension (`lmax`+1)
:   The power spectral density of the function.

# Description

`SHPowerSpectrumDensityC` will calculate the power spectral density of a complex function expressed in complex spherical harmonics. For a given spherical harmonic degree `l`, this is calculated as:

`pspectrum(l) = Sum_{i=1}^2 Sum_{m=0}^l | cilm(i, l+1, m+1) |**2 / (2l + 1)`.

# See also

[shpowerlc](shpowerlc.html), [shpowerdensitylc](shpowerdensitylc.html), [shcrosspowerlc](shcrosspowerlc.html), [shcrosspowerdensitylc](shcrosspowerdensitylc.html), [shpowerspectrumc](shpowerspectrumc.html), [shcrosspowerspectrumc](shcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](shcrosspowerspectrumdensityc.html)
