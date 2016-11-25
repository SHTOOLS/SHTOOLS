# SHPowerDensityL

Compute the power spectral density of a real function for a single spherical harmonic degree.

# Usage

`psd` = SHPowerDensityL (`cilm`, `l`)

# Parameters

`psd` : output, real\*8
:   Power spectral density of the function for the spherical harmonic degree `l`.

`cilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The function expressed in real spherical harmonics.

`l` : input, integer
:   The spherical harmonic degree. This must be less than or equal to `lmaxin`.

# Description

`SHPowerDensityL` will calculate the power spectral density of a function expressed in real 4-pi normalized spherical harmonics for a single degree `l`. This is explicitly calculated as:

`psd = Sum_{i=1}^2 Sum_{m=0}^l cilm(i, l+1, m+1)**2 / (2l + 1)`.

# See also

[shpowerl](shpowerl.html), [shcrosspowerl](shcrosspowerl.html), [shcrosspowerdensityl](shcrosspowerdensityl.html), [shpowerspectrum](shpowerspectrum.html), [shpowerspectrumdensity](shpowerspectrumdensity.html), [shcrosspowerspectrum](shcrosspowerspectrum.html), [shcrosspowerspectrumdensity](shcrosspowerspectrumdensity.html), [shadmitcorr](shadmitcorr.html)
