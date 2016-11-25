# SHCrossPowerDensityL

Compute the cross-power spectral density of two real functions for a single spherical harmonic degree.

# Usage

`cpsd` = SHCrossPowerDensityL (`cilm1`, `cilm2`, `l`)

# Parameters

`cpsd` : output, real\*8
:   Cross-power spectral density of the two real functions for spherical harmonic degree `l`.

`cilm1` : input, real\*8, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The spherical harmonic coefficients of the first real function.

`cilm2` : input, real\*8, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The spherical harmonic coefficients of the second real function.

`l` : input, integer
:   The spherical harmonic degree. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

# Description

`SHCrossPowerDensityL` will calculate the cross-power spectral density of two real functions expressed in real 4-pi normalized spherical harmonics for a single spherical harmonic degree `l`. This is explicitly calculated as:

`cpsd = Sum_{i=1}^2 Sum_{m=0}^l cilm1(i, l+1, m+1) * cilm2(i, l+1, m+1) / (2l + 1)`.

# See also

[shpowerl](shpowerl.html), [shpowerdensityl](shpowerdensityl.html), [shcrosspowerl](shcrosspowerl.html), [shpowerspectrum](shpowerspectrum.html), [shpowerspectrumdensity](shpowerspectrumdensity.html), [shcrosspowerspectrum](shcrosspowerspectrum.html), [shcrosspowerspectrumdensity](shcrosspowerspectrumdensity.html), [shadmitcorr](shadmitcorr.html)
