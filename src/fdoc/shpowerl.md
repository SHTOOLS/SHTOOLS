# SHPowerL

Compute the power of a real function for a single spherical harmonic degree.

# Usage

`power` = SHPowerL (`cilm`, `l`)

# Parameters

`power` : output, real\*8
:   Power of the function for a single spherical harmonic degree `l`.

`cilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The spherical harmonic coefficients of the function.

`l` : input, integer
:   The spherical harmonic degree. This must be less than or equal to `lmaxin`.

# Description

`SHPowerL` will calculate the power of a function expressed in real 4-pi normalized spherical harmonics for a single degree `l`. This is explicitly calculated as:

`power = Sum_{i=1}^2 Sum_{m=0}^l cilm(i, l+1, m+1)**2`.

# See also

[shpowerdensityl](shpowerdensityl.html), [shcrosspowerl](shcrosspowerl.html), [shcrosspowerdensityl](shcrosspowerdensityl.html), [shpowerspectrum](shpowerspectrum.html), [shpowerspectrumdensity](shpowerspectrumdensity.html), [shcrosspowerspectrum](shcrosspowerspectrum.html), [shcrosspowerspectrumdensity](shcrosspowerspectrumdensity.html), [shadmitcorr](shadmitcorr.html)
