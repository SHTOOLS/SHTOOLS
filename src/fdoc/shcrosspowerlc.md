# SHCrossPowerLC

Compute the cross-power of two complex functions for a single spherical harmonic degree.

# Usage

`cpower` = SHCrossPowerLC (`cilm1`, `cilm2`, `l`)

# Parameters

`cpower` : output, complex\*16
:   Cross power of the two complex functions for spherical harmonic degree `l`.

`cilm1` : input, complex\*16, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The first complex function expressed in complex spherical harmonics.

`cilm2` : input, complex\*16, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The second complex function expressed in complex spherical harmonics.

`l` : input, integer
:   The spherical harmonic degree. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

# Description

`SHCrossPowerLC` will calculate the spectral cross-power of two complex functions expressed in complex 4-pi normalized spherical harmonics for a single spherical harmonic degree `l`. This is calculated as:

`cpower = Sum_{i=1}^2 Sum_{m=0}^l cilm1(i, l+1, m+1) * conjg[cilm2(i, l+1, m+1)]`.

# See also

[shpowerlc](shpowerlc.html), [shpowerdensitylc](shpowerdensitylc.html), [shcrosspowerdensitylc](shcrosspowerdensitylc.html), [shpowerspectrumc](shpowerspectrumc.html), [shpowerspectrumdensityc](shpowerspectrumdensityc.html), [shcrosspowerspectrumc](shcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](shcrosspowerspectrumdensityc.html)
