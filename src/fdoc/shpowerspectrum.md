# SHPowerSpectrum

Compute the power spectrum of a real function.

# Usage

call SHPowerSpectrum (`cilm`, `lmax`, `pspectrum`, `exitstatus`)

# Parameters

`cilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The function expressed in real spherical harmonics.

`lmax` : input, integer
:   The maximum spherical harmonic degree used in calculating the power spectrum. This must be less than or equal to `lmaxin`.

`pspectrum` : output, real\*8, dimension (`lmax`+1)
:   The power spectrum of the function.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHPowerSpectrum` will calculate the power spectrum of a function expressed in real 4-pi normalized spherical harmonics. For a given spherical harmonic degree `l`, this is calculated explicitly as:

`pspectrum(l) = Sum_{i=1}^2 Sum_{m=0}^l cilm(i, l+1, m+1)**2`.

# See also

[shpowerl](shpowerl.html), [shpowerdensityl](shpowerdensityl.html), [shcrosspowerl](shcrosspowerl.html), [shcrosspowerdensityl](shcrosspowerdensityl.html), [shpowerspectrumdensity](shpowerspectrumdensity.html), [shcrosspowerspectrum](shcrosspowerspectrum.html), [shcrosspowerspectrumdensity](shcrosspowerspectrumdensity.html), [shadmitcorr](shadmitcorr.html)
