# DHaj

Compute the latitudinal weights used in the Driscoll and Healy (1994) spherical harmonic transform.

# Usage

call DHaj (`n`, `aj`, `exitstatus`)

# Parameters

`n` : input, integer
:   The number of samples in latitude used in the spherical harmonic transform. This must be even.

`aj` : output, real/*8
:   The latitudinal weights used in the spherical harmonic transform.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`DHaj` will calculate the latitudinal weights used in the spherical harmonic transform of Driscoll and Healy (1994; equation 9). The number of samples `n` must be even, and the transform and its inverse are implemented as `SHExpandDH` and `MakeGridDH`, respectively. It is noted that the first element, corresponding to the north pole, is always zero. The element corresponding to the south pole is not included.

# Reference

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

# See also

[shexpanddh](shexpanddh.html), [makegriddh](makegriddh.html), [shexpanddhc](shexpanddhc.html), [makegriddhc](makegriddhc.html)
