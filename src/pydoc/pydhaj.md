# DHaj

Compute the latitudinal weights used in the Driscoll and Healy (1994) spherical harmonic transform.

# Usage

`aj` = DHaj (`n`)

# Returns

`aj` : float, dimension (`n`)
:   The latitudinal weights used in the spherical harmonic transform.

# Parameters

`n` : integer
:   The number of samples in latitude used in the spherical harmonic transform. This must be even.


# Description

`DHaj` will calculate the latitudinal weights used in the spherical harmonic transform of Driscoll and Healy (1994; equation 9). The number of samples `n` must be even, and the transform and its inverse are implemented as `SHExpandDH` and `MakeGridDH`, respectively. It is noted that the first element, corresponding to the north pole, is always zero. The element corresponding to the south pole is not included.

# Reference

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

# See also

[shexpanddh](pyshexpanddh.html), [makegriddh](pymakegriddh.html), [shexpanddhc](pyshexpanddhc.html), [makegriddhc](pymakegriddhc.html)
