# SHMultiply

Multiply two functions and determine the spherical harmonic coefficients of the resulting function.

# Usage

`shout` = SHMultiply (`sh1`, `sh2`, [`lmax1`, `lmax2`, `norm`, `csphase`])

# Returns

`shout` : float, dimension (2, `lmax1`+`lmax2`+1, `lmax1`+`lmax2`+1)
:   The real spherical harmonic coefficients corresponding to the multiplication of `sh1` and `sh2` in the space domain.

# Parameters

`sh1` : float, dimension (2, `lmax1in`+1, `lmax1in`+1)
:   The spherical harmonic coefficients of the first function.

`sh2` : float, dimension (2, `lmax2in`+1, `lmax2in`+1)
:   The spherical harmonic coefficients of the second function.

`lmax1` : integer, optional, default = `lmax1in`
:   The maximum spherical harmonic degree used in evaluting `sh1`.

`lmax2` : integer, optional, default = `lmax2in`
:   The maximum spherical harmonic degree used in evaluting `sh2`.

`norm` : optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

# Description

`SHMultiply` will take two sets of spherical harmonic coefficients, multiply the functions in the space domain, and expand the resulting field in spherical harmonics using `SHExpandGLQ`. The spherical harmonic bandwidth of the resulting field is `lmax1+lmax2`, where `lmax1` and `lmax2` are the bandwidths of the input fields.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

# See also

[shexpandglq](pyshexpandglq.html), [makegridglq](pymakegridglq.html), [shglq](pyshglq.html)
