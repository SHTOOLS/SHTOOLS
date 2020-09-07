# SHMultiply
 
Multiply two functions and determine the spherical harmonic coefficients of the resulting function.

# Usage

call SHMultiply (`cilmout`, `cilm1`, `lmax1`, `cilm2`, `lmax2`, `precomp`, `norm`, `csphase`, `exitstatus`)

# Parameters

`cilmout` : output, real(dp), dimension (2, `lmax1`+`lmax2`+1, `lmax1`+`lmax2`+1)
:   The real spherical harmonic coefficients corresponding to the multiplication of `cilm1` and `cilm2` in the space domain.

`cilm1` : input, real(dp), dimension (2, `lmax1`+1, `lmax1`+1)
:   The spherical harmonic coefficients of the first function.

`lmax1` : input, integer
:   The maximum spherical harmonic degree used in evaluting `cilm1`.

`cilm2` : input, real(dp), dimension (2, `lmax2`+1, `lmax2`+1)
:   The spherical harmonic coefficients of the second function.

`lmax2` : input, integer
:   The maximum spherical harmonic degree used in evaluting `cilm2`.

`precomp` : input, optional, integer, default = 0
:   If 1, the array of Legendre functions `plx` will be precomputed on the Gauss-Legendre quadrature nodes.

`norm` : input, optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHMultiply` will take two sets of spherical harmonic coefficients, multiply the functions in the space domain, and expand the resulting field in spherical harmonics using `SHExpandGLQ`. The spherical harmonic bandwidth of the resulting field is `lmax1+lmax2`, where `lmax1` and `lmax2` are the bandwidths of the input fields. If the optional parameter `precomp` is set, then the array of Legendre functions `plx` will be precomputed on the Gauss-Legendre quadrature nodes.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

# See also

[shexpandglq](shexpandglq.html), [makegridglq](makegridglq.html), [shglq](shglq.html)
