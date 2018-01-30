# SHGLQ

Precompute weights, nodes, and associated Legendre functions used in the Gauss-Legendre quadrature based spherical harmonics routines.

# Usage

call SHGLQ (`lmax`, `zero`, `w`, `plx`, `norm`, `csphase`, `cnorm`, `exitstatus`)

# Parameters

`lmax` : input, integer
:   The maximum spherical harmonic degree of the coefficients to be calculated in the Gauss-Legendre quadrature based spherical harmonic transform routines.

`zero` : output, real\*8, dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude, determined from a call to `PreGLQ`.

`w` : output, real\*8, dimension (`lmax`+1)
:   The weights used in the Gauss-Legendre quadrature over latitude, determined from a call to `PreGLQ`.

`plx` : output, optional, real\*8, dimension (`lmax`+1, (`lmax`+1)\*(`lmax`+2)/2)
:   An array of the associated Legendre functions calculated at the nodes used in the Gauss-Legendre quadrature. 

`norm` : input, optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`cnorm` : input, optional, integer, default = 0
:   If 0 (default), the real normalization of the associated Legendre functions will be used. If 1, the complex normalization of the associated Legendre functions will be used.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHGLQ` will calculate the weights and zeros used in the Gauss-Legendre quadrature based spherical harmonic routines `SHExpandGLQ`, `MakeGridGLQ`, `SHExpandGLQC`, and `MakeGridGLQC`. Optionally, an array of the associated Legendre functions evaluated on the quadrature nodes can be computed as well. If the complex routines are to be used, the optional parameter `cnorm` must be set equal to 1.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m. 

# See also

[shexpandglq](shexpandglq.html), [makegridglq](makegridglq.html), [shexpandglqc](shexpandglqc.html), [makegridglqc](makegridglqc.html), [glqgridcoord](glqgridcoord.html), [preglq](preglq.html)
