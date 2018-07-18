# SHGLQ

Precompute the weights and nodes used in the Gauss-Legendre quadrature based spherical harmonics routines.

# Usage

`zero`, `w` = SHGLQ (`lmax`)

# Returns

`zero` : float, dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude, determined from a call to `PreGLQ`.

`w` : float, dimension (`lmax`+1)
:   The weights used in the Gauss-Legendre quadrature over latitude, determined from a call to `PreGLQ`.

# Parameters

`lmax` : integer
:   The maximum spherical harmonic degree of the coefficients to be calculated in the Gauss-Legendre quadrature based spherical harmonic transform routines.

# Description

`SHGLQ` will calculate the weights and zeros used in the Gauss-Legendre quadrature based spherical harmonic routines `SHExpandGLQ`, `MakeGridGLQ`, `SHExpandGLQC`, and `MakeGridGLQC`.

# See also

[shexpandglq](pyshexpandglq.html), [makegridglq](pymakegridglq.html), [shexpandglqc](pyshexpandglqc.html), [makegridglqc](pymakegridglqc.html), [glqgridcoord](pyglqgridcoord.html)
