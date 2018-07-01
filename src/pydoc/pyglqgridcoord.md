# GLQGridCoord

Compute the latitude and longitude coordinates used in Gauss-Legendre quadrature grids.

# Usage

`latglq`, `longlq` = GLQGridCoord (`lmax`)

# Returns

`latglq` : float, dimension (`lmax`+1)
:   The latitude coordinates of a grid, corresponding to the indices [:,0], in degrees.

`longlq` : float, dimension (2\*`lmax`+1)
:   The longitude coordinates of a grid, corresponding to the indices [0,:], in degrees. The first node is 0 E.

# Parameters

`lmax` : integer
:   The maximum spherical harmonic degree that will be integrated exactly by Gauss-Legendre quadrature.

# Description

`GLQGridCoord` will compute the latitude and longitude coordinates that are used in Gauss-Legendre quadrature grids for performing spherical harmonic transforms and reconstructions. The latitudinal nodes correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees.

# See also

[shglq](pyshglq.html), [shexpandglq](pyshexpandglq.html), [makegridglq](pymakegridglq.html), [shexpandglqc](pyshexpandglqc.html), [makegridglqc](pymakegridglqc.html)
