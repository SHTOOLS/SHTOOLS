# GLQGridCoord

Compute the latitude and longitude coordinates used in Gauss-Legendre quadrature grids.

# Usage

call GLQGridCoord (`latglq`, `longlq`, `lmax`, `nlat`, `nlong`)

# Parameters

`latglq` : output, real\*8, dimension (`lmax`+1)
:   The latitude coordinates of a grid, corresponding to the indices (:,1), in DEGREES.

`longlq` : output, real\*8, dimension (2\*`lmax`+1)
:   The longitude coordinates of a grid, corresponding to the indices (1,:), in DEGREES. The first node is 0 E.

`lmax` : input, integer
:   The maximum spherical harmonic degree that will be integrated exactly by Gauss-Legendre quadrature.

`nlat` : output, integer
:   The number of samples in latitude.

`nlong` : output, integer
:   The number of samples in longitude.

# Description

`GLQGridCoord` will compute the latitude and longitude coordinates that are used in Gauss-Legendre quadrature grids for performing spherical harmonic transforms and reconstructions. The latitudinal nodes correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees.

# See also

[shglq](shglq.html), [shexpandglq](shexpandglq.html), [makegridglq](makegridglq.html), [shexpandglqc](shexpandglqc.html), [makegridglqc](makegridglqc.html), [preglq](preglq.html)
