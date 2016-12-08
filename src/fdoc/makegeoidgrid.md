# MakeGeoidGrid

Create a global map of the geoid.

# Usage

call MakeGeoidGrid (`geoid`, `cilm`, `lmax`, `r0`, `gm`, `potref`, `omega`, `r`, `gridtype`, `order`, `nlat`, `nlong`, `interval`, `lmax_calc`, `a`, `f`, `exitstatus`)

# Parameters

`geoid` : output, real\*8, dimension(`nlat`, `nlong`)
:   A global grid of the height to the potential `potref` above a sphere of radius `r` (or above a flattened ellipsoid if both `a` and `f` are specified). The number of latitude and longitude points depends upon `gridtype`: (1) `lmax+1` by `2lmax+1`, (2) `2lmax+2` by `2lmax+2`, (3) `2lmax+2` by `4lmax+4`, or (4) `180/interval+1` by `360/interval+1`.

`cilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The real spherical harmonic coefficients (geodesy normalized) of the gravitational potential referenced to a spherical interface of radius `r0pot`.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the gravitational-potential coefficients. For `gridtype`s 1, 2 and 3, this determines the number of latitudinal and longitudinal samples.

`r0` : input, real\*8
:   The reference radius of the spherical harmonic coefficients.

`gm` : input, real\*8
:   The product of the gravitational constant and mass of the planet.

`potref` : input, real\*8
:   The value of the potential on the chosen geoid, in SI units.

`omega` : input, real\*8
:   The angular rotation rate of the planet.

`r` : input, real\*8
:   The radius of the reference sphere that the Taylor expansion of the potential is performed on. If `a` and `f` are not specified, the geoid height will be referenced to this spherical interface.

`gridtype` : input, integer
:   The output grid is (1) a Gauss-Legendre quadrature grid whose grid nodes are determined by `lmax`, (2) an equally sampled `n` by `n` grid used with the Driscoll and Healy (1994) sampling theorem, (3) ar a similar `n` by 2`n` grid that is oversampled in longitude, or (4) a 2D Cartesian grid with latitudinal and longitudinal spacing given by `interval`.

`order` : input, integer
:   The order of the Taylor series expansion of the potential about the reference radius `r`. This can be either 1, 2, or 3.

`nlat` : output, integer
:   The number of latitudinal samples.

`nlong` : output, integer
:   The number of longitudinal samples.

`interval`: optional, input, real\*8
:   The latitudinal and longitudinal spacing of the output grid when `gridtype` is 4.

`lmax_calc` : optional, input, integer
:   The maximum degree used in evaluating the spherical harmonic coefficients. This must be less than or equal to `lmax`.

`a` : optional, input, real\*8, default = `r0`
:   The semi-major axis of the flattened ellipsoid that the output grid `geoid` is referenced to. The optional parameter `f` must also be specified.

`f` : optional, input, real\*8, default = 0
:   The flattening `(R_equator-R_pole)/R_equator` of the reference ellipsoid. The optional parameter `a` (i.e., `R_equator`) must be specified.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`MakeGeoidGrid` will create a global map of the geoid, accurate to either first, second, or third order, using the method described in Wieczorek (2007; equation 19-20). The algorithm expands the potential in a Taylor series on a spherical interface of radius `r`, and computes the height above this interface to the potential `potref` exactly from the linear, quadratic, or cubic equation at each grid point. If the optional parameters `a` and `f` are specified, the geoid height will be referenced to a flattened ellipsoid with semi-major axis `a` and flattening `f`. The pseudo-rotational potential is explicitly accounted for by specifying the angular rotation rate `omega` of the planet.

It should be noted that this geoid calculation is only strictly exact when the radius `r` lies above the maximum radius of the planet. Furthermore, the geoid is only strictly valid when it lies above the surface of the planet (it is necessary to know the density structure of the planet when calculating the potential below the surface).

The geoid can be computed on one of four different grids: (1) a Gauss-Legendre quadrature grid (see `MakeGridGLQ`), (2) A `n` by `n` equally sampled grid (see `MakeGridDH`), (3) an `n` by 2`n` equally spaced grid (see `MakeGridDH`), or (4) A 2D Cartesian grid (see `MakeGrid2D`). This routine uses geodesy 4-pi normalized spherical harmonics that exclude the Condon-Shortley phase. This can not be modified.

# References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Wieczorek, M. A. Gravity and topography of the terrestrial planets, Treatise on Geophysics, 10, 165-206, 2007.

# See also

[makegrid2d](makegrid2d.html), [makegridglq](makegridglq.html), [makegriddh](makegriddh.html), [makegravgriddh](makegravgriddh.html), [makegravgradgriddh](makegravgradgriddh.html)
