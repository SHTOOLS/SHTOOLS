# MakeGeoidGridDH

Create a global map of the geoid.

# Usage

`geoid` = MakeGeoidGridDH (`cilm`, `r0`, `gm`, `potref`, [`lmax`, `omega`, `r`, `order`, `lmax_calc`, `a`, `f`])

# Returns

`geoid` : float, dimension (`2lmax+2`, `sampling`\*2`lmax+2`)
:   A global grid of the height to the potential `potref` above a flattened ellipsoid of equatorial radius `a` and flattening `f`.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real spherical harmonic coefficients (geodesy normalized) of the gravitational potential referenced to a spherical interface of radius `r0`.

`r0` : float
:   The reference radius of the spherical harmonic coefficients.

`gm` : float
:   The product of the gravitational constant and mass of the planet.

`potref` : float
:   The value of the potential on the chosen geoid, in SI units.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the gravitational-potential coefficients. This determines the number of latitudinal and longitudinal samples.

`omega` : optional, float, default = 0
:   The angular rotation rate of the planet.

`r` : optional, float, default = `r0`
:   The radius of the reference sphere that the Taylor expansion of the potential is performed on.

`order` : optional, integer, default = 2
:   The order of the Taylor series expansion of the potential about the reference radius `r`. This can be either 1, 2, or 3.

`lmax_calc` : optional, integer, default = `lmax`
:   The maximum degree used in evaluating the spherical harmonic coefficients. This must be less than or equal to `lmax`.

`a` : optional, float, default = `r0`
:   The semi-major axis of the flattened ellipsoid that the output grid `geoid` is referenced to. The optional parameter `f` must also be specified.

`f` : optional, float, default = 0
:   The flattening `(R_equator-R_pole)/R_equator` of the reference ellipsoid. The optional parameter `a` (i.e., `R_equator`) must be specified.

# Description

`MakeGeoidGrid` will create a global map of the geoid, accurate to either first, second, or third order, using the method described in Wieczorek (2007; equation 19-20). The algorithm expands the potential in a Taylor series on a spherical interface of radius `r`, and computes the height above this interface to the potential `potref` exactly from the linear, quadratic, or cubic equation at each grid point. If the optional parameters `a` and `f` are specified, the geoid height will be referenced to a flattened ellipsoid with semi-major axis `a` and flattening `f`. The pseudo-rotational potential is explicitly accounted for by specifying the angular rotation rate `omega` of the planet. 

It should be noted that this geoid calculation is only strictly exact when the radius `r` lies above the maximum radius of the planet. Furthermore, the geoid is only strictly valid when it lies above the surface of the planet (it is necessary to know the density structure of the planet when calculating the potential below the surface).

The default is to calculate grids for use in the Driscoll and Healy routines that are equally spaced (`n` by `2n`), but this can be changed to calculate equally sampled grids (`n` by `n`) by setting the optional argument `sampling` to 1. This routine uses geodesy 4-pi normalized spherical harmonics that exclude the Condon-Shortley phase.

# References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Wieczorek, M. A. Gravity and topography of the terrestrial planets, Treatise on Geophysics, 10, 165-206, 2007.

# See also

[makegrid2d](pymakegrid2d.html), [makegridglq](pymakegridglq.html), [makegriddh](pymakegriddh.html), [makegravgriddh](pymakegravgriddh.html), [makegravgradgriddh](pymakegravgradgriddh.html)
