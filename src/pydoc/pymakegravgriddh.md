# MakeGravGridDH()

Create 2D cylindrical maps on a flattened and rotating ellipsoid of all three components of the gravity, the gravity disturbance, and the gravity potential.

# Usage

rad, theta, phi, total, pot = MakeGravGridDH (cilm, gm, r0, [a, f, lmax, sampling, lmax_calc, omega, normal_gravity, normal_gravity_gm, extend])

# Returns

rad : float, dimension (nlat, nlong)
:   A 2D map of the radial component of the gravity field that conforms to the sampling theorem of Driscoll and Healy (1994). If sampling is 1, the grid is equally sampled and is dimensioned as (n by n), where n is 2lmax+2. If sampling is 2, the grid is equally spaced and is dimensioned as (n by 2n). The first latitudinal band of the grid corresponds to 90 N, the latitudinal sampling interval is 180/n degrees, and the default behavior is to exclude the latitudinal band for 90 S. The first longitudinal band of the grid is 0 E, by default the longitudinal band for 360 E is not included, and the longitudinal sampling interval is 360/n for an equally sampled and 180/n for an equally spaced grid, respectively. If extend is 1, the longitudinal band for 360 E and the latitudinal band for 90 S will be included, which increases each of the dimensions of the grid by 1.

theta : float, dimension (nlat, nlong)
:   A 2D equally sampled or equally spaced grid of the theta component of the gravity field.

phi : float, dimension (nlat, nlong)
:   A 2D equally sampled or equally spaced grid of the phi component of the gravity field.

total : float, dimension (nlat, nlong)
:   A 2D equally sampled or equally spaced grid of the magnitude of the gravity acceleration.

pot : float, dimension (nlat, nlong)
:   A 2D equally sampled or equally spaced grid of the gravity potential.

# Parameters

cilm : float, dimension (2, lmaxin+1, lmaxin+1)
:   The real 4-pi normalized gravitational potential spherical harmonic coefficients. The coefficients c1lm and c2lm correspond to the cosine and sine coefficients, respectively, with c1lm=cilm[0,l,m] and c2lm=cilm[1,l,m].

gm : float
:   The gravitational constant multiplied by the mass of the gravitational potential model.

r0: float
:   The reference radius of the spherical harmonic coefficients.

a : optional, float, default = r0
:   The semi-major axis of the flattened ellipsoid on which the field is computed.

f : optional, float, default = 0
:   The flattening of the reference ellipsoid: f=(R_equator-R_pole)/R_equator.

lmax : optional, integer, default = lmaxin
:   The maximum spherical harmonic degree of the coefficients cilm. This determines the number of samples of the output grids, n=2lmax+2, and the latitudinal sampling interval, 90/(lmax+1).

sampling : optional, integer, default = 2
:   If 1 (default) the output grids are equally sampled (n by n). If 2, the grids are equally spaced (n by 2n).

lmax_calc : optional, integer, default = lmax
:   The maximum spherical harmonic degree used in evaluating the functions. This must be less than or equal to lmax.

omega : optional, float, default = 0
:   The angular rotation rate of the planet.

normal_gravity : optional, integer, default = 1
:   If 1, the normal gravity (the gravitational acceleration on the rotating ellipsoid) will be subtracted from the total gravity, yielding the "gravity disturbance." This is done using Somigliana's formula (after converting geocentric to geodetic coordinates).

normal_gravity_gm : optional, float, default = gm
:   The GM to use when computing the normal gravity. If not specified, the GM of the gravitational potential model will be used.

extend : input, optional, bool, default = False
:   If True, compute the longitudinal band for 360 E and the latitudinal band for 90 S. This increases each of the dimensions of griddh by 1.

# Description

MakeGravGridDH will compute 2-dimensional cylindrical maps of the three spherical vector components, the vector magnitude, and the potential of either the gravity (gravitational force + centrifugal force) or gravitation from the input spherical harmonic coefficients cilm. The grids are equally sampled (n by n) or equally spaced (n by 2n) in latitude and longitude, where latitude is geocentric spherical. The gravitational potential is given by

`V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm}`,

and the gravitational acceleration is

`B = Grad V`.

The coefficients are referenced to a radius r0, and the function is computed on a flattened ellipsoid with semi-major axis a and flattening f. All grids are output in SI units, and the sign of the radial components is positive when directed upwards.

If the optional angular rotation rate omega is specified, the potential and gravity vectors will be calculated in a body-fixed rotating reference frame and will include the contribution from the centrifugal force. If not specified, the output grids will correspond to the external gravitation without the centrifugal force. If normal_gravity is set to 1, the normal gravity will be removed from the magnitude of the gravity vector, yielding the gravity disturbance. To convert from m/s^2 to mGals, multiply the gravity grids by 10^5.

The calculated values should be considered exact only when the radii on the ellipsoid are greater than the maximum radius of the planet (the potential coefficients are simply downward continued in the spectral domain). The components of the gravity vector are calculated using spherical coordinates whose origin is at the center of the planet, and the components are thus not normal to the reference ellipsoid.

The default is to use an input grid that is equally sampled (n by n), but this can be changed to use an equally spaced grid (n by 2n) by the optional argument sampling. The redundant longitudinal band for 360 E and the latitudinal band for 90 S are excluded by default, but these can be computed by specifying the optional argument extend.

# Reference

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.
