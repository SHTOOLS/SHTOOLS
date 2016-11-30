# MakeGravGridDH

Create 2D cylindrical maps on a flattened and rotating ellipsoid of all three components of the gravity field, the gravity disturbance, and the gravitational potential.

# Usage

call MakeGravGridDH (`cilm`, `lmax`, `gm`, `r0`, `a`, `f`, `rad`, `theta`, `phi`, `total`, `n`, `sampling`, `lmax_calc`, `omega`, `normal_gravity`, `pot_grid`, `exitstatus`)

# Parameters

`cilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The real gravitational potential spherical harmonic coefficients. The coefficients c1lm and c2lm refer to the cosine and sine coefficients, respectively, with `c1lm=cilm(1,l+1,m+1)` and `c2lm=cilm(2,l+1,m+1)`.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the coefficients `cilm`. This determines the number of samples of the output grids, `n=2lmax+2`, and the latitudinal sampling interval, `90/(lmax+1)`.

`gm` : input, real\*8
:   The gravitational constant multiplied by the mass of the planet.

`r0`: input, real\*8
:   The reference radius of the spherical harmonic coefficients.

`a` : input, real\*8
:   The semi-major axis of the flattened ellipsoid on which the field is computed.

`f` : input, real\*8
:   The flattening of the reference ellipsoid: `f=(R_equator-R_pole)/R_equator`.

`rad` : output, real\*8, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled (`n` by `n`) or equally spaced (`n` by 2`n`) grid of the radial component of the gravity field corresponding to the input spherical harmonic coefficients `cilm`. The first latitudinal band corresponds to 90 N, the latitudinal band for 90 S is not included, and the latitudinal sampling interval is 180/`n` degrees. The first longitudinal band is 0 E, the longitudinal band for 360 E is not included, and the longitudinal sampling interval is 360/`n` for an equally sampled and 180/`n` for an equally spaced grid, respectively.

`theta` : output, real\*8, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the theta component of the gravity field.

`phi` : output, real\*8, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the phi component of the gravity field.

`total` : output, real\*8, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the magnitude of the gravity acceleration.

`n` : output, integer
:   The number of samples in latitude of the output grids. This is equal to `2lmax+2`.

`sampling` : optional, input, integer, default = 1
:   If 1 (default) the output grids are equally sampled (`n` by `n`). If 2, the grids are equally spaced (`n` by 2`n`).

`lmax_calc` : optional, input, integer
:   The maximum spherical harmonic degree used in evaluating the functions. This must be less than or equal to `lmax`.

`omega` : optional, input, real\*8
:   The angular rotation rate of the planet.

`normal_gravity` : optional, input, integer
:   If 1, the normal gravity (the gravitational acceleration on the ellipsoid) will be subtracted from the total gravity, yielding the "gravity disturbance." This is done using Somigliana's formula (after converting geocentric to geodetic coordinates).

`pot_grid` : output, real\*8, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equaly spaced grid of the gravitational potential.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`MakeGravGridDH` will create 2-dimensional cylindrical maps from the spherical harmonic coefficients `cilm`, equally sampled (`n` by `n`) or equally spaced (`n` by 2`n`) in latitude and longitude, for the three vector components of the gravity field, the magnitude of the gravity field, and the potential (all using geocentric coordinates). The gravitational potential is given by

`V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm}`,

and the gravitational acceleration is

`B = Grad V`.

The coefficients are referenced to a radius `r0`, and the function is computed on a flattened ellipsoid with semi-major axis `a` (i.e., the mean equatorial radius) and flattening `f`. All grids are output in SI units, and the sign of the radial components is positive when directed upwards. If the optional angular rotation rate `omega` is specified, the potential and radial gravitational acceleration will be calculated in a body-fixed rotating reference frame.

To remove the "normal gravity" (the total gravitational acceleration on the ellipsoid) from the magnitude of the total gravity field (to obtain the "gravity disturbance"), set `normal_gravity` to 1. To convert m/s^2 to mGals, multiply the gravity grids by 10^5.

The calculated values should be considered exact only when the radii on the ellipsoid are less than the maximum radius of the planet (the potential coefficients are simply downward continued in the spectral domain). The components of the gravity vector are calculated using spherical coordinates whose origin is at the center of the planet, and the components are thus not normal to the reference ellipsoid.

The default is to calculate grids for use in the Driscoll and Healy routines that are equally sampled (`n` by `n`), but this can be changed to calculate equally spaced grids (`n` by 2`n`) by setting the optional argument `sampling` to 2. The input value of `lmax` determines the number of samples, `n=2lmax+2`, and the latitudinal sampling interval, 90/(`lmax`+1). The first latitudinal band of the grid corresponds to 90 N, the latitudinal band for 90 S is not calculated, and the latitudinal sampling interval is 180/`n` degrees. The first longitudinal band is 0 E, the longitudinal band for 360 E is not calculated, and the longitudinal sampling interval is 360/`n` for equally sampled and 180/`n` for equally spaced grids, respectively.

# References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

# See also

[makegeoidgrid](makegeoidgrid.html), [makegravgradgriddh](makegravgradgriddh.html), [normalgravity](normalgravity.html), [makegriddh](makegriddh.html)
