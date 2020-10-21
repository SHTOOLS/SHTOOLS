# MakeGravGridPoint

Determine the three components of the gravity vector at a single point.

# Usage

`value` = MakeGravGridPoint (`cilm`, `lmax`, `gm`, `r0`, `r`, `lat`, `lon`, `omega`, `dealloc`)

# Parameters

`value` : output, real(dp), dimension (3)
:   Vector components (r, theta, phi) of the gravity at (`r`, `lat`, `lon`).

`cilm` : input, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The real 4-pi normalized spherical harmonic coefficients of the gravitational potential. The coefficients `C1lm` and `C2lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm(1,l+1,m+1)` and `Slm=cilm(2,l+1,m+1)`.

`lmax` : input, integer(int32)
:   The maximum spherical harmonic degree used in evaluating the function.

`gm` : input, real(dp)
:   The gravitational constant multiplied by the mass of the planet.

`r0`: input, real(dp)
:   The reference radius of the spherical harmonic coefficients.

`r`: input, real(dp)
:   The radius to evaluate the gravity field.

`lat` : input, real(dp)
:   The latitude of the point in DEGREES.

`lon` : input, real(dp)
:   The longitude of the point in DEGREES.

`omega` : optional, input, real(dp)
:   The angular rotation rate of the planet.

`dealloc` : input, optional, integer(int32), default = 0
:   0 (default) = Save variables used in the external Legendre function calls. (1) Deallocate this memory at the end of the funcion call.

# Description

`MakeGravGridPoint` will compute the three components of the gravity vector (gravitational force + centrifugal force) at a single point. The input latitude and longitude are in degrees, and the output components of the gravity are in spherical coordinates (r, theta, phi). The gravitational potential is given by

`V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm}`,

and the gravitational acceleration is

`B = Grad V`.

The coefficients are referenced to a radius `r0`, and the output accelerations are in m/s^2. To convert m/s^2 to mGals, multiply by 10^5. If the optional angular rotation rate `omega` is specified, the gravity vector will be calculated in a body-fixed rotating reference frame and will include the contribution of the centrifugal force.

# See also

[makegeoidgrid](makegeoidgrid.html), [makegravgradgriddh](makegravgradgriddh.html), [normalgravity](normalgravity.html), [makegriddh](makegriddh.html)