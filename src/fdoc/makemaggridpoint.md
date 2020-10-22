# MakeMagGridPoint

Determine the three components of the magnetic field vector at a single point.

# Usage

`value` = MakeMagGridPoint (`cilm`, `lmax`, `a`, `r`, `lat`, `lon`, `dealloc`)

# Parameters

`value` : output, real(dp), dimension (3)
:   Vector components (r, theta, phi) of the magnetic field at (`r`, `lat`, `lon`).

`cilm` : input, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The real Schmidt semi-normalized spherical harmonic coefficients of the magnetic potential. The coefficients `C1lm` and `C2lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm(1,l+1,m+1)` and `Slm=cilm(2,l+1,m+1)`.

`lmax` : input, integer(int32)
:   The maximum spherical harmonic degree used in evaluating the function.

`a`: input, real(dp)
:   The reference radius of the spherical harmonic coefficients.

`r`: input, real(dp)
:   The radius to evaluate the magnetic field field.

`lat` : input, real(dp)
:   The latitude of the point in degrees.

`lon` : input, real(dp)
:   The longitude of the point in degrees.

`dealloc` : input, optional, integer(int32), default = 0
:   0 (default) = Save variables used in the external Legendre function calls. (1) Deallocate this memory at the end of the funcion call.

# Description

`MakeMagGridPoint` will compute the three components of the magnetic field vector at a single point. The input latitude and longitude are in degrees, and the output components are in spherical coordinates (r, theta, phi). The magnetic potential is given by

`V = a Sum_{l=0}^lmax (a/r)^(l+1) Sum_{m=-l}^l C_{lm} Y_{lm}`,

and the vector magnetic field is

`B = - Grad V`.

The coefficients are referenced to a radius `a`, and the output values have the same units as the input spherical harmonic coefficients.

# See also

[makemaggradgriddh](makemaggradgriddh.html)