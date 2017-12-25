# MakeGridPoint

Evaluate a real function expressed in real spherical harmonics at a single point.

# Usage

`value` = MakeGridPoint (`cilm`, `lmax`, `lat`, `lon`, `norm`, `csphase`, `dealloc`)

# Parameters

`value` : output, real\*8
:   Value of the function at (`lat`, `lon`).

`cilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The real spherical harmonic coefficients of the function. The coefficients `C1lm` and `C2lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm(1,l+1,m+1)` and `Slm=cilm(2,l+1,m+1)`.

`lmax` : input, integer
:   The maximum spherical harmonic degree used in evaluating the function.

`lat` : input, real\*8
:   The latitude of the point in DEGREES.

`lon` : input, real\*8
:   The longitude of the point in DEGREES.

`norm` : input, optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`dealloc` : input, optional, integer, default = 0
:   0 (default) = Save variables used in the external Legendre function calls. (1) Deallocate this memory at the end of the funcion call.

# Description

`MakeGridPoint` will expand a function expressed in spherical harmonics at a single point. The input latitude and longitude are in degrees. The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

# See also

[makegridpointc](makegridpointc.html), [makegriddh](makegriddh.html), [makegriddhc](makegriddhc.html), [makegridglq](makegridglq.html), [makegridglqc](makegridglqc.html)
