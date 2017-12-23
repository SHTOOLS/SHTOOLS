# PlON_d1

Compute all the orthonormalized Legendre polynomials and first derivatives.

# Usage

`p`, `dp` = PlON_d1 (`lmax`, `z`)

# Returns

`p` : float, dimension (`lmax`+1)
:   An array of orthonormalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l`.

`dp` : float, dimension (`lmax`+1)
:   An array of the first derivatives of the orthonormalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l`.

# Parameters

`lmax` : integer
:   The maximum degree of the Legendre polynomials to be computed.

`z` : float
:   The argument of the Legendre polynomial.

# Description

`PlON_d1` will calculate all of the orthonormalized Legendre polynomials and first derivatives up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula, and the integral of the orthonormalized Legendre polynomials over the interval [-1, 1] is `2/(4pi)`. Note that the derivative of the Legendre polynomials is calculated with respect to its arguement `z`, and not latitude or colatitude. If `z=cos(theta)`, where `theta` is the colatitude, then it is only necessary to multiply `dp` by `-sin(theta)` to obtain the derivative with respect to `theta`.

# See also

[plbar](pyplbar.html), [plbar_d1](pyplbar_d1.html), [plmbar](pyplmbar.html), [plmbar_d1](pyplmbar_d1.html), [plon](pyplon.html), [plmon](pyplmon.html), [plmon_d1](pyplmon_d1.html), [plschmidt](pyplschmidt.html), [plschmidt_d1](pyplschmidt_d1.html), [plmschmidt](pyplmschmidt.html), [plmschmidt_d1](pyplmschmidt_d1.html), [plegendre](pyplegendre.html), [plegendre_d1](pyplegendre_d1.html), [plegendrea](pyplegendrea.html), [plegendrea_d1](pyplegendrea_d1.html)
