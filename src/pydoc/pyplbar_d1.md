# PlBar_d1

Compute all the 4-pi (geodesy) normalized Legendre polynomials and first derivatives.

# Usage

`p`, `dp` = PlBar_d1 (`lmax`, `z`)

# Returns

`p` : float, dimension (`lmax`+1)
:   An array of 4-pi (geodesy) normalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l`.

`dp` : float, dimension (`lmax`+1)
:   An array of the first derivatives of the 4-pi (geodesy) normalized Legendre polynomials up to degree `lmax`.

# Parameters

`lmax` : integer
:   The maximum degree of the Legendre polynomials to be computed.

`z` : float
:   The argument of the Legendre polynomial.

# Description

`PlBar_d1` will calculate all of the 4-pi (geodesy) normalized Legendre polynomials and first derivatives up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula, and the integral of the geodesy-normalized Legendre polynomials over the interval [-1, 1] is 2. Note that the derivative of the Legendre polynomials is calculated with respect to its arguement `z`, and not latitude or colatitude. If `z=cos(theta)`, where theta is the colatitude, then it is only necessary to multiply `dp` by `-sin(theta)` to obtain the derivative with respect to theta.

# See also

[plbar](pyplbar.html), [plmbar](pyplmbar.html), [plmbar_d1](pyplmbar_d1.html), [plon](pyplon.html), [plon_d1](pyplon_d1.html), [plmon](pyplmon.html), [plmon_d1](pyplmon_d1.html), [plschmidt](pyplschmidt.html), [plschmidt_d1](pyplschmidt_d1.html), [plmschmidt](pyplmschmidt.html), [plmschmidt_d1](pyplmschmidt_d1.html), [plegendre](pyplegendre.html), [plegendre_d1](pyplegendre_d1.html), [plegendrea](pyplegendrea.html), [plegendrea_d1](pyplegendrea_d1.html)
