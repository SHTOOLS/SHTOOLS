# PLegendre

Compute all the unnormalized Legendre polynomials.

# Usage

`p` = PLegendre (`lmax`, `z`)

# Returns

`p` : float, dimension (`lmax`+1)
:   An array of unnormalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l`.

# Parameters

`lmax` : integer
:   The maximum degree of the Legendre polynomials to be computed.

`z` : float
:   The argument of the Legendre polynomial.

# Description

`PLegendre` will calculate all of the unnormalized Legendre polynomials up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula. The integral of the Legendre polynomials over the interval [-1, 1] is `2/(2l+1)`.

# See also

[plbar](pyplbar.html), [plbar_d1](pyplbar_d1.html), [plmbar](pyplmbar.html), [plmbar_d1](pyplmbar_d1.html), [plon](pyplon.html), [plon_d1](pyplon_d1.html), [plmon](pyplmon.html), [plmon_d1](pyplmon_d1.html), [plschmidt](pyplschmidt.html), [plschmidt_d1](pyplschmidt_d1.html), [plmschmidt](pyplmschmidt.html), [plmschmidt_d1](pyplmschmidt_d1.html), [plegendre_d1](pyplegendre_d1.html), [plegendrea](pyplegendrea.html), [plegendrea_d1](pyplegendrea_d1.html)
