# PlmSchmidt_d1

Compute all the Schmidt semi-normalized associated Legendre functions and first derivatives.

# Usage

`p`, `dp` = PlmSchmidt_d1 (`lmax`, `z`, [`csphase`, `cnorm`])

# Returns

`p` : float, dimension ((`lmax`+1)\*(`lmax`+2)/2)
:   An array of Schmidt-normalized associated Legendre functions up to degree `lmax`. The index corresponds to `l*(l+1)/2+m`.

`dp` : float, dimension ((`lmax`+1)\*(`lmax`+2)/2)
:   An array of the first derivatives of the Schmidt-normalized associated Legendre functions up to degree `lmax`. The index corresponds to `l*(l+1)/2+m`. 

# Parameters

`lmax` : integer
:   The maximum degree of the associated Legendre functions to be computed.

`z` : float
:   The argument of the associated Legendre functions.

`csphase` : optional, integer, default = 1
:   If 1 (default), the Condon-Shortley phase will be excluded. If -1, the Condon-Shortley phase of (-1)^m will be appended to the associated Legendre functions.

`cnorm` : optional, integer, default = 0
:   If 1, the complex normalization of the associated Legendre functions will be used. The default is to use the real normalization.

# Description

`PlmSchmidt_d1` will calculate all of the Schmidt-normalized associated Legendre functions and first derivatives up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula, and in order to prevent overflows, the scaling approach of Holmes and Featherstone (2002) is utilized. These functions are accurate to about degree 2800. The index of the array corresponding to a given degree `l` and angular order `m` corresponds to `l*(l+1)/2+m`.

The integral of the squared Legendre functions over the interval [-1, 1] is `2*(2-delta(0,m))/(2l+1)`, where delta is the Kronecker delta function. If the optional parameter `cnorm` is set equal to 1, the complex normalization will be used where the integral of the squared Legendre functions over the interval [-1, 1] is `2/(2l+1)`. The default is to exclude the Condon-Shortley phase, but this can be modified by setting the optional argument `csphase` to -1. Note that the derivative of the Legendre functions is calculated with respect to its arguement `z`, and not latitude or colatitude. If `z=cos(theta)`, where `theta` is the colatitude, then it is only necessary to multiply `dp` by `-sin(theta)` to obtain the derivative with respect to `theta`.

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and
order normalised associated Legendre functions, J. Geodesy, 76, 279-
299, 2002.

# See also

[plbar](pyplbar.html), [plbar_d1](pyplbar_d1.html), [plmbar](pyplmbar.html), [plmbar_d1](pyplmbar_d1.html), [plon](pyplon.html), [plon_d1](pyplon_d1.html), [plmon](pyplmon.html), [plmon_d1](pyplmon_d1.html), [plschmidt](pyplschmidt.html), [plmschmidt](pyplmschmidt.html), [plmschmidt_d1](pyplmschmidt_d1.html), [plegendre](pyplegendre.html), [plegendre_d1](pyplegendre_d1.html), [plegendrea](pyplegendrea.html), [plegendrea_d1](pyplegendrea_d1.html)
