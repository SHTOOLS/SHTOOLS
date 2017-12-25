# PlmBar

Compute all the 4-pi (geodesy) normalized associated Legendre functions.

# Usage

call PlmBar (`p`, `lmax`, `z`, `csphase`, `cnorm`, `exitstatus`)

# Parameters

`p` : output, real\*8, dimension ((`lmax`+1)\*(`lmax`+2)/2)
:   An array of 4-pi (geodesy) normalized associated Legendre functions up to degree `lmax`. The index corresponds to `l*(l+1)/2+m+1`, which can be calculated by a call to `PlmIndex`.

`lmax` : input, integer
:   The maximum degree of the associated Legendre functions to be computed. If `lmax = -1`, allocated memory will be deallocated.

`z` : input, real*8
:   The argument of the associated Legendre functions.

`csphase` : input, optional, integer
:   If 1 (default), the Condon-Shortley phase will be excluded. If -1, the Condon-Shortley phase of (-1)^m will be appended to the associated Legendre functions.

`cnorm` : input, optional, integer
:   If 1, the complex normalization of the associated Legendre functions will be used. The default is to use the real normalization.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`PlmBar` will calculate all of the 4-pi (geodesy) normalized associated Legendre functions up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula, and in order to prevent overflows, the scaling approach of Holmes and Featherstone (2002) is utilized. These functions are accurate to about degree 2800. The index of the array corresponding to a given degree `l` and angular order `m` corresponds to `l*(l+1)/2+m+1`, which can be computed by a call to `PlmIndex`. 

The integral of the squared Legendre functions over the interval [-1, 1] is `2*(2-delta(0,m))`, where `delta` is the Kronecker delta function. If the optional parameter `cnorm` is set equal to 1, the complex normalization will be used where the integral of the squared Legendre functions over the interval [-1, 1] is `2`. The default is to exclude the Condon-Shortley phase, but this can be modified by setting the optional argument `csphase` to -1.

This routine saves the three-term recursion factors and square roots of the integers the first time being called. If subsequent calls possess the same value of `lmax`, these will not be recomputed. If you wish to deallocate this memory, which is an array of length `(lmax+1)*(lmax+2)`, recall this routine with `lmax = -1`.

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and
order normalised associated Legendre functions, J. Geodesy, 76, 279-
299, 2002.

# See also

[plbar](plbar.html), [plbar_d1](plbar_d1.html), [plmbar_d1](plmbar_d1.html), [plon](plon.html), [plon_d1](plon_d1.html), [plmon](plmon.html), [plmon_d1](plmon_d1.html), [plschmidt](plschmidt.html), [plschmidt_d1](plschmidt_d1.html), [plmschmidt](plmschmidt.html), [plmschmidt_d1](plmschmidt_d1.html), [plegendre](plegendre.html), [plegendre_d1](plegendre_d1.html), [plegendrea](plegendrea.html), [plegendrea_d1](plegendrea_d1.html), [plmindex](plmindex.html)
