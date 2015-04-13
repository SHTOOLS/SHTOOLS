# PLegendreA_d1

Compute all the unnormalized associated Legendre functions and first derivatives.

# Usage

call PLegendreA_d1 (`p`, `dp`, `lmax`, `z`, `csphase`)

# Parameters

`p` : output, real\*8, dimension ((`lmax`+1)\*(`lmax`+2)/2)
:   An array of unnormalized associated Legendre functions up to degree `lmax`. The index corresponds to `l*(l+1)/2+m+1`, which can be calculated by a call to `PlmIndex`.
	
`dp` : output, real\*8, dimension ((`lmax`+1)\*(`lmax`+2)/2)
:   An array of the first derivatives of the unnormalized associated Legendre functions up to degree `lmax`. The index corresponds to `l*(l+1)/2+m+1`, which can be calculated by a call to `PlmIndex`.

`lmax` : input, integer
:   The maximum degree of the associated Legendre functions to be computed.

`z` : input, real\*8
:   The argument of the associated Legendre functions.

`csphase` : input, integer, optional, default = 1
:   If 1 (default), the Condon-Shortley phase will be excluded. If -1, the Condon-Shortley phase of (-1)^m will be appended to the associated Legendre functions.

# Description

`PLegendreA_d1` will calculate all of the unnormalized associated Legendre functions and first derivatives up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula and hence will overflow for moderate values of `l` and `m`. The index of the array corresponding to a given degree `l` and angular order `m` corresponds to `l*(l+1)/2+m+1`, which can be calculated by a call to `PlmIndex`. The integral of the associated Legendre functions over the interval [-1, 1] is `2*(l+m)!/(l-m)!/(2l+1)`. The default is to exclude the Condon-Shortley phase, but this can be modified by setting the optional argument `csphase` to -1. Note that the derivative of the Legendre polynomials is calculated with respect to its arguement `z`, and not latitude or colatitude. If `z=cos(theta)`, where `theta` is the colatitude, then it is only necessary to multiply `dp` by `-sin(theta)` to obtain the derivative with respect to `theta`.

# See also

[`plbar`](plbar.html), [`plbar_d1`](plbar_d1.html), [`plmbar`](plmbar.html), [`plmbar_d1`](plmbar_d1.html), [`plon`](plon.html), [`plon_d1`](plon_d1.html), [`plmon`](plmon.html), [`plmon_d1`](plmon_d1.html), [`plschmidt`](plschmidt.html), [`plschmidt_d1`](plschmidt_d1.html), [`plmschmidt`](plmschmidt.html), [`plmschmidt_d1`](plmschmidt_d1.html), [`plegendre`](plegendre.html), [`plegendre_d1`](plegendre_d1.html), [`plegendrea`](plegendrea.html), [`plmindex`](plmindex.html)
