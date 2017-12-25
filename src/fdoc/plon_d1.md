# PlON_d1

Compute all the orthonormalized Legendre polynomials and first derivatives.

# Usage

call PlON_d1 (`p`, `dp`, `lmax`, `z`, `exitstatus`)

# Parameters

`p` : output, real\*8, dimension (`lmax`+1)
:   An array of orthonormalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l+1`.

`dp` : output, real\*8, dimension (`lmax`+1)
:   An array of the first derivatives of the orthonormalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l+1`.

`lmax` : input, integer
:   The maximum degree of the Legendre polynomials to be computed.

`z` : input, real\*8
:   The argument of the Legendre polynomial.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`PlON_d1` will calculate all of the orthonormalized Legendre polynomials and first derivatives up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula, and the integral of the orthonormalized Legendre polynomials over the interval [-1, 1] is `2/(4pi)`. Note that the derivative of the Legendre polynomials is calculated with respect to its arguement `z`, and not latitude or colatitude. If `z=cos(theta)`, where `theta` is the colatitude, then it is only necessary to multiply `dp` by `-sin(theta)` to obtain the derivative with respect to `theta`.

# See also

[plbar](plbar.html), [plbar_d1](plbar_d1.html), [plmbar](plmbar.html), [plmbar_d1](plmbar_d1.html), [plon](plon.html), [plmon](plmon.html), [plmon_d1](plmon_d1.html), [plschmidt](plschmidt.html), [plschmidt_d1](plschmidt_d1.html), [plmschmidt](plmschmidt.html), [plmschmidt_d1](plmschmidt_d1.html), [plegendre](plegendre.html), [plegendre_d1](plegendre_d1.html), [plegendrea](plegendrea.html), [plegendrea_d1](plegendrea_d1.html)
