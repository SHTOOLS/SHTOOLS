# PlBar

Compute all the 4-pi (geodesy) normalized Legendre polynomials.

# Usage

call PlBar (`p`, `lmax`, `z`, `exitstatus`)

# Parameters

`p` : output, real/*8, dimension (`lmax`+1)
:   An array of geodesy-normalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l+1`.

`lmax` : input, integer
:   The maximum degree of the Legendre polynomials to be computed.

`z` : input, real/*8
:   The argument of the Legendre polynomial.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`PlBar` will calculate all of the 4-pi (geodesy) normalized Legendre polynomials up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula. The integral of the geodesy-normalized Legendre polynomials over the interval [-1, 1] is 2.

# See also

[plbar_d1](plbar_d1.html), [plmbar](plmbar.html), [plmbar_d1](plmbar_d1.html), [plon](plon.html), [plon_d1](plon_d1.html), [plmon](plmon.html), [plmon_d1](plmon_d1.html), [plschmidt](plschmidt.html), [plschmidt_d1](plschmidt_d1.html), [plmschmidt](plmschmidt.html), [plmschmidt_d1](plmschmidt_d1.html), [plegendre](plegendre.html), [plegendre_d1](plegendre_d1.html), [plegendrea](plegendrea.html), [plegendrea_d1](plegendrea_d1.html)

