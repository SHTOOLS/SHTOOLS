# SHctor

Convert complex spherical harmonics to real form.

# Usage

call SHctor (`ccilm`, `rcilm`, `degmax`, `convention`, `switchcs`, `exitstatus`)

# Parameters

`ccilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input complex spherical harmonic coefficients. `ccilm(1,:,:)` and `ccilm(2,:,:)` correspond to the real and complex part of the coefficients, respectively. Only the positive angular orders are input; the negative orders are assumed to satisfy the relation `C_{l-m}=(-1)^m C_{lm}^*`.

`rcilm` : output, real\*8, dimension (2, `lmaxout`+1, `lmaxout`+1)
:   The output real spherical harmonic coefficients. `rcilm(1,:,:)` and `rcilm(2,:,:)` correspond to the cosine and sine terms, respectively.

`degmax` : input, optional, integer, default = min(`lmaxin`, `lmaxout`)
:   The maximum degree of the output coefficients.

`convention` : input, optional, integer, default = 1
:   If 1 (default), the input and output coefficients will have the same normalization. If 2, orthonormalized coefficients will be converted to real geodesy 4-pi form.

`swtichcs` : input, optional, integer, default = 0
:   If 0 (default), the input and output coefficients will possess the same Condon-Shortley phase convention. If 1, the input coefficients will first be multiplied by (-1)^m.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHctor` will convert complex spherical harmonics of a real function to real form. By default, the dimension of the output array is the minimum of `rcilm(1,:,:)` and `ccilm(1,:,:)`, though this can be changed by specifying the optional parameter `degmax`. The normalization of the input and output coefficients are by default the same, but if the optional argument `convention` is set to 2, this routine will convert from geodesy 4-pi normalized coefficients to orthonormalized coefficients. The Condon-Shortley phase convention between the input an output coefficients can be modified by the optional argument `switchcs`.

# See also

[shrtoc](shrtoc.html)
