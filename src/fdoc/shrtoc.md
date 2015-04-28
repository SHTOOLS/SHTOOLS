# SHrtoc

Convert real spherical harmonics to complex form.

# Usage

call SHrtoc (`rcilm`, `ccilm`, `degmax`, `convention`, `switchcs`)

# Parameters

`rcilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input real spherical harmonic coefficients. `rcilm(1,:,:)` and `rcilm(2,:,:)` correspond to the cosine and sine terms, respectively.
	
`ccilm` : output, real\*8, dimension (2, `lmaxout`+1, `lmaxout`+1)
:   The output complex spherical harmonic coefficients. `ccilm(1,:,:)` and `ccilm(2,:,:)` correspond to the real and complex part of the coefficients, respectively. Only the positive angular orders are output; the negative orders can be calculated from the relation `C_{l-m}=(-1)^m C_{lm}^*`.

`degmax` : input, optional, integer, default = min(`lmaxin`, `lmaxout`)
:   The maximum degree of the output coefficients.

`convention` : input, optional, integer, default = 1
:   If 1 (default), the input and output coefficients will have the same normalization. If 2, real geodesy 4-pi coefficients will be converted to complex orthonormal form.

`swtichcs` : input, optional, integer, default = 0
:   If 0 (default), the input and output coefficients will possess the same Condon-Shortley phase convention. If 1, the input coefficients will first be multiplied by (-1)^m.

# Description

`SHrtoc` will convert real spherical harmonics to complex form. By default, the dimension of the output array is the minimum of `rcilm` and `ccilm`, though this can be changed by specifying the optional parameter `degmax`. The normalization of the input and output coefficients are by default the same, but if the optional argument `convention` is set to 2, this routine will convert from geodesy 4-pi normalized coefficients to orthonormalized coefficients. The Condon-Shortley phase convention between the input an output coefficients can be modified by the optional argument `switchs`.

# See also

[shctor](shctor.html)
