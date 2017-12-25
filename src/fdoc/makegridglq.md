# MakeGridGLQ

Create a 2D map from a set of spherical harmonic coefficients sampled on the Gauss-Legendre quadrature nodes.

# Usage

call MakeGridGLQ (`gridglq`, `cilm`, `lmax`, `plx`, `zero`, `norm`, `csphase`, `lmax_calc`, `exitstatus`)

# Parameters

`gridglq` : output, real\*8, dimension (`lmax`+1, 2\*`lmax`+1)
:   A 2D map of the function sampled on the Gauss-Legendre quadrature nodes.

`cilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real spherical harmonic coefficients of the function. When evaluating the function, the maximum spherical harmonic degree considered is the minimum of `lmax`, `lmaxin`, or `lmax_calc` (if specified). The coefficients `C1lm` and `C2lm` refer to the "cosine" (`Clm`) and "sine" (`Slm`) coefficients, respectively, with `Clm=cilm(1,l+1,m+1)` and `Slm=cilm(2,l+1,m+1)`.

`lmax` : input, integer
:   The maximum spherical harmonic bandwidth of the function. This determines the sampling nodes and dimensions of the output grid.

`plx` : input, optional, real\*8, dimension (`lmax`+1, (`lmax`+1)\*(`lmax`+2)/2)
:   An array of the associated Legendre functions calculated at the Gauss-Legendre quadrature nodes. These are determined from a call to `SHGLQ`. Either `plx` or `zero` must be present, but not both.

`zero` : input, optional, real\*8, dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude, calculated by a call to `SHGLQ`.  Either `plx` or `zero` must be present, but not both.

`norm` : input, optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`lmax_calc` : input, optional, integer, default = `lmax`
:   The maximum spherical harmonic degree used in evaluating the function. This must be less than or equal to `lmax`.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.
# Description

`MakeGridGLQ` will create a 2-dimensional map from a set of input spherical harmonic coefficients sampled on the Gauss-Legendre quadrature nodes. This is the inverse of the routine `SHExpandGLQ`. The latitudinal nodes correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees. When evaluating the function, the maximum spherical harmonic degree that is considered is the minimum of `lmax`, the size of `lmaxin`, or `lmax_calc` (if specified).

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m. The normalized legendre functions are calculated using the scaling algorithm of Holmes and Featherstone (2002), which are accurate to about degree 2800. The unnormalized functions are accurate only to about degree 15.

The reconstruction of the spherical harmonic function may be speeded up by precomputing the Legendre functions on the Gauss-Legendre quadrature nodes in the routine `SHGLQ`. However, given that this array contains on the order of `lmax`**3 entries, this is only feasible for moderate values of `lmax`.

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and
order normalised associated Legendre functions, J. Geodesy, 76, 279-
299, 2002.

# See also

[shexpandglq](shexpandglq.html), [shexpandglqc](shexpandglqc.html), [makegridglqc](makegridglqc.html), [shexpanddh](shexpanddh.html), [makegriddh](makegriddh.html), [shexpanddhc](shexpanddhc.html), [makegriddhc](makegriddhc.html), [shexpandlsq](shexpandlsq.html), [glqgridcoord](glqgridcoord.html), [shglq](shglq.html)
