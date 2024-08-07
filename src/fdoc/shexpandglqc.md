# SHExpandGLQC

Expand a 2D grid sampled on the Gauss-Legendre quadrature nodes into spherical harmonics.

# Usage

call SHExpandGLQC (`cilm`, `lmax`, `gridglq`, `w`, `plx`, `zero`, `norm`, `csphase`, `lmax_calc`, `exitstatus`)

# Parameters

`cilm` : output, complex(dp), dimension (2, `lmax`+1, `lmax`+1) or (2, `lmax_calc`+1, `lmax_calc`+1)
:   The complex spherical harmonic coefficients of the complex function. The first index specifies the coefficient corresponding to the positive and negative order of `m`, respectively, with `Clm=cilm(1,l+1,m+1)` and `Cl,-m =cilm(2,l+1,m+1)`.

`lmax` : input, integer(int32)
:   The spherical harmonic bandwidth of the grid. If `lmax_calc` is not specified, this also corresponds to the maximum spherical harmonic degree of the coefficients `cilm`.

`gridglq` : input, complex(dp), dimension(`lmax`+1, 2\*`lmax`+1)
:   A 2D grid of complex data sampled on the Gauss-Legendre quadrature nodes. The latitudinal nodes correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees. See also `GLQGridCoord>`.

`w` : input, real(dp), dimension (`lmax`+1)
:   The Gauss-Legendre quadrature weights used in the integration over latitude. These are obtained from a call to `SHGLQ`.

`plx` : input, optional, real(dp), dimension (`lmax`+1, (`lmax`+1)*(`lmax`+2)/2)
:   An array of the associated Legendre functions calculated at the Gauss-Legendre quadrature nodes. These are determined from a call to `SHGLQ` with the option `cnorm=1`. Either `plx` or `zero` must be present, but not both.

`zero` : input, optional, real(dp), dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude, calculated by a call to `SHGLQ`.  Either `plx` or `zero` must be present, but not both.

`norm` : input, optional, integer(int32), default = 1
:   1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer(int32), default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`lmax_calc` : input, optional, integer(int32), default = `lmax`
:   The maximum spherical harmonic degree calculated in the spherical harmonic expansion.

`exitstatus` : output, optional, integer(int32)
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHExpandGLQC` will expand a 2-dimensional grid of complex data sampled on the Gauss-Legendre quadrature nodes into complex spherical harmonics. This is the inverse of the routine `MakeGridGLQC`. The latitudinal nodes of the input grid correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees. It is implicitly assumed that the function is bandlimited to degree `lmax`. If the optional parameter `lmax_calc` is specified, the spherical harmonic coefficients will be calculated up to this degree, instead of `lmax`.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m. The normalized legendre functions are calculated in this routine using the scaling algorithm of Holmes and Featherstone (2002), which are accurate to about degree 2800. The unnormalized functions are accurate only to about degree 15. 

The spherical harmonic transformation may be speeded up by precomputing the Legendre functions on the Gauss-Legendre quadrature nodes in the routine `SHGLQ` with the optional parameter `cnorm` set to 1. However, given that this array contains on the order of `lmax`**3 entries, this is only feasible for moderate values of `lmax`.

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and
order normalised associated Legendre functions, J. Geodesy, 76, 279-
299, 2002.

# See also

[makegridglqc](makegridglqc.html), [shexpandglq](shexpandglq.html), [makegridglq](makegridglq.html), [shexpanddh](shexpanddh.html), [makegriddh](makegriddh.html), [shexpanddhc](shexpanddhc.html), [makegriddhc](makegriddhc.html), [shexpandlsq](shexpandlsq.html), [glqgridcoord](glqgridcoord.html), [shglq](shglq.html)
