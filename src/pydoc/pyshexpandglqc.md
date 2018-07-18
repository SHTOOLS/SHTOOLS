# SHExpandGLQC

Expand a 2D grid sampled on the Gauss-Legendre quadrature nodes into spherical harmonics.

# Usage

`cilm` = SHExpandGLQC (`gridglq`, `w`, `zero`, [`norm`, `csphase`, `lmax_calc`])

# Returns

`cilm` : complex, dimension (2, `lmax`+1, `lmax`+1) or (2, `lmax_calc`+1, `lmax_calc`+1)
:   The complex spherical harmonic coefficients of the complex function. The first index specifies the coefficient corresponding to the positive and negative order of `m`, respectively, with `Clm=cilm[0,l,m]` and `Cl,-m =cilm[1,l,m]`.

# Parameters

`gridglq` : complex, dimension (`lmax`+1, 2\*`lmax`+1)
:   A 2D grid of complex data sampled on the Gauss-Legendre quadrature nodes. The latitudinal nodes correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees. See also `GLQGridCoord`.

`w` : float, dimension (`lmax`+1)
:   The Gauss-Legendre quadrature weights used in the integration over latitude. These are obtained from a call to `SHGLQ`.

`zero` : float, dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude, calculated by a call to `SHGLQ`.

`norm` : optional, integer, default = 1
:   1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`lmax_calc` : optional, integer, default = `lmax`
:   The maximum spherical harmonic degree calculated in the spherical harmonic expansion.

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

[makegridglqc](pymakegridglqc.html), [shexpandglq](pyshexpandglq.html), [makegridglq](pymakegridglq.html), [shexpanddh](pyshexpanddh.html), [makegriddh](pymakegriddh.html), [shexpanddhc](pyshexpanddhc.html), [makegriddhc](pymakegriddhc.html), [shexpandlsq](pyshexpandlsq.html), [glqgridcoord](pyglqgridcoord.html), [shglq](pyshglq.html)
