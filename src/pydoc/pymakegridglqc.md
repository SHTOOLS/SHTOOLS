# MakeGridGLQC

Create a 2D complex map from a set of complex spherical harmonic coefficients sampled on the Gauss-Legendre quadrature nodes.

# Usage

`gridglq` = MakeGridGLQC (`cilm`, `zero`, [`lmax`, `norm`, `csphase`, `lmax_calc`])

# Returns

`gridglq` : complex, dimension (`lmax`+1, 2\*`lmax`+1)
:   A 2D complex map of the function sampled on the Gauss-Legendre quadrature nodes.

# Parameters

`cilm` : complex, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex spherical harmonic coefficients of the function. When evaluating the function, the maximum spherical harmonic degree considered is the minimum of `lmax`, `lmaxin`, or `lmax_calc` (if specified). The first index specifies the coefficient corresponding to the positive and negative order of `m`, respectively, with `Clm=cilm[0,l,m+]` and `Cl,-m=cilm[1,l,m]`.

`zero` : float, dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude, calculated by a call to `SHGLQ`.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic bandwidth of the function. This determines the sampling nodes and dimensions of the output grid.

`norm` : optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`lmax_calc` : optional, integer, default = `lmax`
:   The maximum spherical harmonic degree used in evaluating the function. This must be less than or equal to `lmax`.

# Description

`MakeGridGLQC` will create a 2-dimensional complex map from a set of input complex spherical harmonic coefficients sampled on the Gauss-Legendre quadrature nodes. This is the inverse of the routine `SHExpandGLQC`. The latitudinal nodes correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees. When evaluating the function, the maximum spherical harmonic degree that is considered is the minimum of `lmax`, `lmaxin`, or `lmax_calc` (if specified).

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m. The normalized legendre functions are calculated using the scaling algorithm of Holmes and Featherstone (2002), which are accurate to about degree 2800. The unnormalized functions are accurate only to about degree 15.

The reconstruction of the spherical harmonic function may be speeded up by precomputing the Legendre functions on the Gauss-Legendre quadrature nodes in the routine `SHGLQ` with the optional parameter `cnorm` set to 1. However, given that this array contains on the order of `lmax`**3 entries, this is only feasible for moderate values of `lmax`.

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and order normalised associated Legendre functions, J. Geodesy, 76, 279-299, 2002.

# See also

[shexpandglqc](pyshexpandglqc.html), [shexpandglq](pyshexpandglq.html), [makegridglq](pymakegridglq.html), [shexpanddh](pyshexpanddh.html), [makegriddh](pymakegriddh.html), [shexpanddhc](pyshexpanddhc.html), [makegriddhc](pymakegriddhc.html), [shexpandlsq](pyshexpandlsq.html), [glqgridcoord](pyglqgridcoord.html), [shglq](pyshglq.html)
