# MakeGridDH

Create a 2D map from a set of spherical harmonic coefficients using the Driscoll and Healy (1994) sampling theorem.

# Usage

`griddh` = MakeGridDH (`cilm`, [`lmax`, `norm`, `sampling`, `csphase`, `lmax_calc`])

# Returns

`griddh` : float, dimension (2\*`lmax`+2, `sampling`\*(2\*`lmax`+2))
:   A 2D equally sampled (default) or equally spaced map in degrees of the spherical harmonic coefficients `cilm` that conforms to the sampling theorem of Driscoll and Healy (1994). The first latitudinal band corresponds to 90 N, the latitudinal band for 90 S is not included, and the latitudinal sampling interval is 180/`n` degrees, where `n` is 2\*`lmax`+2. The first longitudinal band is 0 E, the longitudinal band for 360 E is not included, and the longitudinal sampling interval is 360/`n` for an equally sampled and 180/`n` for an equally spaced grid, respectively.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real spherical harmonic coefficients of the function. The coefficients `cilm[0,l,m]` and `cilm[1,l,m]` refer to the "cosine" (`Clm`) and "sine" (`Slm`) coefficients, respectively.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the function, which determines the sampling of the output grid.

`norm` : optional, integer, default = 1
:   1 = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics;  4 = orthonormal harmonics.

`sampling` : optional, integer, default = 1
:   If 1, the output grid contains the same number of samples in latitude as in longitude. If 2, the grid is equally spaced in degrees, having twice as many samples in longitude as latitude.

`csphase` : optional, integer, default = 1
:   1 = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`lmax_calc` : optional, integer, default = `lmax`
:   The maximum spherical harmonic degree used in evaluating the  function. This must be less than or equal to `lmax`, and does not affect the number of samples of the output grid.

# Description

`MakeGridDH` will create a 2-dimensional map equally sampled or equally spaced in latitude and longitude from a set of input spherical harmonic coefficients. This grid conforms with the sampling theorem of Driscoll and Healy (1994) and this routine is the inverse of SHExpandDH. The function is evaluated at each longitudinal band by inverse Fourier transforming the sin and cos terms for each degree `l`, and then summing over all degrees. When evaluating the function, the maximum spherical harmonic degree that is considered is the minimum of `lmaxin`, `lmax`, and `lmax_calc` (if specified).

The default is to use an input grid that is equally sampled (`n` by `n`), but this can be changed to use an equally spaced grid (`n` by 2`n`) by use of the optional parameter `sampling`. The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

The normalized legendre functions are calculated using the scaling algorithm of Holmes and Featherstone (2002), which are accurate to about degree 2800. The unnormalized functions are accurate only to about degree 15.

# References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, J. Geodesy, 76, 279- 299, 2002.

# See also

[shexpanddh](pyshexpanddh.html), [makegriddhc](pymakegriddhc.html), [shexpanddhc](pyshexpanddhc.html), [makegridglq](pymakegridglq.html), [shexpandglq](pyshexpandglq.html), [makegridglqc](pymakegridglqc.html), [shexpandglqc](pyshexpandglqc.html), [makegrid2d](pymakegrid2d.html)
