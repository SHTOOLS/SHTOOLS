# MakeGrid2D

Create a 2D cylindrical map of arbitrary grid spacing from a set of spherical harmonic coefficients.

# Usage

`grid` = MakeGrid2D (`cilm`, `interval`, [`lmax`, `norm`, `csphase`, `f`, `a`, `north`, `south`, `east`, `west`, `dealloc`])

# Returns

`grid` : float, dimension ((`north`-`south`)/`interval`+1, (`east`-`west`)/`interval`+1)
:   A 2D equally spaced map of the input spherical harmonic coefficients `cilm`. The  array is in raster format with upper-left and lower-right coordinates of (90 N, 0 E) and (90 S, 360 E), respectively.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real spherical harmonic coefficients to be expanded in the space domain. The coefficients `C1lm` and `C2lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm[0,l,m]` and `Slm=cilm[1,l,m]`. 

`interval` : float
:   The latitudinal and longitudinal spacing of `grid`.

`lmax` :optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the coefficients `cilm` used when calculating the grid.

`norm` : optional, integer, default = 1
:   1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`f` : optional, float, default = 0
:   The flattening of the reference ellipoid that is subtracted from the function. This is given by (`R_equator-R_pole)/R_equator`.

`a` : optional, float, default = 0
:   The semi-major axis of the reference ellispoid that is subtracted from the function.

`north` : optional, float, default = 90
:   The maximum latitude of the output raster grid, in degrees. The default is 90 degrees.

`south` : optional, float, default = -90
:   The minimum latitude of the output raster grid, in degrees. The default is -90 degrees.

`east` : optional, float, default = 360
:   The maximum longitude of the output raster grid, in degrees. The default is 360 degrees.

`west` : optional, float, default = 0
:   The minimum longitude of the output raster grid, in degrees. The default is 0 degrees.

`dealloc` : optional, integer, default = 0
:   0 (default) = Save variables used in the external Legendre function calls. (1) Deallocate this memory at the end of the funcion call.

# Description

`MakeGrid2D` will create a 2-dimensional cylindrical map, equally spaced in (geocentric) latitude and longitude, from a set of input spherical harmonic coefficients. The output grid is in raster format possessing upper-left and lower-right coordinates of (90 N, 0 E) and (90 S, 360 E), repsectively. If the optional parameters `north`, `south`, `east` and `west` are specified, then the output grid will possess upper-left and lower-right coordinates of (`north`, `west`) and (`south`, `east`), repsectively. The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

If the optional arguments `f` and `a` are specified, the output function will be referenced to an ellipsoid with flattening `f` and semimajor axis `a`. The normalized legendre functions are calculated in this routine using the scaling algorithm of Holmes and Featherstone (2002), which are accurate to about degree 2800. The unnormalized functions are accurate only to about degree 15. 

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, J. Geodesy, 76, 279-299, 2002.

# See also

[makegriddh](pymakegriddh.html), [makegriddhc](pymakegriddhc.html), [makegridglq](pymakegridglq.html), [makegridglqc](pymakegridglqc.html), [makegravgriddh](pymakegravgriddh.html), [makemaggriddh](pymakemaggriddh.html)
