# spharm

Compute all the spherical harmonic functions up to a maximum degree and order.

# Usage

`ylm` = spharm (`lmax`, `theta`, `phi`, [`normalization`, `kind`, `csphase`, `packed`, `degrees`])

# Returns

`ylm` : float or complex, dimension (2, `lmax`+1, `lmax`+1) or (2, (`lmax`+1)\*(`lmax`+2)/2)
:   An array of spherical harmonic functions, ylm[i, l, m], where `l` and `m` are the spherical harmonic degree and (positive) order, respectively. The index `i` provides the positive (0) and negative (1) order. If `packed` is True, the array is 2-dimensional with the index of the second column corresponding to `l\*(l+1)/2+m`.

# Parameters

`lmax` : integer
:   The maximum degree of the spherical harmonic functions to be computed.

`theta` : float
:   The colatitude in degrees.

`phi` : float
:   The longitude in degrees.

`normalization` : str, optional, default = '4pi'
:   '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized spherical harmonic functions, respectively.

`kind` : str, optional, default = 'real'
:   'real' or 'complex' spherical harmonic coefficients.

`csphase` : optional, integer, default = 1
:   If 1 (default), the Condon-Shortley phase will be excluded. If -1, the Condon-Shortley phase of (-1)^m will be appended to the spherical harmonic functions.

`packed` : optional, bool, default = False
:   If True, return a 2-dimensional packed array where the index of the second column corresponds to `l\*(l+1)/2+m`, where `l` and `m` are respectively the degree and order.

`degrees` : optional, bool, default = True
:   If True, `colat` and `phi` are expressed in degrees.

# Description

`spharm` will calculate all of the spherical harmonic functions up to degree `lmax` for a given colatitude `theta` and longitude `phi`. Three parameters determine how the spherical harmonic functions are defined. `normalization` can be either '4pi' (default), 'ortho', 'schmidt', or 'unnorm' for 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized spherical harmonic functions, respectively. `kind` can be either 'real' or 'complex', and `csphase` determines whether to include or exclude (default) the Condon-Shortley phase factor.

By default, the routine will return a 3-dimensional array, ylm[i, l, m], where `l` and `m` are the spherical harmonic degree and (positive) order, respectively. The index `i=0` corresponds to the positive orders, whereas `i=1` corresponds to the negative orders. If the optional parameter `packed` is set to True, the output will instead be a 2-dimensional array where the indices of the second column correspond to `l\*(l+1)/2+m`.

The spherical harmonic functions are calculated using the standard three-term recursion formula, and in order to prevent overflows, the scaling approach of Holmes and Featherstone (2002) is utilized. The resulting functions are accurate to about degree 2800. See Wieczorek and Meschede (2018) for exact definitions on how the spherical harmonic functions are defined.

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, J. Geodesy, 76, 279-299, doi:10.1007/s00190-002-0216-2, 2002.

Wieczorek, M. A., and M. Meschede. SHTools — Tools for working with spherical harmonics, Geochem., Geophys., Geosyst., 19, 2574-2592, doi:10.1029/2018GC007529, 2018.

# See also

[pyspharm_lm](pyspharm_lm.html), [pylegendre_lm](pylegendre_lm.html), [pylegendre](pylegendre.html)

