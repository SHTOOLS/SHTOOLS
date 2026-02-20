# spharm_lm()

Compute the spherical harmonic function for specific degrees and orders.

# Usage

ylm = spharm_lm (l, m, theta, phi, [normalization, kind, csphase, degrees])

# Returns

ylm : float or complex, ndarray
The spherical harmonic function ylm, where l and m are the spherical
harmonic degree and order, respectively.

# Parameters

l : integer, array_like
The spherical harmonic degree.

m : integer, array_like
The spherical harmonic order.

theta : float, array_like
The colatitude in degrees. Use radians if 'degrees' is set to False.

phi : float, array_like
The longitude in degrees. Use radians if 'degrees' is set to False.

normalization : str, array_like, optional, default = '4pi'
'4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
orthonormalized, Schmidt semi-normalized, or unnormalized spherical
harmonic functions, respectively.

kind : str, array_like, optional, default = 'real'
'real' or 'complex' spherical harmonic coefficients.

csphase : integer, array_like, optional, default = 1
If 1 (default), the Condon-Shortley phase will be excluded. If -1, the
Condon-Shortley phase of (-1)^m will be appended to the spherical
harmonic functions.

degrees : bool, array_like, optional, default = True
If True, `theta` and `phi` are expressed in degrees.

# Notes

spharm_lm will calculate the spherical harmonic function for specific
degrees l, orders m, colatitudes theta and longitudes phi. Three parameters
determine how the spherical harmonic functions are defined. normalization
can be either '4pi' (default), 'ortho', 'schmidt', or 'unnorm' for 4pi
normalized, orthonormalized, Schmidt semi-normalized, or unnormalized
spherical harmonic functions, respectively. kind can be either 'real' or
'complex', and csphase determines whether to include or exclude (default)
the Condon-Shortley phase factor.

The spherical harmonic functions are calculated using the standard
three-term recursion formula, and in order to prevent overflows, the
scaling approach of Holmes and Featherstone (2002) is utilized.
The resulting functions are accurate to about degree 2800. See Wieczorek
and Meschede (2018) for exact definitions on how the spherical harmonic
functions are defined.

# References

Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
summation and the recursive computation of very high degree and order
normalised associated Legendre functions, J. Geodesy, 76, 279-299,
doi:10.1007/s00190-002-0216-2, 2002.

Wieczorek, M. A., and M. Meschede. SHTools — Tools for working with
spherical harmonics, Geochem., Geophys., Geosyst., 19, 2574-2592,
doi:10.1029/2018GC007529, 2018.
