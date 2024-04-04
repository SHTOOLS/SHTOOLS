# LSQ_G()

Compute the matrix G that is used when computing spherical harmonic coefficients by least squares inversion.

# Usage

g = LSQ_G (lat, lon, lmax, [norm,  csphase])

# Returns

g : float, dimension (nmax, (lmax+1)**2)
:   The matrix G.

# Parameters

lat : float, dimension (nmax)
:   The latitude in degrees of the data points.

lon : float, dimension (nmax)
:   The longitude in degrees of the data points.

lmax : integer
:   The maximum spherical harmonic degree of the inversion.

norm : optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

csphase : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

# Description

LSQ_G will compute the matrix G that is used when computing the spherical harmonic coefficients of an irregularly sampled function by least squares inversion, as used with SHExpandLSQ. The matrix G has dimension (nmax, (lmax+1)**2) where nmax is the number of data points and lmax is the maximum spherical harmonic degree of the expansion. Each element in a given row corresponds to the values of the spherical harmonic functions for a given latitude and longitude. The elements in each row are ordered by increasing degree, with all cosine terms for a given degree followed by all sin terms for the same degree (with increasing order).

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments norm and csphase; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.
