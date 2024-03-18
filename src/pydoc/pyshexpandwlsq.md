# SHExpandWLSQ()

Expand a set of irregularly sampled data points into spherical harmonics using a weighted least squares inversion.

# Usage

cilm, chi2 = SHExpandWLSQ (d, w, lat, lon, lmax, [norm,  csphase])

# Returns

cilm : float, dimension (2, lmax+1, lmax+1)
:   The real spherical harmonic coefficients of the function. The coefficients C0lm and C1lm refer to the cosine (Clm) and sine (Slm) coefficients, respectively, with Clm=cilm[0,l,m] and Slm=cilm[1,l,m].

chi2 : float
:   The residual sum of squares misfit for an overdetermined inversion.

# Parameters

d : float, dimension (nmax)
:   The value of the function at the coordinates (lat, lon).

w : float, dimension (nmax)
:   The weights used in the weighted least squares inversion.

lat : float, dimension (nmax)
:   The latitude in DEGREES corresponding to the value in d.

lon : float, dimension (nmax)
:   The longitude in DEGREES corresponding to the value in d.

lmax : integer
:   The maximum spherical harmonic degree of the output coefficients cilm.

norm : optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

csphase : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

# Description

SHExpandWLSQ will expand a set of irregularly sampled data points into spherical harmonics by a weighted least squares inversion. For this problem, there must be more data points than spherical harmonic coefficients (i.e., nmax>(lmax+1)**2). It is assumed that each measurement is statistically independent (i.e., the weighting matrix is diagonal), and the inversion is performed using the LAPACK routine DGGGLM.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments norm and csphase; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.
