# SHExpandLSQ

Determine the spherical harmonic coefficients of an irregularly sampled function using a (weighted) least squares inversion.

# Usage

call SHExpandLSQ (`cilm`, `d`, `lat`, `lon`, `nmax`, `lmax`, `norm`, `chi2`, `csphase`, `weights`, `g`, `exitstatus`)

# Parameters

`cilm` : output, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The real spherical harmonic coefficients of the function. The coefficients `C1lm` and `C2lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm(1,l+1,m+1)` and `Slm=cilm(2,l+1,m+1)`.

`d` : input, real(dp), dimension (`nmax`)
:   The value of the function at the coordinates (`lat`, `lon`).

`lat` : input, real(dp), dimension (`nmax`)
:   The latitude in degrees corresponding to the value in `d`.

`lon` : input, real(dp), dimension (`nmax`)
:   The longitude in degrees corresponding to the value in `d`.

`nmax` : input, integer(int32)
:   The number of data points.

`lmax` : input, integer(int32)
:   The maximum spherical harmonic degree of the output coefficients `cilm`.

`norm` : input, optional, integer(int32), default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`chi2` : output, optional, real(dp)
:   The residual (weighted) sum of squares misfit.

`csphase` : input, optional, integer(int32), default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`weights` : input, real(dp), dimension (`nmax`)
:   The weights to be applied in a weighted least squares inversion.

`g` : input, real(dp), dimension (`nmax`, `(lmax+1)**2`)
:   The matrix G that is used when performing the least squares inversion

`exitstatus` : output, optional, integer(int32)
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHExpandLSQ` will determine the spherical harmonic coefficients of an irregularly sampled function using a least squares inversion. When the number of data points is greater or equal to the number of spherical harmonic coefficients (i.e., `nmax>=(lmax+1)**2`), the solution of the overdetermined system will be determined. If there are more coefficients than data points, then the solution of the underdetermined system that minimizes the solution norm will be determined. The inversions are performed using the LAPACK routine DGELS.

A weighted least squares inversion will be performed if the optional vector `weights` is specified. The weights should be set equal to the inverse of the data variance, and it is assumed explicitly that each measurement is statistically independent (i.e., the weighting matrix is diagonal). The weighted least squares inversion must be overdetermined, and the inversion is performed using the LAPACK routine DGELS after scaling the data vector and inversion matrix.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m. If the subroutine is to be called several times using the same lat and lon coordinates, the matrix G can be precomputed using the routine LSQ_G.

# See also

[lsq_g](lsq_g.html), [makegriddh](makegriddh.html), [shexpanddh](shexpanddh.html), [makegriddhc](makegriddhc.html), [shexpanddhc](shexpanddhc.html), [makegridglq](makegridglq.html), [shexpandglq](shexpandglq.html), [makegridglqc](makegridglqc.html), [shexpandglqc](shexpandglqc.html), dgels(1)
