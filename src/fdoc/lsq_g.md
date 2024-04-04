# LSQ_G

Compute the matrix G that is used when computing spherical harmonic coefficients by least squares inversion.

# Usage

call LSQ_G (`g`, `lat`, `lon`, `nmax`, `lmax`, `norm`, `csphase`, `exitstatus`)

# Parameters

`g` : output, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The matrix G.

`lat` : input, real(dp), dimension (`nmax`)
:   The latitude in degrees of the data points.

`lon` : input, real(dp), dimension (`nmax`)
:   The longitude in degrees of the data pointss.

`nmax` : input, integer(int32)
:   The number of data points.

`lmax` : input, integer(int32)
:   The maximum spherical harmonic degree of the inversion.

`norm` : input, optional, integer(int32), default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer(int32), default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`exitstatus` : output, optional, integer(int32)
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`LSQ_G` will compute the matrix G that is used when computing the spherical harmonic coefficients of an irregularly sampled function by least squares inversion, as used with SHExpandLSQ. The matrix G has dimension (nmax, (lmax+1)**2) where nmax is the number of data points and lmax is the maximum spherical harmonic degree of the expansion. Each element in a given row corresponds to the values of the spherical harmonic functions for a given latitude and longitude. The elements in each row are ordered by increasing degree, with all cosine terms for a given degree followed by all sin terms for the same degree (with increasing order).

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

# See also

[shexpandlsq](shexpandlsq.html), [makegriddh](makegriddh.html), [shexpanddh](shexpanddh.html), [makegriddhc](makegriddhc.html), [shexpanddhc](shexpanddhc.html), [makegridglq](makegridglq.html), [shexpandglq](shexpandglq.html), [makegridglqc](makegridglqc.html), [shexpandglqc](shexpandglqc.html), dgels(1)
