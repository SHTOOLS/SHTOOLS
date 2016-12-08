# SHConfidence

Compute the probability that two functions are correlated at a given spherical harmonic degree for a given correlation coefficient.

# Usage

`prob` = SHConfidence (`l`, `corr`)

# Returns

`prob` : float
:   Probability that two functions expressed in spherical coefficients with spectral correlation `corr` are correlated at degree `l`.

# Parameters

`l` :  integer
:   The spherical harmonic degree.

`corr` : float
:   The correlation coefficient of the two data sets at degree `l`.

# Description

`SHConfidence` will calculate the probability (between 0 and 1) that two sets of spherical harmonic coefficients with spectral correlation `corr` are linearly correlated at a given degree. This is calculated using equation A7 from Pauer et al. (2006).

# References

Pauer, M, K. Fleming, and O. Cadek, Modeling the dynamic component of the geoid and topography of Venus, J. Geophys. Res., 111, E11012, doi:10.1029/2005JE002511, 2006.

# See also

[shadmitcorr](pyshadmitcorr.html), [spectrum](spectrum.html), [cross_spectrum](cross_spectrum.html)
