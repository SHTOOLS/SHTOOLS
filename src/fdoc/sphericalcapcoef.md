# SphericalCapCoef

Calculate the spherical harmonic coefficients of a spherical cap.

# Usage

call SphericalCapCoef (`coef`, `theta`, `lmax`)

# Parameters

`coef` : output, real\*8, dimension(`lmaxin`+1)
:   The zonal spherical harmonic coefficients of a spherical cap centered over the north pole.
	
`theta` : input, real\*8
:   The angular radius of the spherical cap in radians.

`lmax` : optional, input, integer, default = `lmaxin`
:   The maximum spherical harmonic degree to calculate the spherical harmonic coefficients.

# Description

`SphericalCapCoef` will calculate the spherical harmonic coefficients of a spherical cap centered over the north pole. The zonal coefficients, returned in the array `coef`, are normalized such that the degree-0 term is 1, and are to be used with either the geodesy 4-pi normalized or orthonormalized spherical harmonics.
