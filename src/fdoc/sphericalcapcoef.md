# SphericalCapCoef

Calculate the spherical harmonic coefficients of a spherical cap.

# Usage

call SphericalCapCoef (`coef`, `theta`, `lmax`, `exitstatus`)

# Parameters

`coef` : output, real\*8, dimension(`lmaxin`+1)
:   The zonal spherical harmonic coefficients of a spherical cap centered over the north pole.

`theta` : input, real\*8
:   The angular radius of the spherical cap in radians.

`lmax` : optional, input, integer, default = `lmaxin`
:   The maximum spherical harmonic degree to calculate the spherical harmonic coefficients.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SphericalCapCoef` will calculate the spherical harmonic coefficients of a spherical cap centered over the north pole. The zonal coefficients, returned in the array `coef`, are normalized such that the degree-0 term is 1, and are to be used with either the geodesy 4-pi normalized or orthonormalized spherical harmonics.
