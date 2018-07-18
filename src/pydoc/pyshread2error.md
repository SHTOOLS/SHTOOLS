# SHRead2Error

Read spherical harmonic coefficients from a CHAMP or GRACE-like ascii-formatted file.

# Usage

`cilm`, `error`, `lmax`, `gm`, `r0_pot`, `dot`, `doystart`, `doyend`, `epoch` = SHRead2Error (`filename`, `lmaxin`)

# Returns

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The spherical harmonic coefficients contained in `filename`.

`error` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The uncertainty of the spherical harmonic coefficients contained in `filename`.

`lmax` : integer
:   The maximum spherical harmonic degree of `cilm`. This is the minimum of the dimension of `cilm` and `lmaxin`.

`gm` : float
:   The mass-gravitational constant.

`r0_pot` : float
:   The reference radius of the potential coefficients.

`dot` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The time derivatives of the spherical harmonic coefficients. The dimension of this array can be smaller than `lmax+1`.

`doystart` : float
:   The starting time of the solution.

`doyend` : float
:   The ending time of the solution

`epoch` : float
:   The epoch time for the time derivates.

# Parameters

`filename` : character(*)
:   The ascii-formatted filename containing the spherical harmonic coefficients.

`lmaxin` : integer
:   This spherical harmonic degree controls the dimension of the output array `cilm`. The coefficients between `lmax+1` and `lmaxin` will be set to zero.

# Description

`SHRead2Error` will read spherical harmonic coefficients from a CHAMP or GRACE-like ascii-formatted file into an array `cilm`. The errors and time derivatives associated with the coefficients will be read if the optional arrays `error` and `dot` are specified, respectively. The staring and ending date of the solution are specified by the optional parameters `doystart` and `doyend`, and the epoch of the time derivates is specified by the optional parameter `epoch`. The maximum spherical harmonic degree is read from the file, and the file does not need to be ordered by angular degree or order.

This routine does not read and output all parameters in the file. Records that are read (and at least partially output) include: `EARTH`, `GGM`, `SHM`, `GRCOF2`, `GRDOTA`, `CALSDV`, `gfc`, `gfct`, and `dot`. Comments specified by the record `CMMNT` will be print out to the screen, as will the record names that are not currently implemented.

Each line of the file starts with a character string describing what follows. 

- `EARTH` or `GGM`: `GM`, `R0_POT`
- `SHM`: Maximum spherical harmonic degree of file.
- `GRCOF2`, `CALSDV`, or `gfc`: spherical harmonic coefficients, formatted as (`l`, `m`, `clm`, `slm`) or (`l`, `m`, `clm`, `slm`, `clm_error`, `slm_error`).

# See also

[shread](pyshread.html), [shread2](pyshread2.html), [shreadjpl](pyshreadjpl.html)
