# SHRead2

Read spherical harmonic coefficients from a CHAMP or GRACE-like ascii-formatted file.

# Usage

call SHRead2 (`filename`, `cilm`, `lmax`, `gm`, `r0_pot`, `error`, `dot`, `doystart`, `doyend`, `epoch`, `exitstatus`)

# Parameters

`filename` : input, character(*)
:   The ascii-formatted filename containing the spherical harmonic coefficients.

`cilm` : output, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients contained in `filename`.

`lmax` : output, integer
:   The maximum spherical harmonic degree of `cilm`.

`gm` : output, real\*8
:   The mass-gravitational constant.

`r0_pot` : output, real\*8
:   The reference radius of the potential coefficients.

`error` : output, optional, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The errors corresponding to the spherical harmonic coefficients `cilm`.

`dot` : output, optional, real\*8, dimension (2, `lmaxout`+1, `lmaxout`+1)
:   The time derivatives of the spherical harmonic coefficients. The dimension of this array can be smaller than `lmax+1`.

`doystart` : output, optional, real\*8
:   The starting time of the solution.

`doyend` : output, optional, real\*8
:   The ending time of the solution

`epoch` : output, optional, real\*8
:   The epoch time for the time derivates.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHRead2` will read spherical harmonic coefficients from a CHAMP or GRACE-like ascii-formatted file into an array `cilm`. The errors and time derivatives associated with the coefficients will be read if the optional arrays `error` and `dot` are specified, respectively. The staring and ending date of the solution are specified by the optional parameters `doystart` and `doyend`, and the epoch of the time derivates is specified by the optional parameter `epoch`. The maximum spherical harmonic degree is read from the file, and the file does not need to be ordered by angular degree or order.

This routine does not read and output all parameters in the file. Records that are read (and at least partially output) include: `EARTH`, `GGM`, `SHM`, `GRCOF2`, `GRDOTA`, `CALSDV`, `gfc`, `gfct`, and `dot`. Comments specified by the record `CMMNT` will be print out to the screen, as will the record names that are not currently implemented.

Each line of the file starts with a character string describing what follows. 

- `EARTH` or `GGM`: `GM`, `R0_POT`
- `SHM`: Maximum spherical harmonic degree of file.
- `GRCOF2`, `CALSDV`, or `gfc`: spherical harmonic coefficients, formatted as (`l`, `m`, `clm`, `slm`) or (`l`, `m`, `clm`, `slm`, `clm_error`, `slm_error`).

# See also

[shread](shread.html), [shreadjpl](shreadjpl.html)
