# SHRead

Read spherical harmonic coefficients from an ascii-formatted file.

# Usage

call SHRead (`filename`, `cilm`, `lmax`, `skip`, `header`, `error`, `exitstatus`)

# Parameters

`filename` : input, character(:)
:   The filename of the ascii file containing the spherical harmonic coefficients.

`cilm` : output, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients contained in `filename`.

`lmax` : output, integer
:   The maximum spherical harmonic degree of `cilm`. This is the minimum of the maximum spherical harmonic degree of `filename` and the dimension of `cilm`-1.

`skip` : input, optional, integer
:   The number of lines to skip before parsing `filename`.

`header` : output, optional, real\*8 dimension (`n`)
:   A vector containing the first `n` numbers in the first line of the file (following any skipped lines). 

`error` : output, optional, real\*8 dimension (2, `lmax`+1, `lmax`+1)
:   The errors corresponding to the spherical harmonic coefficients `cilm`.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHRead` will read spherical harmonic coefficients from an ascii-formatted file into an array `cilm`. The maximum spherical harmonic degree that is read is determined by the minimum of the dimension of the input array `cilm`-1 and the maximum degree of the coefficients in the file. If the optional array `skip` is specified, parsing of the file will commence after the first `skip` lines. If the optional array `header` is specified, then the first `n` elements after the skipped lines will be output, where `n` is the length of the array `header`.

The spherical harmonic coefficients in the file are assumed to be ordered by increasing degree `l` and angular order `m` according to the format

`l, m, cilm(1,l+1,m+1), cilm(2,l+1,m+1)`

The actual delimeters (commas, spaces, or tabs) are unimportant. If the optional array `error` is specified, then the error for each coefficient will be read according to the format

`l, m, cilm(1,l+1,m+1), cilm(2,l+1,m+1), error(1,l+1,m+1), error(2,l+1,m+1)`

The ordering of the file is explcitly given by

`l, 0 / l, 1 / l, 2 /l, ... / l, m / l+1, 0 / l+1, 1 / ...`

The first spherical harmonic degree of the filename does not have to be 0; this is determined from the first element after the `skip` and `header` lines.

# See also

[shread2](shread2.html), [shreadjpl](shreadjpl.html)
