# SHReadH

Read spherical harmonic coefficients from an ascii-formatted file with a header line.

# Usage

`cilm`, `lmax`, `header` = SHRead (`filename`, `lmaxin`, `nhead`, [`skip`])

# Returns

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The spherical harmonic coefficients contained in `filename`. Note that dimensions of the output array are determined by `lmaxin`, and not the output `lmax`.

`lmax` : integer
:   The maximum spherical harmonic degree of `cilm`. This is the minimum of the maximum spherical harmonic degree of `filename` and the dimension of `cilm`-1.

`header` : float, dimension (`nhead`)
:   Array of the first `nhead` values in the file's initial header line.

# Parameters

`filename` : character(:)
:   The filename of the ascii file containing the spherical harmonic coefficients.

`lmaxin` : integer
:   This spherical harmonic degree controls the dimension of the output array `cilm`. The coefficients between `lmax+1` and `lmaxin` will be set to zero.

`nhead` : integer 
:   Number of values to read from the file's initial header line.

`skip` : optional, integer, default = 0
:   The number of lines to skip before parsing `filename`.

# Description

`SHRead` will read spherical harmonic coefficients from an ascii-formatted file into an array `cilm`. The maximum spherical harmonic degree that is read is determined by the minimum of the dimension of the input array `cilm`-1 and the maximum degree of the coefficients in the file. If the optional array `skip` is specified, parsing of the file will commence after the first `skip` lines.

The spherical harmonic coefficients in the file are assumed to be ordered by increasing degree `l` and angular order `m` according to the format

`l, m, cilm[0,l,m], cilm[1,l,m]`

The ordering of the file is explcitly given by

`l, 0 / l, 1 / l, 2 /l, ... / l, m / l+1, 0 / l+1, 1 / ...`

The first spherical harmonic degree of the filename does not have to be 0; this is determined from the first element after the `skip` and `header` lines. 

# See also

[shread](pyshread.html), [shreaderror](pyshreaderror.html), [shreaderrorh](pyshreaderrorh.html), [shread2](pyshread2.html), [shread2error](pyshread2error.html), [shreadjpl](pyshreadjpl.html) [shreadjplerror](pyshreadjplerror.html)
