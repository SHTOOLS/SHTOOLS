# SHReadJPLError

Read spherical harmonic coefficients from a JPL ascii-formatted file.

# Usage

`cilm`, `error`, `lmax`, `gm` = SHReadJPLError (`filename`, `lmaxin`, [`formatstring`])

# Returns

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The spherical harmonic coefficients contained in `filename`.

`error` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The uncertainties associated with the spherical harmonic coefficients contained in `filename`.

`lmax` : integer
:   The maximum spherical harmonic degree of `cilm`.

`gm` : float, dimension(2)
:   The mass-gravitational constant and error.

# Parameters

`filename` : character(*)
:   The filename of the JPL ascii-formatted spherical harmonic coefficients.

`lmaxin` : integer
:   This spherical harmonic degree controls the dimension of the output array `cilm`. The coefficients between `lmax+1` and `lmaxin` will be set to zero.

`formatstring` : character*6, optional, default = "E19.12"
:   The format string used to read the elements of `cilm` and `error`. The default is "E19.12".

# Description

`SHReadJPL` will read spherical harmonic coefficients from a JPL ascii formatted file into an array `cilm`. The real numbers are assumed to be formated with the specifier `E19.12`, but this can be changed by specifiying the optional string `formatstring`.

The JPL ascii-formatted file is organized as follows:

- Comment lines starting with "#".
- `gm` (if a gravitational potential file)
- A list of `J_l`, which is `-cilm[0,l,0]`.
- A list of the cosine and sine terms.
- The errors of the above (starting at step 2).

# See also

[shread](pyshread.html), [shread2](pyshread2.html), [shreadjpl](pyshreadjpl.html)
