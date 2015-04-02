# SHReadJPL

Read spherical harmonic coefficients from a JPL ascii-formatted file.

# Usage

call subroutine SHReadJPL (`filename`, `cilm`, `lmax`, `error`, `gm`, `formatstring`)

# Parameters

`filename` : input, character(*)
:   The filename of the JPL ascii formatted spherical harmonic coefficients.
	
`cilm` : output, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients contained in `filename`.

`lmax` : input, integer
:   The maximum spherical harmonic degree of `cilm`.

`error` : output, real\*8, optional, dimension (2, `lmax`+1, `lmax`+1)
:   The errors corresponding to the spherical harmonic coefficients `cilm`.

`gm` : output, real\*8, optional, dimension(2)
:   The mass-gravitational constant and error.

`formatstring` : input character*6, optional, default = "E19.12"
:   The format string used to read the elements of `cilm` and `cilm_error`. The default is "E19.12".

# Description

`SHReadJPL` will read spherical harmonic coefficients from a JPL ascii formatted file into an array `cilm`. The maximum spherical harmonic degree `lmax` must be known a priori. The errors associated with the coefficients `cilm` will be read if the optional array `error` is specified. The real numbers are assumed to be formated with the specifier `E19.12`, but this can be changed by specifiying the optional string `formatstring`. If the optional parameter `gm` is specified, the mass-gravitational constant and error will be output, if present.

The JPL ascii formatted file is organized as follows:

- Comment lines starting with "#".
- `gm` (if a gravitational potential file)
- A list of `J_l`, which is `-cilm(1,l+1,1)`.
- A list of the cosine and sine terms.
- The errors of the above (starting at step 2).

# See also

[shread](shread.html), [shread2](shread2.html)
