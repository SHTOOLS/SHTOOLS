# shread

Read spherical harmonic coefficients from a text file.

# Usage

`coeffs`, `lmaxout` = shread(`filename`, [`lmax`, `skip`])

`coeffs`, `header`, `lmaxout` = shread(`filename`, `header`=True, [`lmax`, `skip`])

`coeffs`, `errors`, `lmaxout` = shread(`filename`, `error`=True, [`lmax`, `skip`])

`coeffs`, `errors`, `header`, `lmaxout` = shread(`filename`, `error`=True, `header`=True, [`lmax`, `skip`])

# Returns

`coeffs` : ndarray, dimension (2, `lmaxout`+1, `lmaxout`+1)
:   The spherical harmonic coefficients.

`errors` : ndarray, dimension (2, `lmaxout`+1, `lmaxout`+1)
:   The errors associated with the spherical harmonic coefficients.

`header` : list of type str
:   A list of values in the header line found before the start of the spherical harmonic coefficients.

`lmaxout` : int
:   The maximum spherical harmonic degree read from the file.

# Parameters

`filename` : str
:   Filename containing the text-formatted spherical harmonic coefficients.

`lmax` : int, optional, default = None
:   The maximum spherical harmonic degree to read from the file. The default is to read the entire file.

`error` : bool, optional, default = False
:   If True, return the errors associated with the spherical harmonic coefficients as a separate array.

`header` : bool, optional, default = False
:   If True, return a list of values in the header line found before the start of the spherical harmonic coefficients.

`skip` : int, optional, default = 0
:   The number of lines to skip before parsing the file.

# Description

This function will read spherical harmonic coefficients from an ascii-formatted text file. The errors associated with the spherical harmonic coefficients, as well as the values in a single header line, can optionally be read by setting the optional parameters `error` and `header` to True. The optional parameter `skip` specifies how many lines should be skipped before attempting to parse the file, and the optional parameter `lmax` specifies the maximum degree to read from the file. Both real and complex spherical harmonic coefficients are supported.

The spherical harmonic coefficients in the file should be formatted as

`l, m, coeffs[0, l, m], coeffs[1, l, m]`

where l and m are the spherical harmonic degree and order, respectively. The terms coeffs[1, l, 0] can be neglected as they are zero. If the errors are to be read, the line should be formatted as For more information, see `shread`.

If the errors are to be read, the line should be formatted as

`l, m, coeffs[0, l, m], coeffs[1, l, m], errors[0, l, m], errors[1, l, m]`

For each value of increasing l, all the angular orders are listed in inceasing
order, from 0 to l.

If a header line is to be read, it should be located directly after the first lines to be skipped, before the start of the spherical harmonic coefficents. The header values are returned as a list, where each value is formatted as a string.

A valid line must contain at least 3 words, and the first two words must be integers. When reading the file, all other lines will be considered as "comments" and will be ignored.

# See also

[shread2](pyshread2.html), [shread2error](pyshread2error.html), [shreadjpl](pyshreadjpl.html) [shreadjplerror](pyshreadjplerror.html)
