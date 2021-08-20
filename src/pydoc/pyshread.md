# shread()

Read shtools-formatted spherical harmonic coefficients from a text file.

# Usage

```python
coeffs, [errors], lmaxout, [header], [header2] = shread(
    filename, [error=True, header=True, header2=True, lmax, skip,
```
        encoding])

# Returns

**coeffs : ndarray, size(2, lmaxout+1, lmaxout+1)**
:   The spherical harmonic coefficients.

**errors : ndarray, size(2, lmaxout+1, lmaxout+1)**
:   The errors associated with the spherical harmonic coefficients.

**lmaxout : int**
:   The maximum spherical harmonic degree read from the file.

**header : list of type str**
:   A list of values in the header line found before the start of the
        spherical harmonic coefficients.

**header2 : list of type str**
:   A list of values in the second header line found before the start of
        the spherical harmonic coefficients.

# Parameters

**filename : str**
:   File name or URL that contains the text-formatted spherical harmonic
        coefficients. filename will be treated as a URL if it starts with
        'http://', 'https://', or 'ftp://'. If filename ends with '.gz' or
        '.zip', the file will be uncompressed before parsing.

**lmax : int, optional, default = None**
:   The maximum spherical harmonic degree to read from the file. The
        default is to read the entire file.

**error : bool, optional, default = False**
:   If True, return the errors associated with the spherical harmonic
        coefficients as a separate array.

**header : bool, optional, default = False**
:   If True, return a list of values in the header line found before the
        start of the spherical harmonic coefficients.

**header2 : bool, optional, default = False**
:   If True, return a list of values in the second header line found before
        the start of the spherical harmonic coefficients.

**skip : int, optional, default = 0**
:   The number of lines to skip before parsing the file.

**encoding : str, optional, default = None**
:   Encoding of the input file. The default is to use the system default.

# Notes

This function will read spherical harmonic coefficients from an
shtools-formatted text file. The errors associated with the spherical
harmonic coefficients, as well as the values in one or two header lines,
can be read optionally by setting the parameters error, header, and header2
to True. The optional parameter skip specifies how many lines should be
skipped before attempting to parse the file, and the optional parameter
lmax specifies the maximum degree to read from the file. Both real and
complex spherical harmonic coefficients are supported.

The spherical harmonic coefficients in the file should be formatted as

l, m, coeffs[0, l, m], coeffs[1, l, m]

where l and m are the spherical harmonic degree and order, respectively.
The terms coeffs[1, l, 0] can be neglected as they are zero. If the errors
are to be read, the line should be formatted as

l, m, coeffs[0, l, m], coeffs[1, l, m], errors[0, l, m], errors[1, l, m]

For each value of increasing l, all the angular orders are listed in
inceasing order, from 0 to l.

If header lines are to be read, they should be located directly after the
first lines to be skipped, before the start of the spherical harmonic
coefficents. The header values are returned as a list, where each value is
formatted as a string. Comment lines will be ignored, where comments start
with '#' or the line is all whitespace.

If filename starts with 'http://', 'https://', or 'ftp://', the file will
be treated as a URL. In this case, the file will be downloaded in its
entirety before it is parsed.

If the filename ends with '.gz' or '.zip', the file will be automatically
uncompressed before parsing. For zip files, archives with only a single
file are supported. Note that reading '.gz' and '.zip' files will be
extremely slow if lmax is not specified.
    