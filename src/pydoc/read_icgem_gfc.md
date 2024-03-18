# read_icgem_gfc()

Read real spherical harmonic gravity coefficients from an ICGEM formatted
file.

# Usage

cilm, gm, r0, [errors] = read_icgem_gfc(filename,
    [errors, lmax, epoch, encoding,
    quiet)

# Returns

cilm : ndarray, size (2, lmax + 1, lmax + 1)
:   Array of '4pi' normalized spherical harmonic coefficients for the
    given epoch.

gm : float
:   Gravitational constant of the model, in m\*\*3/s\*\*2.

r0 : float
:   Reference radius of the model, in meters.

errors : ndarray, optional, shape (2, lmax + 1, lmax + 1)
:   Array of the spherical harmonic error coefficients for the given epoch.

# Parameters

filename : str or pathlib.Path
:   The filename containing the spherical harmonic ICGEM-formatted
    coefficients. filename will be treated as a URL if it starts with
    'http://', 'https://', or 'ftp://'. If filename ends with '.gz' or
    '.zip' (or if the path contains '/zip/'), the file will be
    uncompressed before parsing.

errors : str, optional, default = None
:   Which errors to read. Can be 'unknown', 'calibrated', 'formal' or None.

lmax : int, optional, default = None
:   Maximum degree to read from the file. If lmax is None, less than 0, or
    greater than lmax_model, the maximum degree of the model will be used.

epoch : str or float, optional, default = None
:   The epoch time to calculate time-variable coefficients in YYYYMMDD.DD
    format. If None then the reference epoch t0 of the model will be used.
    If the format of the file is 'icgem2.0' then the epoch must be
    specified.

encoding : str, optional, default = None
:   Encoding of the input file. The default is to use the system default.
    If the default encoding doesn't work, try 'iso-8859-1'.

quiet : bool, default = False
:   If True, suppress warnings about undefined keywords when reading the
    file.

# Notes

This routine reads ICGEM formatted files of gravitational potential models
and outputs arrays of the gravitational potential coefficients, errors, GM,
and the reference radius. If epoch is specified, the coefficients will make
use of the time variable terms in order to compute and return the potential
coefficients for the specified epoch. Otherwise, the coefficients will be
returned for the reference epoch of the model.

Valid keys in the header section include:
        modelname (not used)
    product_type (only 'gravity_field' is allowed)
    earth_gravity_constant or gravity_constant
    radius
    max_degree
    errors ('unknown', 'formal', 'calibrated' or 'calibrated_and_formal')
    tide_system (not used)
    norm (not used)
    format (either None or 'icgem2.0')

Valid keys in the data section include:
        gfc
    gfct
    trnd or dot
    asin
    acos

Data lines starting with an unknown key are ignored.

