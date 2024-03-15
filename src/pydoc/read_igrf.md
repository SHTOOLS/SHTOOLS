# read_igrf()

Read IGRF real spherical harmonic coefficients, and return the magnetic
potential coefficients for the specified year.

# Usage

read_igrm(filename, [year])

# Returns

clm : ndarray, size (2, 14, 14)
:   Array of Schmidt semi-normalized coefficients.

# Parameters

filename : str or pathlib.Path
:   The filename containing the IGRF formatted spherical harmonic
    coefficients. filename will be treated as a URL if it starts with
    'http://', 'https://', or 'ftp://'. If filename ends with '.gz' or
    '.zip', the file will be uncompressed before parsing.

year : float, optional, default = 2020.
:   The year to compute the coefficients.

encoding : str, optional, default = None
:   Encoding of the input file. The default is to use the system default.

# Notes

The current International Geomagnetic Reference Field (IGRF-13) is a
degree 13 time variable model that is valid between 1900 and 2020.
Coefficients are provided in 5 year intervals, and for a given year, the
values of the coefficients are interpolated linearly between adjacent
entries. For years between 2020 and 2025, the coefficients are extrapolated
using the provided secular variation. The reference radius is 6371.2 km.

This routine can read the models IGRF-11, 12, and 13. Prior models have a
different format.

