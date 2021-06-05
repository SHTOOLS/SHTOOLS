"""
Support for reading real magnetic-potential spherical harmonic coefficients
from IGRM formatted files.
"""
import io
import gzip
import zipfile
import numpy as _np
import requests as _requests


def read_igrf(filename, year=2020., encoding=None):
    """
    Read IGRF real spherical harmonic coefficients, and return the magnetic
    potential coefficients for the specified year.

    Usage
    -----
    read_igrm(filename, [year])

    Returns
    -------
    clm : ndarray, size (2, 14, 14)
        Array of Schmidt semi-normalized coefficients.

    Parameters
    ----------
    filename : str
        The filename containing the IGRF formatted spherical harmonic
        coefficients. filename will be treated as a URL if it starts with
        'http://', 'https://', or 'ftp://'. If filename ends with '.gz' or
        '.zip', the file will be uncompressed before parsing.
    year : float, optional, default = 2020.
        The year to compute the coefficients.
    encoding : str, optional, default = None
        Encoding of the input file. The default is to use the system default.

    Notes
    -----
    The current International Geomagnetic Reference Field (IGRF-13) is a
    degree 13 time variable model that is valid between 1900 and 2020.
    Coefficients are provided in 5 year intervals, and for a given year, the
    values of the coefficients are interpolated linearly between adjacent
    entries. For years between 2020 and 2025, the coefficients are extrapolated
    using the provided secular variation. The reference radius is 6371.2 km.

    This routine can read the models IGRF-11, 12, and 13. Prior models have a
    different format.
    """
    lmax = 13
    coeffs = _np.zeros([2, lmax+1, lmax+1])

    # open filename as a text file
    if _isurl(filename):
        _response = _requests.get(filename)
        if filename[-4:] == '.zip':
            zf = zipfile.ZipFile(io.BytesIO(_response.content))
            if len(zf.namelist()) > 1:
                raise Exception('read_igrf can only process zip archives '
                                'that contain a single file. Archive '
                                'contents:\n{}'.format(zf.namelist()))
            f = io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding=encoding)
        else:
            if encoding is not None:
                _response.encoding = encoding
            f = io.StringIO(_response.text)
    elif filename[-3:] == '.gz':
        f = gzip.open(filename, mode='rt', encoding=encoding)
    elif filename[-4:] == '.zip':
        zf = zipfile.ZipFile(filename, 'r')
        if len(zf.namelist()) > 1:
            raise Exception('read_igrf can only process zip archives '
                            'that contain a single file. Archive contents: \n'
                            '{}'.format(zf.namelist()))
        f = io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding=encoding)
    else:
        f = open(filename, 'r', encoding=encoding)

    with f:
        lines = f.readlines()
        year_list = lines[3].split()[3:]
        year_list[-1] = int(year_list[-1].split(sep='-')[1]) + 2000
        year_list = [float(i) for i in year_list]
        year_start = year_list[0]
        year_end = year_list[-1]

        if year < year_start or year > year_end:
            raise ValueError('Year must be between {:d} and {:d}. Input value '
                             'is {:f}'.format(int(year_start), int(year_end),
                                              year))

        for line in lines[4:]:
            line = line.split()
            data = [float(i) for i in line[3:]]

            if line[0] == 'g':
                index = 0
            elif line[0] == 'h':
                index = 1
            else:
                raise ValueError("Error reading file. Expected 'g' or "
                                 "'h' but read {:s}".format(repr(line[0])))

            degree = int(line[1])
            order = int(line[2])

            if year >= year_list[-2]:
                # use SV in file for secular variation
                value = data[-2] + data[-1] * (year - year_list[-2])
                coeffs[index, degree, order] = value

            else:
                # perform linear interpolation between adjacent years
                for i in range(len(year_list) - 2):
                    if year >= year_list[i] and year < year_list[i+1]:
                        value = data[i] + (data[i+1] - data[i]) / \
                            (year_list[i+1] - year_list[i]) * \
                            (year - year_list[i])
                        coeffs[index, degree, order] = value

    return coeffs


def _isurl(filename):
    """
    Determine if filename is a URL. Valid URLs start with
        'http://'
        'https://'
        'ftp://'
    """
    if filename[0:7].lower() == 'http://':
        return True
    elif filename[0:8].lower() == 'https://':
        return True
    elif filename[0:6].lower() == 'ftp://':
        return True
    else:
        return False
