"""
Functions for reading real spherical harmonic coefficients from binary 'bshc'
files provided by Curtin University.
"""
import io
import struct
import gzip
import zipfile
import numpy as _np
import requests as _requests


def read_bshc(filename, lmax=None):
    """
    Read real spherical harmonic coefficients from a binary bshc file.

    Usage
    -----
    coeffs, lmaxout = read_bshc(filename, [lmax])

    Returns
    -------
    coeffs : ndarray, size(2, lmaxout+1, lmaxout+1)
        The spherical harmonic coefficients.
    lmaxout : int
        The maximum spherical harmonic degree read from the file.

    Parameters
    ----------
    filename : str
        File name or URL that contains the spherical harmonic coefficients.
        filename will be treated as a URL if it starts with 'http://',
        'https://', or 'ftp://'. If filename ends with '.gz' or '.zip', the
        file will be uncompressed before parsing.
    lmax : int, optional, default = None
        The maximum spherical harmonic degree to read from the file. The
        default is to read the entire file.

    Notes
    -----
    This function reads real spherical harmonic coefficients from binary
    'bshc'-formatted files as used at Curtin University. The file is composed
    solely of 8-byte floats, starting with the minimum and maximum degree,
    and followed by the cosine coefficients and then sine coefficients
    (with all orders being listed, one degree at a time). For a 10800 degree
    file, the contents are

    0 10800
    C(0,0), C(1,0), C(1,1), C(2,0), C(2,1), ... C(10800,10799), C(10800,10800)
    S(0,0), S(1,0), S(1,1), S(2,0), S(2,1), ... S(10800,10799), S(10800,10800).

    If filename starts with 'http://', 'https://', or 'ftp://', the file will
    be treated as a URL. In this case, the file will be downloaded in its
    entirety before it is parsed.

    If the filename ends with '.gz' or '.zip', the file will be automatically
    uncompressed before parsing. For zip files, archives with only a single
    file are supported.
    """
    if _isurl(filename):
        _response = _requests.get(filename, stream=True)
        if filename[-4:] == '.zip':
            zf = zipfile.ZipFile(io.BytesIO(_response.content))
            if len(zf.namelist()) > 1:
                raise Exception('read_bshc can only process zip archives '
                                'that contain a single file. Archive '
                                'contents:\n{}'.format(zf.namelist()))
            f = zf.open(zf.namelist()[0])
        f = io.BytesIO(_response.content)
    elif filename[-3:] == '.gz':
        f = gzip.open(filename, mode='rb')
    elif filename[-4:] == '.zip':
        zf = zipfile.ZipFile(filename, 'r')
        if len(zf.namelist()) > 1:
            raise Exception('read_bshc can only process zip archives that '
                            'contain a single file. Archive contents: \n'
                            '{}'.format(zf.namelist()))
        f = zf.open(zf.namelist()[0])
    else:
        f = open(filename, 'rb')

    with f:
        degree_start = int(struct.unpack('<d', f.read(8))[0])
        degree_end = int(struct.unpack('<d', f.read(8))[0])
        if lmax is None:
            lmax = degree_end

        if lmax > degree_end:
            raise RuntimeError('lmax must be less than or equal to '
                               '{:d}. Input degree is {:d}.'
                               .format(degree_end, lmax))

        coeffs = _np.zeros((2, lmax+1, lmax+1))

        for degree in range(degree_start, lmax+1):
            n = degree + 1
            coeffs[0, degree, 0:n] = struct.unpack('<{}d'.format(n),
                                                   f.read(8*n))
        for degree in range(lmax+1, degree_end + 1):
            n = degree + 1
            f.read(8*n)

        for degree in range(degree_start, lmax + 1):
            n = degree + 1
            coeffs[1, degree, 0:n] = struct.unpack('<{}d'.format(n),
                                                   f.read(8*n))

        return coeffs, lmax


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
