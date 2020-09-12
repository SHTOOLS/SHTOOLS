"""
Functions for reading and writing real spherical harmonic coefficients in the
binary 'bshc' format used by Curtin University.
"""
import io as _io
import struct as _struct
import gzip as _gzip
import zipfile as _zipfile
import numpy as _np
import requests as _requests
import shutil as _shutil


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
    (with all orders being listed, one degree at a time). For a 100 degree
    file, the contents are

    0 10800
    C(0,0), C(1,0), C(1,1), C(2,0), C(2,1), ... C(100,99), C(100,100)
    S(0,0), S(1,0), S(1,1), S(2,0), S(2,1), ... S(100,99), S(100,100).

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
            zf = _zipfile.ZipFile(_io.BytesIO(_response.content))
            if len(zf.namelist()) > 1:
                raise Exception('read_bshc can only process zip archives '
                                'that contain a single file. Archive '
                                'contents:\n{}'.format(zf.namelist()))
            f = zf.open(zf.namelist()[0])
        f = _io.BytesIO(_response.content)
    elif filename[-3:] == '.gz':
        f = _gzip.open(filename, mode='rb')
    elif filename[-4:] == '.zip':
        zf = _zipfile.ZipFile(filename, 'r')
        if len(zf.namelist()) > 1:
            raise Exception('read_bshc can only process zip archives that '
                            'contain a single file. Archive contents: \n'
                            '{}'.format(zf.namelist()))
        f = zf.open(zf.namelist()[0])
    else:
        f = open(filename, 'rb')

    with f:
        degree_start = int(_struct.unpack('<d', f.read(8))[0])
        degree_end = int(_struct.unpack('<d', f.read(8))[0])
        if lmax is None:
            lmax = degree_end

        if lmax > degree_end:
            raise RuntimeError('lmax must be less than or equal to '
                               '{:d}. Input degree is {:d}.'
                               .format(degree_end, lmax))

        coeffs = _np.zeros((2, lmax+1, lmax+1))

        for degree in range(degree_start, lmax+1):
            n = degree + 1
            coeffs[0, degree, 0:n] = _struct.unpack('<{}d'.format(n),
                                                    f.read(8*n))
        for degree in range(lmax+1, degree_end + 1):
            n = degree + 1
            f.read(8*n)

        for degree in range(degree_start, lmax + 1):
            n = degree + 1
            coeffs[1, degree, 0:n] = _struct.unpack('<{}d'.format(n),
                                                    f.read(8*n))

        return coeffs, lmax


def write_bshc(filename, coeffs, lmax=None):
    """
    Write real spherical harmonic coefficients to a binary bshc file.

    Usage
    -----
    write_bshc(filename, coeffs, [lmax])

    Parameters
    ----------
    filename : str
        File name of the binary 'bshc'-formatted spherical harmonic
        coefficients. If filename ends with '.gz' the file will be
        automatically compressed with gzip.
    coeffs : ndarray, size(2, lmaxin+1, lmaxin+1)
        The spherical harmonic coefficients.
    lmax : int, optional, default = None
        The maximum spherical harmonic degree to write to the file. The
        default is to write all coefficients.

    Notes
    -----
    This function writes real spherical harmonic coefficients to a binary
    'bshc'-formatted file as used at Curtin University. The file is composed
    solely of 8-byte floats, starting with the minimum and maximum degree,
    and followed by the cosine coefficients and then sine coefficients
    (with all orders being listed, one degree at a time). For a 100 degree
    file, the contents are

    0 100
    C(0,0), C(1,0), C(1,1), C(2,0), C(2,1), ... C(100,99), C(100,100)
    S(0,0), S(1,0), S(1,1), S(2,0), S(2,1), ... S(100,99), S(100,100).

    If the filename ends with '.gz', the file will be automatically
    compressed using gzip.
    """
    if lmax is None:
        lmax = coeffs.shape[1] - 1
    else:
        if lmax > coeffs.shape[1] - 1:
            raise ValueError('lmax is greater than the input coefficients. '
                             'lmax = {:d}, lmax of input coefficients = {:d}.'
                             .format(lmax, coeffs.shape[1] - 1))

    if filename[-3:] == '.gz':
        filebase = filename[:-3]
    else:
        filebase = filename

    with open(filebase, 'wb') as f:
        f.write(_struct.pack('<d', 0.))
        f.write(_struct.pack('<d', lmax))

        for degree in range(0, lmax+1):
            n = degree + 1
            f.write(_struct.pack('<{}d'.format(n), *coeffs[0, degree, 0:n]))

        for degree in range(0, lmax+1):
            n = degree + 1
            f.write(_struct.pack('<{}d'.format(n), *coeffs[1, degree, 0:n]))

    if filename[-3:] == '.gz':
        with open(filebase, 'rb') as f_in:
            with _gzip.open(filename, 'wb') as f_out:
                _shutil.copyfileobj(f_in, f_out)


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
