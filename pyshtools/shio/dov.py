"""
Functions for reading and writing spherical harmonic coefficients from text
files that are formatted as [degree, order, value].
"""
import os as _os
import io as _io
import gzip as _gzip
import zipfile as _zipfile
import numpy as _np
import requests as _requests
import shutil as _shutil


def read_dov(filename, lmax=None, error=False, header=False, header2=False,
             skip=0, encoding=None):
    """
    Read spherical harmonic coefficients from a text file formatted as
    [degree, order, value].

    Usage
    -----
    coeffs, [errors], lmaxout, [header], [header2] = read_dov(
        filename, [error=True, header=True, header2=True, lmax, skip,
        encoding])

    Returns
    -------
    coeffs : ndarray, size(2, lmaxout+1, lmaxout+1)
        The spherical harmonic coefficients.
    errors : ndarray, size(2, lmaxout+1, lmaxout+1)
        The errors associated with the spherical harmonic coefficients.
    lmaxout : int
        The maximum spherical harmonic degree read from the file.
    header : list of type str
        A list of values in the header line found before the start of the
        spherical harmonic coefficients.
    header2 : list of type str
        A list of values in the second header line found before the start of
        the spherical harmonic coefficients.

    Parameters
    ----------
    filename : str
        File name or URL that contains the text-formatted spherical harmonic
        coefficients. filename will be treated as a URL if it starts with
        'http://', 'https://', or 'ftp://'. If filename ends with '.gz' or
        '.zip', the file will be uncompressed before parsing.
    lmax : int, optional, default = None
        The maximum spherical harmonic degree to read from the file. The
        default is to read the entire file.
    error : bool, optional, default = False
        If True, return the errors associated with the spherical harmonic
        coefficients as a separate array.
    header : bool, optional, default = False
        If True, return a list of values in the header line found before the
        start of the spherical harmonic coefficients.
    header2 : bool, optional, default = False
        If True, return a list of values in the second header line found before
        the start of the spherical harmonic coefficients.
    skip : int, optional, default = 0
        The number of lines to skip before parsing the file.
    encoding : str, optional, default = None
        Encoding of the input file. The default is to use the system default.

    Notes
    -----
    This function will read spherical harmonic coefficients from a 'dov'-
    formatted text file. The errors associated with the spherical
    harmonic coefficients, as well as the values in one or two header lines,
    can be read optionally by setting the parameters error, header, and header2
    to True. The optional parameter skip specifies how many lines should be
    skipped before attempting to parse the file, and the optional parameter
    lmax specifies the maximum degree to read from the file. Both real and
    complex spherical harmonic coefficients are supported.

    The spherical harmonic coefficients in the file should be formatted as

    l, m, coeffs[0, l, m]
    l, -m, coeffs[1, l, m]

    where l and m are the spherical harmonic degree and order, respectively.
    If the errors are to be read, the line should be formatted as

    l, m, coeffs[0, l, m], errors[0, l, m]
    l, -m, coeffs[0, l, m], errors[1, l, m]

    For each value of increasing l, all the angular orders are listed in
    pairs with inceasing abs(order), from 0 to l.

    If one or two header lines are to be read, they should be located directly
    after the first lines to be skipped, before the start of the spherical
    harmonic coefficents. The header values are returned as a list, where each
    value is formatted as a string. Comment lines will be ignored, where
    comments start with '#' or the line is all whitespace.

    If filename starts with 'http://', 'https://', or 'ftp://', the file will
    be treated as a URL. In this case, the file will be downloaded in its
    entirety before it is parsed.

    If the filename ends with '.gz' or '.zip', the file will be automatically
    uncompressed before parsing. For zip files, archives with only a single
    file are supported. Note that reading '.gz' and '.zip' files will be
    extremely slow if lmax is not specified.
    """
    if _isurl(filename):
        _response = _requests.get(filename)
        if filename[-4:] == '.zip':
            zf = _zipfile.ZipFile(_io.BytesIO(_response.content))
            if len(zf.namelist()) > 1:
                raise Exception('read_dov can only process zip archives '
                                'that contain a single file. Archive '
                                'contents:\n{}'.format(zf.namelist()))
    elif filename[-4:] == '.zip':
        zf = _zipfile.ZipFile(filename, 'r')
        if len(zf.namelist()) > 1:
            raise Exception('read_dov can only process zip archives that '
                            'contain a single file. Archive contents: \n'
                            '{}'.format(zf.namelist()))

    # If lmax is None, determine lmax by reading last line of the file that
    # is not a comment. (Note that this is very slow for zipped and gzipped
    # files. Consider using indexed_gzip when SEEK_END is supported.)
    if lmax is None:
        if _isurl(filename):
            f = _io.BytesIO(_response.content)
            if filename[-4:] == '.zip':
                f = zf.open(zf.namelist()[0])
        elif filename[-3:] == '.gz':
            f = _gzip.open(filename, mode='rb')
        elif filename[-4:] == '.zip':
            f = zf.open(zf.namelist()[0])
        else:
            f = open(filename, 'rb')

        # determine lmax by reading the last non-comment line of the file
        with f:
            line = ''
            if f.seek(0, _os.SEEK_END) == 0:
                raise RuntimeError('File is empty.')
            else:
                f.seek(-1, _os.SEEK_CUR)

            # read backwards to end of preceding line and then read the line
            while _iscomment(line):
                while f.read(1) != b'\n':
                    try:
                        f.seek(-2, _os.SEEK_CUR)
                    except:
                        f.seek(-1, _os.SEEK_CUR)  # beginning of file
                        break

                if f.tell() <= 1:
                    line = f.readline().decode()
                    line = line.replace(',', ' ')
                    if _iscomment(line):
                        raise RuntimeError('Encountered beginning of file '
                                           'while attempting to determine '
                                           'lmax.')
                    break
                else:
                    line = f.readline().decode()
                    line = line.replace(',', ' ')
                    try:
                        f.seek(-len(line)-2, _os.SEEK_CUR)
                    except:
                        raise RuntimeError('Encountered beginning of file '
                                           'while attempting to determine '
                                           'lmax.')
        lmaxout = int(float(line.split()[0]))

    else:
        lmaxout = lmax

    # open file, skip lines, read header, determine lstart, and then read
    # coefficients one line at a time
    if _isurl(filename):
        if encoding is not None:
            _response.encoding = encoding
        f = _io.StringIO(_response.text)
        if filename[-4:] == '.zip':
            f = _io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding=encoding)
    elif filename[-3:] == '.gz':
        f = _gzip.open(filename, mode='rt', encoding=encoding)
    elif filename[-4:] == '.zip':
        f = _io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding=encoding)
    else:
        f = open(filename, 'r', encoding=encoding)

    with f:
        if skip != 0:
            for i in range(skip):
                line = f.readline()
                if line == '':
                    raise RuntimeError('End of file encountered when '
                                       'skipping lines.')

        # read headers
        if header is True:
            line = f.readline()
            if line == '':
                raise RuntimeError('End of file encountered when '
                                   'reading header line.')
            while _iscomment(line):
                line = f.readline()
                if line == '':
                    raise RuntimeError('End of file encountered when '
                                       'reading header line.')
            line = line.replace(',', ' ')
            header_list = line.split()
        if header2 is True:
            line = f.readline()
            if line == '':
                raise RuntimeError('End of file encountered when '
                                   'reading second header line.')
            while _iscomment(line):
                line = f.readline()
                if line == '':
                    raise RuntimeError('End of file encountered when '
                                       'reading second header line.')
            line = line.replace(',', ' ')
            header2_list = line.split()

        # determine the starting degree
        start_position = f.tell()
        line = f.readline()
        if line == '':
            raise RuntimeError('End of file encountered when determining '
                               'value of lstart.')
        while _iscomment(line):
            line = f.readline()
            if line == '':
                raise RuntimeError('End of file encountered when determining '
                                   'value of lstart.')
        line = line.replace(',', ' ')
        lstart = int(float(line.split()[0]))

        # determine if the coefficients are real or complex
        try:
            num = float(line.split()[2])  # noqa F841
            coeffs = _np.zeros((2, lmaxout+1, lmaxout+1))
            kind = 'real'
            if error is True:
                errors = _np.zeros((2, lmaxout+1, lmaxout+1))
        except ValueError:
            try:
                num = _np.complex128(line.split()[2])  # noqa F841
                coeffs = _np.zeros((2, lmaxout+1, lmaxout+1),
                                   dtype=_np.complex128)
                kind = 'complex'
                if error is True:
                    errors = _np.zeros((2, lmaxout+1, lmaxout+1),
                                       dtype=_np.complex128)
            except ValueError:
                raise ValueError('Coefficients can not be converted to '
                                 'either float or complex. Coefficient '
                                 'is {:s}\n'.format(line.split()[2]) +
                                 'Unformatted string is {:s}'.format(line))

        # rewind one line and read coefficients one line at a time
        f.seek(start_position)

        for degree in range(lstart, lmaxout+1):
            for order in range(degree+1):
                line = f.readline()
                if line == '':
                    raise RuntimeError('End of file encountered at '
                                       'degree and order {:d}, {:d}.'
                                       .format(degree, order))
                while _iscomment(line):
                    line = f.readline()
                    if line == '':
                        raise RuntimeError('End of file encountered at '
                                           'degree and order {:d}, {:d}.'
                                           .format(degree, order))
                line = line.replace(',', ' ')
                l = int(float(line.split()[0]))
                m = int(float(line.split()[1]))
                if degree != l or order != m:
                    raise RuntimeError('Degree and order from file do not '
                                       'correspond to expected values.\n '
                                       'Read {:d}, {:d}. Expected {:d}, '
                                       '{:d}.'.format(l, m, degree, order))

                if order > 0:
                    line2 = f.readline()
                    if line2 == '':
                        raise RuntimeError('End of file encountered at '
                                           'degree and order {:d}, {:d}.'
                                           .format(degree, order))
                    while _iscomment(line):
                        line = f.readline()
                        if line == '':
                            raise RuntimeError('End of file encountered at '
                                               'degree and order {:d}, {:d}.'
                                               .format(degree, order))
                    line2 = line2.replace(',', ' ')
                    l2 = int(float(line2.split()[0]))
                    m2 = int(float(line2.split()[1]))
                    if degree != l2 or -order != m2:
                        raise RuntimeError('Degree and order from file do not '
                                           'correspond to expected values.\n '
                                           'Read {:d}, {:d}. Expected {:d}, '
                                           '{:d}.'.format(l2, m2,
                                                          degree, -order))
                    m2 = abs(m2)

                if kind == 'real':
                    coeffs[0, l, m] = _np.float64(line.split()[2])
                    if order > 0:
                        coeffs[1, l, m] = _np.float64(line2.split()[2])
                else:
                    coeffs[0, l, m] = _np.complex128(line.split()[2])
                    if order > 0:
                        coeffs[1, l, m] = _np.complex128(line2.split()[2])

                if error:
                    if len(line.split()) < 4:
                        raise RuntimeError('When reading errors, '
                                           'each line must '
                                           'contain at least 4 elements. '
                                           'Last line is: {:s}'.format(line))

                    if kind == 'real':
                        errors[0, l, m] = _np.float64(line.split()[3])
                        if order > 0:
                            errors[1, l, m] = _np.float64(line2.split()[3])
                    else:
                        errors[0, l, m] = _np.complex128(line.split()[3])
                        if order > 0:
                            errors[1, l, m] = _np.complex128(line2.split()[3])

    if error is True and header is True:
        if header2:
            return coeffs, errors, lmaxout, header_list, header2_list
        else:
            return coeffs, errors, lmaxout, header_list
    elif error is True and header is False:
        return coeffs, errors, lmaxout
    elif error is False and header is True:
        if header2:
            return coeffs, lmaxout, header_list, header2_list
        else:
            return coeffs, lmaxout, header_list
    else:
        return coeffs, lmaxout


def write_dov(filename, coeffs, errors=None, header=None, header2=None,
              lmax=None, encoding=None):
    """
    Write spherical harmonic coefficients to a text file formatted as
    [degree, order, value].

    Usage
    -----
    write_dov(filename, coeffs, [errors, header, header2, lmax, encoding])

    Parameters
    ----------
    filename : str
        File name of the 'dov'-formatted spherical harmonic coefficients. If
        filename ends with '.gz' the file will be automatically compressed with
        gzip.
    coeffs : ndarray, size(2, lmaxin+1, lmaxin+1)
        The spherical harmonic coefficients.
    errors : ndarray, size(2, lmaxin+1, lmaxin+1), optional, default = None
        The errors associated with the spherical harmonic coefficients.
    header : str, optional default = None
        A string to be written directly before the spherical harmonic
        coefficients.
    header2 : str, optional default = None
        A second string to be written directly before the spherical harmonic
        coefficients.
    lmax : int, optional, default = None
        The maximum spherical harmonic degree to write to the file.
    encoding : str, optional, default = None
        Encoding of the output file. The default is to use the system default.

    Notes
    -----
    This function will write spherical harmonic coefficients (and optionally
    the errors) to a text file formatted as [degree, order, value]. If header
    or header2 are specified, these strings will be written first, directly
    before the spherical harmonic coefficients. Both real and complex spherical
    harmonic coefficients are supported.

    The spherical harmonic coefficients in the file will be formatted as pairs
    of lines as

    l, m, coeffs[0, l, m]
    l, -m, coeffs[1, l, m]

    where l and m are the spherical harmonic degree and order, respectively.
    If the errors are included, each pair of lines will be formatted as

    l, m, coeffs[0, l, m], errors[0, l, m]
    l, -m, coeffs[1, l, m], errors[1, l, m]

    For each value of increasing l, all the angular orders are listed in
    inceasing order, from 0 to l.

    If the filename ends with '.gz', the file will be automatically compressed
    using gzip.
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

    with open(filebase, mode='w', encoding=encoding) as file:
        if header is not None:
            file.write(header + '\n')
        if header2 is not None:
            file.write(header2 + '\n')
        for l in range(lmax+1):
            for m in range(l+1):
                if errors is not None:
                    if m == 0:
                        file.write('{:d}, {:d}, {:.16e}, {:.16e}\n'
                                   .format(l, m, coeffs[0, l, m],
                                           errors[0, l, m]))
                    else:
                        file.write('{:d}, {:d}, {:.16e}, {:.16e}\n'
                                   .format(l, m, coeffs[0, l, m],
                                           errors[0, l, m]))
                        file.write('{:d}, {:d}, {:.16e}, {:.16e}\n'
                                   .format(l, -m, coeffs[1, l, m],
                                           errors[1, l, m]))
                else:
                    if m == 0:
                        file.write('{:d}, {:d}, {:.16e}\n'
                                   .format(l, m, coeffs[0, l, m]))
                    else:
                        file.write('{:d}, {:d}, {:.16e}\n'
                                   .format(l, m, coeffs[0, l, m]))
                        file.write('{:d}, {:d}, {:.16e}\n'
                                   .format(l, -m, coeffs[1, l, m]))

    if filename[-3:] == '.gz':
        with open(filebase, 'rb') as f_in:
            with _gzip.open(filename, 'wb') as f_out:
                _shutil.copyfileobj(f_in, f_out)


def _iscomment(line):
    """
    Comment lines begin with the character '#', '', or all whitespace.
    """
    if line == '':
        return True
    if line.isspace():
        return True
    elif line.strip()[0] == '#':
        return True
    else:
        return False


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
