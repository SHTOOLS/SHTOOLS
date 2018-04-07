"""
Functions for reading spherical harmonic coefficients from files.
"""
import os

import numpy as _np


# ==== shread() ====

def shread(filename, lmax=None, error=False, header=False, skip=0):
    """
    Read spherical harmonic coefficients from a text file.

    Usage
    -----
    coeffs, lmaxout = shread(filename, [lmax, skip])
    coeffs, lmaxout, header = shread(filename, header=True, [lmax, skip])
    coeffs, errors, lmaxout = shread(filename, error=True, [lmax, skip])
    coeffs, errors, lmaxout, header = shread(filename, error=True,
                                             header=True, [lmax, skip])

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

    Parameters
    ----------
    filename : str
        Filename containing the text-formatted spherical harmonic coefficients.
    lmax : int, optional, default = None
        The maximum spherical harmonic degree to read from the file. The
        default is to read the entire file.
    error : bool, optional, default = False
        If True, return the errors associated with the spherical harmonic
        coefficients as a separate array.
    header : bool, optional, default = False
        If True, return a list of values in the header line found before the
        start of the spherical harmonic coefficients.
    skip : int, optional, default = 0
        The number of lines to skip before parsing the file.

    Description
    -----------
    This function will read spherical harmonic coefficients from an
    ascii-formatted text file. The errors associated with the spherical
    harmonic coefficients, as well as the values in a single header line, can
    optionally be read by setting the optional parameters error and header to
    True. The optional parameter skip specifies how many lines should be
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

    If a header line is to be read, it should be located directly after the
    first lines to be skipped, before the start of the spherical harmonic
    coefficents. The header values are returned as a list, where each value is
    formatted as a string.

    A valid line must contain at least 3 words, and the first two words must be
    integers. When reading the file, all other lines will be considered as
    "comments" and will be ignored.
    """

    # determine lmax by reading the last non-comment line of the file
    with open(filename, 'rb') as f:
        line = ''
        if f.seek(0, os.SEEK_END) == 0:
            raise RuntimeError('File is empty.')
        else:
            f.seek(-1, os.SEEK_CUR)

        # read backwards to end of preceding line and then read the line
        while _iscomment(line):
            while f.read(1) != b'\n':
                try:
                    f.seek(-2, os.SEEK_CUR)
                except:
                    f.seek(-1, os.SEEK_CUR)  # beginning of file
                    break

            if f.tell() <= 1:
                line = f.readline().decode()
                line = line.replace(',', ' ')
                if _iscomment(line):
                    raise RuntimeError('Encountered beginning of file ' +
                                       'while attempting to determine lmax.')
                break
            else:
                line = f.readline().decode()
                line = line.replace(',', ' ')
                try:
                    f.seek(-len(line)-2, os.SEEK_CUR)
                except:
                    raise RuntimeError('Encountered beginning of file ' +
                                       'while attempting to determine lmax.')

    lmaxfile = int(line.split()[0])
    if lmax is not None:
        lmaxout = min(lmax, lmaxfile)
    else:
        lmaxout = lmaxfile

    # determine if coefficients are real or complex
    try:
        num = float(line.split()[2])
        coeffs = _np.zeros((2, lmaxout+1, lmaxout+1))
        kind = 'real'
        if error is True:
            errors = _np.zeros((2, lmaxout+1, lmaxout+1))
    except ValueError:
        try:
            num = complex(line.split()[2])
            coeffs = _np.zeros((2, lmaxout+1, lmaxout+1), dtype=complex)
            kind = 'complex'
            if error is True:
                errors = _np.zeros((2, lmaxout+1, lmaxout+1), dtype=complex)
        except ValueError:
            raise ValueError('Coefficients can not be converted to ' +
                             'either float or complex. Coefficient ' +
                             'is {:s}\n'.format(line.split()[2]) +
                             'Unformatted string is {:s}'.format(line))

    # determine lstart and read header
    with open(filename, 'r') as f:
        if skip != 0:
            for i in range(skip):
                line = f.readline()
                if line == '':
                    raise RuntimeError('End of file encountered when ' +
                                       'skipping lines.')

        if header is True:
            line = f.readline()
            if line == '':
                    raise RuntimeError('End of file encountered when ' +
                                       'reading header line.')
            line = line.replace(',', ' ')
            header_list = line.split()

        line = f.readline()
        if line == '':
            raise RuntimeError('End of file encountered when determining ' +
                               'value of lstart.')
        line = line.replace(',', ' ')
        while _iscomment(line):
            line = f.readline()
            if line == '':
                raise RuntimeError('End of file encountered when ' +
                                   'determining value of lstart.')
            line = line.replace(',', ' ')
        lstart = int(line.split()[0])

    # read coefficients one line at a time
    with open(filename, 'r') as f:
        if skip != 0:
            for i in range(skip):
                f.readline()
        if header is True:
            f.readline()

        for degree in range(lstart, lmaxout+1):
            for order in range(degree+1):
                line = f.readline()
                if line == '':
                    raise RuntimeError('End of file encountered at ' +
                                       'degree and order {:d}, {:d}.'
                                       .format(degree, order))
                line = line.replace(',', ' ')
                while _iscomment(line):
                    line = f.readline()
                    if line == '':
                        raise RuntimeError('End of file encountered at ' +
                                           'degree and order {:d}, {:d}.'
                                           .format(degree, order))
                    line = line.replace(',', ' ')

                l = int(line.split()[0])
                m = int(line.split()[1])

                if degree != l or order != m:
                    raise RuntimeError('Degree and order from file do not ' +
                                       'correspond to expected values.\n ' +
                                       'Read {:d}, {:d}. Expected {:d}, {:d}.'
                                       .format(degree, order, l, m))

                if kind == 'real':
                    coeffs[0, l, m] = float(line.split()[2])
                    if m > 0:
                        coeffs[1, l, m] = float(line.split()[3])
                else:
                    coeffs[0, l, m] = complex(line.split()[2])
                    if m > 0:
                        coeffs[1, l, m] = complex(line.split()[3])

                if error:
                    if len(line.split()) < 6:
                        raise RuntimeError('When reading errors, ' +
                                           'each line must ' +
                                           'contain at least 6 elements. ' +
                                           'Last line is: {:s}'.format(line))

                    if kind == 'real':
                        errors[0, l, m] = float(line.split()[4])
                        errors[1, l, m] = float(line.split()[5])
                    else:
                        errors[0, l, m] = complex(line.split()[4])
                        errors[1, l, m] = complex(line.split()[5])

    if error is True and header is True:
        return coeffs, errors, lmaxout, header_list
    elif error is True and header is False:
        return coeffs, errors, lmaxout
    elif error is False and header is True:
        return coeffs, lmaxout, header_list
    else:
        return coeffs, lmaxout


def _iscomment(line):
    """
    Determine if a line is a comment line. A valid line contains at least three
    words, with the first two being integers. Note that Python 2 and 3 deal
    with strings differently.
    """
    if line.isspace():
        return True
    elif len(line.split()) >= 3:
        try:  # python 3 str
            if line.split()[0].isdecimal() and line.split()[1].isdecimal():
                return False
        except:  # python 2 str
            if (line.decode().split()[0].isdecimal() and
                    line.split()[1].decode().isdecimal()):
                return False
        return True
    else:
        return True
