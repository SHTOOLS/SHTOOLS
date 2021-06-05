"""
Support for reading and writing real gravitational-potential coefficients in
ICGEM-GFC formatted files.
"""
import io as _io
import gzip as _gzip
import zipfile as _zipfile
import numpy as _np
import shutil as _shutil
import requests as _requests
from pyshtools.utils.datetime import _yyyymmdd_to_year_fraction


def read_icgem_gfc(filename, errors=None, lmax=None, epoch=None,
                   encoding=None, quiet=False):
    """
    Read real spherical harmonic gravity coefficients from an ICGEM formatted
    file.

    Usage
    -----
    cilm, gm, r0, [errors] = read_icgem_gfc(filename,
                                            [errors, lmax, epoch, encoding,
                                             quiet)

    Returns
    -------
    cilm : ndarray, size (2, lmax + 1, lmax + 1)
        Array of '4pi' normalized spherical harmonic coefficients for the
        given epoch.
    gm : float
        Gravitational constant of the model, in m**3/s**2.
    r0 : float
        Reference radius of the model, in meters.
    errors : ndarray, optional, shape (2, lmax + 1, lmax + 1)
        Array of the spherical harmonic error coefficients for the given epoch.

    Parameters
    ----------
    filename : str
        The filename containing the spherical harmonic ICGEM-formatted
        coefficients. filename will be treated as a URL if it starts with
        'http://', 'https://', or 'ftp://'. If filename ends with '.gz' or
        '.zip' (or if the path contains '/zip/'), the file will be
        uncompressed before parsing.
    errors : str, optional, default = None
        Which errors to read. Can be 'unknown', 'calibrated', 'formal' or None.
    lmax : int, optional, default = None
        Maximum degree to read from the file. If lmax is None, less than 0, or
        greater than lmax_model, the maximum degree of the model will be used.
    epoch : str or float, optional, default = None
        The epoch time to calculate time-variable coefficients in YYYYMMDD.DD
        format. If None then the reference epoch t0 of the model will be used.
        If the format of the file is 'icgem2.0' then the epoch must be
        specified.
    encoding : str, optional, default = None
        Encoding of the input file. The default is to use the system default.
        If the default encoding doesn't work, try 'iso-8859-1'.
    quiet : bool, default = False
        If True, suppress warnings about undefined keywords when reading the
        file.

    Notes
    -----
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
    """
    header = {}
    header_keys = ['modelname', 'product_type', 'earth_gravity_constant',
                   'gravity_constant', 'radius', 'max_degree', 'errors',
                   'tide_system', 'norm', 'format']
    valid_err = ('unknown', 'calibrated', 'formal', 'calibrated_and_formal')

    # open filename as a text file
    if _isurl(filename):
        _response = _requests.get(filename)
        if _iszipfile(filename):
            zf = _zipfile.ZipFile(_io.BytesIO(_response.content))
            if len(zf.namelist()) > 1:
                raise Exception('read_icgem_gfc can only process zip archives '
                                'that contain a single file. Archive '
                                'contents:\n{}'.format(zf.namelist()))
            f = _io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding=encoding)
        else:
            if encoding is not None:
                _response.encoding = encoding
            f = _io.StringIO(_response.text)
    elif filename[-3:] == '.gz':
        f = _gzip.open(filename, mode='rt', encoding=encoding)
    elif filename[-4:] == '.zip':
        zf = _zipfile.ZipFile(filename, 'r')
        if len(zf.namelist()) > 1:
            raise Exception('read_icgem_gfc can only process zip archives '
                            'that contain a single file. Archive contents: \n'
                            '{}'.format(zf.namelist()))
        f = _io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding=encoding)
    else:
        f = open(filename, 'r', encoding=encoding)

    with f:
        for line in f:
            if 'end_of_head' in line:
                break
            for key in header_keys:
                if key in line:
                    header[key] = line.strip().split()[1]

        if header['product_type'] != 'gravity_field':
            raise ValueError(
                'This routine reads only gravity_field data products.')

        is_v2 = False
        if 'format' in header and header['format'] == 'icgem2.0':
            is_v2 = True

        if epoch is None and is_v2:
            raise ValueError(
                'epoch must be specified for the "icgem2.0" format.')
        elif epoch is not None:
            epoch = _yyyymmdd_to_year_fraction(epoch)

        if 'earth_gravity_constant' in header:
            gravity_constant = float(header[
                'earth_gravity_constant'].lower().replace('d', 'e'))
        elif 'gravity_constant' in header:
            gravity_constant = float(header[
                'gravity_constant'].lower().replace('d', 'e'))
        else:
            raise ValueError(
                'No standard gravitational constant in the header.')

        radius = float(header['radius'].lower().replace('d', 'e'))

        lmax_model = int(header['max_degree'])
        if lmax is None or lmax < 0 or lmax > lmax_model:
            lmax = lmax_model

        if errors is not None:
            if header['errors'] == 'no':
                raise ValueError('This model has no errors.')
            elif errors not in valid_err[:-1]:
                raise ValueError(
                    'errors can be either "unknown", "formal", "calibrated" '
                    'or None.')
            elif header['errors'] in valid_err and errors in valid_err[:-1]:
                if (errors, header['errors']) == valid_err[2:]:
                    err_cols = (7, 8)
                elif header['errors'] != errors:
                    raise ValueError(
                        'This model has no {} errors.'.format(errors))
                else:
                    err_cols = (5, 6)

        cilm = _np.tile(_np.zeros((lmax + 1, lmax + 1)), (4, 1, 1))
        ref_epoch = _np.zeros((lmax + 1, lmax + 1))
        trnd = _np.zeros_like(cilm)
        periodic = {}

        # read coefficients
        for line in f:
            line = line.lower().strip().split()
            for i in range(1, len(line)):
                line[i] = line[i].replace('d', 'e')

            l, m = int(line[1]), int(line[2])
            if m > lmax:
                break
            if l > lmax:
                continue

            key = line[0]

            value_cs = [float(line[3]), float(line[4]), 0, 0]
            if errors:
                value_cs[2:] = float(line[err_cols[0]]),\
                    float(line[err_cols[1]])

            if key == 'gfc':
                cilm[:, l, m] = value_cs
            elif key == 'gfct':
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-2])
                    t1i = _yyyymmdd_to_year_fraction(line[-1])
                    if not t0i <= epoch < t1i:
                        raise ValueError('epoch is not in the valid '
                                         'time interval of the model.')
                else:
                    t0i = _yyyymmdd_to_year_fraction(line[-1])

                cilm[:, l, m] = value_cs
                ref_epoch[l, m] = t0i
            elif key in ('trnd', 'dot'):
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-2])
                    t1i = _yyyymmdd_to_year_fraction(line[-1])
                    if not t0i <= epoch < t1i:
                        raise ValueError('epoch is not in the valid '
                                         'time interval of the model.')
                trnd[:, l, m] = value_cs
            elif key in ('acos', 'asin'):
                if is_v2:
                    t0i = _yyyymmdd_to_year_fraction(line[-3])
                    t1i = _yyyymmdd_to_year_fraction(line[-2])
                    if not t0i <= epoch < t1i:
                        raise ValueError('epoch is not in the valid '
                                         'time interval of the model.')

                period = float(line[-1])
                if period not in periodic:
                    arr = _np.zeros_like(cilm)
                    periodic[period] = {'acos': arr,
                                        'asin': arr.copy()}

                periodic[period][key][:, l, m] = value_cs
            else:
                if not quiet:
                    print('Unrecognized keyword in data section of ICGEM '
                          'file will be ignored : {:}'.format(key))

    if epoch is None:
        epoch = ref_epoch

    cilm += _time_variable_part(epoch, ref_epoch, trnd, periodic)

    if errors:
        return cilm[:2], gravity_constant, radius, cilm[2:]
    else:
        return cilm[:2], gravity_constant, radius


def write_icgem_gfc(filename, coeffs, errors=None, header=None, lmax=None,
                    modelname=None, product_type='gravity_field',
                    earth_gm=None, gm=None, r0=None, error_kind=None,
                    tide_system='unknown', normalization='4pi', format=None,
                    encoding=None):
    """
    Write real spherical harmonic gravity coefficients to an ICGEM formatted
    file.

    Usage
    -----
    write_icgem_gfc(filename, coeffs, [errors, header, lmax, modelname, gm, r0,
        product_type, earth_gm, error_kind, tide_system, normalization, format,
        encoding)

    Parameters
    ----------
    filename : str
        The filename to save the spherical harmonic ICGEM-formatted
        coefficients. If filename ends with '.gz' the file will be compressed
        using gzip.
    coeffs : ndarray, size (2, lmax + 1, lmax + 1)
        Array of '4pi' or 'unnorm' normalized spherical harmonic coefficients.
    errors : ndarray, optional, shape (2, lmax + 1, lmax + 1)
        Array of the spherical harmonic error coefficients.
    header : str, optional default = None
        An arbitrary string to be written directly before the ICGEM header.
    lmax : int, optional, default = None
        Maximum degree to write to the file. The default is to write all
        coefficients.
    modelname : str, optional, default = None
        The name of the model for 'icgem' formatted files.
    product_type : str, optional, default = 'gravity_field'
        The type of ICGEM product.
    earth_gm : float
        Gravitational constant of the Earth, in m**3/s**2.
    gm : float
        Gravitational constant of the model, in m**3/s**2.
    r0 : float
        Reference radius of the model, in meters.
    error_kind : str, optional, default = None
        Which errors to write. Can be either 'unknown', 'calibrated', or
        'formal'.
    tide_system : str, optional, default = 'unknown'
        The tide system: 'zero_tide', 'tide_free', or 'unknown'.
    normalization : str, optional, default = '4pi'
        The normalization of the spherical harmonic coefficients: either '4pi'
        or 'unnorm'.
    format : str, optional, default = None
        The format of the ICGEM spherical harmonic coefficients.
    encoding : str, optional, default = None
        Encoding of the output file. The default is to use the system default.
    """
    valid_err = ('unknown', 'calibrated', 'formal', 'calibrated_and_formal')
    valid_tide = ('zero_tide', 'tide_free', 'unknown')
    valid_norm = ('4pi', 'unnorm')

    if lmax is None:
        lmax = coeffs.shape[1] - 1

    if errors is None:
        error_kind = 'no'
    else:
        if error_kind not in valid_err[:-1]:
            error_kind = 'unknown'

    if tide_system not in valid_tide:
        raise ValueError('tide_system can be either "zero_tide", "tide_free" '
                         'or "unknown". Input value is {:s}'
                         .format(repr(tide_system)))

    if normalization not in valid_norm:
        raise ValueError('normalization can be either "4pi", or "unnorm". '
                         'Input value is {:s}'.format(repr(normalization)))

    if filename[-3:] == '.gz':
        filebase = filename[:-3]
    else:
        filebase = filename

    with open(filebase, mode='w', encoding=encoding) as file:
        if header is not None:
            file.write(header + '\n')

        # write gfc header
        if errors is not None:
            file.write('begin_of_head ' + 113*'=' + '\n')
        else:
            file.write('begin_of_head ' + 59*'=' + '\n')

        if modelname is not None:
            file.write('{:<28}{:}\n'.format('modelname', modelname))
        if product_type is not None:
            file.write('{:<28}{:}\n'.format('product_type', product_type))
        if earth_gm is not None:
            file.write('{:<28}{:}\n'.format('earth_gravity_constant',
                                            earth_gm))
        if gm is not None:
            file.write('{:<28}{:}\n'.format('gravity_constant', gm))
        if r0 is not None:
            file.write('{:<28}{:}\n'.format('radius', r0))
        if lmax is not None:
            file.write('{:<28}{:}\n'.format('max_degree', lmax))
        if error_kind != 'no':
            file.write('{:<28}{:}\n'.format('errors', error_kind))
        if tide_system is not None:
            file.write('{:<28}{:}\n'.format('tide_system', tide_system))
        if normalization == '4pi':
            file.write('{:<28}{:}\n'.format('norm', 'fully_normalized'))
        elif normalization == 'unnormalized':
            file.write('{:<28}{:}\n'.format('norm', 'unnormalized'))
        if format is not None:
            file.write('{:<28}{:}\n'.format('format', format))

        if errors is not None:
            file.write('\nkey   {:>5}   {:>5}   {:>24}   {:>24}   {:>24}   '
                       '{:>24}\n'.format('L', 'M', 'C', 'S', 'sigma C',
                                         'sigma S'))
            file.write('end_of_head ' + 115*'=' + '\n')
        else:
            file.write('\nkey   {:>5}   {:>5}   {:>24}   {:>24}\n'
                       .format('L', 'M', 'C', 'S'))
            file.write('end_of_head ' + 61*'=' + '\n')

        # write the coefficients
        for l in range(lmax+1):
            for m in range(l+1):
                if errors is not None:
                    file.write('gfc   {:>5d}   {:>5d}   {:24.16e}   {:24.16e}'
                               '   {:24.16e}   {:24.16e}\n'
                               .format(l, m, coeffs[0, l, m], coeffs[1, l, m],
                                       errors[0, l, m], errors[1, l, m]))
                else:
                    file.write('gfc   {:>5d}   {:>5d}   {:24.16e}'
                               '   {:24.16e}\n'
                               .format(l, m, coeffs[0, l, m], coeffs[1, l, m]))

    if filename[-3:] == '.gz':
        with open(filebase, 'rb') as f_in:
            with _gzip.open(filename, 'wb') as f_out:
                _shutil.copyfileobj(f_in, f_out)


def _time_variable_part(epoch, ref_epoch, trnd, periodic):
    """
    Return sum of the time-variable part of the coefficients

    The formula is:
    G(t) = G(t0) + trnd*(t-t0) +
        asin1*sin(2pi/p1 * (t-t0)) + acos1*cos(2pi/p1 * (t-t0)) +
        asin2*sin(2pi/p2 * (t-t0)) + acos2*cos(2pi/p2 * (t-t0))

    This function computes all terms after G(t0).
    """
    delta_t = epoch - ref_epoch
    trend = trnd * delta_t
    periodic_sum = _np.zeros_like(trnd)
    for period in periodic:
        for trifunc in periodic[period]:
            coeffs = periodic[period][trifunc]
            if trifunc == 'acos':
                periodic_sum += coeffs * _np.cos(2 * _np.pi / period * delta_t)
            elif trifunc == 'asin':
                periodic_sum += coeffs * _np.sin(2 * _np.pi / period * delta_t)
    return trend + periodic_sum


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


def _iszipfile(filename):
    """
    Determine if filename is a zip file. Zip files either
        (1) end with '.zip', or
        (2) are located in a subdirectory '/zip/' for files downloaded from
            the ICGEM web site.
    """
    if '/zip/' in filename:
        return True
    elif filename[-4:] == '.zip':
        return True
    else:
        return False
