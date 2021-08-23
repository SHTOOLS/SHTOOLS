# write_icgem_gfc()

Write real spherical harmonic gravity coefficients to an ICGEM formatted
file.

# Usage

write_icgem_gfc(filename, coeffs, [errors, header, lmax, modelname, gm, r0,
    product_type, earth_gm, error_kind, tide_system, normalization, format,
    encoding)

# Parameters

filename : str
:   The filename to save the spherical harmonic ICGEM-formatted
    coefficients. If filename ends with '.gz' the file will be compressed
    using gzip.

coeffs : ndarray, size (2, lmax + 1, lmax + 1)
:   Array of '4pi' or 'unnorm' normalized spherical harmonic coefficients.

errors : ndarray, optional, shape (2, lmax + 1, lmax + 1)
:   Array of the spherical harmonic error coefficients.

header : str, optional default = None
:   An arbitrary string to be written directly before the ICGEM header.

lmax : int, optional, default = None
:   Maximum degree to write to the file. The default is to write all
    coefficients.

modelname : str, optional, default = None
:   The name of the model for 'icgem' formatted files.

product_type : str, optional, default = 'gravity_field'
:   The type of ICGEM product.

earth_gm : float
:   Gravitational constant of the Earth, in m\*\*3/s\*\*2.

gm : float
:   Gravitational constant of the model, in m\*\*3/s\*\*2.

r0 : float
:   Reference radius of the model, in meters.

error_kind : str, optional, default = None
:   Which errors to write. Can be either 'unknown', 'calibrated', or
    'formal'.

tide_system : str, optional, default = 'unknown'
:   The tide system: 'zero_tide', 'tide_free', or 'unknown'.

normalization : str, optional, default = '4pi'
:   The normalization of the spherical harmonic coefficients: either '4pi'
    or 'unnorm'.

format : str, optional, default = None
:   The format of the ICGEM spherical harmonic coefficients.

encoding : str, optional, default = None
:   Encoding of the output file. The default is to use the system default.

