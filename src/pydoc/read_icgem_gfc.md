# read_icgem_gfc

Read spherical harmonic coefficients from an ICGEM GFC ascii-formatted file.

This function only reads files with the gravity field spherical
harmonic coefficients.

# Usage

`cilm`, `gm`, `r0`, [`errors`] = read_icgem_gfc (`filename`, [`errors`, `lmax`, `epoch`])

# Returns

`cilm` : array
:   Array with the coefficients with the shape (2, `lmax` + 1, `lmax` + 1) for the given epoch.

`gm` : float
:   Standard gravitational constant of the model, in m**3/s**2.

`r0` : float
:   Reference radius of the model, in meters.

`errors` : array, optional
:   Array with the errors of the coefficients with the shape (2, `lmax` + 1, `lmax` + 1) for the given epoch.

# Parameters

`filename` : str
:   The ascii-formatted filename containing the spherical harmonic coefficients.

`errors` : str, optional
:   Which errors to read. Can be either "calibrated", "formal" oNone. Default is None.

`lmax` : int, optional
:   Maximum degree to read from the file. If lmax is None, less than 0, or greater than lmax_model, the maximum degree of the model will be used.

`epoch` : str or float, optional
:   The epoch time to calculate time-variable coefficients in YYYYMMDD.DD format. If None then reference epoch t0 of the model will be used. If format of the file is 'icgem2.0' then epoch must be specified.
