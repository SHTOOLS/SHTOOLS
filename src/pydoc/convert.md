# convert

Convert an array of spherical harmonic coefficients to a different normalization convention.

# Usage

`coeffs_out` = convert(`coeffs_in`, [`normalization_in`, `normalization_out`, `csphase_in`, `csphase_out`, `lmax`])

# Returns

`coeffs_out` : ndarray, size (2, lmax+1, lmax+1)
:   An array of spherical harmonic coefficients with the new normalization convention.

# Parameters

`coeffs_in` : ndarray
:   The array of imput spherical harmonic coefficients.

`normalization_in` : str, optional, default = 'None'
:   '4pi', 'ortho', 'schmidt', or 'unnorm', for geodesy 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized coefficients, respectively.

`normalization_out` : str, optional, default = 'None'
:   '4pi', 'ortho', 'schmidt', or 'unnorm', for geodesy 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized coefficients, respectively.

`csphase_in` : int, optional, default = None
:   Condon-Shortley phase convention of the input coefficients: 1 to exclude the phase factor, or -1 to include it.

`csphase_out` : int, optional, default = None
:   Condon-Shortley phase convention of the input coefficients: 1 to exclude the phase factor, or -1 to include it.

`lmax` : int, optional, default = coeffs_in.shape[1] - 1.
:   Maximum spherical harmonic degree to output. If lmax is larger than that of the input coefficients, the output array will be zero padded.

# Description

This routine will convert an array of spherical harmonic coefficients to a different normalization convention and different Condon-Shortley phase convention. Optionally, a different maximum spherical harmonic degree can be specified. If this degree is smaller than that of the input coefficients, the input coefficients will be truncated. If this degree is larger than the input coefficients, then the output coefficients will be zero padded.
