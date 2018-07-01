# spectrum

Return the spectrum of the spherical harmonic coefficients as a function of spherical harmonic degree.

# Usage

`array` = spectrum (`clm`, [`normalization`, `degrees`, `lmax`, `convention`, `unit`, `base`])

# Returns

`array` : ndarray, shape (len(degrees))
:   1-D ndarray of the spectrum.

# Parameters

`clm` : ndarray, shape (2, `lmax` + 1, `lmax` + 1)
:   ndarray containing the spherical harmonic coefficients.

`normalization` : str, optional, default = '4pi'
:   '4pi', 'ortho', 'schmidt', or 'unnorm', for geodesy 4pi normalized, orthonormalized, Schmidt semi-normalized, or unnormalized coefficients, respectively.

`lmax` : int, optional, default = len(clm[0,:,0]) - 1.
:   Maximum spherical harmonic degree to output.

`degrees` : ndarray, optional, default = range(`lmax`+1)
:   Array containing the spherical harmonic degrees where the spectrum is computed.

`convention` : str, optional, default = 'power'
:   The type of spectrum to return: 'power' for power spectrum, 'energy' for energy spectrum, and 'l2norm' for the l2 norm spectrum.

`unit` : str, optional, default = 'per_l'
:   If 'per_l', return the total contribution to the spectrum for each spherical harmonic degree l. If 'per_lm', return the average contribution to the spectrum for each coefficient at spherical harmonic degree l. If 'per_dlogl', return the spectrum per log interval dlog_a(l).

`base` : float, optional, default = 10.
:    The logarithm base when calculating the 'per_dlogl' spectrum.

# Description

This function returns either the power spectrum, energy spectrum, or l2-norm spectrum. Total power is defined as the integral of the function squared over all space, divided by the area the function spans. If the mean of the function is zero, this is equivalent to the variance of the function. The total energy is the integral of the function squared over all space and is 4pi times the total power. The l2-norm is the sum of the magnitude of the coefficients squared.

The output spectrum can be expresed using one of three units. 'per_l' returns the contribution to the total spectrum from all angular orders at degree l. 'per_lm' returns the average contribution to the total spectrum from a single coefficient at degree l, which is equal to the 'per_l' spectrum divided by (2l+1). 'per_dlogl' returns the contribution to the total spectrum from all angular orders over an infinitessimal logarithmic degree band. The contrubution in the band dlog_a(l) is spectrum(l, 'per_dlogl')\*dlog_a(l), where a is the base, and where spectrum(l, 'per_dlogl') is equal to spectrum(l, 'per_l')\*l\*log(a).

# See also

[cross_spectrum](cross_spectrum.html)
