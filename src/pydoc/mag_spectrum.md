# spectrum

Return the spectrum of the magnetic field as a function of spherical harmonic degree.

# Usage

`array` = mag_spectrum (`clm`, `a`, `r``, [`potential`, `normalization`, `degrees`, `lmax`, `convention`, `unit`, `base`])

# Returns

`array` : ndarray, shape (len(degrees))
:   1-D ndarray of the spectrum.

# Parameters

`clm` : ndarray, shape (2, `lmax` + 1, `lmax` + 1)
:   ndarray containing the spherical harmonic coefficients.

`a` : float
:   The reference radius of the spherical harmonic coefficients.

`r` : float
:   The radius at which the spectrum is evaluated.

`potential` : bool, optional, default = False
:   If True, calculate the spectrum of the magnetic potential. Otherwise, calculate the spectrum of the magnetic intensity (default).

`normalization` : str, optional, default = 'schmidt'
:   '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized, orthonormalized, or Schmidt semi-normalized, or unnormalized coefficients, respectively.

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

This function returns either the power spectrum, energy spectrum, or l2-norm spectrum of a global magnetic field. Total power is defined as the integral of the function squared over all space, divided by the area the function spans. If the mean of the function is zero, this is equivalent to the variance of the function. The total energy is the integral of the function squared over all space and is 4pi times the total power. The l2-norm is simply the sum of the magnitude of the spherical harmonic coefficients squared. The default behaviour of this function is to return the Lowes-Mauersberger spectrum, which is the total power of the magnetic field strength as a function of spherical harmonic degree (see below).

The magnetic potential and magnetic field are defined respectively by the two equations

`U = a**2 \sum_{l, m}^L (a/r)**(l+1) glm Ylm`,

`B = - \Del U.`

Here, a is the reference radius of the spherical harmonic coefficients, r is the radius at which the function is evalulated, and L is the maximum spherical harmonic degree of the expansion. If the optional parameter potential is set to True, the spectrum of the magnetic potential will be calculated. Otherwise, the default behavior is to calculate the spectrum of the magnetic field intensity.

The output spectrum can be expresed using one of three units. 'per_l' returns the contribution to the total spectrum from all angular orders at degree l. 'per_lm' returns the average contribution to the total spectrum from a single coefficient at degree l, and is equal to the 'per_l' spectrum divided by (2l+1). 'per_dlogl' returns the contribution to the total spectrum from all angular orders over an infinitessimal logarithmic degree band. The contrubution in the band dlog_a(l) is spectrum(l, 'per_dlogl')*dlog_a(l), where a is the base, and where spectrum(l, 'per_dlogl) is equal to spectrum(l, 'per_l')*l*log(a).

When no optional parameters are specified, the Lowes-Mauersberger power spectrum is calculated. Explicitly, this corrresponds to convention = 'power', unit = 'per_l', and potential = False.
