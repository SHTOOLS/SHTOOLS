"""
Common routines for I/O
"""

import numpy as _np

def _derivative(epoch, ref_epoch, trnd, periodic):
    """Return sum of the time-variable part of the coefficients

    The formula is:
    G(t)=G(t0) + trnd*(t-t0 ) +
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

