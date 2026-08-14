"""Regression tests for DUCC-backed Clenshaw-Curtis scalar grids."""

from unittest.mock import patch

import numpy as np
import pyshtools as pysh

from pyshtools.backends import ducc0_wrapper


_CC_BACKEND_ERROR = 'Clenshaw-Curtis transforms require the DUCC backend.'


def main():
    test_grid()
    test_ducc_is_required()
    if not ducc0_wrapper.available():
        print('DUCC is unavailable; skipping CC transform tests.')
        return
    test_direct_wrapper()
    test_real_transform()
    test_complex_transform()
    test_complex_unnorm()
    test_full_bandwidth()
    test_conventions()
    test_lmax_calc()


def test_grid():
    """Check CC dimensions, coordinates, poles, extension, and xarray I/O."""
    lmax = 5
    for extend in (False, True):
        grid = pysh.SHGrid.from_zeros(lmax, grid='CC', extend=extend)
        expected_nlon = 2 * lmax + 2 + extend
        if grid.grid != 'CC' or grid.data.shape != (lmax + 2, expected_nlon):
            raise AssertionError('CC dimensions or grid type are incorrect.')
        if hasattr(grid, 'sampling'):
            raise AssertionError('CC grids must not define DH sampling.')
        np.testing.assert_allclose(grid.lats(),
                                   np.linspace(90.0, -90.0, lmax + 2))
        np.testing.assert_allclose(
            grid.lons(),
            np.linspace(0.0, 360.0 - 360.0 / expected_nlon, expected_nlon)
            if not extend else np.linspace(0.0, 360.0, expected_nlon),
        )
        restored = pysh.SHGrid.from_xarray(grid.to_xarray(), grid='CC')
        if restored.grid != 'CC' or restored.extend != extend:
            raise AssertionError('CC xarray construction changed the grid.')
        np.testing.assert_allclose(restored.data, grid.data)


def test_ducc_is_required():
    """CC transforms must reject the native backend and missing DUCC."""
    coeffs = pysh.SHCoeffs.from_array(_real_coefficients(3))
    grid = pysh.SHGrid.from_zeros(3, grid='CC', extend=False)
    _assert_ducc_error(lambda: coeffs.expand(grid='CC', backend='shtools'))
    _assert_ducc_error(lambda: grid.expand(backend='shtools'))
    with patch.object(ducc0_wrapper, 'available', return_value=False):
        cc_array = pysh.SHGrid.from_array(np.zeros((5, 8)), grid='CC')
        if cc_array.grid != 'CC':
            raise AssertionError(
                'CC arrays must not require DUCC to construct.'
            )
        _assert_ducc_error(lambda: coeffs.expand(grid='CC', backend='ducc'))
        _assert_ducc_error(lambda: grid.expand(backend='ducc'))


def test_direct_wrapper():
    """Exercise all four public DUCC CC wrapper functions directly."""
    real = _real_coefficients(6)
    real_grid = ducc0_wrapper.MakeGridCC(real, extend=False)
    real_recovered = pysh.SHCoeffs.from_array(
        ducc0_wrapper.SHExpandCC(real_grid)
    )
    _assert_coefficients(real_recovered, pysh.SHCoeffs.from_array(real))

    complex_ = _complex_coefficients(6)
    complex_grid = ducc0_wrapper.MakeGridCCC(complex_, extend=False)
    complex_recovered = pysh.SHCoeffs.from_array(
        ducc0_wrapper.SHExpandCCC(complex_grid)
    )
    _assert_coefficients(complex_recovered, pysh.SHCoeffs.from_array(complex_))


def test_real_transform():
    """Round-trip real coefficients through the public CC grid API."""
    coeffs = pysh.SHCoeffs.from_array(_real_coefficients(8))
    default_grid = coeffs.expand(grid='CC', extend=False)
    _assert_coefficients(default_grid.expand(), coeffs)
    for extend in (False, True):
        grid = coeffs.expand(grid='CC', backend='ducc', extend=extend)
        if extend:
            np.testing.assert_allclose(grid.data[:, 0], grid.data[:, -1])
        _assert_coefficients(grid.expand(backend='ducc'), coeffs)


def test_complex_transform():
    """Round-trip complex coefficients through the public CC grid API."""
    coeffs = pysh.SHCoeffs.from_array(_complex_coefficients(8))
    grid = coeffs.expand(grid='CC', backend='ducc', extend=False)
    _assert_coefficients(grid.expand(backend='ducc'), coeffs)


def test_complex_unnorm():
    """Retain the complex unnormalized order scaling used by DUCC."""
    coefficients = _complex_coefficients(6)
    for csphase in (-1, 1):
        coeffs = pysh.SHCoeffs.from_array(
            coefficients, normalization='unnorm', csphase=csphase
        )
        grid = coeffs.expand(grid='CC', backend='ducc', extend=False)
        _assert_coefficients(grid.expand(backend='ducc'), coeffs, atol=2e-8)


def test_full_bandwidth():
    """Retain sectoral and near-sectoral modes that detect a DH fallback."""
    lmax = 8
    for complex_ in (False, True):
        coefficients = np.zeros(
            (2, lmax + 1, lmax + 1),
            dtype=np.complex128 if complex_ else np.float64,
        )
        coefficients[0, lmax, lmax] = 1.0 + 0.25j if complex_ else 1.0
        coefficients[1, lmax, lmax - 1] = -0.5 + 0.75j if complex_ else -0.5
        coeffs = pysh.SHCoeffs.from_array(coefficients)
        grid = coeffs.expand(grid='CC', backend='ducc', extend=False)
        _assert_coefficients(grid.expand(backend='ducc'), coeffs, atol=3e-10)


def test_conventions():
    """Exercise CC round trips for normalization and phase conventions."""
    coefficients = _real_coefficients(6)
    for normalization in ('4pi', 'schmidt', 'unnorm', 'ortho'):
        atol = 2e-8 if normalization == 'unnorm' else 2e-10
        for csphase in (-1, 1):
            coeffs = pysh.SHCoeffs.from_array(
                coefficients, normalization=normalization, csphase=csphase)
            grid = coeffs.expand(grid='CC', backend='ducc', extend=False)
            _assert_coefficients(grid.expand(backend='ducc'), coeffs,
                                 atol=atol)


def test_lmax_calc():
    """Analyze the full CC geometry before truncating to lmax_calc."""
    coefficients = _real_coefficients(8)
    coeffs = pysh.SHCoeffs.from_array(coefficients)
    grid = coeffs.expand(grid='CC', lmax=8, lmax_calc=3, backend='ducc',
                         extend=False)
    recovered = grid.expand(lmax_calc=3, backend='ducc')
    np.testing.assert_allclose(
        recovered.to_array(normalization='4pi', csphase=1),
        coefficients[:, :4, :4], atol=2e-11)


def _assert_ducc_error(function):
    try:
        function()
    except RuntimeError as error:
        if str(error) != _CC_BACKEND_ERROR:
            raise AssertionError('CC failed with an unclear error.') from error
    else:
        raise AssertionError(
            'CC transform unexpectedly used a non-DUCC backend.'
        )


def _real_coefficients(lmax, seed=4):
    rng = np.random.default_rng(seed)
    coefficients = rng.standard_normal((2, lmax + 1, lmax + 1))
    for degree in range(lmax + 1):
        coefficients[:, degree, degree + 1:] = 0.0
    coefficients[1, :, 0] = 0.0
    return coefficients


def _complex_coefficients(lmax, seed=5):
    rng = np.random.default_rng(seed)
    coefficients = _real_coefficients(lmax, seed=seed).astype(np.complex128)
    coefficients += 1j * rng.standard_normal(coefficients.shape)
    for degree in range(lmax + 1):
        coefficients[:, degree, degree + 1:] = 0.0
    coefficients[1, :, 0] = 0.0
    return coefficients


def _assert_coefficients(recovered, expected, atol=2e-10):
    np.testing.assert_allclose(
        recovered.to_array(normalization=expected.normalization,
                           csphase=expected.csphase),
        expected.to_array(normalization=expected.normalization,
                          csphase=expected.csphase), atol=atol)


if __name__ == '__main__':
    main()
