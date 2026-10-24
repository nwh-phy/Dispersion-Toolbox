"""Local-fit physics parity and noiseless recovery at an independent q bin."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'case_studies/bisb2026/exploration/composition'))
from models import Operator
from dispersion_models import encode, columns
from qbin_fit import local_columns, fit_local


def test_local_columns_equal_frozen_constant_dispersion_operator():
    native = np.arange(200., 3400., 4.)
    d = dict(energy_native=native, energy=native.reshape(-1, 4).mean(1),
        q=np.array([.01]), width=.0025, h=.0066, probe=.0004,
        beam_kV=30., resolution_fwhm=22.)
    op = Operator(d)
    for model in ['A', 'S']:
        p = encode([.85]*2, [1.45]*2, model); p[[5, 7, 12]] = 0.
        a, names = local_columns(op, p, model)
        b, frozen_names = columns(op, p, model)
        assert names == frozen_names
        np.testing.assert_allclose(a, b, atol=1e-12)


def test_independent_local_fit_recovers_two_components():
    native = np.arange(200., 3400., 4.)
    d = dict(energy_native=native, energy=native.reshape(-1, 4).mean(1),
        q=np.array([.015]), width=.0025, h=.0066, probe=.0004,
        beam_kV=30., resolution_fwhm=22.)
    fine = Operator(d, n_radial=512, n_quad=45, n_probe=11)
    op = Operator(d)
    p = encode([.88]*2, [1.48]*2, 'A'); p[[4, 6]] = [.55, .75]; p[[5, 7, 12]]=0.
    col, _ = local_columns(fine, p, 'A')
    y = np.einsum('bek,bk->be', col, np.array([[1500, 2300, 500, 600, 40, 8]]))
    sigma = np.sqrt(y); mask = op.energy >= 500
    single = fit_local(op, y, sigma, 'S', mask, starts=3)
    double = fit_local(op, y, sigma, 'A', mask, starts=4, single=single)
    assert double['success'] and double['cost'] < single['cost']
    assert abs(double['E0_low_eV']-.88) < .005
    assert abs(double['E0_high_eV']-1.48) < .005
    assert double['train_msr'] < .001
