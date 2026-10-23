"""Ordering, nesting and independent-operator recovery of signed dispersion."""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'case_studies/bisb2026/exploration/composition'))
from models import Operator
from dispersion_models import endpoints, encode, transfer, curves, columns, fit_ordered


def data():
    native = np.arange(200., 3400., 4.)
    return dict(energy_native=native, energy=native.reshape(-1, 4).mean(1),
        q=np.array([.0025, .005, .010, .015]), width=.0015,
        h=.0066, probe=.0004, beam_kV=30., resolution_fwhm=22.)


@pytest.mark.parametrize('model', ['A', 'B', 'C', 'D'])
def test_no_crossing_and_both_slope_signs(model):
    rng = np.random.default_rng(43)
    slopes = []
    for _ in range(100):
        p = encode([.8, .8], [1.4, 1.4], model)
        p[:4] = rng.uniform(0, 1, 4)
        el, eh, _, _ = curves(p, model, np.linspace(0, .035, 1001))
        assert np.min(eh-el) >= .025-1e-14
        assert min(el) >= .3 and max(eh) <= 1.8
        slopes.append([el[-1]-el[0], eh[-1]-eh[0]])
    slopes = np.array(slopes)
    for j, varies in enumerate([model in ('B', 'D'), model in ('C', 'D')]):
        if varies:
            assert slopes[:, j].min() < 0 < slopes[:, j].max()
        else:
            assert np.max(abs(slopes[:, j])) == 0


def test_nested_models_have_identical_predictions():
    op = Operator(data())
    p = encode([.8, .8], [1.4, 1.4], 'A')
    ca, names = columns(op, p, 'A')
    for target in ['B', 'C', 'D']:
        ct, nt = columns(op, transfer(p, 'A', target), target)
        np.testing.assert_allclose(ca, ct, atol=1e-14)
        assert nt == names


@pytest.mark.parametrize('change', [.2, -.2])
def test_recovers_signed_low_slope_through_independent_fine_operator(change):
    d = data(); op = Operator(d)
    fine = Operator(d, n_radial=512, n_quad=45, n_probe=11)
    truth = encode([.8, .8+change], [1.43, 1.43], 'B')
    truth[[4, 6]] = [.45, .55]
    truth[[5, 7]] = [0, 0]
    col, _ = columns(fine, truth, 'B')
    coeff = np.broadcast_to([2500, 1800, 350, 500, 30, 8], (4, 6))
    y = np.einsum('bek,bk->be', col, coeff)
    mask = op.energy >= 500
    result = fit_ordered(op, y, np.sqrt(np.maximum(y, 1)), 'B', mask, n_starts=3)
    assert result['success']
    assert result['train_msr'] < .001
    assert result['delta_eV'][0] == pytest.approx(change*.0125/.035, abs=.005)


def test_block_resampling_keeps_joint_frame_indices():
    from dispersion_run import block_indices
    idx = block_indices(250, 5, np.random.default_rng(7))
    assert len(idx) == 250
    assert idx.min() >= 0 and idx.max() < 250
    np.testing.assert_array_equal(np.diff(idx.reshape(-1, 5), axis=1) % 250, 1)
