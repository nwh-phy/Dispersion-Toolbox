"""Physical controls and weighted-division identities for the K audit."""
import sys
from pathlib import Path
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'case_studies/bisb2026/exploration/composition'))
from kinematic_ablation import AblationOperator, corrected_signal, moving_mean
from models import acceptance_nodes, kinematic_2d


@pytest.mark.parametrize('mode', ['none', 'energy_flat'])
def test_integrates_constant_and_linear_radial_response(mode):
    e = np.arange(200., 1800., 4.)
    d = dict(energy_native=e, energy=e.reshape(-1, 4).mean(1), q=np.array([.001, .01]),
             width=.0015, h=.0066, probe=.0004, beam_kV=30., resolution_fwhm=22.)
    op = AblationOperator(d, mode)
    np.testing.assert_allclose(op.H.sum(-1), 1., atol=2e-14)
    for i, center in enumerate(d['q']):
        q, w = acceptance_nodes(center, dq=d['width'], h=d['h'], sigma_probe=d['probe'])
        k = np.ones_like(q) if mode == 'none' else kinematic_2d(q, 1500., 30.)
        expected = (k*w)@q/(k@w)
        np.testing.assert_allclose(op.H[i]@op.radial, expected, atol=1e-15)
        np.testing.assert_allclose(op.KQ[i], expected/.005, atol=1e-12)
    np.testing.assert_allclose(op.H[:, 0], op.H[:, -1])


def test_division_preserves_weighted_residual_when_model_is_divided_too():
    y = np.array([4., 8., 7.]); pred = np.array([3., 9., 5.])
    sigma = np.array([1., 2., 1.]); k = np.array([.1, 3., 4.])
    z, err = corrected_signal(y, sigma, k)
    np.testing.assert_allclose((z-pred/k)/err, (y-pred)/sigma)
    with pytest.raises(ValueError):
        corrected_signal(y, sigma, np.array([0., 1., 1.]))


def test_moving_mean_shrinks_at_edges_without_moving_peak():
    y = np.array([0., 1., 8., 1., 0.])
    np.testing.assert_allclose(moving_mean(y, 3), [.5, 3., 10/3, 3., .5])
