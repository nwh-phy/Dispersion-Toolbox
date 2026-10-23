import importlib.util
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.constants import hbar, e, epsilon_0, m_e

PATH = Path(__file__).resolve().parents[2]/'case_studies/bisb2026/exploration/two_mode_theory/models.py'
spec = importlib.util.spec_from_file_location('two_mode_theory', PATH)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_bilayer_analytic_limits():
    q = np.geomspace(1e-8, .1, 60)
    en, _ = m.bilayer_modes(q, 100, 100, 15)
    np.testing.assert_allclose(en[:, 0]**2, 100*q*(-np.expm1(-q*15)), rtol=1e-9)
    np.testing.assert_allclose(en[:, 1]**2, 100*q*(1+np.exp(-q*15)), rtol=1e-12)
    assert abs(en[0, 0]/q[0]/np.sqrt(1500)-1) < 1e-6
    far, _ = m.bilayer_modes(.005, 100, 200, 1e5)
    np.testing.assert_allclose(far, np.sqrt([.5, 1.]))


def test_bilayer_distance_bound_including_unequal_drude_weights():
    q, low, high = .005, .835, 1.417
    d = m.minimum_distance(q, low, high)
    alpha = (low*low+high*high)/(2*q)
    en, _ = m.bilayer_modes(q, alpha, alpha, d)
    np.testing.assert_allclose(en, [low, high])
    for ratio in [.05, .2, .5, 2, 5, 20]:
        en, _ = m.bilayer_modes(q, alpha, alpha*ratio, d)
        assert en[1]/en[0] >= high/low


def test_bilayer_dark_channel_and_direct_matrix_inverse():
    en, weight = m.bilayer_residues(.005, 100, 100, 0)
    assert weight[0] < 1e-28
    energy, q, a1, a2, d = 1.1, .005, 100, 180, 40
    off = np.sqrt(a1*a2)*q*np.exp(-q*d)
    mat = np.array([[energy**2-a1*q+1j*.4*energy, -off],
                    [-off, energy**2-a2*q+1j*.7*energy]])
    f = np.sqrt([a1*q, a2*q]).astype(complex)
    f[1] *= np.exp(1j*m.longitudinal_q(energy*1000)*d)
    expected = -np.imag(f.conj()@np.linalg.solve(mat, f))
    np.testing.assert_allclose(m.bilayer_loss([energy], [q], a1, a2, d, .4, .7), expected)


def test_hybrid_targets_epsilon_zeros_and_modal_identity():
    pars = m.hybrid_from_targets(.005, .835, 1.417, 1.05, 500)
    energies, _, _ = m.hybrid_modes(.005, *pars, 1.05, 500)
    np.testing.assert_allclose(energies, [.835, 1.417])
    eps = m.hybrid_epsilon(energies, [.005], *pars, 1.05, 500, 1e-10, 1e-10)
    assert np.max(abs(eps)) < 2e-8
    egrid, qs = np.linspace(.03, 4, 600), [.001, .005, .015]
    whole = m.hybrid_loss(egrid, qs, *pars, 1.05, 500, .65, .65)
    parts = m.hybrid_components(egrid, qs, *pars, 1.05, 500, .65)
    np.testing.assert_allclose(whole, parts.sum(2), atol=1e-13, rtol=1e-11)


def test_passivity_and_retarded_poles():
    egrid, qs = np.linspace(.001, 5, 1001), np.geomspace(1e-6, .03, 20)
    for g1, g2 in [(.01, .9), (.7, .2), (.6, .6)]:
        a = m.bilayer_loss(egrid, qs, 200, 350, 30, g1, g2)
        b = m.hybrid_loss(egrid, qs, 200, 80, 1.05, 100, g1, g2)
        assert a.min() >= -1e-12 and b.min() >= -1e-12
    energies, _, _ = m.hybrid_modes(.005, 200, 80, 1.05, 100)
    for E in energies:
        assert np.all(np.roots([1, .6j, -E**2]).imag < 0)


def test_sum_rule_and_zero_lorentz_strength():
    q, a, c, r = .005, 200, 80, 100
    value, _ = quad(lambda E: E*m.hybrid_loss([E], [q], a, c, 1.05, r, .4, .7)[0, 0],
                    0, np.inf, epsabs=1e-8)
    expected = np.pi/2*(a+c)*q/(1+r*q)**2
    assert abs(value/expected-1) < 1e-7
    energies, residues, _ = m.hybrid_modes(q, a, 0, 1.05, r)
    assert np.count_nonzero(residues > 1e-15) == 1
    assert np.any(np.isclose(energies, np.sqrt(a*q/(1+r*q))))


def test_small_q_and_coulomb_off_limits():
    q = np.array([1e-8, 4e-8])
    en, _, _ = m.hybrid_modes(q, 200, 80, 1.05, 100)
    assert abs(en[1, 0]/en[0, 0]-2) < 1e-5
    np.testing.assert_allclose(en[:, 1], 1.05, rtol=2e-6)
    en, weight, _ = m.hybrid_modes(.005, 0, 0, 1.05, 100)
    np.testing.assert_allclose(en, [0, 1.05])
    np.testing.assert_allclose(weight, [0, 0])


def test_stiffness_density_units():
    n_m2, mass, kappa = 1e18, .2*m_e, 4
    alpha_si = hbar**2*e**2*n_m2/(2*epsilon_0*kappa*mass)
    alpha = alpha_si/(e**2*1e-10)
    expected = n_m2/1e4/(kappa*.2)
    np.testing.assert_allclose(m.density_mass_combination(alpha), expected)
