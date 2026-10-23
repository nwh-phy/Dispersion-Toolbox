"""Conditional two-channel electrostatic response; E in eV, q in A^-1.

No crystal structure, density, mass, thickness, or dielectric constant is assigned
to BiSb. Alpha is a charge stiffness, with units eV^2 A. The separated-sheet
model uses a homogeneous nondispersive environment; the coplanar model includes
a real background b(q)=1+r0*q and a dynamic Lorentz polarization explicitly.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy.constants import e, epsilon_0, hbar, m_e

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src/python'))
from qeels_response import longitudinal_q


def bilayer_modes(q, alpha1, alpha2, distance_A):
    """Eigenenergies and normalized charge-coordinate eigenvectors, low first."""
    q = np.asarray(q, dtype=float)
    if np.any(q < 0) or min(alpha1, alpha2) <= 0 or distance_A < 0:
        raise ValueError('Positive stiffnesses and nonnegative q/d are required')
    a, c = alpha1*q, alpha2*q
    off = np.sqrt(a*c)*np.exp(-q*distance_A)
    matrix = np.zeros(q.shape+(2, 2))
    matrix[..., 0, 0] = a
    matrix[..., 1, 1] = c
    matrix[..., 0, 1] = matrix[..., 1, 0] = off
    vals, vecs = np.linalg.eigh(matrix)
    return np.sqrt(np.maximum(vals, 0)), vecs


def bilayer_loss(energy_eV, q, alpha1, alpha2, distance_A,
                 gamma1_eV, gamma2_eV, beam_kV=30.):
    """Positive projected loss with the normal-incidence electron phase.

    Returns shape (n_energy,n_q). Up to a constant, the probability per d^2q
    is q/(q^2+qE^2)^2 times this loss. No azimuthal integration is implied.
    """
    energy = np.atleast_1d(energy_eV)[:, None]
    q = np.atleast_1d(q)[None, :]
    if np.any(energy <= 0) or min(gamma1_eV, gamma2_eV) <= 0:
        raise ValueError('Positive energy and damping are required')
    bilayer_modes(q.ravel(), alpha1, alpha2, distance_A)
    a, c = alpha1*q, alpha2*q
    off = np.sqrt(a*c)*np.exp(-q*distance_A)
    d1 = energy**2-a+1j*energy*gamma1_eV
    d2 = energy**2-c+1j*energy*gamma2_eV
    phase = longitudinal_q(energy*1000, beam_kV)*distance_A
    numerator = a*d2+c*d1+2*np.sqrt(a*c)*off*np.cos(phase)
    return -np.imag(numerator/(d1*d2-off**2))


def bilayer_residues(q, alpha1, alpha2, distance_A, beam_kV=30.):
    """Undamped residues evaluated at each pole; not full damped peak areas."""
    energy, vec = bilayer_modes(q, alpha1, alpha2, distance_A)
    q = np.asarray(q)
    phase = longitudinal_q(energy*1000, beam_kV)*distance_A
    drive1 = np.sqrt(alpha1*q)[..., None]
    drive2 = np.sqrt(alpha2*q)[..., None]*np.exp(1j*phase)
    return energy, abs(drive1*vec[..., 0, :]+drive2*vec[..., 1, :])**2


def minimum_distance(qref, low_eV, high_eV):
    """Smallest d allowing this pole ratio for any two positive Drude weights."""
    if not (qref > 0 and 0 < low_eV < high_eV):
        raise ValueError('Need qref>0 and 0<low<high')
    contrast = (high_eV**2-low_eV**2)/(high_eV**2+low_eV**2)
    return -np.log(contrast)/qref


def hybrid_from_targets(qref, low_eV, high_eV, transition_eV, r0_A):
    """Calibrate two pole positions; this is explicitly not a prediction of them."""
    if not (qref > 0 and 0 < low_eV < transition_eV < high_eV and r0_A >= 0):
        raise ValueError('Bare transition must lie between target pole energies')
    low2, high2, delta2 = low_eV**2, high_eV**2, transition_eV**2
    drude = low2*high2/delta2
    lorentz = low2+high2-delta2-drude
    factor = (1+r0_A*qref)/qref
    return drude*factor, lorentz*factor


def hybrid_modes(q, alpha_drude, alpha_lorentz, transition_eV, r0_A=0.):
    """Longitudinal poles and residues of inverse epsilon, low first."""
    q = np.asarray(q, dtype=float)
    if np.any(q < 0) or min(alpha_drude, alpha_lorentz, r0_A) < 0 or transition_eV <= 0:
        raise ValueError('Nonnegative stiffness/q/r0 and positive transition required')
    b = 1+r0_A*q
    a, c = alpha_drude*q/b, alpha_lorentz*q/b
    off = np.sqrt(a*c)
    matrix = np.zeros(q.shape+(2, 2))
    matrix[..., 0, 0] = a
    matrix[..., 1, 1] = c+transition_eV**2
    matrix[..., 0, 1] = matrix[..., 1, 0] = off
    vals, vecs = np.linalg.eigh(matrix)
    g = np.stack([np.sqrt(a), np.sqrt(c)], axis=-1)
    residues = np.einsum('...i,...ij->...j', g, vecs)**2/b[..., None]
    return np.sqrt(np.maximum(vals, 0)), residues, vecs


def hybrid_epsilon(energy_eV, q, alpha_drude, alpha_lorentz, transition_eV,
                   r0_A, gamma_drude_eV, gamma_lorentz_eV):
    energy = np.atleast_1d(energy_eV)[:, None]
    q = np.atleast_1d(q)[None, :]
    if np.any(energy <= 0) or min(gamma_drude_eV, gamma_lorentz_eV) <= 0:
        raise ValueError('Positive energy and damping are required')
    hybrid_modes(q.ravel(), alpha_drude, alpha_lorentz, transition_eV, r0_A)
    return (1+r0_A*q-alpha_drude*q/(energy*(energy+1j*gamma_drude_eV))
            -alpha_lorentz*q/(energy**2-transition_eV**2+1j*energy*gamma_lorentz_eV))


def hybrid_loss(energy_eV, q, alpha_drude, alpha_lorentz, transition_eV,
                r0_A, gamma_drude_eV, gamma_lorentz_eV):
    return -np.imag(1/hybrid_epsilon(energy_eV, q, alpha_drude, alpha_lorentz,
                                    transition_eV, r0_A, gamma_drude_eV, gamma_lorentz_eV))


def hybrid_components(energy_eV, q, alpha_drude, alpha_lorentz, transition_eV,
                      r0_A, gamma_eV):
    """Exact positive mode decomposition ONLY for equal damping of both channels."""
    if gamma_eV <= 0:
        raise ValueError('Positive damping required')
    poles, residues, _ = hybrid_modes(np.atleast_1d(q), alpha_drude, alpha_lorentz,
                                      transition_eV, r0_A)
    energy = np.atleast_1d(energy_eV)[:, None, None]
    return residues[None, :, :]*gamma_eV*energy / (
        (energy**2-poles[None, :, :]**2)**2+(gamma_eV*energy)**2)


def density_mass_combination(alpha_eV2_A):
    """n2D/[kappa*(m*/me)] in cm^-2; only valid for a parabolic Drude channel.

    alpha_SI=hbar^2 e^2 n/(2 epsilon0 kappa m). Convert eV^2 A to J^2 m.
    This does not return a carrier density without kappa and m*.
    """
    alpha_si = np.asarray(alpha_eV2_A)*e**2*1e-10
    return alpha_si*2*epsilon_0*m_e/(hbar**2*e**2)/1e4
