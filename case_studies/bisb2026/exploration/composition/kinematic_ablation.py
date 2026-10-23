"""Explicit K=1 and energy-flat K controls of the frozen measurement operator.

These operators fit the original counts. They do NOT divide measured counts by
K, and their recovered oscillator parameters are apparent parameters under a
deliberately incomplete measurement model. Frozen models.py is not modified.
"""
import numpy as np
from models import Operator, acceptance_nodes, kinematic_2d


class AblationOperator(Operator):
    def __init__(self, data, mode='none', reference_meV=1500., **kwargs):
        if mode not in ('none', 'energy_flat'):
            raise ValueError(mode)
        super().__init__(data, **kwargs)
        self.mode = mode
        self.H.fill(0); self.K.fill(0); self.KQ.fill(0)
        for b, center in enumerate(self.q):
            if kwargs.get('geometry', 'full') == 'point':
                q, w = np.array([abs(center)]), np.ones(1)
            else:
                q, w = acceptance_nodes(center, dq=float(data['width']), h=self.h,
                    sigma_probe=float(data['probe']), n=kwargs.get('n_quad', 25),
                    n_probe=kwargs.get('n_probe', 7))
            k = (np.ones_like(q) if mode == 'none'
                 else kinematic_2d(q, reference_meV, float(data['beam_kV'])))
            weight = k*w/(k@w)
            self.K[b] = weight.sum()
            self.KQ[b] = weight@(q/.005)
            j = np.maximum(np.searchsorted(self.radial, q, side='right')-1, 0)
            frac = (q-self.radial[j])/(self.radial[j+1]-self.radial[j])
            row = np.bincount(j, weight*(1-frac), minlength=len(self.radial))
            row += np.bincount(j+1, weight*frac, minlength=len(self.radial))
            self.H[b] = row[None, :]


def moving_mean(y, n):
    """MATLAB odd-width movmean with shrinking windows at both edges."""
    if n % 2 != 1:
        raise ValueError('Use an odd smoothing window')
    from scipy.ndimage import uniform_filter1d
    return (uniform_filter1d(np.asarray(y, float), n, mode='constant') /
            uniform_filter1d(np.ones_like(y, dtype=float), n, mode='constant'))


def corrected_signal(signal, sigma, k):
    """Diagnostic division; its errors must be transformed by the same K."""
    if np.any(np.asarray(k) <= 0):
        raise ValueError('K must be positive')
    return signal/k, sigma/k
