"""Ordered, signed B1 dispersion hypotheses sharing the frozen v3 operator.

E0 is the oscillator denominator parameter, never the measured/loss apex.
Both B1 components have identical Q numerator and linewidth freedoms. Energies
are linear in radial Q, a local empirical diagnostic rather than a material
dispersion law extrapolated to Q=0. Endpoint coordinates enforce ordering over
the entire integration support without forcing either slope to be positive.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import least_squares
from models import oscillator, solve_amplitudes

MODELS = ('A', 'B', 'C', 'D', 'S')
NAMES = ['u_start', 'v_start', 'u_end', 'v_end', 'Gamma_low_eV',
         'Gamma_low_exponent', 'Gamma_high_eV', 'Gamma_high_exponent',
         'E_neighbor_eV', 'Gamma_neighbor_eV', 'E_upper_neighbor_eV',
         'Gamma_upper_neighbor_eV', 'upper_slope_eV_per_qref']
LOWER = np.array([0, 0, 0, 0, .15, -1., .15, -1., 1.7, .03, 2.8, .2, 0.])
UPPER = np.array([1, 1, 1, 1, 3., 1., 3., 1., 2.4, .8, 4., 3., .6])
EMIN, EMAX, GAP, QEND = .3, 1.8, .025, .035


def endpoints(p, model):
    if model not in MODELS:
        raise ValueError(model)
    if model == 'S':
        lo = EMIN+(EMAX-EMIN)*np.array([p[0], p[2]])
        return lo, np.full(2, np.nan)
    l0 = EMIN+(EMAX-EMIN-GAP)*p[0]
    h0 = l0+GAP+(EMAX-l0-GAP)*p[1]
    l1, h1 = l0, h0
    if model == 'B':
        l1 = EMIN+(h0-GAP-EMIN)*p[2]
    elif model == 'C':
        h1 = l0+GAP+(EMAX-l0-GAP)*p[3]
    elif model == 'D':
        l1 = EMIN+(EMAX-EMIN-GAP)*p[2]
        h1 = l1+GAP+(EMAX-l1-GAP)*p[3]
    return np.array([l0, l1]), np.array([h0, h1])


def encode(lo, hi, model, rest=None):
    p = np.array([.35, .65, .35, .65, .8, -.2, 1., -.1,
                  2.06, .53, 3.19, 1.14, .15])
    if rest is not None:
        p[4:] = np.asarray(rest)[4:]
    if model == 'S':
        p[[0, 2]] = (np.asarray(lo)-EMIN)/(EMAX-EMIN)
    else:
        p[0] = (lo[0]-EMIN)/(EMAX-EMIN-GAP)
        p[1] = (hi[0]-lo[0]-GAP)/(EMAX-lo[0]-GAP)
        p[2] = ((lo[1]-EMIN)/(hi[0]-GAP-EMIN) if model == 'B'
                else (lo[1]-EMIN)/(EMAX-EMIN-GAP))
        base = lo[0] if model == 'C' else lo[1]
        p[3] = (hi[1]-base-GAP)/(EMAX-base-GAP)
    return np.clip(p, LOWER+1e-8, UPPER-1e-8)


def transfer(p, source, target):
    lo, hi = endpoints(p, source)
    if source == 'S':
        hi = np.minimum(np.maximum(lo+.3, 1.5), EMAX-.001)
        lo = np.minimum(lo, hi-GAP-.001)
    return encode(lo, hi, target, p)


def curves(p, model, q):
    lo, hi = endpoints(p, model)
    t = np.asarray(q)/QEND
    el = lo[0]+(lo[1]-lo[0])*t
    eh = hi[0]+(hi[1]-hi[0])*t
    gl = p[4]*(np.maximum(q, 1e-7)/.005)**p[5]
    gh = p[6]*(np.maximum(q, 1e-7)/.005)**p[7]
    return el, eh, gl, gh


def loss_apex(e0, gamma):
    c = 2*e0**2-gamma**2
    return np.sqrt((c+np.sqrt(c*c+12*e0**4))/6)


def columns(op, p, model, include_constant=True, numerator='q'):
    el, eh, gl, gh = curves(p, model, op.radial)
    cols = []
    for e0, gamma in [(el, gl)] + ([] if model == 'S' else [(eh, gh)]):
        weight = op.radial/.005 if numerator == 'q' else e0**2
        loss = oscillator(op.e[None, :], e0[:, None], gamma[:, None])*weight[:, None]
        cols.append(op.project(np.einsum('beq,qe->be', op.H, loss, optimize=False)))
    b1 = np.stack(cols, -1)
    b1 /= np.maximum(b1.max(axis=1, keepdims=True), 1e-100)
    # Reuse frozen v3 nuisance columns, including its dispersive upper neighbor.
    nuisance, names = op.columns([.9, .9, .8, 0., *p[8:]], 'single',
        include_constant=include_constant, include_high_tail=True, upper_dispersion=True)
    return np.concatenate([b1, nuisance[:, :, 1:]], -1), (
        ['low']+([] if model == 'S' else ['high'])+names[1:])


def active_indices(model, fixed_width=False):
    fixed = {'A': [2, 3], 'B': [3], 'C': [2], 'D': [], 'S': [1, 3, 6, 7]}[model]
    if fixed_width:
        fixed = fixed+[5, 7]
    return np.array([i for i in range(13) if i not in fixed])


def fit_ordered(op, y, sigma, model, mask, n_starts=10, seed=20261023,
                warms=(), include_constant=True, fixed_width=False,
                numerator='q', max_nfev=220, detailed=True):
    active = active_indices(model, fixed_width)
    base = encode([.85, .85], [1.43, 1.43], model)
    if fixed_width:
        base[[5, 7]] = 0.
    def unpack(x):
        p = base.copy(); p[active] = x
        return p
    def evaluate(x, full=False):
        col, names = columns(op, unpack(x), model, include_constant, numerator)
        coeff, pred = solve_amplitudes(col, y, sigma, mask)
        if full:
            return col, names, coeff, pred
        return ((pred-y)/sigma)[:, mask].ravel()
    rng = np.random.default_rng(seed)
    starts = [np.asarray(w).copy() for w in warms]+[base]
    while len(starts) < n_starts:
        lo = rng.uniform(.55, 1.1, 2); hi = rng.uniform(1.2, 1.7, 2)
        p = encode(lo, hi, model)
        p[[4, 6]] = rng.uniform(.4, 1.5, 2)
        p[[5, 7]] = rng.uniform(-.5, .5, 2)
        starts.append(p)
    candidates = []
    best = None
    for p in starts[:n_starts]:
        p = np.clip(p, LOWER+1e-8, UPPER-1e-8)
        r = least_squares(evaluate, p[active], bounds=(LOWER[active], UPPER[active]),
                          x_scale='jac', max_nfev=max_nfev,
                          ftol=2e-7, xtol=2e-7, gtol=2e-7)
        cost = float(r.fun@r.fun)
        candidates.append(dict(cost=cost, success=bool(r.success), nfev=int(r.nfev),
                               parameters=unpack(r.x).tolist()))
        if best is None or cost < best[0]:
            best = cost, r
    cost, r = best
    p = unpack(r.x); col, names, coeff, pred = evaluate(r.x, True)
    lo, hi = endpoints(p, model)
    ratio = (p-LOWER)/(UPPER-LOWER)
    boundary = [NAMES[i] for i in active if ratio[i] < .002 or ratio[i] > .998]
    delta = [(lo[1]-lo[0])*.0125/QEND, (hi[1]-hi[0])*.0125/QEND]
    result = dict(model=model, parameters=p, parameter_names=NAMES, active=active,
        endpoints_low_eV=lo, endpoints_high_eV=hi, delta_eV=delta,
        cost=cost, train_msr=cost/(op.nq*mask.sum()), success=bool(r.success),
        boundary=boundary, candidates=candidates, parameter_count=len(active)+coeff.size,
        include_constant=include_constant, fixed_width=fixed_width, numerator=numerator)
    if detailed:
        result.update(prediction=pred, components=col*coeff[:, None, :],
                      component_names=names, coefficients=coeff)
    return result
