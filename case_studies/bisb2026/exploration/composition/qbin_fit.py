"""Independent local B1 fits using the frozen aperture/energy-response operator.

No cross-q dispersion is imposed. Parameters are constant inside each local
acceptance window; linewidths include unresolved within-window dispersion.
Both hypotheses use identical near-2/3-eV nuisance bands and background.
"""
import numpy as np
from scipy.optimize import least_squares
from dispersion_models import encode, endpoints, LOWER, UPPER, NAMES, loss_apex
from models import oscillator, solve_amplitudes


def local_columns(op, p, model):
    low, high = endpoints(p, model)
    pairs = [(low[0], p[4])]
    names = ['low']
    if model == 'A':
        pairs += [(high[0], p[6])]; names += ['high']
    pairs += [(p[8], p[9]), (p[10], p[11])]
    names += ['near_2eV', 'near_3eV', 'elastic_tail', 'constant']
    cols = [op.project(op.KQ*oscillator(op.e, e0, g)) for e0, g in pairs]
    cols += [np.broadcast_to(op.background, cols[0].shape), np.ones_like(cols[0])]
    col = np.stack(cols, -1)
    return col/np.maximum(col.max(axis=1, keepdims=True), 1e-100), names


def fit_local(op, y, sigma, model, mask, starts=12, seed=20261024, single=None,
              warms=(), max_nfev=350):
    if model not in ('S', 'A') or op.nq != 1:
        raise ValueError('One q bin and S/A model required')
    active = np.array([0, 4, 8, 9, 10, 11] if model == 'S'
                      else [0, 1, 4, 6, 8, 9, 10, 11])
    base = encode([.85, .85], [1.43, 1.43], model)
    base[[5, 7, 12]] = 0.
    def unpack(x):
        p = base.copy(); p[active] = x
        p[2:4] = p[:2]
        return p
    def evaluate(x, full=False):
        col, names = local_columns(op, unpack(x), model)
        coeff, pred = solve_amplitudes(col, y, sigma, mask)
        if full:
            return col, names, coeff, pred
        return ((pred-y)/sigma)[:, mask].ravel()
    init = [np.asarray(p).copy() for p in warms]+[base]
    if single is not None:
        sp = np.asarray(single['parameters']); e0 = endpoints(sp, 'S')[0][0]
        # A contains the single oscillator as either the low or high component.
        if e0 < 1.774:
            p = encode([e0]*2, [min(e0+.3, 1.799)]*2, 'A', sp)
            init.append(p)
        if e0 > .326:
            p = encode([max(.301, e0-.4)]*2, [e0]*2, 'A', sp)
            p[6] = sp[4]; init.append(p)
    rng = np.random.default_rng(seed)
    while len(init) < starts:
        p = encode([rng.uniform(.5, 1.1)]*2, [rng.uniform(1.2, 1.75)]*2, model)
        p[[4, 6]] = rng.uniform(.3, 1.8, 2); init.append(p)
    best = None; candidates = []
    for p in init[:starts]:
        x = np.clip(p[active], LOWER[active]+1e-9, UPPER[active]-1e-9)
        r = least_squares(evaluate, x, bounds=(LOWER[active], UPPER[active]),
            x_scale='jac', max_nfev=max_nfev, ftol=2e-7, xtol=2e-7, gtol=2e-7)
        cost = float(r.fun@r.fun)
        candidates.append(dict(cost=cost, success=bool(r.success), parameters=unpack(r.x)))
        if best is None or cost < best[0]:
            best = cost, r
    cost, r = best; p = unpack(r.x)
    col, names, coeff, pred = evaluate(r.x, True)
    lo, hi = endpoints(p, model)
    ratio = (p-LOWER)/(UPPER-LOWER)
    return dict(model=model, parameters=p, cost=cost, success=bool(r.success),
        train_msr=cost/(op.nq*mask.sum()), prediction=pred,
        components=col*coeff[:, None, :], coefficients=coeff, component_names=names,
        E0_low_eV=lo[0], E0_high_eV=hi[0] if model == 'A' else None,
        low_apex_eV=loss_apex(lo[0], p[4]), gamma_low_eV=p[4],
        boundary=[NAMES[i] for i in active if ratio[i]<.002 or ratio[i]>.998],
        candidates=candidates)
