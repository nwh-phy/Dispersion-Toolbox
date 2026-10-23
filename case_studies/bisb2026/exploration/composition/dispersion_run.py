"""Resumable real-data comparison and refitted paired frame-block bootstrap."""
from __future__ import annotations
import argparse
import json
import time
import tomllib
from pathlib import Path
import numpy as np
from models import Operator
from dispersion_models import fit_ordered, transfer, curves, loss_apex
from io_utils import digest, dump


def block_indices(n, length, rng):
    start = rng.integers(0, n, int(np.ceil(n/length)))
    return ((start[:, None]+np.arange(length)) % n).ravel()[:n]


def scores(d, fit, indices=None, test_indices=None):
    tr = np.arange(250) if indices is None else indices
    te = np.arange(250, 300) if test_indices is None else test_indices
    scale = float(d['beam'][te].sum()/d['beam'][tr].sum())
    residual = (d['frames'][:, te].sum(1)-scale*np.asarray(fit['prediction']))/d['sigma_test']
    energy = d['energy']
    mask = (energy >= 500) & (energy <= 1800)
    return dict(test_b1_msr=float(np.mean(residual[:, mask]**2)),
        test_full_msr=float(np.mean(residual[:, energy >= 500]**2)), test_scale=scale,
        test_b1_by_q=np.mean(residual[:, mask]**2, axis=1))


def main(source, output, sid, protocol, stage):
    cfg = tomllib.loads(protocol.read_text())
    dest = output/sid; dest.mkdir(parents=True, exist_ok=True)
    signature = {name: digest(Path(__file__).with_name(name)) for name in
                 ['models.py', 'dispersion_models.py', 'dispersion_run.py', 'io_utils.py']}
    signature.update(protocol=digest(protocol), data=digest(source/f'{sid}.npz'))
    stamp = dest/'signature.json'
    if stamp.exists() and json.loads(stamp.read_text()) != signature:
        raise ValueError('Stale cache: select a fresh output directory')
    dump(stamp, signature)
    d = dict(np.load(source/f'{sid}.npz')); y = d['frames'][:, :250].sum(1)
    op = Operator(d)
    if stage in ('main', 'sensitivity'):
        for condition in cfg['conditions']:
            name = condition['name']
            if (stage == 'main') != (name == 'main'):
                continue
            path = dest/name; path.mkdir(exist_ok=True)
            window = condition.get('window', cfg['energy_window_meV'])
            mask = (op.energy >= window[0]) & (op.energy <= window[1])
            exclude = condition.get('exclude')
            if exclude:
                mask &= ~((op.energy >= exclude[0]) & (op.energy <= exclude[1]))
            operator = op if 'h_scale' not in condition else Operator(d, h_scale=condition['h_scale'])
            opts = {k: condition[k] for k in ['include_constant', 'fixed_width', 'numerator'] if k in condition}
            for model in cfg['models']:
                target = path/f'{model}.json'
                if target.exists():
                    continue
                warms = []
                own = dest/'main'/f'{model}.json'
                if own.exists():
                    warms.append(json.loads(own.read_text())['parameters'])
                # Include nested solutions, so a larger model cannot lose just
                # because its random initialization misses the simpler fit.
                for smaller in {'B': ['A'], 'C': ['A'], 'D': ['B', 'C', 'A']}.get(model, []):
                    f = path/f'{smaller}.json'
                    if f.exists():
                        p = json.loads(f.read_text())['parameters']
                        warms.append(transfer(p, smaller, model))
                t = time.monotonic()
                f = fit_ordered(operator, y, d['sigma_train'], model, mask, warms=warms,
                    n_starts=cfg['main_starts' if name == 'main' else 'sensitivity_starts'],
                    seed=cfg['seed'], **opts)
                f.update(scores(d, f)); f.update(session=sid, condition=condition,
                    q_Ainv=d['q'], energy_meV=d['energy'], seconds=time.monotonic()-t)
                if exclude:
                    em = (op.energy >= exclude[0]) & (op.energy <= exclude[1])
                    r = (d['frames'][:, 250:].sum(1)-f['test_scale']*f['prediction'])/d['sigma_test']
                    f['test_energy_holdout_msr'] = float(np.mean(r[:, em]**2))
                dump(target, f)
                print(sid, name, model, 'train', round(f['train_msr'], 3),
                      'testB1', round(f['test_b1_msr'], 3), 'delta', np.round(f['delta_eV'], 4),
                      'boundary', f['boundary'], 'sec', round(f['seconds'], 1), flush=True)
    elif stage == 'bootstrap':
        directory = dest/'bootstrap'; directory.mkdir(exist_ok=True)
        reference = {m: json.loads((dest/'main'/f'{m}.json').read_text()) for m in cfg['bootstrap_models']}
        q = np.linspace(.0025, .015, 51)
        mask = (op.energy >= 500) & (op.energy <= 3400)
        for rep in range(cfg['bootstrap_replicates']):
            rng = np.random.default_rng(cfg['seed']+rep*104729)
            tr = block_indices(250, cfg['block_length'], rng)
            te = 250+block_indices(50, cfg['block_length'], rng)
            yy = d['frames'][:, tr].sum(1)
            for model in cfg['bootstrap_models']:
                target = directory/f'{rep:04d}_{model}.json'
                if target.exists():
                    continue
                ref = reference[model]; p = np.array(ref['parameters'])
                jitter = p.copy()
                jitter[:4] += rng.normal(0, .06, 4)
                jitter[4:8] += rng.normal(0, .1, 4)
                f = fit_ordered(op, yy, d['sigma_train'], model, mask, warms=[p, jitter],
                    n_starts=cfg['bootstrap_starts'], seed=cfg['seed']+rep)
                f.update(scores(d, f, tr, te))
                el, eh, gl, gh = curves(np.array(f['parameters']), model, q)
                f.update(q_curve=q, energy_low_eV=el, energy_high_eV=eh,
                    apex_low_eV=loss_apex(el, gl), apex_high_eV=loss_apex(eh, gh),
                    replicate=rep, train_indices=tr, test_indices=te)
                for key in ['components', 'prediction', 'coefficients']:
                    f.pop(key)
                dump(target, f)
            if rep % 10 == 0:
                print(sid, 'bootstrap', rep+1, '/', cfg['bootstrap_replicates'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--session', choices=['590', 'n0', '20w'], required=True)
    parser.add_argument('--stage', choices=['main', 'sensitivity', 'bootstrap'], required=True)
    parser.add_argument('--protocol', type=Path, default=Path(__file__).with_name('dispersion_protocol.toml'))
    args = parser.parse_args()
    main(args.source, args.output, args.session, args.protocol, args.stage)
