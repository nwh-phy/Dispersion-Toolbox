"""Audit fit nesting, provenance, bootstrap coverage and independent quadrature."""
from __future__ import annotations
import argparse
import json
import tomllib
from pathlib import Path
import numpy as np
from models import Operator
from dispersion_models import columns, curves, fit_ordered
from io_utils import digest, dump


def verify(source, output, require_complete=True, block_output=None):
    protocol = Path(__file__).with_name('dispersion_protocol.toml')
    cfg = tomllib.loads(protocol.read_text())
    result = dict(signature_valid=True, nesting=[], convergence=[], coverage=[], temporal=[], inputs=[])
    upstream = json.loads((source/'manifest.json').read_text())
    for row in upstream['sessions']:
        assert digest(row['source']) == row['source_sha256']
        assert digest(source/f"{row['session']}.npz") == row['prepared_sha256']
        result['inputs'].append(dict(session=row['session'], raw_unchanged=True, prepared_unchanged=True,
            raw_sha256=row['source_sha256'], prepared_sha256=row['prepared_sha256']))
    for sid in cfg['sessions']:
        signature = json.loads((output/sid/'signature.json').read_text())
        for name, value in signature.items():
            path = (source/f'{sid}.npz' if name == 'data' else protocol if name == 'protocol'
                    else Path(__file__).with_name(name))
            assert digest(path) == value, (sid, name, 'hash changed')
        d = dict(np.load(source/f'{sid}.npz'))
        op = Operator(d); fine = Operator(d, n_radial=512, n_quad=45, n_probe=11)
        mask = (op.energy >= 500) & (op.energy <= 3400)
        for condition in cfg['conditions']:
            folder = output/sid/condition['name']
            files = [folder/f'{m}.json' for m in cfg['models']]
            if not all(p.exists() for p in files):
                if require_complete:
                    raise ValueError(f'Missing fits: {folder}')
                continue
            fits = {m: json.loads((folder/f'{m}.json').read_text()) for m in cfg['models']}
            for simple, larger in [('A','B'), ('A','C'), ('B','D'), ('C','D')]:
                diff = fits[larger]['cost']-fits[simple]['cost']
                assert diff < .05, (sid, condition, simple, larger, diff)
                result['nesting'].append(dict(session=sid, condition=condition['name'],
                    simple=simple, larger=larger, cost_difference=diff))
            for m, f in fits.items():
                el, eh, _, _ = curves(np.array(f['parameters']), m, op.radial)
                assert np.all(el > 0)
                if m != 'S':
                    assert np.min(eh-el) >= .025-1e-12
                assert np.all(np.asarray(f['coefficients']) >= 0)
                assert np.all(np.isfinite(f['prediction']))
        for model in cfg['models']:
            f = json.loads((output/sid/'main'/f'{model}.json').read_text())
            a, _ = columns(fine, np.array(f['parameters']), model)
            prediction = np.einsum('bek,bk->be', a, np.asarray(f['coefficients']))
            residual = (prediction-np.asarray(f['prediction']))/d['sigma_train']
            row = dict(session=sid, model=model,
                rms_sigma=float(np.sqrt(np.mean(residual[:,mask]**2))),
                max_abs_sigma=float(np.max(abs(residual[:,mask]))))
            # A small residual-scale convergence error is the relevant threshold,
            # not visually similar normalized curves.
            assert row['rms_sigma'] < .05, row
            result['convergence'].append(row)
        for model in cfg['bootstrap_models']:
            fits = [json.loads(p.read_text()) for p in sorted((output/sid/'bootstrap').glob(f'*_{model}.json'))]
            if require_complete:
                assert len(fits) == cfg['bootstrap_replicates']
            result['coverage'].append(dict(session=sid, model=model, fits=len(fits),
                success=sum(f['success'] for f in fits), boundary=sum(bool(f['boundary']) for f in fits),
                multistart_cost_gap_gt1=sum(abs(f['candidates'][0]['cost']-f['candidates'][1]['cost'])>1 for f in fits)))
        # Inspect serial correlation without pretending a 5-frame block captures
        # all calibration/systematic uncertainty or long-term drift.
        band = (op.energy >= 500) & (op.energy <= 1800)
        signal = d['frames'][:,:250,band].sum((0,2))/d['beam'][:250]
        acf = [float(np.corrcoef(signal[:-lag], signal[lag:])[0,1]) for lag in [1,5,10,20]]
        result['temporal'].append(dict(session=sid, normalized_B1_acf_lags_1_5_10_20=acf))
    result['complete'] = require_complete
    if block_output is not None:
        result['block10'] = []
        bprotocol = Path(__file__).with_name('dispersion_block10.toml')
        for sid in cfg['sessions']:
            signature = json.loads((block_output/sid/'signature.json').read_text())
            for name, value in signature.items():
                path = (source/f'{sid}.npz' if name == 'data' else bprotocol if name == 'protocol'
                        else Path(__file__).with_name(name))
                assert digest(path) == value
            references = json.loads((block_output/sid/'reference_signature.json').read_text())
            for name, value in references.items():
                assert digest(output/sid/'main'/name) == value
            fits = [json.loads(p.read_text()) for p in sorted((block_output/sid/'bootstrap').glob('*_D.json'))]
            assert len(fits) == 60
            values = 1000*np.array([f['delta_eV'] for f in fits])
            result['block10'].append(dict(session=sid, fits=60,
                success=sum(f['success'] for f in fits),boundary=sum(bool(f['boundary']) for f in fits),
                delta_E0_percentiles_meV=np.percentile(values,[2.5,50,97.5],axis=0)))
    dump(output/'verification.json', result)
    print((output/'verification.json').read_text())


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--source', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path); p.add_argument('--partial', action='store_true')
    p.add_argument('--block-output',type=Path)
    a = p.parse_args(); verify(a.source, a.output, not a.partial, a.block_output)
