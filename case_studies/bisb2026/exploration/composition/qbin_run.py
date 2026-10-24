"""Count-preserving 3/5 q-channel comparison with independent single spectra."""
from __future__ import annotations
import argparse
import json
import time
import tomllib
from pathlib import Path
import numpy as np
from prepare import prepare
from models import Operator, ROOT
from qbin_fit import fit_local
from io_utils import dump, digest


def read(path):
    return json.loads(path.read_text())


def run(previous, reference, output, protocol):
    cfg = tomllib.loads(protocol.read_text())
    output.mkdir(parents=True, exist_ok=True)
    signature = {n: digest(Path(__file__).with_name(n)) for n in
        ['prepare.py', 'qbin_fit.py', 'qbin_run.py', 'models.py', 'dispersion_models.py', 'io_utils.py']}
    signature.update({n: digest(ROOT/'src/python'/n) for n in
        ['qeels_frames.py', 'qeels_rebin.py', 'qeels_response.py']})
    signature.update(protocol=digest(protocol), previous=digest(previous/'manifest.json'),
                     reference=digest(reference/'manifest.json'))
    stamp = output/'signature.json'
    if stamp.exists() and read(stamp) != signature:
        raise ValueError('Input/code changed: use a new output directory')
    dump(stamp, signature)
    parity = []
    reference_records = {s['session']: s for s in read(reference/'manifest.json')['sessions']}
    for n in cfg['bin_channels']:
        folder = output/f'bin{n}'
        if not (folder/'manifest.json').exists():
            prepare(previous, folder, protocol=protocol, bin_channels=n)
        records = {s['session']: s for s in read(folder/'manifest.json')['sessions']}
        for sid in ['590', 'n0', '20w']:
            assert digest(folder/f'{sid}.npz') == records[sid]['prepared_sha256']
            assert records[sid]['source_sha256'] == reference_records[sid]['source_sha256']
            assert digest(reference/f'{sid}.npz') == reference_records[sid]['prepared_sha256']
            data = dict(np.load(folder/f'{sid}.npz')); old = dict(np.load(reference/f'{sid}.npz'))
            np.testing.assert_allclose(data['energy'], old['energy'], atol=1e-10)
            np.testing.assert_allclose(data['width'], n*data['dq'], atol=1e-15)
            if np.isclose(float(data['width']), float(old['width'])):
                idx = [int(np.argmin(abs(old['q']-q))) for q in data['q']]
                for key in ['frames', 'sigma_train', 'sigma_test', 'sigma_all']:
                    np.testing.assert_allclose(data[key], old[key][idx], rtol=1e-12, atol=1e-7)
                parity.append(dict(session=sid, channels=n, reproduced_old_frames_and_noise=True))
            nt = cfg['train_frames']
            scale = float(data['beam'][nt:].sum()/data['beam'][:nt].sum())
            mask = (data['energy'] >= cfg['fit_window_meV'][0]) & (data['energy'] <= cfg['fit_window_meV'][1])
            score = (data['energy'] >= cfg['score_window_meV'][0]) & (data['energy'] <= cfg['score_window_meV'][1])
            for qi, center in enumerate(data['q']):
                local = {**data, 'q': data['q'][qi:qi+1]}
                op = Operator(local)
                y = data['frames'][qi:qi+1, :nt].sum(1)
                test = data['frames'][qi:qi+1, nt:].sum(1)
                sigma = data['sigma_train'][qi:qi+1]
                dest = folder/'fits'/sid/f'{center:+.4f}'; dest.mkdir(parents=True, exist_ok=True)
                fits = {}
                for model in ['S', 'A']:
                    target = dest/f'{model}.json'
                    if target.exists():
                        f = read(target)
                    else:
                        t = time.monotonic()
                        f = fit_local(op, y, sigma, model, mask, starts=cfg['starts'],
                            seed=cfg['seed']+qi, single=fits.get('S'))
                        r = (test-scale*f['prediction'])/data['sigma_test'][qi:qi+1]
                        f.update(session=sid, channels=n, q_Ainv=float(center), width_Ainv=float(data['width']),
                            test_b1_msr=float(np.mean(r[:, score]**2)), test_scale=scale,
                            seconds=time.monotonic()-t, train_frames=nt)
                        dump(target, f)
                    fits[model] = f
                    assert f['success'], target
                    if model == 'A':
                        assert f['cost'] <= fits['S']['cost']+1e-3, 'Nested model optimized poorly'
                    print(sid, n, center, model, 'test', round(f['test_b1_msr'], 3),
                          'E0', round(f['E0_low_eV'], 3), f['E0_high_eV'],
                          'boundary', f['boundary'], flush=True)
    dump(output/'parity.json', parity)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for key in ['previous', 'reference', 'output']:
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--protocol', type=Path, default=Path(__file__).with_name('qbin_protocol.toml'))
    a = p.parse_args(); run(a.previous, a.reference, a.output, a.protocol)
