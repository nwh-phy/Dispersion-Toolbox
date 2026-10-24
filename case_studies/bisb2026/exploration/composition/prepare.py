"""Extend the verified frame preparation, without editing its archived outputs."""
from __future__ import annotations
import argparse
import json
import sys
import tomllib
from pathlib import Path
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src/python'))
from qeels_frames import momentum_axis, noise_from_frames, sha256_file
from qeels_rebin import fractional_bin_weights


def adaptive_q_bins(dq, native_max, q_max, high_channels):
    """Symmetric, contiguous native-grid bins; omit incomplete outer groups.

    Widths describe integration in q, not errors of the bin centres. Framewise
    drift translation can fractionally share detector pixels at bin edges.
    """
    if high_channels not in (3, 5) or not 0 < native_max < q_max:
        raise ValueError('Require 3/5 high-q channels and 0 < native_max < q_max')
    low = int(np.floor(native_max/dq+1e-8))
    last = int(np.floor(q_max/dq+1e-8))
    groups = [(i, 1) for i in range(1, low+1)]
    groups += [(i+(high_channels-1)/2, high_channels)
               for i in range(low+1, last-high_channels+2, high_channels)]
    q = np.array([-i*dq for i, n in groups[::-1]]+[0.]+[i*dq for i, n in groups])
    counts = np.array([n for i, n in groups[::-1]]+[1]+[n for i, n in groups])
    return q, counts*dq


def prepare(previous: Path, output: Path, protocol: Path | None = None,
            bin_channels: int | None = None):
    if previous.resolve() == output.resolve():
        raise ValueError('Use a separate output directory')
    output.mkdir(parents=True, exist_ok=True)
    if bin_channels is not None and bin_channels not in (1, 3, 5):
        raise ValueError('Only native, 3-channel and 5-channel bins are supported')
    protocol = protocol or Path(__file__).with_name('protocol.toml')
    cfg = tomllib.loads(protocol.read_text())
    old = json.loads((previous / 'manifest.json').read_text())
    manifest = {'protocol': cfg, 'previous_manifest_sha256': sha256_file(previous/'manifest.json'),
                'protocol_sha256': sha256_file(protocol), 'sessions': []}
    if bin_channels is not None:
        manifest['bin_channels'] = bin_channels
        manifest['preparation_sha256'] = sha256_file(Path(__file__))
    for s in old['sessions']:
        sid = s['session']
        d = loadmat(previous/sid/'prepared.mat', simplify_cells=True)
        drift = loadmat(previous/sid/'qdrift.mat', simplify_cells=True)['drift_Ainv']
        raw = np.load(s['source'], mmap_mode='r')
        source_hash = sha256_file(s['source'])
        if source_hash != s['sha256']:
            raise ValueError('Source data differ from verified previous run')
        q = momentum_axis(raw.shape[1], s['q_zero_column'], d['dq'])
        support = d['align_support'].astype(int) - 1
        offsets = d['align_offsets'].astype(int)
        keep = np.flatnonzero(d['align_valid'])
        energy = d['E_raw'][support] - (s['final_energy_origin_meV']-s['metadata_elastic_peak_meV'])
        lo, hi = cfg['prepared_energy_range_meV']
        ix = np.flatnonzero((energy >= lo-1e-8) & (energy <= hi+1e-8))
        n = cfg['energy_sum_channels']; ix = ix[:len(ix)//n*n]
        eb = energy[ix].reshape(-1, n).mean(1)
        width = (s['physical_bin_width_Ainv'] if bin_channels is None
                 else bin_channels*float(d['dq']))
        if 'native_q_max_Ainv' in cfg:
            targets, widths = adaptive_q_bins(float(d['dq']), cfg['native_q_max_Ainv'],
                                              cfg['q_max_Ainv'], bin_channels)
            width = widths
        else:
            targets = np.array(cfg['q_targets_Ainv'])
            widths = np.full(len(targets), width)
        blocks = []; weight_error = 0
        for target, local_width in zip(targets, widths):
            frames = []
            for it, t in enumerate(keep):
                weights = fractional_bin_weights(q, target, local_width, drift[it])
                members = np.flatnonzero(weights > 1e-12)
                if set(members+1) & set(s['invalid_columns_1based']):
                    raise ValueError('Bad detector column in translated bin')
                weight_error = max(weight_error, abs(weights.sum() - local_width/d['dq']))
                y = weights[members] @ raw[t][np.ix_(members, support[ix]+offsets[t])].astype(float)
                frames.append(y.reshape(-1, n).sum(1))
            blocks.append(frames)
        frames = np.array(blocks, dtype=float)  # q, time, summed-energy channel
        beam = d['beam']; nt = cfg['train_frames']
        kwargs = {'smooth_bins': 5}
        sigma_train = np.stack([noise_from_frames(f[:nt], beam[:nt], **kwargs)['sigma'] for f in frames])
        sigma_test = np.stack([noise_from_frames(f[nt:], beam[nt:], **kwargs)['sigma'] for f in frames])
        sigma_all = np.stack([noise_from_frames(f, beam, **kwargs)['sigma'] for f in frames])
        np.savez_compressed(output/f'{sid}.npz', energy=eb, energy_native=energy[ix],
            frames=frames.astype(np.float32), beam=beam, q=targets, width=width, dq=d['dq'],
            h=d['h'], probe=d['sigma_probe'], beam_kV=d['beam_kV'], sigma_train=sigma_train,
            sigma_test=sigma_test, sigma_all=sigma_all, drift_Ainv=drift, offsets=offsets,
            resolution_fwhm=s['zlp_fwhm_after_meV'])
        manifest['sessions'].append(dict(session=sid, source=s['source'], source_sha256=source_hash,
            prepared_sha256=sha256_file(output/f'{sid}.npz'), n_frames=len(keep), n_q=len(targets),
            n_energy=len(eb), width_Ainv=np.asarray(width).tolist(), max_weight_sum_error=weight_error,
            initial_source_sha256=s['sha256'],previous_artifact_sha256={
                name:sha256_file(previous/sid/name) for name in ['prepared.mat','qdrift.mat']}))
        print(f'PREPARED {sid}: {frames.shape}, source hash verified', flush=True)
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--previous', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a=p.parse_args(); prepare(a.previous, a.output)
