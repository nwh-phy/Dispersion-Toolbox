"""Adaptive native/3/5-channel observed B1 ridge, with repeatability checks."""
import argparse
import json
import tomllib
from pathlib import Path
import numpy as np
from prepare import prepare
from ridge_extract import extract, bootstrap_peaks
from io_utils import digest, dump
from kinematic_audit import csv_out


def run(previous, output, protocol):
    cfg = tomllib.loads(protocol.read_text())
    output.mkdir(parents=True, exist_ok=True)
    signature = dict(protocol=cfg, previous_manifest=digest(previous/'manifest.json'),
        code={n: digest(Path(__file__).with_name(n)) for n in
              ['prepare.py','ridge_extract.py','ridge_run.py']})
    sigpath = output/'signature.json'
    if sigpath.exists() and json.loads(sigpath.read_text()) != signature:
        raise ValueError('Code/config changed; use a fresh output directory')
    dump(sigpath, signature)
    rows, candidates, metrics, checks = [], [], [], []
    for n in cfg['bin_channels']:
        base = output/f'bin{n}'
        if not (base/'manifest.json').exists():
            prepare(previous, base, protocol, n)
        manifest = json.loads((base/'manifest.json').read_text())
        assert manifest['protocol_sha256'] == digest(protocol)
        assert manifest['previous_manifest_sha256'] == signature['previous_manifest']
        assert manifest['preparation_sha256'] == digest(Path(__file__).with_name('prepare.py'))
        for session_index, sid in enumerate(['590','n0','20w']):
            record = next(s for s in manifest['sessions'] if s['session']==sid)
            assert digest(base/f'{sid}.npz') == record['prepared_sha256']
            d = dict(np.load(base/f'{sid}.npz')); e = d['energy']; frames = d['frames']
            y = frames.astype(float).sum(1)
            # Identical draws in N=3/N=5 retain covariance in the native region.
            boot = bootstrap_peaks(e, frames, cfg, cfg['seed']+session_index)
            boot[:, abs(d['q']) < 1e-10] = np.nan
            np.savez_compressed(base/f'{sid}_peaks.npz', q=d['q'], bootstrap_peaks_meV=boot)
            local_rows = []
            for qi, q in enumerate(d['q']):
                main = extract(e, y[qi], cfg['smooth_fwhm_meV'], cfg['peak_window_meV'], cfg['prominence_fraction'])
                for c in main['candidates']:
                    candidates.append(dict(session=sid, high_channels=n, q_Ainv=float(q), **c))
                peak = main['peak_meV'] if abs(q)>1e-10 else np.nan
                bp = boot[:, qi]; valid = np.isfinite(bp)
                lo, hi = np.nanpercentile(bp, [2.5,97.5]) if valid.any() else [np.nan,np.nan]
                nt=cfg['train_frames']
                first = extract(e, frames[qi,:nt].astype(float).sum(0), cfg['smooth_fwhm_meV'],
                                cfg['peak_window_meV'], cfg['prominence_fraction'])['peak_meV']
                last = extract(e, frames[qi,nt:].astype(float).sum(0), cfg['smooth_fwhm_meV'],
                               cfg['peak_window_meV'], cfg['prominence_fraction'])['peak_meV']
                alternatives = [extract(e,y[qi],width, cfg['peak_window_meV'], cfg['prominence_fraction'])['peak_meV']
                                for width in cfg['smooth_checks_meV']]
                smooth_change = float(np.max(np.abs(np.array(alternatives)-peak)))
                flags = []
                if not np.isfinite(peak): flags.append('no_interior_peak' if abs(q)>1e-10 else 'elastic_centre')
                if valid.mean()<.9: flags.append('bootstrap_detection_below_90pct')
                if not np.isfinite(hi-lo) or hi-lo>200: flags.append('wide_bootstrap_interval')
                if not np.isfinite(first-last) or abs(first-last)>150: flags.append('temporal_split_shift')
                if not np.isfinite(smooth_change) or smooth_change>100: flags.append('smoothing_sensitive')
                row = dict(session=sid, high_channels=n, q_Ainv=float(q),
                    native_channels=int(round(d['width'][qi]/d['dq'])), width_Ainv=float(d['width'][qi]),
                    peak_meV=float(peak), ci_lo_meV=float(lo), ci_hi_meV=float(hi),
                    bootstrap_detection=float(valid.mean()),
                    bootstrap_within_100meV=float(np.mean(valid & (abs(bp-peak)<100))),
                    first250_peak_meV=float(first), last50_peak_meV=float(last),
                    smooth64_peak_meV=float(alternatives[0]), smooth144_peak_meV=float(alternatives[1]),
                    smoothing_max_shift_meV=smooth_change, candidates=len(main['candidates']),
                    stable=not flags, flags=';'.join(flags))
                rows.append(row); local_rows.append(row)
            for sign in [-1,1]:
                indices = np.flatnonzero(d['q']*sign>0)
                a = indices[abs(abs(d['q'][indices])-.001).argmin()]
                b = indices[abs(abs(d['q'][indices])-.015).argmin()]
                delta = boot[:,b]-boot[:,a]; good=np.isfinite(delta)
                interval = np.percentile(delta[good],[2.5,97.5]) if good.any() else [np.nan,np.nan]
                metric = dict(session=sid,high_channels=n,sign=sign,
                    q_start_Ainv=float(d['q'][a]),q_end_Ainv=float(d['q'][b]),
                    start_peak_meV=local_rows[a]['peak_meV'], end_peak_meV=local_rows[b]['peak_meV'],
                    rise_meV=local_rows[b]['peak_meV']-local_rows[a]['peak_meV'],
                    rise_ci_lo_meV=float(interval[0]), rise_ci_hi_meV=float(interval[1]),
                    bootstrap_pairs=int(good.sum()),end_stable=local_rows[b]['stable'])
                metrics.append(metric)
                print(metric,flush=True)
            checks.append(dict(session=sid,high_channels=n,n_q=len(d['q']),
                n_stable=sum(r['stable'] for r in local_rows),source_verified=True,
                max_weight_sum_error=record['max_weight_sum_error']))
    csv_out(output/'peaks.csv',rows); csv_out(output/'candidates.csv',candidates)
    csv_out(output/'rise.csv',metrics); dump(output/'verification.json',checks)
    dump(output/'rise.json',metrics)
    # Identical single-channel region across the two strategies is a data and
    # algorithm parity check, not a statistical replication.
    parity=[]
    for sid in ['590','n0','20w']:
        a=dict(np.load(output/f'bin3/{sid}.npz')); b=dict(np.load(output/f'bin5/{sid}.npz'))
        ia=abs(a['q'])<=cfg['native_q_max_Ainv']+1e-10
        ib=abs(b['q'])<=cfg['native_q_max_Ainv']+1e-10
        np.testing.assert_array_equal(a['q'][ia],b['q'][ib])
        np.testing.assert_array_equal(a['frames'][ia],b['frames'][ib])
        x=np.load(output/f'bin3/{sid}_peaks.npz')['bootstrap_peaks_meV'][:,ia]
        z=np.load(output/f'bin5/{sid}_peaks.npz')['bootstrap_peaks_meV'][:,ib]
        np.testing.assert_allclose(x,z,equal_nan=True,rtol=0,atol=0)
        parity.append(dict(session=sid,native_q_bins=int(ia.sum()),data_and_bootstrap_exact=True))
    dump(output/'native_parity.json',parity)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--previous',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--protocol',type=Path,default=Path(__file__).with_name('ridge_protocol.toml'))
    a=p.parse_args();run(a.previous,a.output,a.protocol)
