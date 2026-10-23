"""Audit the user's historical 590 map and refit raw-count K ablations.

Run with --historical, --source, --baseline, --output directories. Historical
fitted backgrounds are held fixed for division plots. No image pixels or raw
files are edited. Division is an approximate visualization, not an inversion
of the aperture integral or energy convolution.
"""
from __future__ import annotations
import argparse
import csv
import json
import time
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from models import Operator, acceptance_nodes, kinematic_2d
from dispersion_models import fit_ordered, curves, loss_apex, transfer
from dispersion_run import scores
from kinematic_ablation import AblationOperator, moving_mean, corrected_signal
from io_utils import digest, dump


def read(path):
    return json.loads(path.read_text())


def save(fig, output, name):
    for ext in ['png', 'pdf']:
        fig.savefig(output/f'{name}.{ext}', dpi=180, facecolor='white')
    plt.close(fig)


def csv_out(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def historical_maps(source, output):
    rows, bins, provenance = [], [], []
    for path in sorted((source/'stage2_2d_h5/590_PL2_10w').glob('bin_*.mat')):
        r = loadmat(path, simplify_cells=True)['R']
        f = next(x for x in r['fits'] if x['main'])
        if not f['success'] or abs(r['q']) > .0245:
            continue
        c = f['curves']; e = c['E']; sig = c['observed']-c['background']
        n = int(round(100/np.median(np.diff(e))))
        mask = (e >= 300) & (e <= 1800)
        kold = c['prefactor']/np.interp(1500., e, c['prefactor'])
        # Compare a wide-aperture K on exactly the same old signal/background.
        # This is explicitly not a new joint-background or intrinsic-loss fit.
        qs, ws = acceptance_nodes(r['q'], dq=r['N']*.0005, h=.0066,
                                  sigma_probe=r['sigma_probe'], n=45, n_probe=11)
        kr = kinematic_2d(qs[None, :], e[mask, None], 30.)@ws
        # Only evaluate the positive B1 range, with a 100 meV margin for smoothing.
        use = (e >= 200) & (e <= 1900)
        kw = kinematic_2d(qs[None, :], e[use, None], 30.)@ws
        kw /= np.interp(1500., e[use], kw)
        wide = moving_mean(corrected_signal(sig[use], c['sigma'][use], kw)[0], n)
        wide = wide[(e[use] >= 300)&(e[use] <= 1800)]
        raw = moving_mean(sig, n)[mask]
        div = moving_mean(corrected_signal(sig, c['sigma'], kold)[0], n)[mask]
        energy = e[mask]/1000
        ridge = [float(energy[x.argmax()]) for x in [raw, div, wide]]
        # Count components are the old independently fitted spectra, not new data.
        peaks = np.asarray(c['peaks'])
        if peaks.shape[0] != len(e):
            peaks = peaks.T
        low = np.argmin(f['parameters'][:, 0])
        high = np.argmax(f['parameters'][:, 0])
        predlow = moving_mean(peaks[:, low], n)[mask]
        predhigh = moving_mean(peaks[:, high], n)[mask]
        heights = find_peaks(raw, prominence=.02*raw.max())[0]
        row = dict(q_Ainv=float(r['q']), N=int(r['N']),
            measured_ridge_eV=ridge[0], divided_old_h_ridge_eV=ridge[1],
            divided_actual_h_ridge_eV=ridge[2],
            old_low_E0_eV=float(f['parameters'][low, 0]/1000),
            old_high_E0_eV=float(f['parameters'][high, 0]/1000),
            old_low_count_apex_eV=float(energy[predlow.argmax()]),
            old_high_count_apex_eV=float(energy[predhigh.argmax()]),
            low_fraction_at_total_ridge=float(predlow[raw.argmax()]/
                max(predlow[raw.argmax()]+predhigh[raw.argmax()], 1e-100)),
            local_maxima_eV=';'.join(f'{energy[i]:.3f}' for i in heights),
            divided_old_peak_on_window_edge=bool(np.argmax(div) in (0, len(div)-1)))
        rows.append(row)
        bins.append(dict(q=r['q'], width=r['N']*.0005, energy=energy,
                         signals=[raw, div, wide], low=predlow, high=predhigh))
        provenance.append(dict(path=str(path), sha256=digest(path)))
    if not rows:
        raise ValueError('No historical 590 fit bins')
    rows.sort(key=lambda x: x['q_Ainv']); bins.sort(key=lambda x: x['q'])
    csv_out(output/'historical_ridge.csv', rows)
    dump(output/'historical_inputs.json', provenance)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout='constrained')
    grid = np.arange(-.024375, .0245, .00025)
    titles = ['原图 b：扣背景计数', '原图 c：除以旧 K，h=0.0005',
              '诊断对照：除以 K，h=0.0066']
    for k, ax in enumerate(axes):
        m = np.full((len(energy), len(grid)), np.nan)
        for b in bins:
            ix = abs(grid-b['q']) < b['width']/2-1e-10
            y = b['signals'][k]
            m[:, ix] = (y/y.max())[:, None]
        im = ax.pcolormesh(grid, energy, m, cmap='turbo', vmin=.3, vmax=1,
                           shading='auto', rasterized=True)
        ax.plot([b['q'] for b in bins], [b['energy'][b['signals'][k].argmax()] for b in bins],
                'wo', ms=3.2, mec='#202020', mew=.5)
        ax.set(title=titles[k], xlabel='q / Å⁻¹', ylabel='能量损失 / eV', ylim=(.3, 1.8))
    fig.colorbar(im, ax=axes, shrink=.85, label='每列最大值归一化')
    fig.suptitle('590 真实数据：白点为每列总谱最大值，不能预先认作同一条低能分支\n'
                 '三个面板沿用同一份旧背景；100 meV 滑动均值；后两图为近似除法诊断', fontsize=13)
    save(fig, output, '01_historical_map_audit')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    positive = [r for r in rows if r['q_Ainv'] > 0]
    q = [r['q_Ainv'] for r in positive]
    for key, label, color in [('measured_ridge_eV', '扣背景实测总谱顶', '#333333'),
            ('divided_old_h_ridge_eV', '除以旧 K 后总谱顶', '#ce8847'),
            ('divided_actual_h_ridge_eV', '除以实际狭缝 K 后总谱顶', '#227e96')]:
        axes[0].plot(q, [r[key] for r in positive], 'o-', ms=4, label=label, color=color)
    axes[0].axvspan(.0025, .015, color='#aaaaaa', alpha=.12, label='上一轮拟合 q 范围')
    axes[0].set(xlabel='正 q / Å⁻¹', ylabel='总谱最大值能量 / eV', ylim=(.3, 1.85))
    axes[0].legend(fontsize=9)
    for target in [.001, .0025, .005, .01]:
        b = min(bins, key=lambda b: abs(b['q']-target))
        axes[1].plot(b['energy'], b['signals'][0]/b['signals'][0].max(),
                     label=f'q={b["q"]:.4f}')
    axes[1].set(xlabel='能量损失 / eV', ylabel='归一化扣背景计数', ylim=(0, 1.05))
    axes[1].legend(fontsize=10)
    fig.suptitle('亮脊的上升确实存在；随 q 改变的是整个宽谱的最大值位置')
    save(fig, output, '02_historical_ridge_and_spectra')
    print('Historical map audited:', len(rows), 'bins', flush=True)
    return rows


def ablations(source, baseline, output, starts):
    summary = []
    for sid in ['590', 'n0', '20w']:
        d = dict(np.load(source/f'{sid}.npz'))
        expected = next(x for x in read(source/'manifest.json')['sessions'] if x['session'] == sid)
        assert digest(source/f'{sid}.npz') == expected['prepared_sha256']
        y = d['frames'][:, :250].sum(1)
        mask = (d['energy'] >= 500)&(d['energy'] <= 3400)
        original = {m: read(baseline/sid/'main'/f'{m}.json') for m in ['A', 'D']}
        fits = {'full': original}
        for mode in ['none', 'energy_flat']:
            op = AblationOperator(d, mode)
            fits[mode] = {}
            for model in ['A', 'D']:
                dest = output/sid/mode; dest.mkdir(parents=True, exist_ok=True)
                target = dest/f'{model}.json'
                if target.exists():
                    f = read(target)
                else:
                    warms = [original[model]['parameters']]
                    if model == 'D':
                        warms.append(transfer(fits[mode]['A']['parameters'], 'A', 'D'))
                    t = time.monotonic()
                    f = fit_ordered(op, y, d['sigma_train'], model, mask,
                                    n_starts=starts, warms=warms, seed=20261024)
                    f.update(scores(d, f), seconds=time.monotonic()-t,
                             kinematic_mode=mode, session=sid)
                    dump(target, f)
                fits[mode][model] = f
                print(sid, mode, model, 'delta', np.round(f['delta_eV'], 4),
                      'test', round(f['test_b1_msr'], 4), 'boundary', f['boundary'], flush=True)
        for mode, mf in fits.items():
            for model, f in mf.items():
                el, eh, gl, gh = curves(np.array(f['parameters']), model, np.array([.0025, .015]))
                summary.append(dict(session=sid, kinematic_mode=mode, model=model,
                    delta_low_E0_meV=float(f['delta_eV'][0]*1000),
                    low_E0_start_meV=el[0]*1000, low_E0_end_meV=el[1]*1000,
                    delta_low_apex_meV=float(np.diff(loss_apex(el, gl))[0]*1000),
                    train_msr=f['train_msr'], test_b1_msr=f['test_b1_msr'],
                    success=f['success'], boundary=';'.join(f['boundary'])))
    csv_out(output/'ablation_summary.csv', summary)
    dump(output/'ablation_summary.json', summary)
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), layout='constrained')
    q = np.linspace(.0025, .015, 100)
    for ax, sid in zip(axes, ['590', 'n0', '20w']):
        for mode, label, color in [('full', '完整运动学 K', '#227e96'),
            ('none', '忽略运动学：K=1', '#c4594b'), ('energy_flat', 'K 保留角度权重，冻结能量依赖', '#ac8944')]:
            path = (baseline/sid/'main/D.json' if mode == 'full' else output/sid/mode/'D.json')
            f = read(path); el, eh, gl, gh = curves(np.array(f['parameters']), 'D', q)
            ax.plot(q, el, color=color, lw=2, label=label)
        ax.set(title=sid, xlabel='Q / Å⁻¹', ylabel='低能分量 E₀ / eV', ylim=(.4, 1.4))
        ax.grid(alpha=.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=3, fontsize=10)
    fig.suptitle('同一批原始计数重新拟合：关闭 K 与除以 K 是两种不同操作\n'
                 '固定旧模型、q 点、拟合能窗和噪声；K=1 得到的是表观振子参数', fontsize=13)
    save(fig, output, '03_same_counts_refit')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for key in ['historical', 'source', 'baseline', 'output']:
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--starts', type=int, default=12)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    font = Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'pdf.fonttype': 42})
    signature = {name: digest(Path(__file__).with_name(name)) for name in
                 ['kinematic_audit.py', 'kinematic_ablation.py', 'models.py',
                  'dispersion_models.py', 'dispersion_run.py', 'io_utils.py']}
    signature.update(starts=args.starts, sources={k: str(getattr(args, k).resolve())
                     for k in ['historical', 'source', 'baseline']})
    for sid in ['590', 'n0', '20w']:
        signature[sid] = digest(args.source/f'{sid}.npz')
        for model in ['A', 'D']:
            signature[f'{sid}_{model}'] = digest(args.baseline/sid/'main'/f'{model}.json')
    stamp = args.output/'signature.json'
    if stamp.exists() and read(stamp) != signature:
        raise ValueError('Inputs changed: choose a fresh output directory')
    dump(stamp, signature)
    historical_maps(args.historical, args.output)
    ablations(args.source, args.baseline, args.output, args.starts)
