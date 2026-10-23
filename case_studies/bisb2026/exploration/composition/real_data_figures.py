"""Figures from the three measured B1 acquisitions and frozen v3 fits.

Rebuilds the positive-q maps from read-only raw frames using the previously
verified alignment. Does not fit new models or divide an aperture-integrated
count spectrum by a nominal-q kinematic factor.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import zipfile
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

from io_utils import digest, dump
from models import ROOT
from summarize import apex, sorted_apices
from qeels_frames import momentum_axis
from qeels_rebin import fractional_bin_weights

SESSIONS = ['590', 'n0', '20w']
COLORS = ['#008c77', '#ca5947']


def read(path):
    return json.loads(path.read_text())


def save(fig, output, name):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(output/f'{name}.{ext}', dpi=190, facecolor='white')
    plt.close(fig)


def verify_file(path, expected):
    actual = digest(path)
    if actual != expected:
        raise ValueError(f'Changed input: {path}')
    return dict(path=str(path), sha256=actual)


def measured_map(previous, session, record, data):
    """No model or background subtraction enters this count reconstruction."""
    raw = np.load(record['source'], mmap_mode='r')
    calibration = loadmat(previous/session/'prepared.mat', simplify_cells=True)
    drift = loadmat(previous/session/'qdrift.mat', simplify_cells=True)['drift_Ainv']
    qraw = momentum_axis(raw.shape[1], record['q_zero_column'], calibration['dq'])
    q = np.linspace(.0025, .015, int(round(.0125/calibration['dq']))+1)
    support = calibration['align_support'].astype(int)-1
    offset = calibration['align_offsets'].astype(int)
    keep = np.flatnonzero(calibration['align_valid'])
    energy = calibration['E_raw'][support]-(record['final_energy_origin_meV']-
                                            record['metadata_elastic_peak_meV'])
    ix = np.flatnonzero((energy >= data['energy_native'][0]-1e-8)&
                       (energy <= data['energy_native'][-1]+1e-8))
    np.testing.assert_allclose(energy[ix], data['energy_native'], atol=1e-7)
    np.testing.assert_allclose(drift, data['drift_Ainv'])
    np.testing.assert_array_equal(offset, data['offsets'])
    groups = len(ix)//len(data['energy'])
    counts = np.zeros((len(q), len(data['energy'])))
    train = np.zeros_like(counts)
    bad = np.asarray(record['invalid_columns_1based'])-1
    for k, t in enumerate(keep):
        weights = np.stack([fractional_bin_weights(qraw, center, float(data['width']), drift[k])
                            for center in q])
        members = np.flatnonzero(weights.sum(0) > 1e-12)
        if np.intersect1d(members, bad).size:
            raise ValueError('Bad detector column enters displayed map')
        np.testing.assert_allclose(weights.sum(1), float(data['width']/data['dq']), atol=1e-10)
        aligned = raw[t][np.ix_(members, support[ix]+offset[t])].astype(float)
        frame = (weights[:, members]@aligned).reshape(len(q), -1, groups).sum(-1)
        counts += frame
        if k < 250:
            train += frame
    # A separate direct raw-data extraction must reproduce every existing +q bin.
    errors = []
    for j in np.flatnonzero(data['q'] > 0):
        n = int(np.argmin(abs(q-data['q'][j])))
        np.testing.assert_allclose(q[n], data['q'][j], atol=1e-12)
        expected = data['frames'][j].astype(float).sum(0)
        np.testing.assert_allclose(counts[n], expected, rtol=2e-7, atol=.005)
        np.testing.assert_allclose(train[n], data['frames'][j, :250].astype(float).sum(0),
                                   rtol=2e-7, atol=.005)
        errors.append(float(np.max(abs(counts[n]-expected)/np.maximum(expected, 1))))
    return dict(q=q, energy=data['energy']/1000, counts=counts, train_counts=train,
                test_counts=counts-train, parity_max_relative_error=max(errors))


def run(previous, composition, output):
    output.mkdir(parents=True, exist_ok=True)
    old = read(previous/'manifest.json')
    manifest = read(composition/'manifest.json')
    inputs = [verify_file(previous/'manifest.json', manifest['previous_manifest_sha256'])]
    records = {s['session']: s for s in old['sessions']}
    frozen = {s['session']: s for s in manifest['sessions']}
    font = Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'path'})
    datasets, maps, fits, checks = {}, {}, {}, []
    for sid in SESSIONS:
        s = frozen[sid]
        inputs.append(verify_file(Path(s['source']), s['source_sha256']))
        inputs.append(verify_file(composition/f'{sid}.npz', s['prepared_sha256']))
        for name, expected in s['previous_artifact_sha256'].items():
            inputs.append(verify_file(previous/sid/name, expected))
        signature = read(composition/'fits_v3'/sid/'signature.json')
        for name in ['models.py', 'compare.py', 'io_utils.py']:
            verify_file(Path(__file__).with_name(name), signature[name])
        verify_file(Path(__file__).with_name('protocol_v3.toml'), signature['protocol'])
        assert signature['data'] == s['prepared_sha256']
        d = dict(np.load(composition/f'{sid}.npz'))
        fpath = composition/'fits_v3'/sid/'main_double.json'
        inputs.append(dict(path=str(fpath), sha256=digest(fpath)))
        f = read(fpath)
        assert f['success'] and f['model'] == 'double' and not f['flat']
        np.testing.assert_allclose(np.asarray(f['components']).sum(-1), f['prediction'], rtol=1e-12)
        np.testing.assert_allclose(f['q_Ainv'], d['q'])
        np.testing.assert_allclose(f['energy_meV'], d['energy'])
        assert d['frames'].shape[1] == 300
        maps[sid] = measured_map(previous, sid, records[sid], d)
        datasets[sid], fits[sid] = d, f
        checks.append(dict(session=sid, frames=300, raw_hash_verified=True,
                           dense_q_bins=len(maps[sid]['q']),
                           parity_max_relative_error=maps[sid]['parity_max_relative_error'],
                           fit_boundary_parameters=[n for n, b in zip(f['parameter_names'], f['boundary']) if b]))
        np.savez_compressed(output/f'{sid}_measured_map.npz', **maps[sid])
        print('REBUILT', sid, checks[-1], flush=True)

    # 1. Real measured maps: actual detector bins, no model-generated pixels.
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.1), layout='constrained')
    for ax, sid in zip(axes, SESSIONS):
        m, d = maps[sid], datasets[sid]
        band = (m['energy'] >= .5)&(m['energy'] <= 1.8)
        z = m['counts'][:, band]
        z = z/z.max(1, keepdims=True)
        im = ax.pcolormesh(m['q'], m['energy'][band], z.T, cmap='magma',
                          shading='auto', vmin=0, vmax=1, rasterized=True)
        ax.set(title=f'{sid} | Δq = {float(d["width"]):.5f} Å⁻¹',
               xlabel=r'$q_x$ / Å$^{-1}$', ylabel='能量损失 / eV')
        ax.set_xticks([.0025, .005, .01, .015], ['0.0025', '0.005', '0.010', '0.015'])
    fig.colorbar(im, ax=axes, shrink=.84, label='每个 q 的 B1 能窗内归一化计数')
    fig.suptitle('真实数据：590、n0、20w 的 B1 谱带\n'
                 '各 300 帧；沿用能量 / q 漂移校正，16 meV 分箱；未扣背景、未分峰、未插值补 q 点', fontsize=14)
    save(fig, output, '01_measured_qe_maps')

    # 2. Held-out real counts with the already frozen training-only prediction.
    fig, axes = plt.subplots(3, 3, figsize=(14.5, 10.4), layout='constrained')
    scores = []
    for row, sid in enumerate(SESSIONS):
        d, f = datasets[sid], fits[sid]
        e = d['energy']/1000
        band = (e >= .5)&(e <= 1.8)
        test = d['frames'][:, 250:].astype(float).sum(1)
        scale = d['beam'][250:].sum()/d['beam'][:250].sum()
        pred = np.asarray(f['prediction'])*scale
        parts = np.asarray(f['components'])*scale
        ap, order = sorted_apices(f)
        score = float(np.mean(((test-pred)/d['sigma_test'])[:, band]**2))
        scores.append(dict(session=sid, heldout_B1_mean_squared_standardized_residual=score))
        for col, qi in enumerate([4, 5, 7]):
            ax = axes[row, col]
            ax.errorbar(e[band], test[qi, band], yerr=d['sigma_test'][qi, band],
                        fmt='o', ms=2.5, color='#30343b', elinewidth=.6, alpha=.72,
                        label='实测：最后 50 帧')
            ax.plot(e[band], pred[qi, band], color='#24648c', lw=2, label='训练模型预测总谱')
            for rank, idx in enumerate(order):
                ax.plot(e[band], parts[qi, band, idx], color=COLORS[rank], lw=1.8,
                        label=['低能 B1 拟合贡献', '高能 B1 拟合贡献'][rank])
            ax.plot(e[band], parts[qi, band, 2:].sum(1), color='#898989', ls=':', lw=1.6,
                    label='邻近带与背景')
            ax.set(title=f'{sid} | q = {d["q"][qi]:.4f} Å⁻¹', xlim=(.5, 1.8),
                   ylim=(0, 1.16*np.max(test[qi, band]+d['sigma_test'][qi, band])))
            if col == 0:
                ax.set_ylabel('计数 / 16 meV / 50 帧')
            if row == 2:
                ax.set_xlabel('能量损失 / eV')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=5, fontsize=10)
    fig.suptitle('真实留出谱上的双分量拟合：黑点来自实验，彩线来自已有模型\n'
                 '前 250 帧拟合，后 50 帧按束流比例预测；两条 B1 分量均已包含运动学、接收积分与能量展宽', fontsize=14)
    save(fig, output, '02_measured_spectra_components')

    # 3. Separate measured maxima from model-dependent intrinsic loss apices.
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.4), layout='constrained')
    ridge_rows = []
    for ax, sid in zip(axes, SESSIONS):
        d, f = datasets[sid], fits[sid]
        e = d['energy']/1000
        band = (e >= .5)&(e <= 1.8)
        positive = np.flatnonzero(d['q'] > 0)
        qs = d['q'][positive]
        observed = gaussian_filter1d(d['frames'].astype(float).sum(1), 2, axis=-1)
        observed_peak = e[band][observed[positive][:, band].argmax(1)]
        prediction = np.asarray(f['prediction'])
        predicted_peak = e[band][prediction[positive][:, band].argmax(1)]
        p = f['parameters']
        denseq = np.linspace(.0025, .015, 301)
        x = denseq/.005
        e0 = p[0]*np.sqrt(x/(1+p[1]*(x-1)))
        g = p[2]*x**p[3]
        intrinsic = np.stack([apex(e0, g), np.full_like(denseq, apex(p[4], p[5]))])
        _, order = sorted_apices(f)
        for rank, idx in enumerate(order):
            # Dashed = energy/width constant by model design, not discovered.
            ax.plot(denseq, intrinsic[idx], color=COLORS[rank], lw=2,
                    ls='--' if idx == 1 else '-',
                    label=['低能裸损失峰顶（拟合）', '高能裸损失峰顶（拟合）'][rank])
        ax.plot(qs, predicted_peak, 's-', color='#24648c', lw=1.3, ms=5,
                label='模型预测的计数谱顶')
        ax.plot(qs, observed_peak, 'ko', ms=5, label='实测计数谱顶')
        ax.set(title=sid, xlabel=r'$q_x$ / Å$^{-1}$', ylabel='能量 / eV', ylim=(.45, 1.8))
        ax.set_xticks([.0025, .005, .01, .015], ['0.0025', '0.005', '0.010', '0.015'])
        ax.grid(alpha=.15)
        if any(f['boundary']):
            ax.text(.03, .96, '屏蔽参数触及上界', transform=ax.transAxes, va='top', fontsize=9,
                    color='#7c5627')
        for qi, qm, obs, pred in zip(positive, qs, observed_peak, predicted_peak):
            xx = qm/.005
            ap = np.array([apex(p[0]*np.sqrt(xx/(1+p[1]*(xx-1))), p[2]*xx**p[3]),
                           apex(p[4], p[5])])
            ridge_rows.append(dict(session=sid, q_Ainv=qm, measured_smoothed_peak_eV=obs,
                predicted_count_peak_eV=pred, fitted_low_loss_apex_eV=ap[order[0]],
                fitted_high_loss_apex_eV=ap[order[1]],
                fixed_energy_component_rank=int(np.flatnonzero(order == 1)[0]+1)))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4, fontsize=10)
    fig.suptitle('真实谱顶与条件分解：上升谱脊不必等于一条模式的本征色散\n'
                 '彩色裸峰顶由既有拟合反推；虚线分量的能量与宽度被设为常数；黑点仅为取峰作 32 meV 高斯平滑', fontsize=13)
    save(fig, output, '03_measured_ridge_and_fitted_components')

    with (output/'ridge_values.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(ridge_rows[0])); w.writeheader(); w.writerows(ridge_rows)
    source_files = [Path(__file__), Path(__file__).with_name('REAL_DATA_FIGURES.md'),
        Path(__file__).with_name('models.py'), Path(__file__).with_name('summarize.py'),
        Path(__file__).with_name('io_utils.py'), ROOT/'src/python/qeels_frames.py',
        ROOT/'src/python/qeels_rebin.py', composition/'manifest.json']
    inputs += [dict(path=str(p), sha256=digest(p)) for p in source_files]
    dump(output/'manifest.json', dict(kind='real measured data and frozen conditional fit figures',
        checks=checks, scores=scores, input_files=inputs,
        map_processing='read-only raw frames, frozen E/q drift, measured q bins, 16 meV sums, column display normalization',
        model_limit='One B1 component has fixed energy/width; the other allows screened dispersion. No new fit or direct deconvolution.',
        display_only_ridge_smoothing_sigma_meV=32, python=sys.version, numpy=np.__version__))
    (output/'READ_ME.md').write_text(Path(__file__).with_name('REAL_DATA_FIGURES.md').read_text())
    with zipfile.ZipFile(output/'b1_real_data_figures.zip', 'w', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(output.iterdir()):
            if p.suffix in {'.png', '.pdf', '.svg', '.md', '.csv', '.json'}:
                z.write(p, p.name)
    print('DONE', output, json.dumps(scores), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, default=ROOT/'paper_results/qeels_exploration_20261023')
    parser.add_argument('--composition', type=Path, default=ROOT/'paper_results/b1_composition_20261023')
    parser.add_argument('--output', type=Path, default=ROOT/'paper_results/b1_real_data_figures_20261023')
    args = parser.parse_args()
    run(args.previous, args.composition, args.output)
