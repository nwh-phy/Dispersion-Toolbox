"""Controlled figures separating material modes from q-EELS visibility.

Uses the frozen illustrative parameters in config.toml, never raw acquisitions.
All spectra are computed before plotting; each displayed column is normalized
within the stated energy window. No ridges are hand drawn or fitted to images.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from models import ROOT, hybrid_components, hybrid_from_targets, hybrid_modes
from qeels_response import acceptance_nodes, kinematic_2d, longitudinal_q


def aperture_matrix(q_centers, radial, geometry, n_quad=31, n_probe=7):
    """Interpolate a radial count spectrum at the existing area-quadrature nodes.

    The two weights at neighboring radial grid points sum to the original
    area weight. There is no extra radial Jacobian in the Cartesian aperture.
    """
    rows = []
    for center in q_centers:
        q, w = acceptance_nodes(center, dq=geometry['bin_width_Ainv'],
            h=geometry['slot_halfwidth_Ainv'],
            sigma_probe=geometry['probe_sigma_Ainv'], n=n_quad, n_probe=n_probe)
        if q.min() < radial[0] or q.max() > radial[-1]:
            raise ValueError('Radial support does not cover the aperture')
        ix = np.searchsorted(radial, q, side='right')-1
        ix = np.minimum(ix, len(radial)-2)
        t = (q-radial[ix])/(radial[ix+1]-radial[ix])
        row = np.bincount(ix, weights=w*(1-t), minlength=len(radial))
        row += np.bincount(ix+1, weights=w*t, minlength=len(radial))
        rows.append(row)
    result = np.asarray(rows)
    np.testing.assert_allclose(result.sum(1), 1, atol=1e-13)
    return result


def single_loss(energy, q, e_ref, q_ref, gamma):
    ep2 = e_ref**2*np.asarray(q)/q_ref
    en = np.asarray(energy)[:, None]
    return ep2[None, :]*gamma*en/((en**2-ep2[None, :])**2+(gamma*en)**2)


def normalize(a, mask):
    return a/np.maximum(a[mask].max(axis=0, keepdims=True), 1e-100)


def instrument(counts, energy, fwhm_meV):
    sigma = fwhm_meV/1000/np.sqrt(8*np.log(2))/np.diff(energy).mean()
    return gaussian_filter1d(counts, sigma, axis=0, mode='constant', truncate=8)


def peak_counts(spectra, energy, mask):
    return [len(find_peaks(y, prominence=.03*y.max())[0]) for y in spectra[mask].T]


def save_figure(fig, output, name):
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(output/f'{name}.{extension}', dpi=190, facecolor='white')
    plt.close(fig)


def heatmap(ax, q, energy, spectrum, mask, poles, title, show_ridge=False):
    z = normalize(spectrum, mask)
    im = ax.pcolormesh(q, energy[mask], z[mask], shading='auto', cmap='magma',
                      vmin=0, vmax=1, rasterized=True)
    for j, p in enumerate(np.atleast_2d(poles)):
        ax.plot(q, p, '--', color='#76e5e9', lw=1.3,
                label='本征能量（无阻尼）' if j == 0 else None)
    if show_ridge:
        ax.plot(q, energy[mask][np.argmax(spectrum[mask], axis=0)],
                color='white', lw=1.4, label='计数 / 损失谱顶')
    ax.set(title=title, xlabel=r'$q_x$ / Å$^{-1}$', ylabel='能量 / eV',
           ylim=(energy[mask][0], energy[mask][-1]))
    ax.set_xticks([.0025, .005, .010, .015], ['0.0025', '0.005', '0.010', '0.015'])
    return im


def run(config_path, output):
    output.mkdir(parents=True, exist_ok=True)
    cfg = tomllib.loads(config_path.read_text())
    geometry = cfg['standalone']
    font = Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.labelsize': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'path', 'pdf.fonttype': 42})
    energy = np.arange(.004, 4.504, .004)
    q = np.linspace(.0025, .015, 101)
    radial = np.geomspace(1e-8, .04, 2400)
    mask = (energy >= .2)&(energy <= 2.2)
    beam, qref = cfg['beam_kV'], cfg['q_reference_Ainv']
    lo, hi = geometry['low_target_eV'], geometry['high_target_eV']
    delta, r0 = cfg['transition_demo_eV'], cfg['screening_demo_A']
    gamma, sharp = cfg['gamma_demo_eV'], cfg['gamma_resolved_eV']
    ad, al = hybrid_from_targets(qref, lo, hi, delta, r0)
    aperture = aperture_matrix(q, radial, geometry)
    kp = kinematic_2d(q[None, :], energy[:, None]*1000, beam)
    kr = kinematic_2d(radial[None, :], energy[:, None]*1000, beam)
    fwhm = geometry['energy_fwhm_meV']

    one = single_loss(energy, q, lo, qref, gamma)
    one_k = one*kp
    one_full = instrument((single_loss(energy, radial, lo, qref, gamma)*kr)@aperture.T,
                          energy, fwhm)
    single_pole = lo*np.sqrt(q/qref)
    single_counts = [peak_counts(x, energy, mask) for x in [one, one_k, one_full]]
    print('Single-mode maximum number of spectral maxima:', [max(x) for x in single_counts], flush=True)

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 5.0), layout='constrained')
    for ax, y, title in zip(axes, [one, one_k, one_full],
        ['A  材料损失函数 L(q, E)', 'B  只乘运动学因子 K(q, E)', 'C  再加入接收范围与能量展宽']):
        im = heatmap(ax, q, energy, y, mask, single_pole, title, show_ridge=True)
    axes[0].legend(fontsize=9, loc='upper left', facecolor='#222222', labelcolor='white')
    fig.colorbar(im, ax=axes, shrink=.86, label='每个 q 的谱内归一化强度')
    fig.suptitle('单模对照：测量过程未把这支色散模变成两条水平谱带\n'
                 rf'示意模型；$E_p\propto\sqrt{{q}}$，Γ = {gamma:g} eV；{beam:g} kV；能量分辨率 {fwhm:.1f} meV', fontsize=15)
    save_figure(fig, output, '01_single_mode_control')

    qwide = np.linspace(0, .015, 501)
    unscreened, _, _ = hybrid_modes(qwide, ad, al, delta, 0)
    screened, _, _ = hybrid_modes(qwide, ad, al, delta, r0)
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), layout='constrained')
    axes[0].plot(qwide, np.sqrt(ad*qwide), color='#277da8', lw=2.5, label='自由载流子 plasmon')
    axes[0].axhline(delta, color='#8b5b2b', ls=':', label='待耦合的电子共振 Δ')
    for ax, energies, title in zip(axes[1:], [unscreened, screened],
        ['B  加入动态电子极化：两支混合模', 'C  再加入非局域屏蔽：两支趋平']):
        ax.plot(qwide, energies[:, 0], color='#00897b', lw=2.5, label='低能纵向模')
        ax.plot(qwide, energies[:, 1], color='#c45b44', lw=2.5, label='高能纵向模')
        ax.axhline(delta, color='#8b5b2b', ls=':', label='裸电子共振 Δ')
        ax.set_title(title)
    axes[0].set_title('A  从一支 √q 色散开始')
    for ax in axes:
        ax.set(xlabel=r'$Q$ / Å$^{-1}$', ylabel='本征能量 / eV', xlim=(0, .015), ylim=(0, 4.4))
        ax.set_xticks([0, .005, .01, .015])
        ax.axvspan(.0025, .015, color='#e6edf2', alpha=.5, zorder=0)
        ax.grid(alpha=.15)
        ax.legend(fontsize=9, loc='upper left')
    axes[2].annotate('B1 参考尺度\n约 0.8 与 1.4–1.5 eV', xy=(.009, 1.48), xytext=(.006, 2.4),
                     arrowprops={'arrowstyle': '->', 'color': '#555555'}, fontsize=10)
    fig.suptitle('双模与变平发生在材料响应中：这三幅图均未使用运动学因子\n'
                 f'αD = {ad:.1f}、αL = {al:.1f} eV² Å 保持不变；Δ = {delta:g} eV；C 中 r₀ = {r0/10:g} nm（示意参数）', fontsize=14)
    save_figure(fig, output, '02_material_origin')

    poles = hybrid_modes(q, ad, al, delta, r0)[0].T
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.2), layout='constrained')
    hybrids = {}
    for row, gg in enumerate([sharp, gamma]):
        loss = hybrid_components(energy, q, ad, al, delta, r0, gg).sum(-1)
        radial_parts = hybrid_components(energy, radial, ad, al, delta, r0, gg)
        parts = np.stack([instrument((radial_parts[:, :, j]*kr)@aperture.T, energy, fwhm)
                          for j in range(2)], axis=-1)
        full = parts.sum(-1)
        hybrids[gg] = dict(loss=loss, point=loss*kp, full=full, parts=parts)
        for col, (yy, title) in enumerate(zip([loss, loss*kp, full],
            ['材料损失函数', '× 运动学因子', '+ 接收积分与能量展宽'])):
            im = heatmap(axes[row, col], q, energy, yy, mask, poles,
                         f'{chr(65+3*row+col)}  {title}')
            axes[row, col].set_ylabel(f'能量 / eV\nΓ = {gg:g} eV')
    fig.colorbar(im, ax=axes, shrink=.9, label='每个 q 的谱内归一化强度')
    fig.suptitle('同一个双模响应：运动学改变可见强度，较大阻尼使两模重叠\n'
                 '青色虚线为同一组无阻尼模式能量；右列以标称 q 作参考；逐 q 归一化不代表绝对谱权重', fontsize=14)
    save_figure(fig, output, '03_two_mode_visibility')

    # Validate the new radial interpolation and quadrature, including a direct
    # node evaluation independent of the interpolation construction.
    check_q = np.array([.0025, .005, .01, .015])
    coarse = aperture_matrix(check_q, radial, geometry)
    fine_r = np.geomspace(1e-9, .04, 4800)
    fine = aperture_matrix(check_q, fine_r, geometry, n_quad=51, n_probe=9)
    fine_k = kinematic_2d(fine_r[None, :], energy[:, None]*1000, beam)
    errors = {}
    response_cases = {'single': lambda en, mom: single_loss(en, mom, lo, qref, gamma)}
    for gg in [sharp, gamma]:
        response_cases[f'hybrid_{gg}'] = lambda en, mom, gg=gg: hybrid_components(
            en, mom, ad, al, delta, r0, gg).sum(-1)
    for name, response in response_cases.items():
        a = instrument((response(energy, radial)*kr)@coarse.T, energy, fwhm)
        b = instrument((response(energy, fine_r)*fine_k)@fine.T, energy, fwhm)
        err = float(np.max(abs(normalize(a, mask)[mask]-normalize(b, mask)[mask])))
        errors[name] = err
        if err > .005:
            raise RuntimeError(f'{name} aperture refinement error {err:.5g} exceeds 0.5%')
    direct_q, direct_w = acceptance_nodes(.005, dq=geometry['bin_width_Ainv'],
        h=geometry['slot_halfwidth_Ainv'], sigma_probe=geometry['probe_sigma_Ainv'], n=31, n_probe=7)
    ee = np.linspace(.2, 2.2, 301)
    response = response_cases[f'hybrid_{sharp}']
    direct = (response(ee, direct_q)*kinematic_2d(direct_q[None, :], ee[:, None]*1000, beam))@direct_w
    compressed = (response(ee, radial)*kinematic_2d(radial[None, :], ee[:, None]*1000, beam))@coarse[1]
    errors['direct_vs_radial_at_q005'] = float(np.max(abs(direct-compressed))/direct.max())
    assert errors['direct_vs_radial_at_q005'] < .002
    # A constant spectrum is invariant under the area-average projection.
    np.testing.assert_allclose(np.ones((2, len(radial)))@coarse.T, 1, atol=1e-13)
    hv = float(1/longitudinal_q(1000, beam))
    rows = []
    for center in check_q:
        j = int(np.argmin(abs(q-center)))
        row = dict(q_Ainv=float(q[j]), single_pole_eV=float(single_pole[j]))
        for name, spectrum in [('loss', one), ('kinematic', one_k), ('aperture', one_full)]:
            row[f'single_{name}_peak_eV'] = float(energy[mask][np.argmax(spectrum[mask, j])])
        row['hybrid_low_pole_eV'], row['hybrid_high_pole_eV'] = map(float, poles[:, j])
        rows.append(row)
    with (output/'ridge_values.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    source_files = [Path(__file__), Path(__file__).with_name('models.py'), config_path,
                    Path(__file__).with_name('KINEMATIC_FIGURES.md'),
                    ROOT/'src/python/qeels_response.py']
    summary = dict(status='illustration and controlled mechanism test; not a fit to experimental data',
        q_range_Ainv=[float(q.min()), float(q.max())], energy_display_eV=[.2, 2.2],
        parameters=dict(alpha_D_eV2_A=ad, alpha_L_eV2_A=al, Delta_eV=delta, r0_A=r0,
                        gamma_eV=[sharp, gamma], beam_kV=beam, hbar_v_eV_A=hv, geometry=geometry),
        single_mode_peak_count_range=[dict(stage=name, minimum=min(x), maximum=max(x))
            for name, x in zip(['loss', 'kinematic', 'aperture'], single_counts)],
        displayed_peak_prominence_fraction=.03,
        refinement_errors=errors, selected_ridges=rows,
        files=[dict(path=str(p.resolve().relative_to(ROOT)) if p.resolve().is_relative_to(ROOT)
                    else str(p.resolve()), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
               for p in source_files], python=sys.version, numpy=np.__version__)
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')
    np.savez_compressed(output/'spectra.npz', q=q, energy_eV=energy,
        single_poles_eV=single_pole, hybrid_poles_eV=poles,
        single_loss=one, single_point_counts=one_k, single_full_counts=one_full,
        hybrid_sharp_loss=hybrids[sharp]['loss'], hybrid_sharp_counts=hybrids[sharp]['full'],
        hybrid_broad_loss=hybrids[gamma]['loss'], hybrid_broad_counts=hybrids[gamma]['full'])
    (output/'READ_ME.md').write_text(Path(__file__).with_name('KINEMATIC_FIGURES.md').read_text())
    print(json.dumps(dict(refinement_errors=errors, selected_ridges=rows), indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('config.toml'))
    parser.add_argument('--output', type=Path, default=ROOT/'paper_results/b1_kinematic_figures_20261023')
    args = parser.parse_args()
    run(args.config, args.output)
