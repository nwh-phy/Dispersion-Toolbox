"""Run the small, reproducible theory experiment; original acquisitions are not read."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.signal import find_peaks

from models import (ROOT, bilayer_loss, bilayer_modes, bilayer_residues,
                    density_mass_combination, hybrid_components, hybrid_from_targets,
                    hybrid_loss, hybrid_modes, minimum_distance)
from qeels_response import kinematic_2d, longitudinal_q

spec = importlib.util.spec_from_file_location('composition_forward_operator',
    Path(__file__).resolve().parent.parent/'composition/models.py')
forward = importlib.util.module_from_spec(spec)
spec.loader.exec_module(forward)


def save_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def save_csv(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def save_fig(fig, out, name):
    for ax in fig.axes:
        ax.title.set_fontsize(11)
        if ax.get_xscale() == 'linear':
            ax.xaxis.set_major_locator(MaxNLocator(5))
    fig.savefig(out/f'{name}.png', dpi=180)
    fig.savefig(out/f'{name}.pdf')
    plt.close(fig)


def projected_components(op, alpha_d, alpha_l, delta, r0, gamma):
    loss = hybrid_components(op.e, op.radial, alpha_d, alpha_l, delta, r0, gamma)
    # Reuse the existing full-Q acceptance operator and Gaussian energy kernel.
    return np.stack([op.project(np.einsum('bkr,kr->bk', op.H, loss[:, :, j]))
                     for j in range(2)], axis=-1)


def integrate_low_fraction(components, energy_meV):
    mask = (energy_meV >= 500)&(energy_meV <= 1800)
    areas = components[:, mask].sum(1)
    return areas[:, 0]/areas.sum(1)


def run(config_path, composition, out, standalone=False):
    out.mkdir(parents=True, exist_ok=True)
    cfg = tomllib.loads(config_path.read_text())
    targets, input_paths = [], [config_path]
    for sid in ([] if standalone else ['590', 'n0', '20w']):
        path = composition/'fits_v3'/sid/'main_double.json'
        f = json.loads(path.read_text())
        energy = sorted([f['parameters'][0], f['parameters'][4]])
        targets.append(dict(session=sid, low_E0_eV=energy[0], high_E0_eV=energy[1],
                            minimum_distance_A=float(minimum_distance(cfg['q_reference_Ainv'], *energy))))
        input_paths.append(path)
    if standalone:
        ref = cfg['standalone']
        lo, hi = ref['low_target_eV'], ref['high_target_eV']
        targets.append(dict(session='frozen_reference', low_E0_eV=lo, high_E0_eV=hi,
                            minimum_distance_A=float(minimum_distance(cfg['q_reference_Ainv'], lo, hi))))
    low, high = np.mean([[r['low_E0_eV'], r['high_E0_eV']] for r in targets], axis=0)
    qref = cfg['q_reference_Ainv']
    dmin = float(minimum_distance(qref, low, high))
    alpha = (low**2+high**2)/(2*qref)
    beam = cfg['beam_kV']
    ep, residue = bilayer_residues(qref, alpha, alpha, dmin, beam)
    area = residue/ep*kinematic_2d(qref, ep*1000, beam)
    dark_fraction = float(area[0]/area.sum())
    q = np.linspace(.0001, .02, 240)
    energy = np.linspace(.03, 2.8, 1000)
    distance_rows = []
    for d in cfg['distance_scan_A']+[dmin]:
        modes, weight = bilayer_residues(qref, alpha, alpha, d, beam)
        ar = weight/modes*kinematic_2d(qref, modes*1000, beam)
        distance_rows.append(dict(distance_A=d, low_eV=float(modes[0]), high_eV=float(modes[1]),
                                  low_pole_area_fraction=float(ar[0]/ar.sum())))
    save_csv(out/'bilayer_distance.csv', distance_rows)
    save_csv(out/'empirical_targets.csv', targets)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), layout='constrained')
    for d in [3, 30, dmin]:
        modes, _ = bilayer_modes(q, alpha, alpha, d)
        color = axes[0].plot(q, modes[:, 1], label=f'd={d/10:.1f} nm')[0].get_color()
        axes[0].plot(q, modes[:, 0], ls='--', color=color)
    axes[0].scatter([qref, qref], [low, high], c='black', marker='x', zorder=5)
    axes[0].set(xlabel='Q (Å$^{-1}$)', ylabel='Pole energy (eV)', ylim=(0, 2.3),
                title='Two Drude sheets; fixed total stiffness')
    axes[0].legend(fontsize=8)
    for gamma in [cfg['gamma_resolved_eV'], cfg['gamma_demo_eV']]:
        y = bilayer_loss(energy, [qref], alpha, alpha, dmin, gamma, gamma, beam)[:, 0]
        y *= kinematic_2d(qref, energy*1000, beam)
        axes[1].plot(energy, y/y.max(), label=f'gamma={gamma:g} eV')
    for E in [low, high]:
        axes[1].axvline(E, color='gray', ls=':', lw=.8)
    axes[1].set(xlabel='Energy (eV)', ylabel='Relative count response', xlim=(.3, 2),
                title=f'Normal incidence: weak lower mode\nd={dmin/10:.1f} nm, Q={qref:g} Å$^{{-1}}$')
    axes[1].legend(fontsize=8)
    ratios = np.geomspace(.03, 30, 300)
    for d in [3, 30, dmin, 300]:
        vals = []
        for r in ratios:
            modes, _ = bilayer_modes(qref, 2*alpha/(1+r), 2*alpha*r/(1+r), d)
            vals.append(modes[1]/modes[0])
        axes[2].plot(ratios, vals, label=f'd={d/10:.1f} nm')
    axes[2].axhline(high/low, color='black', ls=':', label='Empirical E0 ratio')
    axes[2].set(xscale='log', yscale='log', xlabel='Drude stiffness ratio', ylabel='E+ / E-',
                title='Unequal Drude weights\nincrease the energy ratio')
    axes[2].legend(fontsize=7)
    fig.suptitle('Conditional mechanism test; target E0 values are not measured plasmon eigenenergies', fontsize=10)
    save_fig(fig, out, 'bilayer_constraints')

    delta, r0, gamma = cfg['transition_demo_eV'], cfg['screening_demo_A'], cfg['gamma_demo_eV']
    ad, al = hybrid_from_targets(qref, low, high, delta, r0)
    qmap = np.linspace(.0005, .02, 110)
    emap = np.linspace(.1, 2.2, 701)
    poles, residues, _ = hybrid_modes(qmap, ad, al, delta, r0)
    loss = hybrid_loss(emap, qmap, ad, al, delta, r0, gamma, gamma)
    counts = loss*kinematic_2d(qmap[None, :], emap[:, None]*1000, beam)
    counts /= counts.max(0, keepdims=True)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), layout='constrained')
    axes[0].pcolormesh(qmap, emap, counts, cmap='magma', shading='auto', rasterized=True)
    axes[0].plot(qmap, poles[:, 0], 'c--', lw=1, label='Longitudinal poles')
    axes[0].plot(qmap, poles[:, 1], 'c--', lw=1)
    axes[0].axhline(delta, color='white', ls=':', lw=1, label='Bare transition Delta')
    axes[0].set(xlabel='Q (Å$^{-1}$)', ylabel='Energy (eV)', title='Drude + dynamic polarization: count map')
    axes[0].legend(fontsize=7)
    for g in [cfg['gamma_resolved_eV'], gamma]:
        components = hybrid_components(energy, [qref], ad, al, delta, r0, g)[:, 0]
        yy = components.sum(1)*kinematic_2d(qref, energy*1000, beam)
        axes[1].plot(energy, yy/yy.max(), label=f'gamma={g:g} eV')
    axes[1].set(xlabel='Energy (eV)', ylabel='Relative count response', xlim=(.3, 2),
                title='Same poles, different damping')
    axes[1].legend(fontsize=8)
    coupling = np.linspace(0, 1, 201)
    # Scale the static background polarization together with the dynamic ones.
    curves = np.array([hybrid_modes(qref, ad*x, al*x, delta, r0*x)[0] for x in coupling])
    axes[2].plot(coupling, curves[:, 0], label='Lower longitudinal mode')
    axes[2].plot(coupling, curves[:, 1], label='Upper longitudinal mode')
    axes[2].axhline(delta, color='gray', ls=':', label='Bare electronic transition')
    axes[2].set(xlabel='Coulomb coupling scale lambda', ylabel='Pole energy (eV)',
                title='Track the origin of both modes')
    axes[2].legend(fontsize=7)
    fig.suptitle(f'Illustrative: Delta={delta:g} eV, r0={r0/10:g} nm; two E0 targets matched by construction', fontsize=10)
    save_fig(fig, out, 'hybrid_mechanism')
    save_csv(out/'hybrid_dispersion.csv', [dict(q_Ainv=float(x), low_eV=float(p[0]), high_eV=float(p[1]),
        low_residue=float(s[0]), high_residue=float(s[1])) for x, p, s in zip(qmap, poles, residues)])

    # Full-Q spectra use the previously frozen 590 geometry; no new raw processing.
    if standalone:
        native = np.arange(200., 3400., 4.)
        qs = np.array(cfg['momentum_scan_Ainv'])
        data = dict(energy=native.reshape(-1, 4).mean(1), energy_native=native,
                    q=np.r_[-qs[::-1], qs], h=ref['slot_halfwidth_Ainv'], width=ref['bin_width_Ainv'],
                    probe=ref['probe_sigma_Ainv'], beam_kV=beam, resolution_fwhm=ref['energy_fwhm_meV'])
    else:
        data_path = composition/'590.npz'
        data = dict(np.load(data_path))
        input_paths.append(data_path)
    op = forward.Operator(data)
    cp = projected_components(op, ad, al, delta, r0, gamma)
    fraction = integrate_low_fraction(cp, op.energy)
    empirical = []
    if not standalone:
        area_path = composition/'summary/component_areas.csv'
        input_paths.append(area_path)
        empirical = list(csv.DictReader(area_path.open()))
    selected = [i for i, x in enumerate(data['q']) if x > 0]
    scan, fraction_rows = [], []
    for rr in cfg['screening_scan_A']:
        for fraction_delta in cfg['transition_fraction_scan']:
            dd = low+(high-low)*fraction_delta
            aa, cc = hybrid_from_targets(qref, low, high, dd, rr)
            for gg in cfg['gamma_scan_eV']:
                pp = projected_components(op, aa, cc, dd, rr, gg)
                ff = integrate_low_fraction(pp, op.energy)
                pp_sum = pp.sum(-1)
                mask = (op.energy >= 500)&(op.energy <= 1800)
                peaks, _ = find_peaks(pp_sum[5, mask], prominence=.03*pp_sum[5, mask].max())
                eends = hybrid_modes(np.array([.0025, .015]), aa, cc, dd, rr)[0]
                row = dict(r0_A=rr, transition_eV=dd, gamma_eV=gg, alpha_D_eV2_A=aa,
                           alpha_L_eV2_A=cc, density_over_kappa_mass_cm2=float(density_mass_combination(aa)),
                           count_peaks_q005=len(peaks), low_energy_change_eV=float(eends[1, 0]-eends[0, 0]),
                           high_energy_change_eV=float(eends[1, 1]-eends[0, 1]))
                for i in selected:
                    row[f'fraction_q{data["q"][i]:.4f}'] = float(ff[i])
                scan.append(row)
    save_csv(out/'hybrid_parameter_scan.csv', scan)
    print(f'Completed {len(scan)} hybrid parameter combinations', flush=True)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), layout='constrained')
    for i in selected:
        yy = cp[i].sum(-1)
        axes[0].plot(op.energy/1000, yy/yy.max(), label=f'q={data["q"][i]:g}')
        fraction_rows.append(dict(q_Ainv=float(data['q'][i]), low_fraction=float(fraction[i])))
    axes[0].set(xlabel='Energy (eV)', ylabel='Relative count response', xlim=(.4, 2),
                title='One model through the 590 aperture')
    axes[0].legend(fontsize=7)
    qs = data['q'][selected]
    for sid, marker in [('590', 'o'), ('n0', 's'), ('20w', '^')]:
        rr = [r for r in empirical if r['session'] == sid and float(r['q_Ainv']) > 0]
        if not rr:
            continue
        axes[1].plot([float(x['q_Ainv']) for x in rr], [float(x['low_fraction_of_B1']) for x in rr],
                     marker=marker, ls=':', alpha=.7, label=f'{sid}: empirical decomposition')
    axes[1].plot(qs, fraction[selected], 'k-', lw=2, label='Illustrative hybrid; 590 geometry')
    allfractions = np.array([[r[f'fraction_q{x:.4f}'] for x in qs] for r in scan])
    axes[1].fill_between(qs, allfractions.min(0), allfractions.max(0), color='gray', alpha=.15,
                         label='All scanned hypotheses (not CI)')
    axes[1].set(xlabel='Nominal q (Å$^{-1}$)', ylabel='Lower component / B1 components', ylim=(0, 1),
                title='Mode weights\nwithout independent q amplitudes')
    axes[1].legend(fontsize=6)
    for dd in [low+.2*(high-low), low+.5*(high-low), low+.8*(high-low)]:
        combos = [hybrid_from_targets(qref, low, high, dd, rr)[0] for rr in cfg['screening_scan_A']]
        axes[2].plot(np.array(cfg['screening_scan_A'])/10, density_mass_combination(combos)/1e14,
                     'o-', label=f'Delta={dd:.2f} eV')
    axes[2].set(xlabel='Background screening length r0 (nm)', ylabel='n / [kappa (m*/me)]  ($10^{14}$ cm$^{-2}$)',
                title='Required stiffness; not a carrier density')
    axes[2].legend(fontsize=7)
    fig.suptitle('Comparison to model-derived targets; these curves have not been fitted to raw spectra', fontsize=10)
    save_fig(fig, out, 'observable_constraints')
    save_csv(out/'hybrid_demo_fractions.csv', fraction_rows)

    fine = forward.Operator(data, n_radial=640, n_quad=45, n_probe=9)
    fine_cp = projected_components(fine, ad, al, delta, r0, gamma)
    mask = (op.energy >= 500)&(op.energy <= 1800)
    # Both normalizations are internal conditioning; compare shapes, not raw H norms.
    shape = cp.sum(-1); fine_shape = fine_cp.sum(-1)
    shape /= shape[:, mask].max(1)[:, None]
    fine_shape /= fine_shape[:, mask].max(1)[:, None]
    max_shape_error = float(np.max(abs(shape[:, mask]-fine_shape[:, mask])))
    max_fraction_error = float(np.max(abs(fraction-integrate_low_fraction(fine_cp, fine.energy))))
    if max_shape_error > .01 or max_fraction_error > .003:
        raise RuntimeError('Forward projection has not converged to the declared tolerance')
    save_json(out/'summary.json', dict(
        status='conditional theory exploration; no material prediction and no raw-data mechanism fit',
        standalone=standalone,
        empirical_E0_targets_eV=[float(low), float(high)], qref_Ainv=qref,
        minimum_drude_sheet_distance_A=dmin,
        bilayer_weak_damping_low_area_fraction_at_qref=dark_fraction,
        hbar_v_eV_A=float(1/longitudinal_q(1000, beam)),
        hybrid_demo=dict(Delta_eV=delta, r0_A=r0, gamma_eV=gamma,
                         alpha_D_eV2_A=float(ad), alpha_L_eV2_A=float(al),
                         density_over_kappa_mass_cm2=float(density_mass_combination(ad))),
        scan_count=len(scan), scan_decreasing_low_fraction_count=sum(
            r['fraction_q0.0150'] < r['fraction_q0.0025'] for r in scan),
        scan_two_count_peak_count=sum(r['count_peaks_q005'] == 2 for r in scan),
        forward_max_normalized_shape_error=max_shape_error,
        forward_max_fraction_error=max_fraction_error,
        demo_fractions=fraction_rows))
    input_paths += sorted(Path(__file__).resolve().parent.glob('*.py'))
    input_paths += [Path(forward.__file__), ROOT/'src/python/qeels_response.py']
    save_json(out/'manifest.json', dict(config=cfg, files=[dict(path=str(p.resolve()),
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in input_paths],
        numpy=np.__version__, python=sys.version))
    print(json.dumps(json.loads((out/'summary.json').read_text()), indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('config.toml'))
    parser.add_argument('--composition', type=Path, default=ROOT/'paper_results/b1_composition_20261023')
    parser.add_argument('--output', type=Path, default=ROOT/'paper_results/b1_two_mode_theory_20261023')
    parser.add_argument('--standalone', action='store_true', help='Use the frozen reference scales; no data files required')
    args = parser.parse_args()
    run(args.config, args.composition, args.output, args.standalone)
