"""Independent quadrature checks and peak-definition comparison for the audit."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from kinematic_audit import read, save, csv_out
from kinematic_ablation import AblationOperator, moving_mean
from models import acceptance_nodes, kinematic_2d
from dispersion_models import columns, curves, loss_apex
from io_utils import digest, dump


def run(source, baseline, output):
    font = Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size': 11, 'pdf.fonttype': 42})
    signature = read(output/'signature.json')
    for name, sha in signature.items():
        if name.endswith('.py'):
            assert digest(Path(__file__).with_name(name)) == sha
    checks, rows = [], []
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout='constrained')
    for ax, sid in zip(axes, ['590', 'n0', '20w']):
        assert digest(source/f'{sid}.npz') == signature[sid]
        d = dict(np.load(source/f'{sid}.npz'))
        for mode in ['none', 'energy_flat']:
            op = AblationOperator(d, mode)
            fine = AblationOperator(d, mode, n_radial=512, n_quad=45, n_probe=11)
            for model in ['A', 'D']:
                f = read(output/sid/mode/f'{model}.json')
                p = np.array(f['parameters']); c = np.array(f['coefficients'])
                col, _ = columns(op, p, model)
                pred = np.einsum('bek,bk->be', col, c)
                np.testing.assert_allclose(pred, f['prediction'], atol=1e-8, rtol=1e-12)
                mask = d['energy'] >= 500
                residual = ((pred-d['frames'][:, :250].sum(1))/d['sigma_train'])[:, mask]
                np.testing.assert_allclose((residual**2).sum(), f['cost'], rtol=1e-10)
                fine_col, _ = columns(fine, p, model)
                fine_pred = np.einsum('bek,bk->be', fine_col, c)
                rms = float(np.sqrt(np.mean(((fine_pred-pred)/d['sigma_train'])[:, mask]**2)))
                assert rms < .03, (sid, mode, model, rms)
                assert f['success']
                best = min(x['cost'] for x in f['candidates'])
                check = dict(session=sid, mode=mode, model=model, fine_operator_rms_sigma=rms,
                    converged_starts=sum(x['success'] for x in f['candidates']),
                    starts=len(f['candidates']), boundary=f['boundary'],
                    near_best_starts=sum(x['cost']-best < 1e-3 for x in f['candidates']))
                if model == 'D':
                    a = read(output/sid/mode/'A.json')
                    assert f['cost'] <= a['cost']+1e-7
                checks.append(check)
        f = read(baseline/sid/'main/D.json')
        e = d['energy']/1000; band = (e >= .5)&(e <= 1.8)
        qi = np.flatnonzero(d['q'] > 0); q = d['q'][qi]
        comp = np.array(f['components'])
        signal = d['frames'][:, :250].sum(1)-comp[:, :, 2:].sum(-1)
        # Six 16-meV channels approximate the old 100-meV plotting window;
        # use seven (112 meV) so the centered window has no half-bin shift.
        observed = np.array([moving_mean(y, 7) for y in signal])
        ridge = e[band][observed[qi][:, band].argmax(1)]
        pred_b1 = comp[:, :, :2].sum(-1)
        predicted = e[band][pred_b1[qi][:, band].argmax(1)]
        low_count = e[band][comp[qi, :, 0][:, band].argmax(1)]
        el, eh, gl, gh = curves(np.array(f['parameters']), 'D', q)
        apex = loss_apex(el, gl)
        for y, label, style, color in [(ridge, '实测总谱顶（扣当前拟合背景）', 'o-', '#333333'),
                (predicted, '模型两部分相加后的计数谱顶', 's--', '#9b74b1'),
                (el, '低能振子参数 E₀', '-', '#227e96'),
                (apex, '低能裸损失峰顶', '--', '#008b74'),
                (low_count, '低能分量的计数峰顶', ':', '#c7833b')]:
            ax.plot(q, y, style, color=color, lw=1.8, ms=5, label=label)
        for i in range(len(q)):
            rows.append(dict(session=sid, q_Ainv=q[i], measured_background_subtracted_peak_eV=ridge[i],
                predicted_B1_count_peak_eV=predicted[i], low_E0_eV=el[i],
                low_loss_apex_eV=apex[i], low_count_peak_eV=low_count[i]))
        ax.set(title=sid, xlabel='正 q / Å⁻¹', ylabel='能量 / eV', ylim=(.45, 1.65))
        ax.grid(alpha=.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=3, fontsize=10)
    fig.suptitle('同一个拟合中的不同“峰位”：总谱顶、分量峰顶、E₀ 必须分开\n'
                 '前 250 帧；完整 K 模型 D；谱顶仅在 0.5–1.8 eV 内寻找；背景含邻近谱带', fontsize=13)
    save(fig, output, '04_peak_definitions')
    csv_out(output/'peak_definitions.csv', rows)
    # Independent quadrature convergence for the corrected historical map.
    wide = []
    for center in [.0005, .001, .0025, .005]:
        ks = []
        for n, p in [(45, 11), (65, 15)]:
            q, w = acceptance_nodes(center, dq=.0005, h=.0066,
                sigma_probe=.00041245201526075604, n=n, n_probe=p)
            k = kinematic_2d(q[None, :], np.array([300., 700., 1200., 1500., 1800.])[:, None], 30.)@w
            ks.append(k/k[3])
        err = float(np.max(abs(ks[0]/ks[1]-1)))
        assert err < .005
        wide.append(dict(q_Ainv=center, relative_K_shape_error=err))
    for row in read(output/'historical_inputs.json'):
        assert digest(row['path']) == row['sha256']
    report = dict(fits=checks, wide_aperture_convergence=wide,
        prepared_data_and_historical_inputs_hash_unchanged=True,
        code_sha256=digest(Path(__file__)))
    dump(output/'verification.json', report)
    print('VERIFIED', len(checks), 'fits; max fine operator RMS / sigma',
          max(x['fine_operator_rms_sigma'] for x in checks), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for key in ['source', 'baseline', 'output']:
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    run(a.source, a.baseline, a.output)
