"""Read-only summaries and numerical checks of the requested 3/5-q-bin run."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from qbin_run import read
from qbin_fit import local_columns
from models import Operator
from dispersion_models import columns as dispersion_columns
from io_utils import digest, dump
from kinematic_audit import csv_out, save


def summarize(output, baseline):
    font = Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    rows, changes, checks, mix = [], [], [], []
    data = {}; fits = {}
    for sid in ['590', 'n0', '20w']:
        for n in [3, 5]:
            d = dict(np.load(output/f'bin{n}/{sid}.npz')); data[sid, n] = d
            record = next(s for s in read(output/f'bin{n}/manifest.json')['sessions'] if s['session']==sid)
            assert digest(output/f'bin{n}/{sid}.npz') == record['prepared_sha256']
            band = (d['energy'] >= 500)&(d['energy'] <= 1800)
            for qi, center in enumerate(d['q']):
                local = {**d, 'q':d['q'][qi:qi+1]}
                op = Operator(local)
                for model in ['S', 'A']:
                    f = read(output/f'bin{n}/fits/{sid}/{center:+.4f}/{model}.json')
                    fits[sid, n, qi, model] = f
                    col, _ = local_columns(op, np.array(f['parameters']), model)
                    pred = np.einsum('bek,bk->be', col, np.array(f['coefficients']))
                    np.testing.assert_allclose(pred, f['prediction'], atol=1e-8, rtol=1e-12)
                    r = (pred-d['frames'][qi:qi+1,:250].sum(1))/d['sigma_train'][qi:qi+1]
                    np.testing.assert_allclose(np.sum(r[:,d['energy']>=500]**2), f['cost'], rtol=1e-10)
                    if model == 'A':
                        assert f['cost'] <= fits[sid,n,qi,'S']['cost']+1e-3
                    assert f['success']
                    best = min(c['cost'] for c in f['candidates'])
                    near = [c for c in f['candidates'] if c['cost']-best < 1e-3]
                    checks.append(dict(session=sid, channels=n, q=float(center), model=model,
                        starts=len(f['candidates']), converged_starts=sum(c['success'] for c in f['candidates']),
                        near_best_starts=len(near), boundary=f['boundary']))
                    comp = np.array(f['components'])[0]
                    frac = (comp[band,0].sum()/comp[band,:2].sum()) if model=='A' else 1.
                    y = d['frames'][qi].astype(float).sum(0)
                    snr = float(np.median(y[band]/d['sigma_all'][qi,band]))
                    rows.append(dict(session=sid, channels=n, q_Ainv=float(center),
                        width_Ainv=float(d['width']), model=model, E0_low_meV=f['E0_low_eV']*1000,
                        E0_high_meV=f['E0_high_eV']*1000 if model=='A' else '',
                        gamma_low_meV=f['gamma_low_eV']*1000, low_apex_meV=f['low_apex_eV']*1000,
                        low_B1_area_fraction=float(frac), median_count_snr=snr,
                        test_b1_msr=f['test_b1_msr'], boundary=';'.join(f['boundary'])))
                # A fine-operator check at +q=.01, bin5 for each session/model.
                if n==5 and abs(center-.01)<1e-12:
                    fine = Operator(local,n_radial=512,n_quad=45,n_probe=11)
                    for model in ['S','A']:
                        f = fits[sid,n,qi,model]
                        col,_ = local_columns(fine,np.array(f['parameters']),model)
                        pred=np.einsum('bek,bk->be',col,np.array(f['coefficients']))
                        rms=float(np.sqrt(np.mean(((pred-np.array(f['prediction']))/
                            d['sigma_train'][qi:qi+1])[:,band]**2)))
                        assert rms < .01
                        checks.append(dict(session=sid,channels=n,q=float(center),model=model,
                                           fine_operator_rms_sigma=rms))
        for qi, center in enumerate(data[sid,3]['q']):
            a=next(r for r in rows if r['session']==sid and r['channels']==3 and r['q_Ainv']==center and r['model']=='A')
            b=next(r for r in rows if r['session']==sid and r['channels']==5 and r['q_Ainv']==center and r['model']=='A')
            changes.append(dict(session=sid,q_Ainv=float(center),snr_ratio_5_over_3=b['median_count_snr']/a['median_count_snr'],
                E0_low_change_meV=b['E0_low_meV']-a['E0_low_meV'],
                gamma_low_change_meV=b['gamma_low_meV']-a['gamma_low_meV'],
                low_fraction_3=a['low_B1_area_fraction'],low_fraction_5=b['low_B1_area_fraction'],
                improvement_3=1-fits[sid,3,qi,'A']['test_b1_msr']/fits[sid,3,qi,'S']['test_b1_msr'],
                improvement_5=1-fits[sid,5,qi,'A']['test_b1_msr']/fits[sid,5,qi,'S']['test_b1_msr']))
        # Conditional estimate of width-induced mixing using the frozen global
        # dispersive fit, without refitting it or fixing local energies to it.
        f=read(baseline/sid/'main/D.json'); centers=[]
        for n in [3,5]:
            d=data[sid,n]; op=Operator(d)
            col,_=dispersion_columns(op,np.array(f['parameters']),'D')
            centers.append(np.sum(col[:,band,:2]*d['energy'][None,band,None],axis=1)/col[:,band,:2].sum(1))
        for qi,q in enumerate(d['q']):
            mix.append(dict(session=sid,q_Ainv=float(q),
                low_centroid_change_meV=float(centers[1][qi,0]-centers[0][qi,0]),
                high_centroid_change_meV=float(centers[1][qi,1]-centers[0][qi,1])))
    csv_out(output/'fit_summary.csv',rows); csv_out(output/'bin_changes.csv',changes)
    csv_out(output/'conditional_dispersion_mixing.csv',mix)
    dump(output/'verification.json',dict(checks=checks,summary_code_sha256=digest(Path(__file__))))

    fig,axes=plt.subplots(3,4,figsize=(17,10),layout='constrained')
    for row,sid in enumerate(['590','n0','20w']):
        for qi,ax in enumerate(axes[row]):
            for n,color in [(3,'#ba6847'),(5,'#237895')]:
                d=data[sid,n]; e=d['energy']/1000; band=(e>=.5)&(e<=1.8)
                y=d['frames'][qi].astype(float).sum(0)/n/1000
                err=d['sigma_all'][qi]/n/1000
                ax.plot(e[band],y[band],color=color,lw=1.2,label=f'{n} q 通道')
                ax.fill_between(e[band],(y-err)[band],(y+err)[band],color=color,alpha=.15)
            change=next(x for x in changes if x['session']==sid and x['q_Ainv']==d['q'][qi])
            ax.set(title=f'{sid}  q={d["q"][qi]:+.3f} | SNR ×{change["snr_ratio_5_over_3"]:.2f}',
                   xlabel='能量损失 / eV',ylabel='千计数 / q 通道 / 16 meV',xlim=(.5,1.8))
            ax.grid(alpha=.12)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=2)
    fig.suptitle('高 q 真实单谱：仅比较 3、5 个相邻 q 通道\n'
                 '300 帧；先校正漂移再加计数；除以通道数便于比较；阴影为计数噪声估计；未扣背景、未平滑',fontsize=14)
    save(fig,output,'01_measured_3_vs_5')

    for sign,name in [(1,'positive'),(-1,'negative')]:
        fig,axes=plt.subplots(3,2,figsize=(12,10),layout='constrained')
        for row,sid in enumerate(['590','n0','20w']):
            d=data[sid,5]; e=d['energy']/1000; band=(e>=.5)&(e<=1.8)
            indices=np.flatnonzero(np.sign(d['q'])==sign)
            for ax,qi in zip(axes[row],indices):
                f=fits[sid,5,qi,'A']; single=fits[sid,5,qi,'S']; scale=f['test_scale']/5
                test=d['frames'][qi,250:].sum(0)/5; err=d['sigma_test'][qi]/5
                comp=np.array(f['components'])[0]*scale
                ax.errorbar(e[band],test[band],yerr=err[band],fmt='.',ms=3,color='#333333',
                            elinewidth=.6,alpha=.7,label='留出 50 帧实测')
                ax.plot(e[band],np.array(single['prediction'])[0,band]*scale,'--',color='#969696',label='单 B1 成分')
                ax.plot(e[band],np.array(f['prediction'])[0,band]*scale,color='#237895',lw=1.8,label='双 B1 总谱')
                for j,color,label in [(0,'#c47b3c','低能 B1'),(1,'#258b71','高能 B1')]:
                    ax.plot(e[band],comp[band,j],color=color,lw=1.5,label=label)
                ax.plot(e[band],comp[band,2:].sum(1),':',color='#777777',label='邻近带与背景')
                ax.set(title=f'{sid}  q={d["q"][qi]:+.3f}',xlabel='能量损失 / eV',
                       ylabel='计数 / q 通道 / 16 meV',xlim=(.5,1.8))
        handles,labels=axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=10)
        fig.suptitle('5 q 通道：独立单谱拟合对未参与训练的最后 50 帧的预测\n'
                     '前 250 帧拟合；含运动学与接收积分；低、高成分身份仍是经验分解',fontsize=13)
        save(fig,output,f'02_bin5_heldout_{name}')
    fig,axes=plt.subplots(1,3,figsize=(13,4.7),layout='constrained')
    for ax,sid in zip(axes,['590','n0','20w']):
        for n,color in [(3,'#ba6847'),(5,'#237895')]:
            rr=[r for r in rows if r['session']==sid and r['channels']==n and r['model']=='A']
            ax.plot(np.arange(4),[r['E0_low_meV'] for r in rr],'o',color=color,ms=7,label=f'{n} q 通道')
        ax.set_xticks(np.arange(4),['−.015','−.010','+.010','+.015'])
        ax.set(title=sid,xlabel='q / Å⁻¹',ylabel='局部拟合低能 E₀ / meV')
        ax.grid(alpha=.15)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=2)
    fig.suptitle('改变 q-bin 后的低能分量位置：各点独立拟合，未施加跨 q 色散约束\n'
                 '这里只比较最优值；未给出 bootstrap 参数置信区间',fontsize=13)
    save(fig,output,'03_local_low_energy_stability')
    print('Summary complete. Fits',len(rows),'all raw spectra and heldout comparisons saved.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    a=p.parse_args();summarize(a.output,a.baseline)
