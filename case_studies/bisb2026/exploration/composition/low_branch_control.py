"""Fixed-response projection control, with parameters taken from one real fit.

This is a conditional shape demonstration, not a fit of all spectra to a flat
mode. E0 and Gamma come from the +q=0.0025 fit and are never tuned to match an
observed rise. The full aperture and q-binning pattern still enter K(Q,E).
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from models import Operator,oscillator
from low_branch_run import read,count_apex
from kinematic_audit import csv_out,save
from io_utils import digest,dump


def run(output):
    font=Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    mappings=read(output/'mappings.json');rows=[];metrics=[]
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
    for sid,ax in zip(['590','n0','20w'],axes):
        d=dict(np.load(output/f'bin5/{sid}.npz'))
        group=sorted([m for m in mappings if m['session']==sid and m['high_channels']==5 and m['q_Ainv']>0],
                     key=lambda m:m['q_Ainv'])
        ref=min(group,key=lambda m:abs(m['q_Ainv']-.0025))
        f=read(output/'fits'/ref['key']/'main.json')['A']
        e0=f['E0_low_eV'];gamma=f['gamma_low_eV'];local=[]
        for m in group:
            qi=abs(d['q']-m['q_Ainv']).argmin()
            op=Operator({**d,'q':d['q'][qi:qi+1],'width':float(d['width'][qi])})
            curve=op.project(op.KQ*oscillator(op.e,e0,gamma))[0]
            apex,edge=count_apex(op.energy,curve)
            actual=read(output/'fits'/m['key']/'main.json')['A']['low_count_apex_eV']
            r=dict(session=sid,q_Ainv=m['q_Ainv'],fixed_E0_eV=e0,fixed_gamma_eV=gamma,
                   fixed_response_count_apex_eV=apex,actual_local_low_count_apex_eV=actual)
            rows.append(r);local.append(r)
            if m['key']==ref['key']:np.testing.assert_allclose(apex,f['low_count_apex_eV'],atol=1e-10)
        a=min(local,key=lambda r:abs(r['q_Ainv']-.001));b=min(local,key=lambda r:abs(r['q_Ainv']-.005))
        metrics.append(dict(session=sid,fixed_E0_eV=e0,fixed_gamma_eV=gamma,
            actual_rise_eV=b['actual_local_low_count_apex_eV']-a['actual_local_low_count_apex_eV'],
            fixed_response_rise_eV=b['fixed_response_count_apex_eV']-a['fixed_response_count_apex_eV']))
        q=[r['q_Ainv'] for r in local]
        ax.plot(q,[r['fixed_response_count_apex_eV'] for r in local],color='#247b99',lw=2,label='固定响应：只改变测量权重')
        ax.plot(q,[r['actual_local_low_count_apex_eV'] for r in local],'o',color='#be7535',ms=3.5,
                label='真实数据逐谱分解的最佳值')
        ax.axvspan(.006,.02,color='#888888',alpha=.08)
        ax.set(title=f'{sid} | 固定 E₀={e0:.3f}，Γ={gamma:.3f} eV',
               xlabel='正 q / Å⁻¹',ylabel='低能成分计数峰顶 / eV',ylim=(.3,1.3),xlim=(0,.020))
        ax.grid(alpha=.15)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=2)
    fig.suptitle('固定 E₀ 和 Γ 的响应，也能产生上升的计数峰顶\n'
                 '参数直接取自 q=0.0025 的真实拟合；此图仅为投影对照，未把所有谱重新拟合成平坦模式',fontsize=13)
    save(fig,output,'05_fixed_response_control')
    csv_out(output/'fixed_response_control.csv',rows)
    dump(output/'fixed_response_control.json',dict(metrics=metrics,code_sha256=digest(Path(__file__))))
    print(metrics)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output)
