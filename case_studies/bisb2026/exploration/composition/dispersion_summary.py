"""Plots and machine-readable tables from real-data ordered-model fits only."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from dispersion_models import curves, loss_apex
from io_utils import dump

COLORS = ['#008c77', '#ca5947']
LABELS = dict(A='A 两支恒定', B='B 低能色散', C='C 高能色散', D='D 两支色散', S='S 单成分')
CONDITIONS = dict(main='主设置', no_constant='去掉常数背景', half_h='接受角减半',
    lower350='下限 350 meV', upper3000='上限 3000 meV', fixed_width='线宽不随 q 变',
    ep2_numerator='改变谱权重 Q 依赖', holdout1100_1500='挖去 1.1–1.5 eV')


def read(path):
    return json.loads(path.read_text())


def save(fig, dest, name):
    fig.savefig(dest/f'{name}.png', dpi=180, bbox_inches='tight')
    fig.savefig(dest/f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def csv_table(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def summarize(source, output):
    font = Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(str(font)); plt.rcParams['font.family'] = 'Arial Unicode MS'
    plt.rcParams.update({'axes.spines.top':False, 'axes.spines.right':False,
        'font.size':11, 'axes.unicode_minus':False})
    dest = output/'summary'; dest.mkdir(exist_ok=True)
    sessions = ['590','n0','20w']; rows=[]; bootrows=[]; intervals={}
    main = {}; boot = {}; q = np.linspace(.0025,.015,51)
    for sid in sessions:
        main[sid] = {m:read(output/sid/'main'/f'{m}.json') for m in LABELS}
        for name in CONDITIONS:
            for m in LABELS:
                f = read(output/sid/name/f'{m}.json')
                p = np.array(f['parameters']); el,eh,gl,gh = curves(p,m,q)
                rows.append(dict(session=sid,condition=name,model=m,
                    delta_E0_low_meV=1000*f['delta_eV'][0],delta_E0_high_meV=1000*f['delta_eV'][1],
                    delta_apex_low_meV=1000*(loss_apex(el,gl)[-1]-loss_apex(el,gl)[0]),
                    delta_apex_high_meV=1000*(loss_apex(eh,gh)[-1]-loss_apex(eh,gh)[0]),
                    train_msr=f['train_msr'],test_b1_msr=f['test_b1_msr'],
                    test_energy_holdout_msr=f.get('test_energy_holdout_msr',float('nan')),
                    success=f['success'],boundary=';'.join(f['boundary'])))
        for m in ['B','D']:
            b = [read(p) for p in sorted((output/sid/'bootstrap').glob(f'*_{m}.json'))]
            if len(b) != 120:
                raise ValueError(f'Incomplete bootstrap {sid}/{m}: {len(b)}')
            boot[sid,m] = b
            for quantity in ['E0_low','E0_high','apex_low','apex_high']:
                key = quantity.replace('E0','energy')+'_eV'
                delta = 1000*np.array([x[key][-1]-x[key][0] for x in b])
                lo,med,hi = np.percentile(delta,[2.5,50,97.5])
                row = dict(session=sid,model=m,quantity=quantity,median_meV=med,
                    lower_meV=lo,upper_meV=hi,positive_fraction=float(np.mean(delta>0)),
                    success=sum(x['success'] for x in b),boundary=sum(bool(x['boundary']) for x in b))
                bootrows.append(row); intervals[sid,m,quantity]=row
    csv_table(dest/'model_scores.csv',rows);csv_table(dest/'bootstrap_intervals.csv',bootrows)
    paired=[]
    for sid in sessions:
        diff=np.array([d['test_b1_msr']-b['test_b1_msr'] for b,d in zip(boot[sid,'B'],boot[sid,'D'])])
        paired.append(dict(session=sid,D_minus_B_test_msr_percentiles=np.percentile(diff,[2.5,50,97.5]),
            D_better_fraction=float(np.mean(diff<0))))
    dump(dest/'paired_prediction.json',paired)

    fig,axes=plt.subplots(2,3,figsize=(13.5,8),sharex=True,sharey='row',layout='constrained')
    for j,sid in enumerate(sessions):
        f=main[sid]['D']; el,eh,gl,gh=curves(np.array(f['parameters']),'D',q)
        for i,(key,pairs) in enumerate([('energy',[(el,gl),(eh,gh)]),('apex',[(el,gl),(eh,gh)])]):
            ax=axes[i,j]
            for k,(branch,(e0,g)) in enumerate(zip(['low','high'],pairs)):
                samples=np.array([b[f'{key}_{branch}_eV'] for b in boot[sid,'D']])
                a,z=np.percentile(samples,[2.5,97.5],axis=0)
                value=e0 if i==0 else loss_apex(e0,g)
                ax.fill_between(q*1000,a,z,color=COLORS[k],alpha=.16)
                ax.plot(q*1000,value,color=COLORS[k],lw=2,label=['低能成分','高能成分'][k])
            ax.grid(alpha=.15); ax.set_title(f'{sid}  |  '+['振子能量 $E_0$','模型损失谱的峰顶'][i])
            if i==1:ax.set_xlabel(r'$|q|$ ($10^{-3}$ Å$^{-1}$)')
            if j==0:ax.set_ylabel('能量 (eV)')
    axes[0,0].legend(loc='upper left',fontsize=9)
    fig.suptitle('真实 B1 数据：四个 |q|、正负两侧的全局拟合\n实线为双色散模型；阴影为 120 次重新拟合的逐点 95% 范围',fontsize=14)
    save(fig,dest,'01_energy_and_apex')

    fig,axes=plt.subplots(1,3,figsize=(13.5,4.5),sharey=True,layout='constrained')
    for ax,sid in zip(axes,sessions):
        for k,m in enumerate(['B','D']):
            row=intervals[sid,m,'E0_low'];value=1000*main[sid][m]['delta_eV'][0]
            ax.plot([row['lower_meV'],row['upper_meV']],[k,k],lw=4,color=['#999999',COLORS[0]][k])
            ax.plot(value,k,'o',color=['#666666',COLORS[0]][k]);ax.plot(row['median_meV'],k,'|',color='black')
        ax.axvline(0,color='#333',lw=.8,ls='--');ax.set_title(sid)
        ax.set_yticks([0,1],['高能固定','高能也放开']);ax.set_ylim(-.5,1.5)
        ax.set_xlabel(r'低能 Δ$E_0$ (meV)');ax.grid(axis='x',alpha=.15)
    fig.suptitle('高能成分的约束会改变低能色散判断\n区间：帧块重采样 95%；圆点：原始数据；短线：重采样中位数',fontsize=14)
    save(fig,dest,'02_low_dispersion_intervals')

    fig,axes=plt.subplots(1,3,figsize=(13.5,4.5),sharey=True,layout='constrained')
    for ax,sid in zip(axes,sessions):
        vals=[main[sid][m]['test_b1_msr'] for m in LABELS]
        ax.bar(np.arange(5),vals,color=['#999','#008c77','#ca5947','#24648c','#cab485'])
        for k,v in enumerate(vals):ax.text(k,v+.012,f'{v:.3f}',ha='center',fontsize=10)
        ax.set_xticks(np.arange(5),['A\n恒定','B\n低色散','C\n高色散','D\n双色散','S\n单成分'])
        ax.set_title(sid);ax.set_ylim(0,1.8);ax.grid(axis='y',alpha=.15)
    axes[0].set_ylabel('留出 B1 标准化残差平方均值（越低越好）')
    fig.suptitle('前 250 帧拟合，预测最后 50 帧：两成分优势与色散优势要分开评估',fontsize=14)
    save(fig,dest,'03_heldout_scores')

    fig,axes=plt.subplots(3,3,figsize=(13.5,10),sharex=True,layout='constrained')
    for row,sid in enumerate(sessions):
        d=dict(np.load(source/f'{sid}.npz'));f=main[sid]['D'];ee=d['energy']/1000
        measured=d['frames'][:,250:].sum(1);components=np.array(f['components'])*f['test_scale']
        for col,qvalue in enumerate([.0025,.005,.015]):
            j=int(np.argmin(abs(d['q']-qvalue)));ax=axes[row,col]
            scale=measured[j,(ee>=.5)&(ee<=1.8)].max();y=measured[j]/scale
            ax.plot(ee,y,'o',ms=2,color='#222',label='留出实测计数')
            ax.plot(ee,components[j].sum(1)/scale,color='#24648c',lw=1.6,label='总预测')
            for k in range(2):ax.plot(ee,components[j,:,k]/scale,color=COLORS[k],lw=1.5,label=['低能分量','高能分量'][k])
            ax.plot(ee,components[j,:,2:].sum(1)/scale,color='#999',ls=':',lw=1,label='背景及邻近响应')
            ax.set_xlim(.5,1.8);ax.set_ylim(-.03,1.15)
            ax.set_title(f'{sid} | q = {qvalue:g} Å'+r'$^{-1}$',fontsize=11)
            if row==2:ax.set_xlabel('能量损失 (eV)')
            if col==0:ax.set_ylabel('计数 / 本谱 B1 最大值')
    axes[0,0].legend(fontsize=8,loc='upper right')
    fig.suptitle('真实留出谱与双色散分解：仅为显示按同一系数缩放，拟合使用原始计数',fontsize=14)
    save(fig,dest,'04_heldout_real_spectra')

    fig,axes=plt.subplots(1,3,figsize=(14,5.5),sharey=True,layout='constrained')
    for ax,sid in zip(axes,sessions):
        for k,m in enumerate(['B','D']):
            vals=[next(r['delta_E0_low_meV'] for r in rows if r['session']==sid and r['condition']==c and r['model']==m) for c in CONDITIONS]
            ax.plot(vals,np.arange(len(CONDITIONS))+(k-.5)*.15,'o',color=['#888',COLORS[0]][k],label=LABELS[m])
        ax.axvline(0,color='#555',lw=.8,ls='--');ax.set_title(sid)
        ax.set_yticks(np.arange(len(CONDITIONS)),list(CONDITIONS.values()));ax.invert_yaxis()
        ax.set_xlabel(r'低能 Δ$E_0$ (meV)');ax.grid(axis='x',alpha=.15)
    axes[0].legend(fontsize=9,loc='lower left')
    fig.suptitle('系统敏感性：每个条件都重新拟合；模型设定变化不包含在重采样阴影内',fontsize=14)
    save(fig,dest,'05_sensitivity')

    fig,axes=plt.subplots(1,3,figsize=(13.5,4.5),layout='constrained')
    decomposition=[]
    for ax,sid in zip(axes,sessions):
        b=boot[sid,'D'];delta=np.array([f['delta_eV'][0] for f in b])*1000
        eta=np.array([f['parameters'][5] for f in b])
        rho=float(np.corrcoef(delta,eta)[0,1])
        ax.scatter(delta,eta,s=18,alpha=.45,color=COLORS[0],edgecolors='none')
        ax.axvline(0,color='#555',ls='--',lw=.8)
        ax.set_title(f'{sid} | 相关系数 {rho:.2f}')
        ax.set_xlabel(r'低能 $\Delta E_0$ (meV)');ax.set_ylabel(r'线宽指数 $\eta_L$')
        for m in ['A','B','C','D']:
            f=main[sid][m];el,eh,gl,gh=curves(np.array(f['parameters']),m,np.array([.0025,.015]))
            a00=loss_apex(el[0],gl[0]);a10=loss_apex(el[1],gl[0])
            a01=loss_apex(el[0],gl[1]);a11=loss_apex(el[1],gl[1])
            decomposition.append(dict(session=sid,model=m,
                E0_change_meV=1000*(el[1]-el[0]),Gamma_change_meV=1000*(gl[1]-gl[0]),
                apex_change_meV=1000*(a11-a00),
                energy_path_contribution_meV=500*((a10-a00)+(a11-a01)),
                width_path_contribution_meV=500*((a01-a00)+(a11-a10)),
                D_bootstrap_delta_eta_correlation=rho))
    fig.suptitle('两种解释会相互补偿：振子能量上移与阻尼变窄\n每点是一次重新拟合；线宽按 Γ(q) ∝ q 的 η 次幂变化',fontsize=14)
    save(fig,dest,'06_energy_width_tradeoff')
    csv_table(dest/'apex_decomposition.csv',decomposition)
    dump(dest/'summary.json',dict(intervals=bootrows,paired_prediction=paired,
        real_data_sessions=sessions,model_scores=rows))


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();summarize(a.source,a.output)
