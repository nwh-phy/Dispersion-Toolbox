"""Build audit tables and figures from checkpoints; never filter silently."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SESSIONS=['590','n0','20w']
VARIANTS=['one_peak','small_h','no_acceptance','no_kinematic','poisson_weights','unaligned','symmetric','upper_2100']


def read_csv(path):
    with Path(path).open() as f:
        rows=list(csv.DictReader(f))
    for r in rows:
        for k,v in r.items():
            if k in {'session','status','variant','nominal_status','variant_status'}:continue
            try:r[k]=float(v)
            except (ValueError,TypeError):pass
    return rows


def write_csv(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def valid(r):
    return r['success']==1 and r['status']=='converged' and not r['boundary1'] and not r['boundary2']


def savefig(fig,folder,name):
    fig.savefig(folder/(name+'.png'),dpi=180)
    fig.savefig(folder/(name+'.pdf'));plt.close(fig)


def summarize(folder):
    dest=folder/'summary';dest.mkdir(exist_ok=True)
    rows=[];summaries=[];bs=[];sensitivity=[];prediction=[]
    for sid in SESSIONS:
        session_rows=read_csv(folder/f'fits_{sid}.csv')
        if len(session_rows)!=120 or any(r['session']!=sid for r in session_rows):
            raise ValueError(f'{sid}: expected 120 correctly identified fit rows')
        rows.extend(session_rows)
        summaries.append(json.loads((folder/sid/'summary.json').read_text()))
        bp=folder/sid/'bootstrap.csv'
        if bp.exists():bs.extend(read_csv(bp))
    write_csv(dest/'all_fits.csv',rows)
    for sid in SESSIONS:
        rr=[r for r in rows if r['session']==sid and not r['native']]
        base={int(r['bin']):r for r in rr if r['variant']=='nominal'}
        for r in rr:
            if r['variant'] in VARIANTS:
                b=base[int(r['bin'])]
                sensitivity.append(dict(session=sid,q_Ainv=r['q_Ainv'],variant=r['variant'],
                    nominal_status=b['status'],variant_status=r['status'],both_unflagged=valid(b) and valid(r),
                    delta_E01=r['E01']-b['E01'],delta_apex1=r['apex1']-b['apex1'],
                    delta_E02=r['E02']-b['E02'],delta_apex2=r['apex2']-b['apex2'],
                    chi2_nominal=b['chi2'],chi2_variant=r['chi2']))
        for b,nom in base.items():
            a=next((r for r in rr if r['bin']==b and r['variant']=='train_one'),None)
            z=next((r for r in rr if r['bin']==b and r['variant']=='train_two'),None)
            if a and z:
                prediction.append(dict(session=sid,q_Ainv=nom['q_Ainv'],one_status=a['status'],two_status=z['status'],
                    one_chi2=a['predictive_chi2'],two_chi2=z['predictive_chi2'],
                    difference_one_minus_two=a['predictive_chi2']-z['predictive_chi2'],
                    two_predictive_bias=z['predictive_bias']))
    write_csv(dest/'sensitivity.csv',sensitivity);write_csv(dest/'temporal_prediction.csv',prediction)
    nominal=[r for r in rows if not r['native'] and r['variant']=='nominal']
    write_csv(dest/'nominal.csv',nominal)
    # All selected q are shown, including boundary fits (crosses).
    fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained',sharex=True)
    for j,sid in enumerate(SESSIONS):
        for sign,ls in [(-1,'--'),(1,'-')]:
            rr=sorted([r for r in nominal if r['session']==sid and np.sign(r['q_Ainv'])==sign],key=lambda x:abs(x['q_Ainv']))
            q=np.array([abs(r['q_Ainv']) for r in rr])
            for p,col in [(1,'tab:blue'),(2,'tab:orange')]:
                axs[0,j].plot(q,[r[f'E0{p}'] for r in rr],ls,color=col,label=f'Component {p}, {sign:+d}q')
                axs[1,j].plot(q,[r[f'apex{p}'] for r in rr],ls,color=col)
                for ax,key in [(axs[0,j],f'E0{p}'),(axs[1,j],f'apex{p}')]:
                    for x,r in zip(q,rr):ax.scatter(x,r[key],marker='o' if valid(r) else 'x',color=col,s=25)
        axs[0,j].set(title=sid,ylabel='Oscillator E0 (meV)')
        axs[1,j].set(xlabel='|q| (Å$^{-1}$)',ylabel='Intrinsic component apex (meV)')
    axs[0,0].legend(fontsize=7)
    fig.suptitle('Fixed two-component model; crosses mark flagged fits; lines are guides')
    savefig(fig,dest,'nominal_energy_definitions')
    # Compare binned and native fits without relabeling either as truth.
    fig,axs=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
    binning=[]
    for j,sid in enumerate(SESSIONS):
        rr=[r for r in nominal if r['session']==sid]
        for r in rr:
            n=next(x for x in rows if x['session']==sid and x['native'] and x['variant']=='nominal' and abs(x['q_Ainv']-r['q_Ainv'])<1e-9)
            binning.append(dict(session=sid,q_Ainv=r['q_Ainv'],native_status=n['status'],binned_status=r['status'],
                delta_E01_binned_minus_native=r['E01']-n['E01'],delta_apex1_binned_minus_native=r['apex1']-n['apex1']))
            axs[j].scatter(r['q_Ainv'],r['E01']-n['E01'],marker='o' if valid(r) and valid(n) else 'x',color='tab:blue')
        axs[j].axhline(0,color='k',lw=.7);axs[j].set(title=sid,xlabel='Signed q (Å$^{-1}$)',ylabel='E01 binned - native (meV)')
    savefig(fig,dest,'binning_sensitivity');write_csv(dest/'binning.csv',binning)
    fig,axs=plt.subplots(2,3,figsize=(12,6.5),layout='constrained')
    # A single broad component has no one-to-one parameter correspondence
    # with the first component of a two-component decomposition.
    plot_variants=[v for v in VARIANTS if v!='one_peak']
    for j,sid in enumerate(SESSIONS):
        rr=[r for r in sensitivity if r['session']==sid and abs(r['q_Ainv']-.005)<1e-9]
        for ax,key in [(axs[0,j],'delta_E01'),(axs[1,j],'delta_apex1')]:
            values=[next(r[key] for r in rr if r['variant']==v) for v in plot_variants]
            colors=['tab:blue' if next(r['both_unflagged'] for r in rr if r['variant']==v) else 'tab:red' for v in plot_variants]
            ax.barh(plot_variants,values,color=colors);ax.axvline(0,color='k',lw=.7);ax.tick_params(axis='y',labelsize=7)
            ax.set(xlabel=('E01' if key=='delta_E01' else 'Apex1')+' change (meV)',title=sid+' at q=+0.005 Å$^{-1}$')
    fig.suptitle('Sensitivity to explicit alternatives; red includes flagged/non-comparable fits')
    savefig(fig,dest,'model_sensitivity')
    fig,axs=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
    for j,sid in enumerate(SESSIONS):
        rr=sorted([r for r in prediction if r['session']==sid],key=lambda r:r['q_Ainv'])
        axs[j].plot([r['q_Ainv'] for r in rr],[r['one_chi2'] for r in rr],'o-',label='One component')
        axs[j].plot([r['q_Ainv'] for r in rr],[r['two_chi2'] for r in rr],'s-',label='Two components')
        axs[j].axhline(1,color='.5',ls=':');axs[j].set(title=sid,xlabel='Signed q (Å$^{-1}$)',ylabel='Held-out mean squared standard residual')
    axs[0].legend(fontsize=8);savefig(fig,dest,'temporal_prediction')
    # Raw-count overlays retain the elastic/background term and residuals.
    fig,axs=plt.subplots(2,3,figsize=(12,6.5),layout='constrained')
    for j,sid in enumerate(SESSIONS):
        r=next(x for x in nominal if x['session']==sid and abs(x['q_Ainv']-.005)<1e-9)
        b=int(r['bin']);f=loadmat(folder/sid/'fit'/f'bin{b:02d}_nominal.mat',simplify_cells=True)['result']
        one=loadmat(folder/sid/'fit'/f'bin{b:02d}_one_peak.mat',simplify_cells=True)['result']
        e=f['energy_meV'];mask=(e>=300)&(e<=1800)
        ax=axs[0,j];ax.plot(e[mask],f['observed'][mask],color='.6',lw=.7,label='Raw aligned sum')
        ax.plot(e[mask],f['prediction'][mask],color='k',lw=1,label='Two components')
        ax.plot(e[mask],one['prediction'][mask],color='tab:red',lw=1,label='One component')
        ax.plot(e[mask],f['background'][mask],ls=':',label='ZLP + auxiliary background')
        for k in range(2):ax.plot(e[mask],f['peaks'][mask,k],ls='--',lw=1,label=f'Component {k+1}')
        ax.set(title=sid+', q=+0.005 Å$^{-1}$',ylabel='Counts')
        axs[1,j].plot(e[mask],(f['residual']/f['sigma'])[mask],lw=.7,color='k')
        axs[1,j].axhline(0,color='.5',lw=.7)
        axs[1,j].set(xlabel='Energy loss (meV)',ylabel='Residual / empirical sigma')
    axs[0,0].legend(fontsize=7)
    savefig(fig,dest,'representative_fits')
    # Bootstrap percentiles are conditional on convergence and an unflagged joint fit.
    # Every rejected draw remains in the raw table and denominator.
    bootstrap=[];paired=[]
    if bs:
        write_csv(dest/'all_bootstrap.csv',bs)
        for sid in SESSIONS:
            for q in [-.005,.005,.015]:
                for block in [1,10]:
                    rr=[r for r in bs if r['session']==sid and abs(r['q_Ainv']-q)<1e-9 and r['block_length']==block]
                    for key in ['E01','E02','apex1','apex2']:
                        good=[r for r in rr if valid(r) and np.isfinite(r[key])];v=np.array([r[key] for r in good])
                        lo,med,hi=np.percentile(v,[2.5,50,97.5]) if len(v) else [np.nan]*3
                        half=np.array([r[key] for r in good if r['replicate']<=64])
                        hlo,hhi=np.percentile(half,[2.5,97.5]) if len(half)>3 else [np.nan]*2
                        all_finite=np.array([r[key] for r in rr if r['success']==1 and np.isfinite(r[key])])
                        alo,ahi=np.percentile(all_finite,[2.5,97.5]) if len(all_finite) else [np.nan]*2
                        bootstrap.append(dict(session=sid,q_Ainv=q,block_length=block,parameter=key,
                            total=len(rr),unflagged=len(good),valid_fraction=len(good)/max(len(rr),1),
                            median=med,lo95=lo,hi95=hi,width95=hi-lo,first64_width95=hhi-hlo,
                            all_success_lo95=alo,all_success_hi95=ahi,
                            sd=np.std(v,ddof=1) if len(v)>1 else np.nan))
            for block in [1,10]:
                for key in ['E01','apex1']:
                    negative={r['replicate']:r for r in bs if r['session']==sid and abs(r['q_Ainv']+.005)<1e-9 and r['block_length']==block}
                    positive={r['replicate']:r for r in bs if r['session']==sid and abs(r['q_Ainv']-.005)<1e-9 and r['block_length']==block}
                    delta=[positive[k][key]-negative[k][key] for k in negative.keys()&positive.keys() if valid(positive[k]) and valid(negative[k])]
                    if delta:
                        lo,med,hi=np.percentile(delta,[2.5,50,97.5]);paired.append(dict(session=sid,block_length=block,parameter=key,
                            n_valid=len(delta),positive_minus_negative_median=med,lo95=lo,hi95=hi))
        write_csv(dest/'bootstrap_intervals.csv',bootstrap);write_csv(dest/'paired_q_asymmetry.csv',paired)
        fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
        for j,sid in enumerate(SESSIONS):
            validity=[]
            for block,offset,col in [(1,-.00015,'tab:blue'),(10,.00015,'tab:orange')]:
                rr=[r for r in bootstrap if r['session']==sid and r['block_length']==block and r['parameter']=='E01']
                x=np.array([r['q_Ainv']+offset for r in rr]);med=np.array([r['median'] for r in rr])
                errors=np.array([[r['median']-r['lo95'] for r in rr],[r['hi95']-r['median'] for r in rr]])
                axs[j].errorbar(x,med,yerr=errors,fmt='o',ms=4,color=col,capsize=3,label=f'Block length {block}')
                validity.append(('IID' if block==1 else 'Block')+': '+', '.join(f"{r['valid_fraction']:.0%}" for r in rr))
            axs[j].set(title=sid+'\nUnflagged draws (-q, +q, high q):\n'+'; '.join(validity),xlabel='Signed q (Å$^{-1}$)',ylabel='E01 and conditional 95% interval (meV)')
            axs[j].title.set_fontsize(9)
        axs[0].legend(fontsize=8);savefig(fig,dest,'bootstrap_E01')
    synthetic=read_csv(folder/'synthetic'/'synthetic_fits.csv')
    noiseless=[r for r in synthetic if r['replicate']==0 and r['n_true']==r['n_fit']]
    synthetic_errors=[abs(r['E01']-700) for r in noiseless]+[abs(r['E02']-1600) for r in noiseless if r['n_true']==2]
    false_double=[r for r in synthetic if r['replicate']>0 and r['n_true']==1 and r['n_fit']==2]
    sstats=dict(total_fits=len(synthetic),noiseless_correct_model_max_E0_error_meV=max(synthetic_errors),
        noisy_single_fit_as_double=len(false_double),noisy_single_double_unflagged=sum(r['status']=='converged' and not r['boundary'] for r in false_double),
        caveat='Only four noisy draws per condition; this is a falsification example, not a calibrated false-positive rate.')
    write_csv(dest/'synthetic_fit_results.csv',synthetic)
    fig,axs=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
    for j,level in enumerate([80,5]):
        for nc,m,col in [(1,'o','tab:blue'),(2,'s','tab:orange')]:
            rr=[r for r in synthetic if r['height']==level and r['n_fit']==nc and r['n_true']==nc and r['replicate']>0]
            axs[j].scatter([r['q_Ainv'] for r in rr],[r['E01']-700 for r in rr],marker=m,color=col,label=f'{nc}-component truth/model')
        axs[j].axhline(0,color='k',lw=.7);axs[j].set(title=f'Injected height {level} counts',xlabel='q (Å$^{-1}$)',ylabel='Recovered E01 - truth (meV)')
    axs[0].legend(fontsize=8);savefig(fig,dest,'synthetic_recovery')
    audit=dict(expected_main_fits=360,actual_main_fits=len(rows),nominal_binned=len(nominal),
        nominal_unflagged=sum(valid(r) for r in nominal),nominal_flagged=sum(not valid(r) for r in nominal),
        numerical_failures=sum(not r['success'] for r in rows),
        status_counts={s:sum(r['status']==s for r in rows) for s in sorted({r['status'] for r in rows})},
        bootstrap_draw_fits=len(bs),bootstrap_expected=2304,synthetic=sstats,
        temporal_two_better=sum(r['difference_one_minus_two']>0 for r in prediction),temporal_total=len(prediction),
        raw_frame_summaries=summaries)
    (dest/'audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='raw_frame_summaries'},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);summarize(p.parse_args().output)
