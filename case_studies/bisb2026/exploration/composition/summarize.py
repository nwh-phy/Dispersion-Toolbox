"""Evidence tables and figures; no silent exclusion of failed/boundary fits."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io_utils import dump

MODELS=['single','skew_single','double','continuum']
LABELS=['One DL','One skew band','Two DL','DL + spectral continuum']
COLORS=['#8b8e96','#9c79b8','#176b87','#d6842d']


def read(path):return json.loads(path.read_text())


def csv_write(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def apex(e,g):
    c=2*e*e-g*g
    return np.sqrt((c+np.sqrt(c*c+12*e**4))/6)


def sorted_apices(f):
    p=f['parameters'];m=f['model']
    a=p[0] if m=='skew_single' else apex(p[0],p[2])
    if m=='double':a=np.array([a,apex(p[4],p[5])])
    elif m=='continuum':a=np.array([a,p[4]])
    else:a=np.array([a])
    return a,np.argsort(a)


def summarize(out):
    dest=out/'summary';dest.mkdir(exist_ok=True)
    rows=[];components=[];paired=[];ridges=[]
    fig,axes=plt.subplots(3,3,figsize=(12,9),sharex=True)
    residual_fig,resaxes=plt.subplots(3,1,figsize=(10,8),sharex=True)
    score_fig,scoreaxes=plt.subplots(1,3,figsize=(12,3.6),sharey=True)
    for si,sid in enumerate(['590','n0','20w']):
        d=dict(np.load(out/f'{sid}.npz'));e=d['energy'];band=(e>=500)&(e<=1800)
        y=d['frames'][:,:250].sum(1);yt=d['frames'][:,250:].sum(1)
        scale=d['beam'][250:].sum()/d['beam'][:250].sum()
        folder=out/'fits_v3'/sid
        fits={m:read(folder/f'main_{m}.json') for m in MODELS}
        for path in sorted(folder.glob('*.json')):
            if path.name=='signature.json':continue
            f=read(path);pred=np.array(f['prediction']);a,order=sorted_apices(f)
            res=(yt-scale*pred)/d['sigma_test']
            row=dict(session=sid,condition=f['condition'],model=f['model'],success=f['success'],
                boundary_count=sum(f['boundary']),train_score=f['mean_squared_residual'],
                test_fitwindow_score=f['test_mean_squared_residual'],test_B1_score=float(np.mean(res[:,band]**2)),
                low_apex_eV=float(a.min()),high_apex_eV=float(a.max()) if len(a)>1 else np.nan,
                parameter_count=f['parameter_count'],seconds=f['seconds'])
            rows.append(row)
        for mi,m in enumerate(MODELS):
            f=fits[m];pred=np.array(f['prediction']);res=(yt-scale*pred)/d['sigma_test']
            scoreaxes[si].bar(mi,np.mean(res[:,band]**2),color=COLORS[mi])
            resaxes[si].plot(e/1000,res[5],color=COLORS[mi],label=LABELS[mi],lw=1)
        scoreaxes[si].set(title=sid,xticks=range(4),xticklabels=['DL','Skew','2 DL','DL+C'])
        scoreaxes[si].axhline(1,color='k',ls=':',lw=1)
        resaxes[si].axhline(0,color='k',lw=.6);resaxes[si].set(ylim=(-5,5),ylabel=f'{sid}: residual / sigma')
        # Plot counts held out in time; model components use only training data.
        f=fits['double'];a,order=sorted_apices(f)
        cc=np.array(f['components']);cp=np.array(fits['continuum']['prediction'])
        for col,qi in enumerate([4,5,7]):
            ax=axes[si,col]
            ax.errorbar(e[::2]/1000,yt[qi,::2],yerr=d['sigma_test'][qi,::2],fmt='.',ms=2,
                        color='#555b63',alpha=.65,elinewidth=.4,label='Last 50 frames')
            ax.plot(e/1000,scale*np.array(f['prediction'])[qi],color=COLORS[2],label='Two DL prediction')
            ax.plot(e/1000,scale*cp[qi],color=COLORS[3],ls='--',lw=1,label='DL + continuum')
            for j,idx in enumerate(order):
                ax.plot(e/1000,scale*cc[qi,:,idx],color=['#2d9c80','#b15d56'][j],lw=1.2,
                        label=f'B1 shape {j+1} ({a[idx]:.2f} eV apex)')
            ax.plot(e/1000,scale*cc[qi,:,2:].sum(1),color='#a8a8a8',ls=':',label='Neighbors + background')
            ax.set(xlim=(.5,2.55),ylim=(0,max(yt[qi,band])*1.25),title=f'{sid}, q = {d["q"][qi]:.4f} Å$^{{-1}}$')
        for qi,q in enumerate(d['q']):
            areas=cc[qi,band,:].sum(0)
            components.append(dict(session=sid,q_Ainv=q,model='double',
                low_fraction_of_B1=float(areas[order[0]]/areas[:2].sum()),
                neighbor_background_fraction=float(areas[2:].sum()/areas.sum()),
                low_apex_eV=float(a.min()),high_apex_eV=float(a.max())))
        # Conditional block bootstrap of the last 50 measured frames. Preserve
        # all q and energy correlations by drawing whole frames together.
        rng=np.random.default_rng(20261024);test=d['frames'][:,250:];beam=d['beam'][250:]
        draws=[]
        for _ in range(2000):
            starts=rng.integers(0,50,10)
            draws.append(((starts[:,None]+np.arange(5))%50).ravel())
        pairs=[('single','double'),('skew_single','double'),('double','continuum')]
        for ma,mb in pairs:
            pa=np.array(fits[ma]['prediction']);pb=np.array(fits[mb]['prediction']);deltas=[]
            for draw in draws:
                yy=test[:,draw].sum(1);ss=beam[draw].sum()/d['beam'][:250].sum()
                ra=(yy-ss*pa)/d['sigma_test'];rb=(yy-ss*pb)/d['sigma_test']
                deltas.append(float(np.mean((ra*ra-rb*rb)[:,band])))
            ra=(yt-scale*pa)/d['sigma_test'];rb=(yt-scale*pb)/d['sigma_test']
            paired.append(dict(session=sid,model_a=ma,model_b=mb,
                delta_score_a_minus_b=float(np.mean((ra*ra-rb*rb)[:,band])),
                ci_low=float(np.quantile(deltas,.025)),ci_high=float(np.quantile(deltas,.975)),
                block_frames=5,bootstrap_repeats=2000,
                meaning='Conditional on fitted models, variance and preprocessing; not training/calibration uncertainty'))
        flat=read(folder/'flat_width_double.json');ap,od=sorted_apices(flat)
        for qi in [4,5,6,7]:
            observed=gaussian_filter1d(d['frames'][qi].sum(0),2)
            pred=np.array(flat['prediction'])[qi];valid=(e>=500)&(e<=1800)
            er=e[valid]
            ridges.append(dict(session=sid,q_Ainv=d['q'][qi],
                observed_smoothed_ridge_meV=float(er[np.argmax(observed[valid])]),
                fixed_intrinsic_model_count_ridge_meV=float(er[np.argmax(pred[valid])]),
                intrinsic_low_apex_meV=float(ap.min()*1000),intrinsic_high_apex_meV=float(ap.max()*1000)))
    axes[0,0].legend(fontsize=6.5,loc='upper right')
    for ax in axes[-1]:ax.set_xlabel('Energy loss (eV)')
    for ax in axes[:,0]:ax.set_ylabel('Counts / 16 meV / 50 frames')
    fig.tight_layout();fig.savefig(dest/'heldout_components.png',dpi=180);fig.savefig(dest/'heldout_components.pdf');plt.close(fig)
    resaxes[0].legend(ncol=4,fontsize=8);resaxes[-1].set_xlabel('Energy loss (eV)');resaxes[-1].set_xlim(.5,3.4)
    residual_fig.tight_layout();residual_fig.savefig(dest/'heldout_residuals.png',dpi=160);plt.close(residual_fig)
    scoreaxes[0].set_ylabel('Held-out mean squared residual, B1')
    score_fig.tight_layout();score_fig.savefig(dest/'model_comparison.png',dpi=180);plt.close(score_fig)
    csv_write(dest/'model_scores.csv',rows);csv_write(dest/'component_areas.csv',components)
    csv_write(dest/'paired_time_bootstrap.csv',paired);csv_write(dest/'ridge_counterexample.csv',ridges)
    fig,axs=plt.subplots(1,3,figsize=(12,3.5),sharey=True)
    for ax,sid in zip(axs,['590','n0','20w']):
        rr=[r for r in ridges if r['session']==sid];q=[r['q_Ainv'] for r in rr]
        ax.plot(q,[r['observed_smoothed_ridge_meV'] for r in rr],'o',color='#444',label='Observed count ridge')
        ax.plot(q,[r['fixed_intrinsic_model_count_ridge_meV'] for r in rr],'-',color=COLORS[2],label='Flat intrinsic bands: count ridge')
        ax.axhline(rr[0]['intrinsic_low_apex_meV'],color='#2d9c80',ls='--',label='Fixed low loss apex')
        ax.axhline(rr[0]['intrinsic_high_apex_meV'],color='#b15d56',ls='--',label='Fixed high loss apex')
        ax.set(title=sid,xlabel='q (Å$^{-1}$)',ylim=(450,1650))
    axs[0].set_ylabel('Energy (meV)');axs[0].legend(fontsize=7)
    fig.tight_layout();fig.savefig(dest/'ridge_counterexample.png',dpi=180);plt.close(fig)
    synth=[]
    for folder in sorted((out/'synthetic').glob('*')):
        for file in sorted(folder.glob('rep*.json')):
            f=read(file);synth.append(dict(truth=f['truth'],model=f['model'],replicate=f['replicate'],
                noise_free=f['noise_free'],success=f['success'],test_score=f['test_score'],
                fit_score=f['mean_squared_residual'],negative_channels=f['negative_simulated_channels']))
    csv_write(dest/'synthetic_scores.csv',synth)
    dump(dest/'counts.json',dict(real_fits=len(rows),real_success=sum(r['success'] for r in rows),
        synthetic_fits=len(synth),synthetic_success=sum(r['success'] for r in synth)))
    print(json.dumps(paired,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);summarize(p.parse_args().output)
