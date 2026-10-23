"""Separate held-out intensity scaling from spectral-shape mismatch.

The held-out optimal scale is a *post hoc diagnostic*, never a validation score.
This analysis does not change or refit any previously selected model.
"""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src/python'))
from qeels_frames import momentum_axis


def run(out):
    dest=out/'summary';dest.mkdir(exist_ok=True)
    manifest=json.loads((out/'manifest.json').read_text());rows=[];qrows=[]
    fig,axs=plt.subplots(3,3,figsize=(12,9),layout='constrained')
    for j,s in enumerate(manifest['sessions']):
        sid=s['session'];D=loadmat(out/sid/'prepared.mat',simplify_cells=True)
        raw=np.load(s['source'],mmap_mode='r');q=momentum_axis(raw.shape[1],s['q_zero_column'],s['dq_Ainv'])
        inds=np.flatnonzero(abs(q)<=.002+1e-12)
        support=D['align_support'].astype(int)-1;offset=D['align_offsets'].astype(int)
        keep=np.flatnonzero(D['align_valid'])
        # D.E_raw is centered on the summed reference, while the final energy
        # origin is the median frame peak. Recover that exact final origin.
        e0=s['final_energy_origin_meV']-s['metadata_elastic_peak_meV']
        er=D['E_raw'][support]-e0
        energy_indices=support[abs(er)<=16+1e-9]
        profiles=np.array([raw[t][np.ix_(inds,energy_indices+offset[t])].sum(1,dtype=float) for t in keep])
        centers=(profiles*q[inds]).sum(1)/profiles.sum(1)
        widths=np.sqrt((profiles*(q[inds][None,:]-centers[:,None])**2).sum(1)/profiles.sum(1))
        qrows.append(dict(session=sid,elastic_center_first250_Ainv=float(centers[:-50].mean()),elastic_center_last50_Ainv=float(centers[-50:].mean()),
            elastic_center_change_pixels=float((centers[-50:].mean()-centers[:-50].mean())/s['dq_Ainv']),
            width_first250_Ainv=float(widths[:-50].mean()),width_last50_Ainv=float(widths[-50:].mean())))
        axs[0,j].plot(keep+1,centers/s['dq_Ainv'],lw=.8);axs[0,j].axvline(250,color='k',ls=':')
        axs[0,j].set(title=sid,xlabel='Frame',ylabel='Central elastic centroid (q pixels)')
        for qvalue,col in [(-.005,'tab:blue'),(.005,'tab:orange'),(.015,'tab:green')]:
            b=int(np.flatnonzero((abs(D['q']-qvalue)<1e-9)&~D['native'].astype(bool))[0])
            F=D['frames'][:,:,b].astype(float);mask=(D['E']>=300)&(D['E']<=1800)
            trace=F[mask].sum(0)/D['beam']; trace/=trace.mean()
            # Ten-frame means are displayed; native traces remain in prepared.mat.
            axs[1,j].plot(np.arange(5,301,10),trace.reshape(30,10).mean(1),'-o',ms=2,lw=.8,color=col,label=f'q={qvalue:+.3f}')
            f=loadmat(out/sid/'fit'/f'bin{b+1:02d}_train_two.mat',simplify_cells=True)['result']
            e=f['energy_meV'];isig=(e>=300)&(e<=1800)
            idx=np.searchsorted(D['E'],e);y=F[:, -50:].sum(1)[idx]
            sg=D['test_sigma'][:,b][idx];pred=f['prediction']
            alpha_beam=D['beam'][-50:].sum()/D['beam'][:-50].sum()
            alpha_opt=np.sum(y[isig]*pred[isig]/sg[isig]**2)/np.sum(pred[isig]**2/sg[isig]**2)
            before=np.mean(((y-alpha_beam*pred)/sg)[isig]**2)
            after=np.mean(((y-alpha_opt*pred)/sg)[isig]**2)
            ratio=(F[mask,-50:].sum()/D['beam'][-50:].sum())/(F[mask,:-50].sum()/D['beam'][:-50].sum())
            rows.append(dict(session=sid,q_Ainv=qvalue,heldout_band_per_elastic_ratio=ratio,
                strict_predictive_residual2=before,posthoc_scale_relative_to_elastic=alpha_opt/alpha_beam,
                posthoc_scaled_residual2=after,interpretation='Post hoc scale consumes holdout information; not independent validation'))
            axs[2,j].plot(e[isig],((y-alpha_beam*pred)/sg)[isig],color=col,alpha=.6,lw=.5)
        axs[1,j].axvline(250,color='k',ls=':');axs[1,j].set(xlabel='Frame',ylabel='Band / elastic reference, relative')
        axs[2,j].axhline(0,color='k',lw=.5);axs[2,j].set(xlabel='Energy (meV)',ylabel='Strict temporal prediction residual / sigma')
    axs[1,0].legend(fontsize=8)
    fig.savefig(dest/'temporal_diagnostics.png',dpi=180);fig.savefig(dest/'temporal_diagnostics.pdf');plt.close(fig)
    for name,r in [('time_scaling_diagnostic',rows),('elastic_q_drift',qrows)]:
        with (dest/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(r[0]));w.writeheader();w.writerows(r)
    print(json.dumps(dict(temporal_scaling=rows,elastic_q_drift=qrows),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);run(p.parse_args().output)
