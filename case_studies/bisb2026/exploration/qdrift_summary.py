"""Audit centroid-window sensitivity and show the post hoc q-drift comparison."""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'src/python'))
from qeels_frames import momentum_axis


def run(out):
    dest=out/'summary';manifest=json.loads((out/'manifest.json').read_text());robust=[]
    fig,axs=plt.subplots(1,3,figsize=(12,3.4),layout='constrained')
    for j,s in enumerate(manifest['sessions']):
        sid=s['session'];D=loadmat(out/sid/'prepared.mat',simplify_cells=True);Q=loadmat(out/sid/'qdrift.mat',simplify_cells=True)
        raw=np.load(s['source'],mmap_mode='r');q=momentum_axis(raw.shape[1],s['q_zero_column'],s['dq_Ainv'])
        support=D['align_support'].astype(int)-1;offset=D['align_offsets'].astype(int);keep=np.flatnonzero(D['align_valid'])
        ea=D['E_raw'][support]-(s['final_energy_origin_meV']-s['metadata_elastic_peak_meV'])
        for qw in [.0015,.002,.003]:
            inds=np.flatnonzero(abs(q)<=qw+1e-12)
            for ew in [16,24,40]:
                ei=support[abs(ea)<=ew+1e-9]
                p=np.array([raw[t][np.ix_(inds,ei+offset[t])].sum(1,dtype=float) for t in keep])
                centers=(p*q[inds]).sum(1)/p.sum(1);drift=median_filter(centers,size=11,mode='nearest')-np.median(centers)
                robust.append(dict(session=sid,q_reference_halfwidth_Ainv=qw,energy_halfwidth_meV=ew,
                    end_vs_train_change_pixels=float((centers[-50:].mean()-centers[:-50].mean())/s['dq_Ainv']),
                    correlation_with_nominal=float(np.corrcoef(drift,Q['drift_Ainv'])[0,1]),
                    rms_difference_from_nominal_pixels=float(np.sqrt(np.mean((drift-Q['drift_Ainv'])**2))/s['dq_Ainv'])))
                axs[j].plot(keep+1,drift/s['dq_Ainv'],color='.65',lw=.7,alpha=.7)
        axs[j].plot(keep+1,Q['drift_Ainv']/s['dq_Ainv'],color='k',lw=1.2,label='Nominal reference')
        axs[j].set(title=sid,xlabel='Frame',ylabel='Relative elastic-centroid shift (pixels)')
    axs[0].legend(fontsize=8);fig.suptitle('Reference-window sensitivity: 3 q windows × 3 energy windows')
    fig.savefig(dest/'qdrift_window_sensitivity.png',dpi=180);fig.savefig(dest/'qdrift_window_sensitivity.pdf');plt.close(fig)
    with (dest/'qdrift_reference_sensitivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(robust[0]));w.writeheader();w.writerows(robust)
    with (dest/'qdrift_intensity.csv').open() as f:intensity=list(csv.DictReader(f))
    with (dest/'qdrift_fits.csv').open() as f:fits=list(csv.DictReader(f))
    with (dest/'temporal_prediction.csv').open() as f:original=list(csv.DictReader(f))
    fig,axs=plt.subplots(2,3,figsize=(12,6),layout='constrained')
    for j,sid in enumerate(['590','n0','20w']):
        rr=[r for r in intensity if r['session']==sid];x=np.arange(len(rr));labels=[f"{float(r['q_Ainv']):+.3f}" for r in rr]
        axs[0,j].bar(x-.18,[(float(r['original_ratio'])-1)*100 for r in rr],.36,label='Energy alignment only')
        axs[0,j].bar(x+.18,[(float(r['q_corrected_ratio'])-1)*100 for r in rr],.36,label='Also relative q drift')
        axs[0,j].set(title=sid,xticks=x,xticklabels=labels,ylabel='Last50 / first250 band ratio change (%)')
        pred=[]
        for r in rr:
            qv=float(r['q_Ainv'])
            old=next(a for a in original if a['session']==sid and abs(float(a['q_Ainv'])-qv)<1e-9)
            new=next(a for a in fits if a['session']==sid and a['variant']=='train_two' and abs(float(a['q_Ainv'])-qv)<1e-9)
            pred.append((float(old['two_chi2']),float(new['predictive_chi2'])))
        axs[1,j].bar(x-.18,[a for a,b in pred],.36);axs[1,j].bar(x+.18,[b for a,b in pred],.36)
        axs[1,j].axhline(1,color='k',ls=':',lw=.8)
        axs[1,j].set(xticks=x,xticklabels=labels,xlabel='Signed q (Å$^{-1}$)',ylabel='Held-out squared standard residual')
    axs[0,0].legend(fontsize=7);fig.suptitle('Controlled post hoc q-drift check; uncertainty bands remain from the original protocol')
    fig.savefig(dest/'qdrift_improvement.png',dpi=180);fig.savefig(dest/'qdrift_improvement.pdf');plt.close(fig)
    print(json.dumps({sid:{'minimum_trace_correlation':min(r['correlation_with_nominal'] for r in robust if r['session']==sid),
                            'max_rms_pixels':max(r['rms_difference_from_nominal_pixels'] for r in robust if r['session']==sid)} for sid in ['590','n0','20w']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);run(p.parse_args().output)
