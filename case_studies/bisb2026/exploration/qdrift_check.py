"""Post hoc controlled q-drift experiment using elastic-reference centroids."""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat,savemat
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'src/python'))
from qeels_frames import momentum_axis,noise_from_frames
from qeels_rebin import fractional_bin_weights


def run(out):
    manifest=json.loads((out/'manifest.json').read_text());rows=[]
    for s in manifest['sessions']:
        sid=s['session'];folder=out/sid;D=loadmat(folder/'prepared.mat',simplify_cells=True)
        raw=np.load(s['source'],mmap_mode='r');q=momentum_axis(raw.shape[1],s['q_zero_column'],s['dq_Ainv'])
        support=D['align_support'].astype(int)-1;offset=D['align_offsets'].astype(int);keep=np.flatnonzero(D['align_valid'])
        e0=s['final_energy_origin_meV']-s['metadata_elastic_peak_meV'];ea=D['E_raw'][support]-e0
        elastic_indices=support[abs(ea)<=16+1e-9];qq=np.flatnonzero(abs(q)<=.002+1e-12)
        profile=np.array([raw[t][np.ix_(qq,elastic_indices+offset[t])].sum(1,dtype=float) for t in keep])
        centers=(profile*q[qq]).sum(1)/profile.sum(1)
        # Retain the original q reference; translate relative to the sequence
        # median only. Eleven-frame median suppresses centroid shot noise.
        drift=median_filter(centers,size=11,mode='nearest')-np.median(centers)
        ii=np.flatnonzero((ea>=D['E'][0]-1e-7)&(ea<=D['E'][-1]+1e-7));assert len(ii)==len(D['E'])
        base_bins=[int(np.flatnonzero((abs(D['q']-x)<1e-9)&~D['native'].astype(bool))[0]) for x in [-.005,.005,.015]]
        frames=[];sigmas=[];trains=[];tests=[]
        for b in base_bins:
            f=[];errors=[]
            for it,t in enumerate(keep):
                w=fractional_bin_weights(q,D['q'][b],D['width'][b],drift[it]);members=np.flatnonzero(w>1e-10)
                if any(k+1 in s['invalid_columns_1based'] for k in members):raise ValueError('Bad channel in translated bin')
                f.append(w[members]@raw[t][np.ix_(members,support[ii]+offset[t])].astype(float))
                errors.append(abs(w.sum()-D['width'][b]/D['dq']))
            F=np.array(f)
            sg=noise_from_frames(F,D['beam'])['sigma'];train=noise_from_frames(F[:-50],D['beam'][:-50])['sigma'];test=noise_from_frames(F[-50:],D['beam'][-50:])['sigma']
            band=(D['E']>=300)&(D['E']<=1800)
            original=D['frames'][:,:,b].astype(float).T
            def ratio(x):return (x[-50:,band].sum()/D['beam'][-50:].sum())/(x[:-50,band].sum()/D['beam'][:-50].sum())
            rows.append(dict(session=sid,q_Ainv=D['q'][b],original_ratio=ratio(original),q_corrected_ratio=ratio(F),
                max_weight_sum_error=max(errors),drift_span_pixels=float(np.ptp(drift)/D['dq'])))
            frames.append(F.T);sigmas.append(sg);trains.append(train);tests.append(test)
        savemat(folder/'qdrift.mat',dict(E=D['E'][:,None],frames=np.stack(frames,-1),sigma=np.stack(sigmas,1),
            train_sigma=np.stack(trains,1),test_sigma=np.stack(tests,1),beam=D['beam'][:,None],base_bins=np.array(base_bins)+1,
            q=D['q'][base_bins],drift_Ainv=drift,centers_Ainv=centers),do_compression=True)
        np.savetxt(folder/'q_drift.csv',np.c_[keep+1,centers,drift],delimiter=',',header='frame,elastic_centroid_Ainv,applied_relative_offset_Ainv',comments='')
    with (out/'summary/qdrift_intensity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(rows,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);run(p.parse_args().output)
