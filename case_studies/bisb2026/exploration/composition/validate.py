"""Independent-grid recovery and conditional noise calibration of model choice."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from models import Operator,fit
from io_utils import digest,dump

MODELS=['single','skew_single','double','continuum']


def simulate(out,truth,replicates=16):
    d=dict(np.load(out/'590.npz'));dest=out/'synthetic'/truth;dest.mkdir(parents=True,exist_ok=True)
    signature={f:digest(Path(__file__).with_name(f)) for f in ['models.py','validate.py']}
    signature['data']=digest(out/'590.npz');signature['replicates']=replicates
    stamp=dest/'signature.json'
    if stamp.exists() and json.loads(stamp.read_text())!=signature:raise ValueError('Stale simulation cache')
    if not stamp.exists():dump(stamp,signature)
    op=Operator(d);dense=Operator(d,n_radial=512,n_quad=35,n_probe=9)
    pars=[.85,.7,.55,0.]
    if truth=='skew_single':pars[2]=.45
    if truth=='double':pars += [1.45,.9]
    elif truth=='continuum':pars += [1.35,.35]
    pars += [2.05,.3,3.15,1.,.12]
    a,names=dense.columns(pars,truth,include_high_tail=True,upper_dispersion=True)
    weights=[1.]
    if truth in ('double','continuum'):weights += [.8]
    weights += [.18,.4,.002,.02]
    coeff=np.tile(weights,(len(d['q']),1))
    if truth in ('double','continuum'):
        coeff[:,1] *= .4+.6*abs(d['q'])/.005
    mu=np.einsum('bek,bk->be',a,coeff)
    m=(op.energy>=500)&(op.energy<=3400);band=(op.energy>=500)&(op.energy<=1800)
    counts=d['frames'][:,:250].sum(1)
    factor=counts[:,band].mean(1)/mu[:,band].mean(1);mu*=factor[:,None]
    scale=d['beam'][250:].sum()/d['beam'][:250].sum()
    dump(dest/'truth.json',dict(parameters=pars,component_names=names,model=truth,
        coefficients=coeff*factor[:,None],mean_train=mu,scale_test=scale,
        noise='Gaussian with real frame-difference variance, independent train/test draws; conditional on this count level',
        generator_grid=dict(radial=512,quadrature=35,probe=9),fit_grid=dict(radial=160,quadrature=25,probe=7)))
    del dense,a
    for rep in range(-1,replicates):
        rng=np.random.default_rng(20261024+1000*MODELS.index(truth)+rep+1)
        y=mu.copy() if rep<0 else mu+rng.normal(size=mu.shape)*d['sigma_train']
        yt=scale*mu if rep<0 else scale*mu+rng.normal(size=mu.shape)*d['sigma_test']
        for model in MODELS:
            file=dest/f'rep{rep+1:02d}_{model}.json'
            if file.exists():continue
            f=fit(op,y,d['sigma_train'],model,m,n_starts=6,include_high_tail=True,
                  upper_dispersion=True,band_limited=True,seed=20261024+rep+1)
            f.update(truth=truth,replicate=rep,noise_free=rep<0,
                test_score=float(np.mean(((yt-scale*f['prediction'])/d['sigma_test'])[:,band]**2)),
                negative_simulated_channels=int(np.count_nonzero(y<0)))
            # Keep the sufficient fitted outputs, not simulated frame cubes.
            f.pop('components');dump(file,f)
        print('SYNTHETIC',truth,rep+1,'/',replicates,flush=True)


def convergence(out):
    rows=[]
    for sid in ['590','n0','20w']:
        d=dict(np.load(out/f'{sid}.npz'));op=Operator(d)
        fine=Operator(d,n_radial=512,n_quad=45,n_probe=11)
        for model in MODELS:
            f=json.loads((out/'fits_v3'/sid/f'main_{model}.json').read_text())
            a,_=op.columns(f['parameters'],model,include_high_tail=True,upper_dispersion=True)
            b,_=fine.columns(f['parameters'],model,include_high_tail=True,upper_dispersion=True)
            pred=np.einsum('bek,bk->be',b,np.array(f['coefficients']))
            e=(pred-np.array(f['prediction']))/d['sigma_train']
            mask=(op.energy>=500)&(op.energy<=3400)
            rows.append(dict(session=sid,model=model,max_abs_sigma=float(np.max(abs(e[:,mask]))),
                rms_sigma=float(np.sqrt(np.mean(e[:,mask]**2))),
                max_component_shape_difference=float(np.max(abs(a[:,mask]-b[:,mask])))))
        del fine,op
    dump(out/'convergence.json',rows)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path)
    p.add_argument('--truth',choices=MODELS);p.add_argument('--replicates',type=int,default=16)
    p.add_argument('--convergence',action='store_true');a=p.parse_args()
    if a.convergence:convergence(a.output)
    else:simulate(a.output,a.truth,a.replicates)
