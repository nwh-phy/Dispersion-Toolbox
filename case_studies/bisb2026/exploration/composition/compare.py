"""Configured expanded-window comparisons with resumable signed results."""
from __future__ import annotations
import argparse
import json
import time
import tomllib
from pathlib import Path
import numpy as np
from models import Operator,fit
from io_utils import digest,dump

MODELS=['single','skew_single','double','continuum']


def run(out,session,protocol):
    settings=tomllib.loads(protocol.read_text())
    d=dict(np.load(out/f'{session}.npz'));dest=out/settings['fit_directory']/session;dest.mkdir(parents=True,exist_ok=True)
    signature={f:digest(Path(__file__).with_name(f)) for f in ['models.py','compare.py','io_utils.py']}
    signature['protocol']=digest(protocol)
    signature['data']=digest(out/f'{session}.npz');stamp=dest/'signature.json'
    if stamp.exists() and json.loads(stamp.read_text())!=signature:
        raise ValueError('Stale cache: use a new output directory')
    if not stamp.exists():dump(stamp,signature)
    specs=[dict(name='main',window=[500,3400]),dict(name='flat',window=[500,3400],flat=True),
        dict(name='flat_width',window=[500,3400],flat=True,fixed_values={3:0.}),
        dict(name='low350',window=[350,3400]),
        *[dict(name=f'upper{v}',window=[500,v]) for v in [1800,2200,2600,3000]],
        dict(name='half_h',window=[500,3400],h_scale=.5),
        dict(name='point',window=[500,3400],geometry='point'),
        dict(name='no_constant',window=[500,3400],include_constant=False),
        *[dict(name=f'holdout_{l}_{h}',window=[500,3400],exclude=[l,h])
          for l,h in [(600,800),(1100,1500),(1800,2200)]]]
    ops={}
    y=d['frames'][:,:250].sum(1);yt=d['frames'][:,250:].sum(1)
    scale=d['beam'][250:].sum()/d['beam'][:250].sum()
    for spec in specs:
        cfg=dict(spec);name=cfg.pop('name');window=cfg.pop('window');exclude=cfg.pop('exclude',None)
        okw={k:cfg.pop(k) for k in ['h_scale','geometry'] if k in cfg};key=tuple(sorted(okw.items()))
        if key not in ops:ops[key]=Operator(d,**okw)
        op=ops[key];mask=(op.energy>=window[0])&(op.energy<=window[1]);ev=mask.copy()
        if exclude:
            ev=(op.energy>=exclude[0])&(op.energy<=exclude[1]);mask &= ~ev
        for model in settings['models']:
            path=dest/f'{name}_{model}.json'
            if path.exists():continue
            main=dest/f'main_{model}.json';warm=None
            if main.exists():warm=json.loads(main.read_text())['parameters']
            t=time.monotonic()
            model_options={k:settings.get(k,False) for k in ['include_high_tail','band_limited','upper_dispersion']}
            f=fit(op,y,d['sigma_train'],model,mask,n_starts=settings['starts'],warm=warm,**model_options,**cfg)
            f.update(session=session,condition=name,window=window,exclude=exclude,
                seconds=time.monotonic()-t,q_Ainv=d['q'],energy_meV=d['energy'],geometry=okw,
                test_mean_squared_residual=float(np.mean(((yt-scale*f['prediction'])/d['sigma_test'])[:,ev]**2)),
                train_heldout_mean_squared_residual=float(np.mean(((y-f['prediction'])/d['sigma_train'])[:,ev]**2)),
                train_frames=250,test_frames=50)
            dump(path,f)
            print(session,name,model,'fit',round(f['mean_squared_residual'],3),
                'test',round(f['test_mean_squared_residual'],3),'success',f['success'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--session',required=True,choices=['590','n0','20w'])
    p.add_argument('--protocol',type=Path,default=Path(__file__).with_name('protocol_v3.toml'))
    a=p.parse_args();run(a.output,a.session,a.protocol)
