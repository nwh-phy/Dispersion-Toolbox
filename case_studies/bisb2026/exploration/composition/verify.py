"""Artifact completeness, source immutability and numerical acceptance checks."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from io_utils import digest,dump


def verify(out,previous):
    manifest=json.loads((out/'manifest.json').read_text());checks=[]
    def check(name,ok,detail=None):checks.append(dict(name=name,passed=bool(ok),detail=detail))
    check('previous_manifest_unchanged',digest(previous/'manifest.json')==manifest['previous_manifest_sha256'])
    allfits=[]
    for s in manifest['sessions']:
        sid=s['session'];path=out/f'{sid}.npz';d=dict(np.load(path))
        check(f'{sid}:raw_immutable',digest(s['source'])==s['source_sha256'])
        check(f'{sid}:prepared_immutable',digest(path)==s['prepared_sha256'])
        for name,value in s['previous_artifact_sha256'].items():
            check(f'{sid}:upstream_immutable:{name}',digest(previous/sid/name)==value)
        check(f'{sid}:complete_frames',d['frames'].shape==(8,300,200))
        old=loadmat(previous/sid/'qdrift.mat',simplify_cells=True)
        errors=[]
        for j,q in enumerate(old['q']):
            qi=int(np.argmin(abs(d['q']-q)))
            for k in range(len(d['energy'])):
                native=d['energy_native'][4*k:4*k+4]
                if native[-1]>old['E'][-1]+1e-7:continue
                idx=np.array([np.argmin(abs(old['E']-e)) for e in native])
                assert np.max(abs(old['E'][idx]-native))<1e-7
                expected=old['frames'][idx,:,j].sum(0)
                errors.append(float(np.max(abs(d['frames'][qi,:,k]-expected))))
        check(f'{sid}:previous_processing_parity',max(errors)<2e-4,max(errors))
        folder=out/'fits_v3'/sid;stamp=json.loads((folder/'signature.json').read_text())
        for f in ['models.py','compare.py','io_utils.py']:
            check(f'{sid}:source_hash:{f}',digest(Path(__file__).with_name(f))==stamp[f])
        check(f'{sid}:protocol_hash',digest(Path(__file__).with_name('protocol_v3.toml'))==stamp['protocol'])
        fits=[json.loads(p.read_text()) for p in folder.glob('*.json') if p.name!='signature.json']
        check(f'{sid}:56_fits',len(fits)==56,len(fits));allfits+=fits
        check(f'{sid}:predictions_finite_nonnegative',all(np.all(np.isfinite(f['prediction'])) and
              np.min(f['prediction'])>=0 for f in fits))
    conv=json.loads((out/'convergence.json').read_text())
    check('quadrature:12_comparisons',len(conv)==12)
    check('quadrature:max_error_below_0.1_sigma',max(r['max_abs_sigma'] for r in conv)<.1,
          max(r['max_abs_sigma'] for r in conv))
    synthetic=[]
    for model in ['single','skew_single','double','continuum']:
        fs=[json.loads(p.read_text()) for p in (out/'synthetic'/model).glob('rep*.json')]
        check(f'{model}:68_synthetic_fits',len(fs)==68,len(fs));synthetic+=fs
        f=json.loads((out/'synthetic'/model/f'rep00_{model}.json').read_text())
        truth=json.loads((out/'synthetic'/model/'truth.json').read_text())
        check(f'{model}:noiseless_independent_grid_recovery',
            f['mean_squared_residual']<.003 and abs(f['parameters'][0]-truth['parameters'][0])<.015,
            dict(score=f['mean_squared_residual'],energy_error_eV=f['parameters'][0]-truth['parameters'][0]))
    result=dict(checks=checks,passed=sum(c['passed'] for c in checks),total=len(checks),
        real_fits=len(allfits),real_converged=sum(f['success'] for f in allfits),
        real_with_parameter_boundary=sum(any(f['boundary']) for f in allfits),
        synthetic_fits=len(synthetic),synthetic_converged=sum(f['success'] for f in synthetic),
        source_sha256={p.name:digest(p) for p in Path(__file__).parent.glob('*.py')})
    dump(out/'verification.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['checks','source_sha256']},indent=2))
    if not all(c['passed'] for c in checks):raise RuntimeError('Failed artifact verification; see verification.json')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True);a=p.parse_args();verify(a.output,a.previous)
