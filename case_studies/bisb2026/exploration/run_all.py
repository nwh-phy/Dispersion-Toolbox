"""Single-command local reproduction; raw data never enter subprocess code.

Use --resume only with the same config/preparation source/data. MATLAB fit
checkpoints independently verify the prepared arrays and fitter source hashes.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'src/python'))
from qeels_frames import sha256_file


def matlab_string(path):
    return "'"+str(path).replace("'","''")+"'"


def run(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    logs=out/'driver_logs';logs.mkdir(exist_ok=True)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    def command(argv,log):
        print('Running:',Path(log).name,flush=True)
        with Path(log).open('w') as stream:
            subprocess.run(argv,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,check=True)
    def mat(body,log):
        prefix='addpath('+matlab_string(HERE)+');'
        command([args.matlab,'-singleCompThread','-batch',prefix+body],log)
    manifest=out/'manifest.json'
    if args.resume:
        m=json.loads(manifest.read_text())
        if sha256_file(HERE/'datasets.toml')!=m['config_sha256']:
            raise ValueError('Preparation config changed; use a fresh output directory')
        for p in [HERE/'run.py',ROOT/'src/python/qeels_frames.py']:
            if sha256_file(p)!=m['code_sha256'][str(p.relative_to(ROOT))]:
                raise ValueError('Preparation source changed; use a fresh output directory')
        for s in m['sessions']:
            source=Path(s['source'])
            if args.data_root.resolve() not in source.parents or sha256_file(source)!=s['sha256'] or sha256_file(source.with_suffix('.json'))!=s['metadata_sha256']:
                raise ValueError('Raw source, metadata, or selected data root changed')
    elif manifest.exists():
        raise ValueError('Existing output: choose --resume or a fresh directory')
    else:
        command([sys.executable,str(HERE/'run.py'),'prepare','--data-root',str(args.data_root.resolve()),'--output',str(out)],logs/'prepare.log')
    command([sys.executable,'-m','pytest','tests/python','-q'],logs/'python_tests.log')
    syn=out/'synthetic'
    # These independent examples are cheap; synthetic fitting is cached only
    # when the existing synthetic inputs and driver match their saved stamp.
    syn_signature={str(p.relative_to(ROOT)):sha256_file(p) for p in [HERE/'mechanism_demo.py',HERE/'bisb_explore_synthetic.m',ROOT/'src/python/qeels_response.py',ROOT/'src/fitting/qe_zlp_joint_fit.m',ROOT/'src/fitting/qe_kinematic_prefactor.m']}
    stamp=syn/'driver_signature.json'
    if args.resume and (syn/'fits').exists() and not stamp.exists():
        raise ValueError('Existing synthetic fits have no driver signature; audit them or use a fresh output')
    if stamp.exists():
        if json.loads(stamp.read_text())!=syn_signature:raise ValueError('Synthetic source changed; use a fresh output')
    else:
        command([sys.executable,str(HERE/'mechanism_demo.py'),'--output',str(syn)],logs/'mechanism.log')
        stamp.write_text(json.dumps(syn_signature,indent=2))
    command([sys.executable,str(HERE/'screening_demo.py'),'--output',str(syn)],logs/'screening.log')
    mat('bisb_explore_fit('+matlab_string(out)+",'check');",logs/'parity.log')
    def session(sid):
        for stage in ['fit','bootstrap']:
            mat('bisb_explore_fit('+matlab_string(out)+','+matlab_string(stage)+','+matlab_string(sid)+');',logs/f'{sid}_{stage}.log')
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures=[pool.submit(session,sid) for sid in ['590','n0','20w']]
        futures.append(pool.submit(mat,'bisb_explore_synthetic('+matlab_string(syn)+');',logs/'synthetic.log'))
        for f in futures:f.result()
    command([sys.executable,str(HERE/'summarize.py'),'--output',str(out)],logs/'summary.log')
    command([sys.executable,str(HERE/'temporal_diagnostics.py'),'--output',str(out)],logs/'temporal_diagnostics.log')
    command([sys.executable,str(HERE/'qdrift_check.py'),'--output',str(out)],logs/'qdrift_prepare.log')
    mat('bisb_explore_qdrift('+matlab_string(out)+');',logs/'qdrift_fit.log')
    command([sys.executable,str(HERE/'qdrift_summary.py'),'--output',str(out)],logs/'qdrift_summary.log')
    audit=json.loads((out/'summary/audit.json').read_text())
    assert audit['actual_main_fits']==360 and audit['bootstrap_draw_fits']==2304 and audit['synthetic']['total_fits']==80
    command([sys.executable,str(HERE/'verify_run.py'),'--output',str(out)],logs/'verify.log')
    print('Complete:',out/'summary/audit.json',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data-root',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--matlab',default='/Applications/MATLAB_R2026a.app/bin/matlab');p.add_argument('--jobs',type=int,choices=range(1,5),default=3)
    p.add_argument('--resume',action='store_true');run(p.parse_args())
