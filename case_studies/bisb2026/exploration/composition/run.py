"""One-command reproducible B1 comparison, with at most two CPU jobs by default."""
from __future__ import annotations
import argparse
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from io_utils import digest

HERE=Path(__file__).resolve().parent


def command(script,*args):
    subprocess.run([sys.executable,str(HERE/script),*map(str,args)],check=True)


def simulation_complete(output,truth):
    folder=output/'synthetic'/truth;stamp=folder/'signature.json'
    if not stamp.exists():return False
    expected={f:digest(HERE/f) for f in ['models.py','validate.py']}
    expected.update(data=digest(output/'590.npz'),replicates=16)
    if json.loads(stamp.read_text())!=expected:raise ValueError('Stale simulation cache')
    return (folder/'truth.json').exists() and all(
        (folder/f'rep{rep:02d}_{model}.json').exists() for rep in range(17)
        for model in ['single','skew_single','double','continuum'])


def run(previous,output,resume,workers):
    if not 1<=workers<=4:raise ValueError('Use one to four CPU jobs')
    output=output.resolve();previous=previous.resolve()
    if (output/'manifest.json').exists():
        if not resume:raise ValueError('Existing run: use --resume or a new output directory')
        m=json.loads((output/'manifest.json').read_text())
        if m['previous_manifest_sha256']!=digest(previous/'manifest.json'):
            raise ValueError('Previous source manifest changed')
        for s in m['sessions']:
            if digest(output/f'{s["session"]}.npz')!=s['prepared_sha256']:
                raise ValueError('Prepared data changed')
            for name,value in s['previous_artifact_sha256'].items():
                if digest(previous/s['session']/name)!=value:
                    raise ValueError('Upstream calibration/frames changed')
    else:
        command('prepare.py','--previous',previous,'--output',output)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        jobs=[pool.submit(command,'compare.py','--output',output,'--session',s) for s in ['590','n0','20w']]
        for j in jobs:j.result()
        jobs=[pool.submit(command,'validate.py','--output',output,'--truth',t)
              for t in ['single','skew_single','double','continuum']
              if not (resume and simulation_complete(output,t))]
        for j in jobs:j.result()
    if not (resume and (output/'convergence.json').exists()):
        command('validate.py','--output',output,'--convergence')
    command('summarize.py','--output',output)
    command('verify.py','--output',output,'--previous',previous)
    print('COMPLETE:',output,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--resume',action='store_true')
    p.add_argument('--workers',type=int,default=2);a=p.parse_args();run(a.previous,a.output,a.resume,a.workers)
