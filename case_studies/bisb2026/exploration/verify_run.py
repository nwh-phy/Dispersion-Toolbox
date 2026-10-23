"""Audit input immutability, exact fit inputs, coverage and numerical checks."""
import argparse
import csv
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import numpy as np
from scipy.io import loadmat

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'src/python'))
from qeels_frames import sha256_file


def verify(out):
    m=json.loads((out/'manifest.json').read_text());checks=[]
    for s in m['sessions']:
        source=Path(s['source']);folder=out/s['session'];prepared=folder/'prepared.mat'
        checks.extend([dict(check='raw_sha256',session=s['session'],passed=sha256_file(source)==s['sha256']),
                       dict(check='metadata_sha256',session=s['session'],passed=sha256_file(source.with_suffix('.json'))==s['metadata_sha256']),
                       dict(check='retained_count_conservation',session=s['session'],passed=s['max_retained_count_relative_error']<1e-12)])
        for stage in ['fit','bootstrap']:
            stamp=loadmat(folder/stage/'signature.mat',simplify_cells=True)['signature']
            expected=dict(prepared=sha256_file(prepared),runner=sha256_file(HERE/'bisb_explore_fit.m'),
                          fitter=sha256_file(ROOT/'src/fitting/qe_zlp_joint_fit.m'),kinematic=sha256_file(ROOT/'src/fitting/qe_kinematic_prefactor.m'))
            checks.append(dict(check=stage+'_cache_signature',session=s['session'],passed=stamp==expected))
        mat=loadmat(prepared,simplify_cells=True)
        checks.append(dict(check='positive_noise_uniform_energy',session=s['session'],passed=bool(np.all(mat['sigma']>0) and np.allclose(np.diff(mat['E']),4))))
    audit=json.loads((out/'summary/audit.json').read_text())
    checks.extend([dict(check='360_main_fits',passed=audit['actual_main_fits']==360),
                   dict(check='80_synthetic_fits',passed=audit['synthetic']['total_fits']==80),
                   dict(check='2304_bootstrap_fits',passed=audit['bootstrap_draw_fits']==2304),
                   dict(check='noiseless_recovery_below_1meV',passed=audit['synthetic']['noiseless_correct_model_max_E0_error_meV']<1)])
    with (out/'summary/qdrift_fits.csv').open() as f:qr=list(csv.DictReader(f))
    with (out/'summary/qdrift_intensity.csv').open() as f:qi=list(csv.DictReader(f))
    checks.append(dict(check='18_posthoc_qdrift_fits',passed=len(qr)==18 and all(r['success']=='1' for r in qr)))
    checks.append(dict(check='q_bin_overlap_conservation',passed=len(qi)==9 and max(float(r['max_weight_sum_error']) for r in qi)<1e-12))
    paths=list(HERE.glob('*.py'))+list(HERE.glob('*.m'))+list(HERE.glob('*.toml'))+list((ROOT/'src/python').glob('*.py'))
    evidence=dict(time_unix=time.time(),python=sys.version,platform=platform.platform(),
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        worktree_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True).splitlines(),
        source_sha256={str(p.relative_to(ROOT)):sha256_file(p) for p in paths},
        packages=subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True).splitlines(),
        checks=checks,all_passed=all(c['passed'] for c in checks))
    (out/'validation').mkdir(exist_ok=True)
    (out/'validation/verification.json').write_text(json.dumps(evidence,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(dict(all_passed=evidence['all_passed'],n_checks=len(checks)),indent=2))
    if not evidence['all_passed']:raise AssertionError('Run verification failed; inspect validation/verification.json')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);verify(p.parse_args().output.resolve())
