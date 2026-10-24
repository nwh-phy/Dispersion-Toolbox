"""Independent low-B1-component analysis on adaptive q bins.

Stages prepare/main/bootstrap/checks checkpoint separately. Native low-q bins
are fitted once and shared by the 3/5 high-q strategies. Resampling uses the
same frame-block weights across q; train/test draws remain independent.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import time
import tomllib
import numpy as np
from prepare import prepare
from qbin_fit import fit_local, local_columns
from models import Operator, ROOT
from kinematic_ablation import AblationOperator
from dispersion_models import endpoints,loss_apex
from ridge_extract import block_weights
from io_utils import digest,dump


BOOT_NAMES=['low_E0_eV','high_E0_eV','low_loss_apex_eV','low_count_apex_eV',
            'high_count_apex_eV','low_area_fraction','low_gamma_eV','high_gamma_eV',
            'test_gain_fraction','test_delta_cost','success','boundary','apex_inversion']


def read(path):
    return json.loads(Path(path).read_text())


def local_data(task):
    d=dict(np.load(task['data_path']));i=task['qi']
    d['q']=d['q'][i:i+1];d['width']=float(d['width'][i])
    for k in ['frames','sigma_train','sigma_test','sigma_all']:
        d[k]=d[k][i:i+1]
    return d


def count_apex(energy,y):
    use=np.flatnonzero((energy>=200)&(energy<=1800));v=y[use];j=int(v.argmax())
    if j in (0,len(v)-1):return float(energy[use[j]]/1000),True
    a,b,c=v[j-1:j+2];den=a-2*b+c
    dx=np.clip(.5*(a-c)/den,-.5,.5) if den<0 else 0.
    return float((energy[use[j]]+dx*np.median(np.diff(energy)))/1000),False


def quantities(op,fit):
    comp=np.asarray(fit['components'])[0];band=(op.energy>=300)&(op.energy<=1800)
    low,edge=count_apex(op.energy,comp[:,0])
    if fit['model']=='A':
        high,edge_high=count_apex(op.energy,comp[:,1])
        frac=comp[band,0].sum()/max(comp[band,:2].sum(),1e-100)
        gap=fit['E0_high_eV']-fit['E0_low_eV']
    else:
        high=np.nan;edge_high=False;frac=1.;gap=np.nan
    return dict(low_count_apex_eV=low,high_count_apex_eV=high,
        low_count_apex_on_edge=edge,high_count_apex_on_edge=edge_high,
        low_B1_area_fraction=float(frac),E0_separation_eV=float(gap),
        count_apex_order_inversion=bool(np.isfinite(high) and low>=high))


def fit_pair(task,cfg,mode='full',window=None,bootstrap=False):
    d=local_data(task);op=Operator(d) if mode=='full' else AblationOperator(d,'none')
    nt=cfg['train_frames'];y=d['frames'][:,:nt].astype(float).sum(1)
    test=d['frames'][:,nt:].astype(float).sum(1)
    sigma=d['sigma_train'];window=window or cfg['fit_window_meV']
    mask=(op.energy>=window[0])&(op.energy<=window[1])
    score=(op.energy>=cfg['score_window_meV'][0])&(op.energy<=cfg['score_window_meV'][1])
    scale=float(d['beam'][nt:].sum()/d['beam'][:nt].sum());fits={}
    for m in ['S','A']:
        f=fit_local(op,y,sigma,m,mask,starts=cfg['starts'],seed=cfg['seed']+task['seed_index'],single=fits.get('S'))
        residual=(test-scale*np.asarray(f['prediction']))/d['sigma_test']
        f.update(quantities(op,f),test_b1_msr=float(np.mean(residual[:,score]**2)),
            test_scale=scale,kinematic_mode=mode,fit_window_meV=window)
        fits[m]=f
    if fits['A']['cost']>fits['S']['cost']+1e-3:
        raise RuntimeError('Double model failed to recover nested single fit')
    return fits


def worker_main(item):
    task,cfg,output,stage=item;dest=Path(output)/'fits'/task['key'];dest.mkdir(parents=True,exist_ok=True)
    t=time.monotonic()
    if stage=='main':
        path=dest/'main.json'
        if not path.exists():dump(path,fit_pair(task,cfg))
    elif stage=='checks':
        for mode,window,name in [('none',cfg['fit_window_meV'],'noK'),('full',cfg['check_window_meV'],'window400')]:
            path=dest/f'{name}.json'
            if not path.exists():dump(path,fit_pair(task,cfg,mode,window))
    else:
        path=dest/'bootstrap.npz'
        if path.exists():return task['key'],0.
        d=local_data(task);op=Operator(d);nt=cfg['train_frames'];nr=cfg['bootstrap_repeats']
        main=read(dest/'main.json')
        rngseed=cfg['seed']+['590','n0','20w'].index(task['session'])*10000
        wt=block_weights(nt,nr,cfg['bootstrap_block_frames'],rngseed)
        wv=block_weights(len(d['beam'])-nt,nr,cfg['bootstrap_block_frames'],rngseed+4000)
        ys=wt@d['frames'][0,:nt].astype(float)
        tests=wv@d['frames'][0,nt:].astype(float)
        scales=(wv@d['beam'][nt:])/(wt@d['beam'][:nt])
        mask=(op.energy>=cfg['fit_window_meV'][0])&(op.energy<=cfg['fit_window_meV'][1])
        band=(op.energy>=cfg['score_window_meV'][0])&(op.energy<=cfg['score_window_meV'][1])
        result=[]
        for b in range(nr):
            fits={}
            for m in ['S','A']:
                fits[m]=fit_local(op,ys[b:b+1],d['sigma_train'],m,mask,starts=cfg['bootstrap_starts'],
                    seed=rngseed+b,single=fits.get('S'),warms=[main[m]['parameters']])
                r=(tests[b:b+1]-scales[b]*np.asarray(fits[m]['prediction']))/d['sigma_test']
                fits[m]['test_cost']=float(np.sum(r[:,band]**2))
            a=fits['A'];s=fits['S'];v=quantities(op,a)
            result.append([a['E0_low_eV'],a['E0_high_eV'],a['low_apex_eV'],v['low_count_apex_eV'],
                v['high_count_apex_eV'],v['low_B1_area_fraction'],a['gamma_low_eV'],a['parameters'][6],
                1-a['test_cost']/max(s['test_cost'],1e-100),s['test_cost']-a['test_cost'],
                float(a['success'] and s['success'] and a['cost']<=s['cost']+1e-3),
                float(bool(a['boundary'])),float(v['count_apex_order_inversion'])])
        np.savez_compressed(path,values=np.array(result),names=np.array(BOOT_NAMES))
    return task['key'],time.monotonic()-t


def setup(previous,output,protocol):
    cfg=tomllib.loads(protocol.read_text());output.mkdir(parents=True,exist_ok=True)
    signature=dict(protocol=cfg,previous_manifest=digest(previous/'manifest.json'),
        code={n:digest(Path(__file__).with_name(n)) for n in ['low_branch_run.py','prepare.py',
            'qbin_fit.py','models.py','dispersion_models.py','kinematic_ablation.py','ridge_extract.py','io_utils.py']},
        shared={n:digest(ROOT/'src/python'/n) for n in ['qeels_frames.py','qeels_rebin.py','qeels_response.py']})
    stamp=output/'signature.json'
    if stamp.exists() and read(stamp)!=signature:raise ValueError('Signature changed: use a new output directory')
    dump(stamp,signature)
    tasks=[];mappings=[]
    for n in cfg['bin_channels']:
        base=output/f'bin{n}'
        if not (base/'manifest.json').exists():prepare(previous,base,protocol,n)
        manifest=read(base/'manifest.json')
        assert manifest['protocol_sha256']==digest(protocol)
        assert manifest['preparation_sha256']==signature['code']['prepare.py']
        for sid in ['590','n0','20w']:
            rec=next(r for r in manifest['sessions'] if r['session']==sid)
            assert digest(base/f'{sid}.npz')==rec['prepared_sha256']
            d=dict(np.load(base/f'{sid}.npz'))
            for i,q in enumerate(d['q']):
                if abs(q)<cfg['fit_q_min_Ainv']-1e-10:continue
                native=abs(q)<=cfg['native_q_max_Ainv']+1e-10
                key=f'{sid}/{"native" if native else f"bin{n}"}/{q:+.5f}'
                mappings.append(dict(session=sid,high_channels=n,q_Ainv=float(q),width_Ainv=float(d['width'][i]),key=key))
                if native and n==5:continue
                tasks.append(dict(session=sid,key=key,data_path=str((base/f'{sid}.npz').resolve()),qi=i,
                                  seed_index=len(tasks),q_Ainv=float(q),native=native))
    dump(output/'tasks.json',tasks);dump(output/'mappings.json',mappings)
    return cfg,tasks


def run(previous,output,protocol,stage,workers):
    cfg,tasks=setup(previous,output,protocol)
    if stage=='prepare':return
    if stage=='checks':
        selected=[]
        for sid in ['590','n0','20w']:
            for sign in [-1,1]:
                available=[t for t in tasks if t['session']==sid and t['q_Ainv']*sign>0
                           and ('bin5' in t['key'] or t['native'])]
                for q in cfg['anchor_q_Ainv']:
                    selected.append(min(available,key=lambda t:abs(abs(t['q_Ainv'])-q)))
        tasks=list({t['key']:t for t in selected}.values())
    print('STAGE',stage,'TASKS',len(tasks),'WORKERS',workers,flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(worker_main,(t,cfg,str(output),stage)) for t in tasks]
        for i,f in enumerate(as_completed(futures)):
            key,seconds=f.result();print(f'{stage} {i+1}/{len(tasks)} {key} {seconds:.1f}s',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--previous',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--protocol',type=Path,default=Path(__file__).with_name('low_branch_protocol.toml'))
    p.add_argument('--stage',choices=['prepare','main','bootstrap','checks'],required=True)
    p.add_argument('--workers',type=int,default=6)
    a=p.parse_args();run(a.previous,a.output,a.protocol,a.stage,a.workers)
