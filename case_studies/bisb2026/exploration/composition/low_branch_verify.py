"""Readback, source integrity, aperture convergence, and compact count export."""
import argparse
from pathlib import Path
import numpy as np
from low_branch_run import read,local_data,quantities,BOOT_NAMES
from models import Operator
from qbin_fit import local_columns
from io_utils import digest,dump


def verify(output,ridge):
    signature=read(output/'signature.json');cfg=signature['protocol']
    for name,sha in signature['code'].items():assert digest(Path(__file__).with_name(name))==sha
    tasks=read(output/'tasks.json');records=[];fine=[];raw=[];all_counts={};n_models=0
    for n in [3,5]:
        manifest=read(output/f'bin{n}/manifest.json')
        for rec in manifest['sessions']:
            sid=rec['session'];path=output/f'bin{n}/{sid}.npz'
            assert digest(path)==rec['prepared_sha256']
            if n==3:
                sha=digest(rec['source']);assert sha==rec['source_sha256']
                raw.append(dict(session=sid,sha256=sha,unchanged=True))
            d=dict(np.load(path));old=dict(np.load(ridge/f'bin{n}/{sid}.npz'))
            np.testing.assert_array_equal(d['q'],old['q'])
            ie,je=np.where(np.isclose(d['energy'][:,None],old['energy'][None,:],rtol=0,atol=1e-8))
            assert len(ie)>100
            np.testing.assert_array_equal(d['frames'][...,ie],old['frames'][...,je])
            records.append(dict(session=sid,channels=n,common_energy_channels=len(ie),all_frames_exact=True))
            nt=cfg['train_frames'];prefix=f'{sid}_bin{n}_'
            for k in ['q','width','energy','sigma_train','sigma_test','beam']:
                all_counts[prefix+k]=d[k]
            all_counts[prefix+'train_counts']=d['frames'][:,:nt].astype(float).sum(1)
            all_counts[prefix+'test_counts']=d['frames'][:,nt:].astype(float).sum(1)
    convergence_keys=set()
    for sid in ['590','n0','20w']:
        eligible=[t for t in tasks if t['session']==sid and (t['native'] or 'bin5' in t['key'])]
        for target in [.0005,.001,.005,.015]:
            convergence_keys.add(min([t for t in eligible if t['q_Ainv']>0],key=lambda t:abs(t['q_Ainv']-target))['key'])
    booleans=[]
    for t in tasks:
        d=local_data(t);op=Operator(d);dest=output/'fits'/t['key'];f=read(dest/'main.json')
        mask=(op.energy>=cfg['fit_window_meV'][0])&(op.energy<=cfg['fit_window_meV'][1])
        for model,x in f.items():
            col,_=local_columns(op,np.array(x['parameters']),model)
            pred=np.einsum('bek,bk->be',col,np.array(x['coefficients']))
            np.testing.assert_allclose(pred,x['prediction'],rtol=1e-11,atol=1e-8)
            residual=(pred-d['frames'][:,:cfg['train_frames']].astype(float).sum(1))/d['sigma_train']
            np.testing.assert_allclose((residual[:,mask]**2).sum(),x['cost'],rtol=1e-10)
            q=quantities(op,x)
            for k in ['low_count_apex_eV','low_B1_area_fraction']:
                np.testing.assert_allclose(q[k],x[k],atol=1e-12)
            n_models+=1
        assert f['A']['cost']<=f['S']['cost']+1e-3
        z=np.load(dest/'bootstrap.npz');assert z['names'].tolist()==BOOT_NAMES
        assert z['values'].shape==(cfg['bootstrap_repeats'],len(BOOT_NAMES))
        assert np.all(np.isfinite(z['values']))
        booleans.append(dict(key=t['key'],success_fraction=float(z['values'][:,10].mean())))
        if t['key'] in convergence_keys:
            opfine=Operator(d,n_radial=512,n_quad=45,n_probe=11)
            for model,x in f.items():
                col,_=local_columns(opfine,np.array(x['parameters']),model)
                pred=np.einsum('bek,bk->be',col,np.array(x['coefficients']))
                band=(op.energy>=300)&(op.energy<=1800)
                rms=float(np.sqrt(np.mean(((pred-np.asarray(x['prediction']))/d['sigma_train'])[:,band]**2)))
                assert rms<.25
                fine.append(dict(key=t['key'],model=model,rms_sigma=rms))
    for sid in ['590','n0','20w']:
        a=dict(np.load(output/f'bin3/{sid}.npz'));b=dict(np.load(output/f'bin5/{sid}.npz'))
        ia=abs(a['q'])<=.006+1e-10;ib=abs(b['q'])<=.006+1e-10
        np.testing.assert_array_equal(a['frames'][ia],b['frames'][ib])
    np.savez_compressed(output/'measured_spectra.npz',**all_counts)
    dump(output/'verification.json',dict(main_models_reconstructed=n_models,
        independent_q_bins=len(tasks),bootstrap_refits=2*cfg['bootstrap_repeats']*len(tasks),
        energy_and_q_parity=records,raw_sources=raw,convergence=fine,
        bootstrap_convergence=booleans,code_sha256=digest(Path(__file__))))
    print('Verified',len(tasks),'q bins,',n_models,'main fits; source hashes unchanged;',
          'fine operator max RMS sigma',max(r['rms_sigma'] for r in fine))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ridge',type=Path,required=True);a=p.parse_args();verify(a.output,a.ridge)
