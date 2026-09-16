"""Read-only numerical review of the supplied C20 v4 packet.

Usage: python check_delivery.py /path/to/extracted_packet /path/to/new_checks
Requires NumPy, pandas, SciPy. Does NOT run MATLAB or refit experiment spectra.
All MATLAB indices are explicitly converted to zero-based indices.
MATLAB compact MAT storage is loaded using mat_dtype=True; numerical arrays
are additionally cast to float before powers/summations.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat


def load(path: Path) -> dict:
    return loadmat(path, simplify_cells=True, mat_dtype=True)


def peak(E, par, model):
    E = np.asarray(E, dtype=float)
    E0, w, A = np.asarray(par, dtype=float)
    if model == 'lorentz_symmetric':
        return A / np.pi * (w / 2) / ((E - E0) ** 2 + (w / 2) ** 2)
    if model == 'lorentz':
        return A * E * w / ((E ** 2 - E0 ** 2) ** 2 + E ** 2 * w ** 2)
    raise ValueError(model)


def reconstruct(E, p, n, scale, ampunit, model):
    E = np.asarray(E, dtype=float)
    p = np.asarray(p, dtype=float)
    n = int(n)
    pars = p[2:2 + 3*n].reshape(n, 3) * [1000, 1000, scale * ampunit]
    comp = np.column_stack([peak(E, pp, model) for pp in pars])
    bg = scale * p[0] * (E / 1000) ** (-p[1])
    if len(p) > 2 + 3*n:
        bg = bg + scale * p[-1]
    return bg, comp, bg + comp.sum(axis=1), pars


def emit(df, output, name):
    df.to_csv(output / name, index=False, float_format='%.17g')


def main(root: Path, output: Path):
    root = root.resolve()
    output = output.resolve()
    if output == root or output.is_relative_to(root):
        raise ValueError('Audit output must be outside the delivered packet.')
    output.mkdir(parents=True, exist_ok=True)
    summary = {'method': 'Independent Python array reconstruction; no MATLAB execution, no refit',
               'source': str(root)}
    mf = pd.read_csv(root / 'FILE_MANIFEST.csv')
    rows = []
    for r in mf.itertuples():
        p = root / r.path.replace('\\', '/').strip()
        assert p.resolve().is_relative_to(root)
        rows.append([r.path, p.stat().st_size == r.bytes,
                     hashlib.sha256(p.read_bytes()).hexdigest() == r.sha256])
    manifest = pd.DataFrame(rows, columns=['path', 'bytes_match', 'sha256_match'])
    assert manifest.bytes_match.all() and manifest.sha256_match.all()
    emit(manifest, output, 'manifest_checks.csv')
    summary['manifest_files_verified'] = len(manifest)

    b = pd.read_csv(root / 'audit/boundary_by_parameter_native.csv')
    c = pd.read_csv(root / 'audit/candidate_solution_families.csv')
    rows = []
    for (key, n, start), g in b.groupby(['key', 'n_components', 'start'], sort=False):
        n = int(n)
        g = g.sort_values('raw_parameter_index')
        p, lb, ub = (g[x].to_numpy(float) for x in ['raw_value', 'raw_lb', 'raw_ub'])
        lo = np.isfinite(p) & (abs(p-lb) < 1e-5*np.maximum(1, abs(lb)))
        iu = np.isfinite(ub)
        hi = np.zeros(len(p), dtype=bool)
        hi[iu] = np.isfinite(p[iu]) & (abs(p[iu]-ub[iu]) < 1e-5*np.maximum(1, abs(ub[iu])))
        lo[~iu] = np.isfinite(p[~iu]) & (abs(p[~iu]-lb[~iu]) < 1e-5)
        pars = p[2:2+3*n].reshape(n,3)
        order = np.argsort(pars[:,0], kind='stable')
        inverse = np.zeros(n, dtype=int); inverse[order] = np.arange(1,n+1)
        mapped = np.array_equal(np.repeat(inverse,3), g.ordered_component.to_numpy()[2:2+3*n])
        flag = np.array_equal(lo | hi, g.raw_boundary.to_numpy(bool))
        flag &= np.array_equal(lo, g.lower_hit.to_numpy(bool)) and np.array_equal(hi, g.upper_hit.to_numpy(bool))
        native = g.native_value.to_numpy(float)
        err = np.max(abs(native-p*g.native_factor.to_numpy(float))/np.maximum(1,abs(native)))
        assert mapped and flag and err < 5e-12
        rows.append([key,n,start,mapped,flag,err])
    emit(pd.DataFrame(rows, columns=['key','n','start','mapping_ok','flags_ok','native_relative_error']),
         output,'parent_mapping_checks.csv')
    hits = b[(b.selected == 1) & ((b.lower_hit == 1) | (b.upper_hit == 1))]
    emit(hits, output, 'selected_parent_boundary_hits.csv')
    summary['parent_mapping'] = {'candidates':len(rows),'parameter_rows':len(b),
        'selected_rows':int(b.selected.sum()),'boundary_types':hits.groupby('boundary_type').size().to_dict(),
        'raw_order_swapped_candidates':int((c.order == '[2;1]').sum()),
        'all_parent_objectives_independently_recomputed':False,
        'reason':'Packet has all parent parameter rows/mappings but not every original parent observation array.'}

    sets = load(root / '590_PL2_10w/solver_background_comparison.mat')
    fitrows, candrows, witnessrows = [], [], []
    for name in ['legacy','independent','background']:
        for d in sets[name]:
            for f in d['fits']:
                n = int(f['n_components']); E = f['energy_meV']; Y = f['observed']; scale = f['scale']
                for k, cd in enumerate(f['candidates'],1):
                    if not np.isfinite(cd['p']).all():
                        continue
                    _,_,pred,_ = reconstruct(E,cd['p'],n,scale,f['ampunit'],f['peak_model'])
                    Q = np.sum(((Y-pred)/scale)**2)
                    err = abs(Q-cd['objective'])
                    assert err < 1e-10*max(1,Q)
                    candrows.append([name,d['key'],n,k,Q,cd['objective'],err,cd['exitflag']])
                selected = f['candidates'][int(f['selected_start'])-1]
                bg,comp,pred,pars = reconstruct(E,selected['p'],n,scale,f['ampunit'],f['peak_model'])
                order=np.argsort(pars[:,0],kind='stable'); pars=pars[order]
                savedpars=np.asarray(f['parameters']).reshape(n,3)
                savedcomp=np.asarray(f['components']).reshape(len(E),n)
                parerr=np.max(abs(pars-savedpars)/np.maximum(1,abs(savedpars)))
                cerr=np.max(abs(comp[:,order]-savedcomp))/scale
                perr=np.max(abs(pred-f['prediction']))/scale
                rerr=np.max(abs((Y-pred)-f['residual']))/scale
                qerr=abs(np.sum(((Y-pred)/scale)**2)-f['normalized_sse'])
                orderok=np.array_equal(np.atleast_1d(f['component_order']),order+1)
                feasible=[x['objective'] for x in f['candidates'] if x['exitflag']>0 and np.isfinite(x['objective'])]
                selectok=abs(selected['objective']-min(feasible))<1e-12
                assert max(parerr,cerr,perr,rerr,qerr)<1e-9 and orderok and selectok
                fitrows.append([name,d['key'],n,parerr,cerr,perr,rerr,qerr,orderok,selectok,selected['firstorderopt']])
                w=f['witness']
                if w['available']:
                    _,_,yp,_=reconstruct(E,w['p'],n,scale,f['ampunit'],f['peak_model'])
                    qw=np.sum(((Y-yp)/scale)**2)
                    assert abs(qw-w['h0_objective'])<1e-10
                    assert np.all(w['p']>=f['lb']) and np.all(w['p']<=f['ub'])
                    witnessrows.append([name,d['key'],qw,w['h0_objective'],f['normalized_sse']])
    fitdf=pd.DataFrame(fitrows,columns=['set','key','n','parameter_error','component_error','prediction_error','residual_error','objective_error','order_ok','selection_ok','firstorderopt'])
    canddf=pd.DataFrame(candrows,columns=['set','key','n','start','Q_recomputed','Q_saved','absolute_error','exitflag'])
    emit(fitdf,output,'real_fit_checks.csv'); emit(canddf,output,'real_candidate_checks.csv')
    emit(pd.DataFrame(witnessrows,columns=['set','key','witness_Q','H0_Q','H1_optimized_Q']),output,'witness_checks.csv')
    summary['real_reconstruction']={'models':len(fitdf),'components':int(fitdf.n.sum()),'candidates':len(canddf),
        'max_candidate_objective_error':float(canddf.absolute_error.max()),'witnesses':len(witnessrows),
        'max_curve_relative_error':float(fitdf.component_error.max()),'candidate_exitflags':canddf.exitflag.value_counts().to_dict()}

    # Audit a genuine scope-label bug: exported Wref field integrates the whole fit window.
    pp=pd.read_csv(root/'audit/component_parameters.csv'); area_rows=[]
    for r in pp.itertuples():
        hi=int(re.search(r'_W(\d+)_',r.key).group(1))
        Eref=np.arange(300,1801,4,dtype=float); Efit=np.arange(300,hi+1,4,dtype=float)
        par=[r.E0_meV,r.native_width_meV,r.native_A]
        ar=np.trapezoid(peak(Eref,par,r.peak_model),Eref)
        af=np.trapezoid(peak(Efit,par,r.peak_model),Efit)
        assert abs(af-r.area_Wref_300_1800)<1e-10*max(1,abs(af))
        area_rows.append([r.key,r.n_components,r.ordered_component,hi,r.area_Wref_300_1800,ar,af,(r.area_Wref_300_1800/ar-1)*100])
    area=pd.DataFrame(area_rows,columns=['key','n','component','fit_high_meV','saved_Wref_field','correct_Wref_300_1800','fit_window_area','relative_overstatement_percent'])
    emit(area,output,'parent_area_scope_audit.csv')
    summary['area_scope_bug']={'affected_parent_component_rows':int((area.fit_high_meV!=1800).sum()),
       'current_300_1800_fits_affected':False,'max_percent':float(area.relative_overstatement_percent.max())}

    bins=load(root/'590_PL2_10w/centered_binned_spectra.mat')
    fr=load(root/'590_PL2_10w/sequence_block_spectra.mat')['frames']
    supp=np.asarray(fr['support'],int)-1; offs=np.asarray(fr['alignment']['measured_offset_pixels'],int)
    X=np.asarray(fr['member_spectra_A0'],float); XA=np.asarray(fr['member_spectra_A1'],float)
    assert np.array_equal(X.sum(axis=1),bins['member_spectra'])
    max_alignment_error=0.
    for t,off in enumerate(offs):
        max_alignment_error=max(max_alignment_error,float(np.max(abs(X[supp+off,t,:]-XA[:,t,:]))))
    assert max_alignment_error==0
    bin_rows=[]
    for bb in bins['bins']:
        ix=[list(bins['members']).index(x) for x in np.atleast_1d(bb['source_channel'])]
        ss=np.asarray(bins['member_spectra'],float)[:,ix].sum(axis=1)
        err1=np.max(abs(ss-bb['sum'])); err2=np.max(abs(ss/bb['N']-bb['mean']))
        assert err1<1e-9 and err2<1e-9
        bin_rows.append([bb['target_q'],bb['N'],bb['q_Ainv'],bb['q_width'],err1,err2])
    emit(pd.DataFrame(bin_rows,columns=['target_q','N','actual_q','q_width','sum_error','mean_error']),output,'bin_checks.csv')
    block_rows=[]
    for j,memb in enumerate(fr['n3_members']):
        ix=[list(fr['native_members']).index(x) for x in memb]
        s0=X[supp,:,:][:,:,ix].mean(axis=2); s1=XA[:,:,ix].mean(axis=2)
        for block in range(1,7):
            use=fr['block_id']==block
            er0=np.max(abs(s0[:,use].sum(axis=1)-fr['block_sum_A0'][:,block-1,j]))
            er1=np.max(abs(s1[:,use].sum(axis=1)-fr['block_sum_A1'][:,block-1,j]))
            erm=np.max(abs(fr['block_per_frame_mean_A0'][:,block-1,j]*use.sum()-fr['block_sum_A0'][:,block-1,j]))
            assert max(er0,er1,erm)<1e-8
            block_rows.append([j+1,block,er0,er1,erm])
    emit(pd.DataFrame(block_rows,columns=['target','block','A0_error','A1_error','time_mean_scale_error']),output,'block_checks.csv')
    pd.crosstab(fr['block_id'],fr['q_center_channel']).to_csv(output/'q_peak_channel_by_block.csv')
    mask=(fr['E']>=300)&(fr['E']<=1800); E=np.asarray(fr['E'][mask],float); frrows=[]
    for j,q in enumerate(fr['targets']):
        y0=fr['sequence_bin_A0'][mask,:,j].sum(axis=1); y1=fr['sequence_bin_A1'][mask,:,j].sum(axis=1)
        a0=np.trapezoid(y0,E); a1=np.trapezoid(y1,E)
        ab=np.trapezoid(fr['block_per_frame_mean_A0'][mask,:,j],E,axis=0)
        frrows.append([q,(a1/a0-1)*100,np.linalg.norm(y1-y0)/np.linalg.norm(y0),
                       (ab[-1]/ab[0]-1)*100,ab.min(),ab.max()])
    emit(pd.DataFrame(frrows,columns=['q','A1_area_change_percent','A1_A0_relative_L2','block6_vs1_area_percent','block_area_min','block_area_max']),output,'frame_derived_diagnostics.csv')
    summary['frames']={'A1_member_error':max_alignment_error,'L1_selected_members_recomputed':True,
       'common_support_matlab_indices':[int(supp[0]+1),int(supp[-1]+1)],'blocks_recomputed':len(block_rows),
       'q_center_channel_counts':dict(zip(*[list(map(int,x)) for x in np.unique(fr['q_center_channel'],return_counts=True)])),
       'raw_NPY_JSON_not_in_packet':True}

    simrows=[]; recovery=[]; sim_candidates=0; sim_maxerr=0.
    for phase in (1,2):
        for scenario in ['static_single','q_mixed_single','overlapping_double','background_mismatch']:
            data=load(root/f'validation/simulation_phase{phase}_{scenario}.mat')
            for it,t in enumerate(data['trials'],1):
                s=t['simulation']; fs=t['fits']; E=np.asarray(s['E'],float)
                for f in fs:
                    for cd in f['candidates']:
                        if not np.isfinite(cd['p']).all(): continue
                        _,_,pred,_=reconstruct(E,cd['p'],f['n_components'],f['scale'],f['ampunit'],f['peak_model'])
                        Q=np.sum(((f['observed']-pred)/f['scale'])**2)
                        err=abs(Q-cd['objective']); sim_maxerr=max(sim_maxerr,float(err));sim_candidates+=1
                        assert err<1e-9*max(1,Q)
                f=fs[1]; par=np.asarray(f['parameters']).reshape(2,3)
                ar=np.trapezoid(f['components'],E,axis=0); fraction=ar[1]/ar.sum()
                gain=1-f['normalized_sse']/fs[0]['normalized_sse']
                split=gain>.1 and min(fraction,1-fraction)>.05 and par[1,0]-par[0,0]>=4
                simrows.append([phase,scenario,it,s['seed'],gain,split,*par[:,0],*par[:,1],fraction])
                if phase==2 and scenario=='overlapping_double':
                    truth=np.asarray(s['truth_parameters']); tar=np.array([np.trapezoid(peak(E,p,'lorentz_symmetric'),E) for p in truth])
                    recovery.append([it,*((par-truth).flatten()),fraction-tar[1]/tar.sum()])
    sr=pd.DataFrame(simrows,columns=['phase','scenario','trial','seed','gain','split','E1','E2','width1','width2','P2_area_fraction'])
    supplied=pd.read_csv(root/'validation/simulation_trials.csv')
    joined=sr.merge(supplied,on=['phase','scenario','trial','seed'])
    assert len(joined)==440 and (joined.split==joined.engineering_split.astype(bool)).all()
    emit(sr,output,'simulation_checks.csv')
    rr=pd.DataFrame(recovery,columns=['trial','err_E1','err_width1','err_A1','err_E2','err_width2','err_A2','err_P2_fraction'])
    emit(rr,output,'hypothetical_double_recovery_errors.csv')
    summary['simulations']={'trials':len(sr),'candidates_recomputed':sim_candidates,'max_Q_error':sim_maxerr,
       'all_flags_match':True,'engineering_only':True,
       'pilot_split_counts':sr[sr.phase==2].groupby('scenario').split.sum().astype(int).to_dict(),
       'double_95pct_abs_parameter_error':rr.drop(columns='trial').abs().quantile(.95).to_dict()}

    profiles=load(root/'validation/profiles/profile_arrays.mat')['profiles']; pr=[]
    for v in profiles:
        f=v['fit']; E=np.asarray(f['energy_meV'],float)
        for k,t in enumerate(v['trials'],1):
            _,co,pred,pa=reconstruct(E,t['p'],2,f['scale'],f['ampunit'],f['peak_model'])
            Q=np.sum(((f['observed']-pred)/f['scale'])**2)
            if v['variable']=='separation_meV': value=pa[1,0]-pa[0,0]
            elif v['variable']=='P1_width_meV': value=pa[0,1]
            else:
                ar=np.trapezoid(co,E,axis=0); value=ar[1]/ar.sum()
            pr.append([v['q'],v['variable'],v['value'],k,abs(Q-t['Q']),value-v['value'],t['exitflag']])
    pcheck=pd.DataFrame(pr,columns=['q','parameter','fixed','start','Q_error','constraint_error','exitflag'])
    assert pcheck.Q_error.max()<1e-10
    assert abs(pcheck[pcheck.exitflag>0].constraint_error).max()<1e-8
    emit(pcheck,output,'profile_checks.csv')
    summary['profiles']={'fixed_parameter_points':len(profiles),'trials':len(pcheck),
         'failed_trials_retained':int((pcheck.exitflag<=0).sum()),'not_confidence_intervals':True}
    (output/'CHECK_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('packet',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();main(a.packet,a.output)
