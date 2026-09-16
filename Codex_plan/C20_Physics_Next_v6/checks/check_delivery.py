from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.io import loadmat

base=Path('/mnt/data/C20_v5_delivery')
out=Path('/mnt/data/C20_Physics_Next_v6/checks')
summary={'scope':'Targeted independent readback of supplied v5 selected model arrays; no MATLAB re-optimization or access to external NPY.'}
manifest=pd.read_csv(base/'FILE_MANIFEST.csv')
print(manifest.columns.tolist())
errors=[]
for r in manifest.to_dict('records'):
    name=str(r.get('path',r.get('file',''))).replace('\\','/')
    f=base/name
    if not f.is_file(): errors.append({'path':name,'reason':'missing'});continue
    h=hashlib.sha256(f.read_bytes()).hexdigest()
    if h.lower()!=str(r['sha256']).lower():errors.append({'path':name,'reason':'hash'})
summary['manifest']={'records':len(manifest),'mismatches':errors}

def peak(E,c,w,A,m):
    E=np.asarray(E,dtype=np.float64)
    if m=='lorentz_symmetric':return A*(w/2)/np.pi/((E-c)**2+(w/2)**2)
    return A*E*w/((E**2-c**2)**2+E**2*w**2)
rows=[]
A1=loadmat(base/'A0_A1/fit_details.mat',simplify_cells=True)['A1']
for d in A1:
    for f in d['fits']:
        p=np.asarray(f['parameters']).reshape(f['n_components'],3)
        E=f['energy_meV']; Y=f['observed'];
        cs=np.column_stack([peak(E,*z,f['peak_model']) for z in p])
        ref=np.asarray(f['components']).reshape(len(E),-1)
        y=f['background']+cs.sum(axis=1)
        obj=np.sum((Y-y)**2)/f['scale']**2
        rows.append({'key':d['key'],'n':f['n_components'],'component_relative_error':np.max(abs(cs-ref))/max(1,np.max(abs(ref))), 'prediction_relative_error':np.max(abs(y-f['prediction']))/max(1,np.max(abs(Y))), 'objective_abs_error':abs(obj-f['normalized_sse'])})
summary['A1_selected']={'models':len(rows),'max_component_relative_error':max(r['component_relative_error'] for r in rows),'max_objective_abs_error':max(r['objective_abs_error'] for r in rows)}
pd.DataFrame(rows).to_csv(out/'A1_checks.csv',index=False)
rows=[]
M=loadmat(base/'member_models/all_member_fits.mat',simplify_cells=True)['allmembers']
for d in M:
    for f in d['fits']:
        E=f['E'];q=f['q'];n=int(f['n_components']); Q=len(q);p=f['p'];beta=(q-q.min())/np.ptp(q)
        bg=(E/1000)**(-p[0]);B=bg[:,None]*p[1:Q+1][None,:]
        C=np.zeros((len(E),Q,n))
        for j in range(n):
            b=1+Q+j*(3+Q);center=((1-beta)*p[b]+beta*p[b+1])*1000;w=p[b+2]*1000
            for k in range(Q):C[:,k,j]=peak(E,center[k],w,p[b+3+k]*f['ampunit'],f['peak_model'])
        order=np.atleast_1d(f['order']).astype(int)-1;C=C[:,:,order]*f['scale'];Y=B*f['scale']+C.sum(axis=2)
        ref=np.asarray(f['components']).reshape(len(E),Q,n)
        obj=np.sum((Y-f['observed'])**2)/f['scale']**2
        rows.append({'key':d['key'],'n':n,'prediction_relative_error':np.max(abs(Y-f['prediction']))/max(1,np.max(abs(f['observed']))), 'component_relative_error':np.max(abs(C-ref))/max(1,np.max(abs(ref))), 'objective_abs_error':abs(obj-f['objective'])})
summary['member_selected']={'models':len(rows),'max_prediction_relative_error':max(r['prediction_relative_error'] for r in rows),'max_objective_abs_error':max(r['objective_abs_error'] for r in rows)}
pd.DataFrame(rows).to_csv(out/'member_checks.csv',index=False)
area=loadmat(base/'audit/area_arrays.mat',simplify_cells=True)['area_records'];maxerr=0;count=0
for r in area:
    E=r['E'];C=np.asarray(r['components']).reshape(len(E),-1);mask=(E>=300)&(E<=1800)
    val=np.trapezoid(C[mask],E[mask],axis=0);ref=np.atleast_1d(r['areas']['area_reference_window'])
    maxerr=max(maxerr,float(np.max(abs(val-ref)/np.maximum(1,abs(ref)))))
    count+=len(val)
summary['reference_area_check']={'component_rows':count,'max_relative_error':maxerr}
changes=pd.read_csv(base/'A0_A1/parameter_changes.csv');d=changes[changes.n==2]
summary['A1_effect']={'max_abs_native_center_shift_meV':float(d.delta_E0.abs().max()),'max_abs_fraction_change':float(abs(d.fraction_A1-d.fraction_A0).max())}
summary['acquisition_status']='Delivery records direct user confirmation: 300 continuous frames, same region, no deliberate scan/region change/beam adjustment. This does not establish independence/stationarity.'
(out/'CHECK_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False,indent=2))
# Physics-oriented, conditional candidate table on three A1 bins, NOT a full dispersion.
rows=[]
for di,d in enumerate(A1):
    f=d['fits'][1];E=f['energy_meV']; C=np.asarray(f['components']);areas=np.trapezoid(C,E,axis=0)
    for j in range(2):
        r=d['unit'];
        rows.append({'region':di//2+1,'q_Ainv':r['q_Ainv'],'domain':'A1, mean q N3, sum 300 frames','model':f['peak_model'],'component':j+1,'native_E0_meV':f['parameters'][j,0],'model_apex_on_4meV_grid_meV':float(E[np.argmax(C[:,j])]),'native_width_meV':f['parameters'][j,1],'area_W300_1800':areas[j],'fraction_W300_1800':areas[j]/sum(areas),'status':'conditional effective spectral component; no CI or quasiparticle assignment'})
pd.DataFrame(rows).to_csv(out.parent/'CURRENT_PHYSICS_CANDIDATES.csv',index=False)
