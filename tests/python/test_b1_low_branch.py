"""Independent low-component recovery and energy-semantics checks."""
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'case_studies/bisb2026/exploration/composition'))
from models import Operator
from dispersion_models import encode
from qbin_fit import fit_local,local_columns
from low_branch_run import quantities,count_apex


def data(q):
    e=np.arange(200.,3400.,4.)
    return dict(energy_native=e,energy=e.reshape(-1,4).mean(1),q=np.array([q]),
        width=.0005 if q<.006 else .0025,h=.0066,probe=.0004,beam_kV=30.,resolution_fwhm=22.)


def test_local_low_branch_can_rise_or_fall_without_q_prior():
    for q,e0 in [(.001,.72),(.005,1.02),(.015,.82)]:
        d=data(q);fine=Operator(d,n_radial=512,n_quad=45,n_probe=11);op=Operator(d)
        p=encode([e0]*2,[1.5]*2,'A');p[[4,6]]=[.45,.8];p[[5,7,12]]=0.
        col,_=local_columns(fine,p,'A')
        y=np.einsum('bek,bk->be',col,np.array([[1200.,1800.,400.,500.,10.,3.]]))
        mask=op.energy>=300;sigma=np.sqrt(y)
        s=fit_local(op,y,sigma,'S',mask,starts=3)
        a=fit_local(op,y,sigma,'A',mask,starts=6,single=s)
        assert abs(a['E0_low_eV']-e0)<.015
        assert abs(a['E0_high_eV']-1.5)<.015
        assert a['cost']<s['cost'] and a['success']
        warm=fit_local(op,y,sigma,'A',mask,starts=4,single=s,warms=[a['parameters']])
        assert warm['cost']<=a['cost']+1e-3


def test_reported_component_apex_is_not_oscillator_E0():
    op=Operator(data(.001));p=encode([.8]*2,[1.45]*2,'A');p[[4,6]]=[.6,1.];p[[5,7,12]]=0
    col,_=local_columns(op,p,'A');comp=col*np.array([1000.,1500.,500.,600.,10.,2.])[None,None,:]
    f=dict(model='A',components=comp,E0_low_eV=.8,E0_high_eV=1.45)
    v=quantities(op,f)
    assert .2<v['low_count_apex_eV']<.8
    assert 0<v['low_B1_area_fraction']<1
    assert not v['count_apex_order_inversion']


def test_count_apex_marks_censoring_at_energy_boundary():
    e=np.arange(206.,3400.,16.)
    apex,edge=count_apex(e,e**-3)
    assert edge and apex==e[0]/1000
