"""Physical limits and independent quadrature checks for the B1 comparison."""
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'case_studies/bisb2026/exploration/composition'))
from models import Operator,oscillator,lognormal_density,fit


def example_data(q=None):
    native=np.arange(200.,2792.,4.)
    return dict(energy_native=native,energy=native.reshape(-1,4).mean(1),
        q=np.array([-.005,.0025,.005,.010]) if q is None else np.array(q),
        width=.0015,h=.0066,probe=.0004,beam_kV=30.,resolution_fwhm=22.)


@pytest.fixture(scope='module')
def op():return Operator(example_data())


def test_oscillator_energy_sum_rule():
    for e0,gamma in [(.7,.6),(1.4,1.5),(.4,.08)]:
        val,_=quad(lambda e:e*oscillator(e,e0,gamma),0,np.inf,epsabs=1e-10)
        assert val==pytest.approx(np.pi/2,rel=1e-8)


def test_continuum_is_positive_normalized_and_unimodal():
    for mode,width in [(.7,.4),(1.3,.2),(2.,1.2)]:
        val,_=quad(lambda z:lognormal_density(np.exp(z),mode,width)*np.exp(z),-20,20)
        assert val==pytest.approx(1.,rel=1e-7)
        e=np.geomspace(mode/4,mode*4,10001);s=lognormal_density(e,mode,width)
        assert np.all(s>0)
        assert e[np.argmax(s)]==pytest.approx(mode,rel=.001)
        assert np.count_nonzero(np.diff(np.sign(np.diff(s))))==1


def test_projection_conserves_constant_counts_and_positive_response(op):
    z=op.project(np.ones((len(op.q),len(op.e))))
    np.testing.assert_allclose(z,4.,rtol=1e-14)
    assert np.all(op.low([.9,.7,.6,0])>=0)


def test_even_q_geometry(op):
    # ±q are geometrically identical; anisotropy is not silently introduced.
    z=op.low([.9,.7,.6,.1]);np.testing.assert_allclose(z[0],z[2],rtol=1e-12)


def test_point_limit_matches_direct_spectrum():
    d=example_data([.005]);p=Operator(d,geometry='point',n_radial=512)
    e0=.9;g=.6
    direct=p.project(p.K*oscillator(p.e,e0,g)*e0**2)
    actual=p.low([e0,.7,g,0])
    np.testing.assert_allclose(actual/actual.max(),direct/direct.max(),rtol=2e-4)


def test_aperture_and_radial_quadrature_converge(op):
    fine=Operator(example_data(),n_radial=512,n_quad=45,n_probe=11)
    a=op.low([.85,.65,.5,.2]);b=fine.low([.85,.65,.5,.2])
    a/=a.max(1,keepdims=True);b/=b.max(1,keepdims=True)
    assert np.max(abs(a-b))<.003


def test_full_aperture_cannot_be_replaced_by_kinematic_average(op):
    full=op.low([.9,.2,.2,0])
    # A narrow, dispersive resonance is broadened by Q acceptance. Its shape
    # differs from a local pole times the exactly averaged kinematic factor.
    x=abs(op.q)/.005;ep2=.9**2*x/(1+.2*(x-1))
    local=op.project(op.K*oscillator(op.e[None,:],np.sqrt(ep2)[:,None],.2)*ep2[:,None])
    full/=full.max(1,keepdims=True);local/=local.max(1,keepdims=True)
    assert np.max(abs(full[1]-local[1]))>.1


def test_noiseless_recovery_with_independent_fine_operator(op):
    fine=Operator(example_data(),n_radial=512,n_quad=45,n_probe=11)
    pars=np.array([.9,.55,.5,.1,2.05,.25])
    a,_=fine.columns(pars,'single')
    coeff=np.broadcast_to([2000.,150.,5.,2.],(len(op.q),4))
    y=np.einsum('bek,bk->be',a,coeff)
    mask=(op.energy>=500)&(op.energy<=2600)
    r=fit(op,y,np.sqrt(np.maximum(y,1)),'single',mask,n_starts=6)
    assert r['success']
    assert r['mean_squared_residual']<.003
    assert r['parameters'][0]==pytest.approx(.9,abs=.015)
