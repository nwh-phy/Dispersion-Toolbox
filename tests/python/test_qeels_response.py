import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'/'python'))
from qeels_response import longitudinal_q, kinematic_2d, dl_loss, acceptance_nodes, averaged_kinematic, measured_loss
from qeels_response import screened_plasmon_energy, dispersive_aperture_loss


def test_relativistic_energy_to_momentum_conversion():
    assert abs(longitudinal_q(1000,30)-.001543)<.000005
    assert longitudinal_q(1000,80)<longitudinal_q(1000,30)


def test_2d_kinematic_limits_and_symmetry():
    assert kinematic_2d(0.,1000.)==0.
    q=.2
    assert abs(kinematic_2d(q,1000)*q**3-1)<.001
    e=np.array([300.,1000.,1800.])
    np.testing.assert_allclose(averaged_kinematic(.005,e,dq=.0005,h=.0066),averaged_kinematic(-.005,e,dq=.0005,h=.0066),rtol=1e-12)


def test_dl_energy_weighted_sum_rule():
    for e0,g in [(700.,500.),(1200.,100.),(700.,1200.)]:
        value,_=quad(lambda x:x*dl_loss(x,e0,g),0,np.inf,epsabs=1e-7)
        assert abs(value-np.pi/2)<1e-6


def test_acceptance_constant_response_and_quadrature_convergence():
    q,w=acceptance_nodes(.0025,dq=.0005,h=.0066,sigma_probe=.0003)
    assert abs(w.sum()-1)<1e-14 and np.all(w>=0) and np.all(q>=0)
    e=np.linspace(300,1800,101)
    a=averaged_kinematic(.0025,e,dq=.0005,h=.0066,n=31)
    b=averaged_kinematic(.0025,e,dq=.0005,h=.0066,n=61)
    assert np.max(abs(a/b-1))<1e-5


def test_response_positive_linear_and_resolution_area_conserving():
    e=np.arange(-2000,10004,4.)
    args=dict(q_center=.005,components=[(700.,500.,1.)],h=.0066,dq=.0005)
    a=measured_loss(e,**args)
    b=measured_loss(e,**args,resolution_fwhm_meV=24)
    assert np.all(a>=0) and np.all(b>=0)
    assert abs(np.sum(b)/np.sum(a)-1)<1e-5
    doubled=measured_loss(e,.005,[(700.,500.,2.)],h=.0066,dq=.0005)
    np.testing.assert_allclose(doubled,2*a)


def test_screened_drude_limits_and_aperture_quadrature():
    A=1e8;r0=200.
    assert screened_plasmon_energy(0,A,r0)==0
    assert abs(screened_plasmon_energy(1e-8,A,r0)**2/(A*1e-8)-1)<3e-6
    assert abs(screened_plasmon_energy(1e4,A,r0)**2/(A/r0)-1)<1e-6
    np.testing.assert_allclose(screened_plasmon_energy(np.array([.001,.004]),A,0),np.sqrt(A*np.array([.001,.004])))
    e=np.arange(100.,2000.,10.)
    a=dispersive_aperture_loss(e,.0025,A,r0,400.,dq=.0005,h=.0066,n=35)
    b=dispersive_aperture_loss(e,.0025,A,r0,400.,dq=.0005,h=.0066,n=71)
    assert np.all(a>=0)
    assert np.max(abs(a/b-1))<1e-5
