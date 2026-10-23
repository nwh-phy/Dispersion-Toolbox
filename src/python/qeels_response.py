"""Minimal quasistatic 2D measurement model, in meV and inverse angstrom.

This is a mechanism/sensitivity model, not a Bi/BiSb dielectric prediction.
Normal incidence; isotropic in-plane response; no retarded or multiple-scattering terms.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss

HBAR_EV_S = 6.582119569e-16
C_M_S = 299792458.0


def longitudinal_q(energy_meV, beam_kV=30.):
    if beam_kV <= 0: raise ValueError('Positive beam energy required')
    gamma = 1 + beam_kV / 510.99895
    v = C_M_S*np.sqrt(1-gamma**-2)
    return np.asarray(energy_meV)*1e-3/HBAR_EV_S/v*1e-10


def dl_loss(energy_meV, e0_meV, gamma_meV):
    e=np.asarray(energy_meV)
    if e0_meV <= 0 or gamma_meV <= 0: raise ValueError('Positive oscillator parameters required')
    return e*gamma_meV/((e**2-e0_meV**2)**2+e**2*gamma_meV**2)


def kinematic_2d(q_Ainv, energy_meV, beam_kV=30.):
    q=np.asarray(q_Ainv)
    qe=longitudinal_q(energy_meV,beam_kV)
    den=(q*q+qe*qe)**2
    return np.divide(q,den,out=np.zeros(np.broadcast_shapes(q.shape,qe.shape)),where=den>0)


def acceptance_nodes(q_center, dq=0., h=0., sigma_probe=0., n=25, n_probe=7):
    """Area-average weights, normalised to one; absolute collection efficiency omitted."""
    if min(dq,h,sigma_probe) < 0: raise ValueError('Widths must be nonnegative')
    x,w=leggauss(n)
    qx=q_center + dq*x/2 if dq else np.array([q_center])
    wx=w/2 if dq else np.ones(1)
    qy=h*x if h else np.zeros(1);wy=w/2 if h else np.ones(1)
    if sigma_probe:
        g, wg=hermgauss(n_probe);g=g*np.sqrt(2)*sigma_probe;wg=wg/np.sqrt(np.pi)
        qx=(qx[:,None]+g).ravel();wx=(wx[:,None]*wg).ravel()
        qy=(qy[:,None]+g).ravel();wy=(wy[:,None]*wg).ravel()
    q=np.hypot(qx[:,None],qy[None,:]).ravel()
    weights=(wx[:,None]*wy[None,:]).ravel()
    return q,weights/weights.sum()


def averaged_kinematic(q_center, energy_meV, **kwargs):
    beam=kwargs.pop('beam_kV',30.)
    q,w=acceptance_nodes(q_center,**kwargs)
    return kinematic_2d(q[None,:],np.asarray(energy_meV)[:,None],beam)@w


def measured_loss(energy_meV, q_center, components, dq=0., h=0., sigma_probe=0.,
                  resolution_fwhm_meV=0., beam_kV=30., temperature_K=300., n=25):
    """Components are (E0, Gamma, amplitude), all local loss parameters.

    Response is evaluated on the extended axis before convolution. Gain-side
    detailed balance is imposed before instrumental broadening.
    """
    e=np.asarray(energy_meV);de=np.median(np.diff(e))
    if de<=0 or not np.allclose(np.diff(e),de):raise ValueError('Uniform increasing energy required')
    s=resolution_fwhm_meV/np.sqrt(8*np.log(2))
    half=int(np.ceil(8*s/de)) if s else 0
    ee=e[0]+np.arange(-half,len(e)+half)*de
    positive=np.abs(ee)
    loss=sum(a*dl_loss(positive,e0,g) for e0,g,a in components)
    K=averaged_kinematic(q_center,np.maximum(positive,1e-6),dq=dq,h=h,sigma_probe=sigma_probe,n=n,beam_kV=beam_kV)
    result=loss*K*np.where(ee<0,np.exp(-positive/(.08617333262*temperature_K)),1.)
    if not half:return result
    lag=np.arange(-half,half+1)*de;kernel=np.exp(-.5*(lag/s)**2);kernel/=kernel.sum()
    return np.convolve(result,kernel,mode='valid')


def screened_plasmon_energy(q_Ainv, stiffness_meV2_A, screening_length_A):
    """Long-wavelength Drude/RPA pole: E_p^2 = A*q/(1+r0*q).

    A has units meV^2 Å and absorbs n_2D/m* and environmental screening.
    Neither A nor r0 is assigned to any material by this function.
    """
    q=np.asarray(q_Ainv)
    if np.any(q<0) or stiffness_meV2_A<=0 or screening_length_A<0:
        raise ValueError('Nonnegative momentum/screening and positive stiffness required')
    return np.sqrt(stiffness_meV2_A*q/(1+screening_length_A*q))


def dispersive_aperture_loss(energy_meV,q_center,stiffness_meV2_A,screening_length_A,
                            gamma_meV,dq=0.,h=0.,n=35,beam_kV=30.):
    """Integrate K(Q,E)*[-Im 1/epsilon(Q,E)] over the aperture.

    Epsilon = 1-E_p(Q)^2/[E(E+i*gamma)]. This computes the *positive*
    energy response before resolution broadening. Both the pole and its
    spectral weight vary with Q; averaging K alone cannot reproduce this.
    """
    e=np.asarray(energy_meV)
    if np.any(e<=0) or gamma_meV<=0:
        raise ValueError('Positive loss energy and damping required')
    q,w=acceptance_nodes(q_center,dq=dq,h=h,n=n)
    ep2=screened_plasmon_energy(q,stiffness_meV2_A,screening_length_A)**2
    loss=e[:,None]*gamma_meV*ep2/((e[:,None]**2-ep2)**2+(e[:,None]*gamma_meV)**2)
    return (kinematic_2d(q[None,:],e[:,None],beam_kV)*loss)@w
