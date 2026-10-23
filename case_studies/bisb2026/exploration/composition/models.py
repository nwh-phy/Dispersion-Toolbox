"""Conditional positive-response models of B1, projected through the aperture.

All inputs use meV and inverse angstrom; optimizer energy parameters use eV.
The log-normal component is a spectral-density alternative, not a microscopic
calculation of an electron-hole continuum. Absolute oscillator strength is free.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.ndimage import convolve1d
from scipy.optimize import least_squares, nnls

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT/'src/python'))
from qeels_response import acceptance_nodes, kinematic_2d


def oscillator(e, e0, gamma):
    return e*gamma / ((e*e-e0*e0)**2 + (e*gamma)**2)


def lognormal_density(e, mode, width):
    # Parameterize by the mode, not the mean of the log-normal distribution.
    mu = np.log(mode) + width*width
    return np.exp(-.5*((np.log(e)-mu)/width)**2) / (e*width*np.sqrt(2*np.pi))


class Operator:
    def __init__(self, data, h_scale=1., n_radial=160, n_quad=25, n_probe=7,
                 geometry='full', resolution_scale=1.):
        self.energy = np.asarray(data['energy'])
        native = np.asarray(data['energy_native'])
        self.group = len(native)//len(self.energy)
        self.de = np.median(np.diff(native)); self.pad=32
        self.e_meV = native[0] + np.arange(-self.pad, len(native)+self.pad)*self.de
        self.e = self.e_meV/1000
        if np.any(self.e <= 0):
            raise ValueError('This high-energy model needs a strictly positive extended grid')
        self.q = np.asarray(data['q']); self.nq = len(self.q)
        self.radial = np.linspace(1e-7, .035, n_radial)
        self.h = float(data['h'])*h_scale
        self.H = np.zeros((self.nq, len(self.e), n_radial))
        self.KQ = np.zeros((self.nq, len(self.e)))
        self.K = self.KQ.copy()
        for b, center in enumerate(self.q):
            if geometry == 'point':
                q, w = np.array([abs(center)]), np.ones(1)
            else:
                q,w=acceptance_nodes(center, dq=float(data['width']), h=self.h,
                    sigma_probe=float(data['probe']), n=n_quad, n_probe=n_probe)
            if q.max() >= self.radial[-1]:
                raise ValueError('Radial support too small')
            k = kinematic_2d(q[None,:], self.e_meV[:,None], float(data['beam_kV']))
            # Normalize for conditioning only. Per-spectrum amplitudes carry it.
            norm = np.mean(k[len(k)//2]@w)
            k = k*w[None,:]/norm
            self.K[b] = k.sum(1); self.KQ[b] = k@(q/.005)
            j = np.searchsorted(self.radial,q,side='right')-1; j=np.maximum(j,0)
            frac = (q-self.radial[j])/(self.radial[j+1]-self.radial[j])
            # Compress only the intrinsic spectrum onto the radial grid.
            # K is evaluated at the original quadrature nodes, not interpolated.
            for z in np.unique(j):
                mask=j==z
                self.H[b,:,z] += k[:,mask]@(1-frac[mask])
                self.H[b,:,z+1] += k[:,mask]@frac[mask]
        sigma=float(data['resolution_fwhm'])*resolution_scale/np.sqrt(8*np.log(2))
        lags=np.arange(-self.pad,self.pad+1)*self.de
        self.kernel=np.exp(-.5*(lags/max(sigma,1e-9))**2); self.kernel/=self.kernel.sum()
        self.background = (500/self.energy)**3

    def project(self, curve):
        z=convolve1d(curve,self.kernel,axis=-1,mode='constant')
        z=z[...,self.pad:-self.pad]
        return z.reshape(*z.shape[:-1],-1,self.group).sum(-1)

    def low(self, pars, flat=False, skew=False):
        e0, screening_fraction, gamma, exponent=pars
        x=self.radial/.005
        ep2=e0*e0*x/(1+screening_fraction*(x-1))
        if flat:
            ep2=np.full_like(ep2,e0*e0)
        g=gamma*x**exponent
        if skew:
            loss=lognormal_density(self.e[None,:],np.sqrt(ep2)[:,None],g[:,None])
        else:
            loss=oscillator(self.e[None,:],np.sqrt(ep2)[:,None],g[:,None])
        # Drude dielectric loss has numerator Ep^2. In the flat alternative,
        # a finite-frequency polarizability instead has the 2D optical Q weight.
        loss*= x[:,None] if flat else ep2[:,None]
        counts=np.einsum('beq,qe->be',self.H,loss,optimize=False)
        return self.project(counts)

    def columns(self, pars, model, flat=False, include_constant=True, include_high_tail=False,
                upper_dispersion=False):
        cols=[self.low(pars[:4],flat=flat,skew=(model=='skew_single'))]; names=['screened_response']
        i=4
        if model=='double':
            cols.append(self.project(self.KQ*oscillator(self.e,pars[i],pars[i+1])));names.append('additional_resonance');i+=2
        elif model=='continuum':
            cols.append(self.project(self.KQ*lognormal_density(self.e,pars[i],pars[i+1])));names.append('continuum');i+=2
        elif model not in ('single','skew_single'):
            raise ValueError(model)
        cols.append(self.project(self.KQ*oscillator(self.e,pars[i],pars[i+1])));names.append('near_2eV');i+=2
        if include_high_tail:
            if upper_dispersion:
                eq=pars[i]+pars[i+2]*(self.radial/.005-1)
                response=oscillator(self.e[None,:],eq[:,None],pars[i+1])*(self.radial[:,None]/.005)
                cols.append(self.project(np.einsum('beq,qe->be',self.H,response,optimize=False)))
            else:
                cols.append(self.project(self.KQ*oscillator(self.e,pars[i],pars[i+1])))
            names.append('near_3eV')
        cols.append(np.broadcast_to(self.background,(self.nq,len(self.energy))));names.append('elastic_tail')
        if include_constant:
            cols.append(np.ones((self.nq,len(self.energy))));names.append('constant')
        columns=np.stack(cols,-1)
        scale=np.max(columns,axis=1,keepdims=True)
        return columns/np.maximum(scale,1e-100),names


def limits(model, high_tail=False, band_limited=False, upper_dispersion=False):
    # No desired plasmon energy is imposed: the one-mode hypothesis spans B1.
    low=[.3,0,.08,-1.5];high=[2.7,.9999,3.,1.5]
    names=['Eref_eV','screening_fraction','Gamma_ref_eV','Gamma_q_exponent']
    if model=='skew_single':
        low[2]=.10;high[2]=1.5
        names[0]='mode_ref_eV';names[2]='log_width_ref';names[3]='log_width_q_exponent'
    if model=='double':
        low += [.5,.15];high += [3.2,4.]
        names += ['Ehigh_eV','Gamma_high_eV']
    elif model=='continuum':
        low += [.5,.10];high += [3.2,1.5]
        names += ['continuum_mode_eV','continuum_log_width']
    low += [1.7,.03]; high += [2.4,.8]
    names += ['E_neighbor_eV','Gamma_neighbor_eV']
    if high_tail:
        low += [2.8,.2]; high += [4.,3.]
        names += ['E_upper_neighbor_eV','Gamma_upper_neighbor_eV']
        if upper_dispersion:
            low += [0.]; high += [.6]; names += ['upper_slope_eV_per_qref']
    if band_limited:
        high[0]=1.8
        if model in ('double','continuum'):high[4]=1.8
    return np.array(low),np.array(high),names


def solve_amplitudes(columns, y, sigma, mask):
    coefficients=[];prediction=[]
    for a, yy, ss in zip(columns,y,sigma):
        c,_=nnls(a[mask]/ss[mask,None], yy[mask]/ss[mask], maxiter=1000)
        coefficients.append(c);prediction.append(a@c)
    return np.array(coefficients),np.array(prediction)


def fit(op, y, sigma, model, mask, n_starts=6, seed=20261024,
        warm=None, flat=False, include_constant=True, include_high_tail=False,
        max_nfev=250, fixed_values=None, band_limited=False, upper_dispersion=False):
    lo,hi,names=limits(model,include_high_tail,band_limited,upper_dispersion)
    fixed = np.array([False]*len(lo))
    if flat: fixed[1]=True
    fixed_values={} if fixed_values is None else fixed_values
    for k in fixed_values:
        fixed[k]=True
    def unpack(p):
        a=(lo+hi)/2;a[~fixed]=p
        for k,v in fixed_values.items():a[k]=v
        return a
    def evaluate(p,detail=False):
        a=unpack(p)
        col,cnames=op.columns(a,model,flat,include_constant,include_high_tail,upper_dispersion)
        c,pred=solve_amplitudes(col,y,sigma,mask)
        if detail:return pred,c,col,cnames
        return ((pred-y)/sigma)[:,mask].ravel()
    rng=np.random.default_rng(seed)
    initial=[.85,.8,.8,0]
    if model=='double':initial += [1.6,1.5]
    elif model=='continuum':initial += [1.4,.6]
    initial += [2.05,.3]
    if include_high_tail:initial += [3.3,.8]
    if upper_dispersion:initial += [.15]
    starts=[np.array(initial)]
    if warm is not None and len(warm)==len(lo):starts.insert(0,np.array(warm))
    while len(starts)<n_starts:
        v=lo+(hi-lo)*rng.uniform(.05,.95,len(lo))
        v[2]=rng.uniform(.3,1.5);v[3]=rng.uniform(-.6,.6)
        starts.append(v)
    best=None;candidates=[]
    for s in starts[:n_starts]:
        p=np.clip(s,lo+1e-8,hi-1e-8)[~fixed]
        r=least_squares(evaluate,p,bounds=(lo[~fixed],hi[~fixed]),x_scale='jac',
                        max_nfev=max_nfev,ftol=2e-7,xtol=2e-7,gtol=2e-7)
        cost=float(r.fun@r.fun)
        candidates.append(dict(cost=cost,success=bool(r.success),status=int(r.status),nfev=int(r.nfev)))
        if best is None or cost<best[0]:best=(cost,r)
    cost,r=best;pred,c,col,cnames=evaluate(r.x,True);pars=unpack(r.x)
    boundary=((pars-lo)/(hi-lo)<.002)|((hi-pars)/(hi-lo)<.002)
    return dict(model=model,parameters=pars,parameter_names=names,coefficients=c,
        component_names=cnames,prediction=pred,components=col*c[:,None,:],cost=cost,
        mean_squared_residual=float(cost/(op.nq*mask.sum())),success=bool(r.success),
        boundary=boundary,candidates=candidates,flat=flat,include_constant=include_constant,
        include_high_tail=include_high_tail,fixed_values=fixed_values,
        band_limited=band_limited,upper_dispersion=upper_dispersion,
        parameter_count=int((~fixed).sum()+c.size))
