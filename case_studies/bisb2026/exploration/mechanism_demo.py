"""Illustrative forward response and independent synthetic fit inputs.

The numerical oscillator energies/widths are test cases, not material inputs.
"""
import argparse
from pathlib import Path
import sys
import json
import numpy as np
from scipy.io import savemat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src/python'))
from qeels_response import dl_loss, averaged_kinematic, measured_loss, longitudinal_q
from qeels_frames import dl_apex


def run(output):
    output.mkdir(parents=True,exist_ok=True)
    e=np.arange(-300.,2204.,4.); ep=e[e>0]
    qs=np.linspace(.0005,.03,80);e0=700.;gamma=500.
    maps=[];ridges=[]
    for h in [0.,.0066]:
        spectra=np.array([measured_loss(e,q,[(e0,gamma,1.)],dq=.0005,h=h,resolution_fwhm_meV=24.) for q in qs])
        spectra=spectra[:,e>0]
        maps.append(spectra/spectra.max(axis=1)[:,None]);ridges.append(ep[np.argmax(spectra,axis=1)])
    fig,axs=plt.subplots(1,3,figsize=(13,3.6),layout='constrained')
    for k,h in enumerate([0.,.0066]):
        axs[k].pcolormesh(qs,ep,maps[k].T,shading='auto',cmap='magma')
        axs[k].axhline(dl_apex(e0,gamma),color='cyan',ls='--',lw=1,label='Intrinsic loss apex')
        axs[k].plot(qs,ridges[k],color='white',lw=1,label='Count-spectrum maximum')
        axs[k].set(xlabel='Nominal q (Å$^{-1}$)',ylabel='Energy (meV)',ylim=(0,1600),title=f'Synthetic constant E0; h={h:g} Å$^{{-1}}$')
    axs[0].legend(fontsize=7)
    axs[2].plot(qs,ridges[0],label='Channel width only')
    axs[2].plot(qs,ridges[1],label='Finite perpendicular aperture')
    axs[2].axhline(e0,color='k',ls=':',label='Fixed oscillator E0')
    axs[2].axhline(dl_apex(e0,gamma),color='k',ls='--',label='Fixed loss apex')
    axs[2].set(xlabel='Nominal q (Å$^{-1}$)',ylabel='Energy (meV)',title='Kinematics alone shifts the measured ridge')
    axs[2].legend(fontsize=8)
    fig.savefig(output/'mechanism_demo.png',dpi=180);fig.savefig(output/'mechanism_demo.pdf');plt.close(fig)
    np.savetxt(output/'mechanism_ridge.csv',np.c_[qs,ridges[0],ridges[1]],delimiter=',',header='q_Ainv,ridge_h0_meV,ridge_h0066_meV',comments='')
    # Independent exact construction of the joint-fit model, including kernel.
    lag=np.arange(-300.,304.,4.);ee=np.arange(e[0]-300,e[-1]+304,4.)
    pearson=lambda x:1e5*(1+(x/8.)**2/3.)**-3 + 3000*(1+(x/25.)**2)**-1
    kernel=pearson(lag);kernel/=kernel.sum();zlp=pearson(e)
    rng=np.random.default_rng(20261023);ys=[];truth=[];qcase=[];ncase=[];levels=[];replicas=[]
    for q in [.0025,.015]:
        K=averaged_kinematic(q,np.maximum(abs(ee),250),dq=.0005,h=.0066,n=53)
        profiles=[]
        for E0,G in [(700,500),(1600,900)]:
            response=dl_loss(abs(ee),E0,G)*K*np.where(ee<0,np.exp(-abs(ee)/(.08617333262*300)),1.)
            profiles.append(np.convolve(response,kernel,mode='valid'))
        for nc in [1,2]:
            for level in [80.,5.]:
                mu=zlp+sum(profiles[j]/profiles[j].max()*level*(1 if j==0 else .7) for j in range(nc))
                for rep in range(5):
                    y=mu if rep==0 else .82*rng.poisson(mu/.82)
                    ys.append(y);truth.append(mu);qcase.append(q);ncase.append(nc);levels.append(level);replicas.append(rep)
    savemat(output/'synthetic.mat',dict(E=e[:,None],Y=np.array(ys).T,mu=np.array(truth).T,
        q=np.array(qcase),n_true=np.array(ncase),height=np.array(levels),replicate=np.array(replicas),
        E0_true=np.array([700.,1600.]),Gamma_true=np.array([500.,900.])),do_compression=True)
    (output/'mechanism_summary.json').write_text(json.dumps(dict(
        status='illustrative mechanism, not material prediction',beam_kV=30,qE_1eV_Ainv=float(longitudinal_q(1000)),
        E0_meV=e0,Gamma_meV=gamma,intrinsic_apex_meV=float(dl_apex(e0,gamma)),
        ridge_h0_first_last=[float(ridges[0][0]),float(ridges[0][-1])],
        ridge_aperture_first_last=[float(ridges[1][0]),float(ridges[1][-1])],n_synthetic=len(ys)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);run(p.parse_args().output)
