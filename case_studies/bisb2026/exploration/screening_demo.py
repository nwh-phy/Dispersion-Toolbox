"""Dimensioned, illustrative screening/aperture experiment; no BiSb parameters."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src/python'))
from qeels_response import screened_plasmon_energy, dispersive_aperture_loss


def run(dest):
    dest.mkdir(parents=True,exist_ok=True)
    q=np.linspace(.00025,.04,100);e=np.arange(4.,2504.,4.)
    # Illustrative axes deliberately span the observed energy/momentum scale.
    # Stiffness is fixed across screening variants, not tuned to produce a peak.
    A=1e8;rs=[0.,50.,200.,500.];g=400.;curves=[]
    fig,ax=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    for r in rs:
        y=screened_plasmon_energy(q,A,r);curves.append(y)
        ax[0].plot(q,y,label=f'r0={r:g} Å')
    ax[0].set(xlabel='q (Å$^{-1}$)',ylabel='Pole energy (meV)',title='Screening at fixed stiffness A')
    ax[0].legend(fontsize=8)
    for h in [0.,.0033,.0066]:
        rows=np.array([dispersive_aperture_loss(e,x,A,200.,g,dq=.0005,h=h) for x in q])
        ridge=e[rows.argmax(1)]
        ax[1].plot(q,ridge,label=f'h={h:g} Å$^{{-1}}$')
        np.savetxt(dest/f'dispersive_ridge_h{h:.4f}.csv',np.c_[q,ridge],delimiter=',',header='q_Ainv,measured_ridge_meV',comments='')
    ax[1].plot(q,screened_plasmon_energy(q,A,200.),'k--',label='Intrinsic pole')
    ax[1].set(xlabel='Nominal q (Å$^{-1}$)',ylabel='Energy (meV)',title='Integrate dispersive response over Q')
    ax[1].legend(fontsize=8)
    pole=screened_plasmon_energy(q,A,200.)
    for d in [5.,20.]:
        ax[2].plot(q,pole*np.sqrt(1+np.exp(-q*d)),label=f'Symmetric, d={d:g} Å')
        ax[2].plot(q,pole*np.sqrt(1-np.exp(-q*d)),ls='--',label=f'Antisymmetric, d={d:g} Å')
    ax[2].set(xlabel='q (Å$^{-1}$)',ylabel='Conditional pole energy (meV)',title='Two identical conducting sheets: poles only')
    ax[2].legend(fontsize=7)
    fig.suptitle('Illustrative Drude limit of RPA; parameters are not material predictions',fontsize=11)
    fig.savefig(dest/'screening_aperture.png',dpi=180);fig.savefig(dest/'screening_aperture.pdf');plt.close(fig)
    np.savetxt(dest/'screening_dispersion.csv',np.c_[q,*curves],delimiter=',',header='q_Ainv,'+','.join(f'E_r0_{r:g}_meV' for r in rs),comments='')
    (dest/'screening_inputs.json').write_text(json.dumps(dict(stiffness_meV2_A=A,screening_lengths_A=rs,gamma_meV=g,
        coupled_sheet_distances_A=[5,20],status='dimensioned illustrative scan, no measured material parameters'),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);run(p.parse_args().output)
