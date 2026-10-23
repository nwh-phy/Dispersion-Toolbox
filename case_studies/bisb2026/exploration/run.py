#!/usr/bin/env python3
"""Reproducible local driver; source data are opened only for reading."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
import tomllib
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "python"))
from qeels_frames import (energy_axis, momentum_axis, invalid_channels,
                          integer_alignment, apply_integer_alignment,
                          noise_from_frames, sha256_file, elastic_origin)


def write_json(path, value):
    def default(v):
        if isinstance(v, np.ndarray): return v.tolist()
        if isinstance(v, np.generic): return v.item()
        if isinstance(v, Path): return str(v)
        raise TypeError(type(v).__name__)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, default=default)+"\n")


def write_csv(path, rows):
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def fwhm_axis(y, e):
    k = np.argmax(y); half = y[k]/2
    left = np.flatnonzero(y[:k] < half); right = np.flatnonzero(y[k:] < half)
    if not len(left) or not len(right): return np.nan
    l=left[-1]; r=k+right[0]
    xl = e[l]+(half-y[l])/(y[l+1]-y[l])*(e[l+1]-e[l])
    xr = e[r-1]+(half-y[r-1])/(y[r]-y[r-1])*(e[r]-e[r-1])
    return xr-xl


def prepare(args):
    output=args.output.resolve(); data=args.data_root.resolve()
    if output == data or data in output.parents:
        raise ValueError("Output must not be inside the source-data directory")
    output.mkdir(parents=True, exist_ok=True)
    config_path=Path(__file__).with_name("datasets.toml")
    cfg=tomllib.loads(config_path.read_text())
    manifest={"config":cfg,"config_sha256":sha256_file(config_path),
              "code_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
              "code_sha256":{str(f.relative_to(ROOT)):sha256_file(f) for f in
                [Path(__file__),ROOT/'src/python/qeels_frames.py',ROOT/'src/python/qeels_response.py',
                 Path(__file__).with_name('bisb_explore_fit.m')]},
              "sessions":[],"created_unix":time.time()}
    checks=[]
    for session in cfg['sessions']:
        t0=time.monotonic(); sid=session['id']; dest=output/sid; dest.mkdir(exist_ok=True)
        source=data/session['folder']/session['file']; meta_path=source.with_suffix('.json')
        meta=json.loads(meta_path.read_text()); raw=np.load(source,mmap_mode='r')
        if raw.ndim != 3 or not meta.get('is_sequence'):
            raise ValueError(f'{sid}: expected explicitly marked sequence [frame,q,E]')
        e_metadata=energy_axis(meta,raw.shape[-1]); q=momentum_axis(raw.shape[1],session['q_zero_column'],session['dq_Ainv'])
        bad=invalid_channels(raw,cfg['invalid_count_threshold'])
        reference_idx=np.flatnonzero((abs(q)<=cfg['reference_abs_q_Ainv']+1e-12)&~bad)
        reference=np.sum(raw[:,reference_idx,:],axis=1,dtype=float)
        e,origin=elastic_origin(reference,e_metadata,cfg['raw_elastic_search_window_meV'])
        al=integer_alignment(reference,e,cfg['alignment_window_meV'])
        valid=np.flatnonzero(al['valid']); support=al['support']; ea=al['energy_meV']
        ar=apply_integer_alignment(reference,al)
        beam=np.sum(ar[:,abs(ea)<=100],axis=1)
        if not np.all(beam>0): raise ValueError('Non-positive elastic beam reference')
        rawmap=np.sum(raw[valid],axis=0,dtype=float).T; rawmap[:,bad]=np.nan
        # Summing per-frame shifted pixels; no interpolation and no wrap.
        am=np.zeros((len(support),raw.shape[1]))
        for t in valid: am+=np.asarray(raw[t,:,support+al['offsets'][t]],float)
        am[:,bad]=np.nan
        qview=(abs(q)<=.031)&~bad
        ekeep=(ea>=cfg['fit_energy_range_meV'][0])&(ea<=cfg['fit_energy_range_meV'][1])
        energy=ea[ekeep]
        rows=[]; Fs=[]; Fus=[]; noises=[]; unnoises=[]; trainnoises=[]; testnoises=[]; qs=[]; widths=[]; native=[]; members=[]
        conservation=[]
        for width in [1,session['bin_channels']]:
            for target in sorted(np.r_[-np.array(cfg['targets_abs_q_Ainv']),cfg['targets_abs_q_Ainv']]):
                c=int(np.argmin(abs(q-target))); ii=np.arange(c-width//2,c+width//2+1)
                if np.any(ii<0) or np.any(ii>=len(q)) or np.any(bad[ii]):
                    raise ValueError(f'{sid}: invalid planned bin {target}, N={width}')
                f=np.sum(raw[:,ii,:],axis=1,dtype=float)
                aligned=apply_integer_alignment(f,al)
                expected=sum(np.sum(f[t,support+al['offsets'][t]]) for t in valid)
                err=abs(aligned.sum()-expected)/max(abs(expected),1)
                conservation.append(err)
                fu=f[valid][:,support][:,ekeep]; aligned=aligned[:,ekeep]
                noise=noise_from_frames(aligned,beam)
                unnoise=noise_from_frames(fu,beam)
                trainnoise=noise_from_frames(aligned[:-50],beam[:-50])
                testnoise=noise_from_frames(aligned[-50:],beam[-50:])
                y=aligned.sum(axis=0); sig=(energy>=300)&(energy<=1800)
                ridge=energy[sig][np.argmax(gaussian_filter1d(y[sig],3))]
                row=dict(id=len(rows)+1,q_Ainv=float(q[ii].mean()),n_channels=width,
                    q_width_Ainv=width*session['dq_Ainv'],counts_signal=float(y[sig].sum()),
                    raw_band_centroid_meV=float(np.sum(y[sig]*energy[sig])/y[sig].sum()),
                    smoothed_count_ridge_meV=float(ridge),g_median=float(np.median(noise['g'][sig])),
                    raw_variance_ratio=float(np.median(noise['raw_ratio'][sig])),
                    corrected_variance_ratio=float(np.median(noise['corrected_ratio'][sig])),
                    lag1_time=float(np.median(noise['rho_time'][sig])))
                rows.append(row); Fs.append(aligned.T); Fus.append(fu.T)
                noises.append(noise['sigma']);unnoises.append(unnoise['sigma'])
                trainnoises.append(trainnoise['sigma']);testnoises.append(testnoise['sigma'])
                qs.append(q[ii].mean()); widths.append(width*session['dq_Ainv'])
                native.append(width==1); members.append(q[ii])
        frames=np.stack(Fs,axis=-1); unaligned=np.stack(Fus,axis=-1)
        cell=np.empty((1,len(members)),dtype=object)
        for k,m in enumerate(members): cell[0,k]=m
        # Narrow central five-channel angular second moment, used only as a diagnostic/variant.
        elastic=np.sum(am[abs(ea)<=16],axis=0); k0=session['q_zero_column']-1
        kk=np.arange(k0-2,k0+3); weights=elastic[kk]/elastic[kk].sum()
        center=np.sum(weights*q[kk]); probe=np.sqrt(max(np.sum(weights*(q[kk]-center)**2)-session['dq_Ainv']**2/12,0))
        sio.savemat(dest/'prepared.mat',dict(E=energy[:,None],frames=frames.astype('float32'),
            unaligned_frames=unaligned.astype('float32'),sigma=np.stack(noises,axis=1),
            unaligned_sigma=np.stack(unnoises,axis=1),train_sigma=np.stack(trainnoises,axis=1),test_sigma=np.stack(testnoises,axis=1),
            q=np.array(qs)[None,:],width=np.array(widths)[None,:],native=np.array(native)[None,:],
            members=cell,beam=beam[:,None],dq=session['dq_Ainv'],h=session['nominal_h_Ainv'],
            beam_kV=session['beam_kV'],sigma_probe=probe,counts_per_electron=meta['metadata']['hardware_source'].get('counts_per_electron',np.nan),
            reference=reference.T,E_raw=e[:,None],align_offsets=al['offsets'][None,:],
            align_support=(support+1)[:,None],align_valid=al['valid'][None,:]),do_compression=True)
        np.savez_compressed(dest/'maps.npz',energy=ea,q=q,aligned=am,raw=rawmap,
            energy_raw=e-e[al['anchor']],offsets=al['offsets'],beam=beam,valid=valid,bad=bad)
        before=fwhm_axis(np.sum(reference[valid],axis=0),e)
        after=fwhm_axis(ar.sum(axis=0),ea)
        hw=meta['metadata']['hardware_source']
        measured_beam_kV=meta['metadata']['instrument']['high_tension']/1000
        if not np.isclose(measured_beam_kV,session['beam_kV']):
            raise ValueError('Configured beam energy differs from acquisition metadata')
        summary=dict(session=sid,source=str(source),sha256=sha256_file(source),
            metadata_sha256=sha256_file(meta_path),shape=list(raw.shape),dtype=str(raw.dtype),
            energy_step_meV=float(np.median(np.diff(e))),counts_per_electron=hw.get('counts_per_electron'),
            metadata_elastic_peak_meV=origin,final_energy_origin_meV=float(e_metadata[al['anchor']]),beam_kV=measured_beam_kV,
            exposure_seconds=hw.get('exposure'),invalid_columns_1based=(np.flatnonzero(bad)+1),
            frames_used=len(valid),frames_excluded=(np.flatnonzero(~al['valid'])+1),
            energy_drift_span_meV=float(np.ptp(al['offsets'][valid])*np.median(np.diff(e))),
            zlp_fwhm_before_meV=float(before),zlp_fwhm_after_meV=float(after),
            beam_relative_std=float(np.std(beam,ddof=1)/np.mean(beam)),
            central_probe_sigma_diagnostic_Ainv=float(probe),
            max_retained_count_relative_error=max(conservation),
            q_zero_column=session['q_zero_column'],dq_Ainv=session['dq_Ainv'],
            physical_bin_width_Ainv=session['bin_channels']*session['dq_Ainv'],seconds=time.monotonic()-t0)
        manifest['sessions'].append(summary); write_json(dest/'summary.json',summary);write_csv(dest/'spectra.csv',rows)
        checks.append(dict(session=sid,valid_frames=len(valid),all_counts_conserved=max(conservation)<1e-12,
                           all_energy_steps_4meV=bool(np.allclose(np.diff(energy),4))))
        fig,axs=plt.subplots(2,2,figsize=(10,6),layout='constrained')
        axs[0,0].plot(valid+1,al['offsets'][valid]*np.median(np.diff(e)),lw=1);axs[0,0].set(xlabel='Frame',ylabel='Energy offset (meV)',title=f'{sid}: common elastic reference')
        axs[0,1].plot(valid+1,beam/beam.mean(),lw=1);axs[0,1].set(xlabel='Frame',ylabel='Relative elastic intensity')
        axs[1,0].plot(e-e[al['anchor']],reference[valid].sum(axis=0)/reference[valid].sum(axis=0).max(),label=f'Sum first: {before:.1f} meV')
        axs[1,0].plot(ea,ar.sum(axis=0)/ar.sum(axis=0).max(),label=f'Align first: {after:.1f} meV')
        axs[1,0].set(xlim=(-80,80),xlabel='Energy (meV)',ylabel='Normalised elastic reference');axs[1,0].legend(fontsize=8)
        for sign in [-1,1]:
            rr=[r for r in rows if r['n_channels']>1 and np.sign(r['q_Ainv'])==sign]
            axs[1,1].plot([abs(r['q_Ainv']) for r in rr],[r['g_median'] for r in rr],'o-',label=f'{sign:+d}q')
        axs[1,1].axhline(.82,color='k',ls=':',label='metadata: 0.82 counts/electron')
        axs[1,1].set(xlabel='|q| (Å$^{-1}$)',ylabel='Variance / mean counts');axs[1,1].legend(fontsize=8)
        fig.savefig(dest/'frame_diagnostics.png',dpi=180);fig.savefig(dest/'frame_diagnostics.pdf');plt.close(fig)
        fig,ax=plt.subplots(figsize=(7,4),layout='constrained'); ev=(ea>=200)&(ea<=2200)
        z=am[np.ix_(ev,qview)]; z=z/np.maximum(z.sum(axis=0),1)
        ax.pcolormesh(q[qview],ea[ev],z,shading='auto',cmap='magma',vmin=0,vmax=np.nanpercentile(z,99))
        ax.set(xlabel='Signed q (Å$^{-1}$)',ylabel='Energy loss (meV)',title=f'{sid}: aligned count map, column-normalised')
        fig.savefig(dest/'aligned_map.png',dpi=180);plt.close(fig)
        print(json.dumps(summary,default=lambda a:a.tolist() if isinstance(a,np.ndarray) else str(a)),flush=True)
        del raw,frames,unaligned,Fs,Fus
    write_json(output/'manifest.json',manifest);write_json(output/'preparation_checks.json',checks)
    if not all(c['all_counts_conserved'] for c in checks): raise AssertionError('Count conservation failed')


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare'])
    p.add_argument('--data-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();prepare(args)


if __name__=='__main__': main()
