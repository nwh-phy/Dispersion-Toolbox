"""Nonempty parity checks, source integrity, and published-apex readback."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from io_utils import digest,dump
from ridge_extract import extract


def verify(output, reference):
    rows=list(csv.DictReader((output/'peaks.csv').open()))
    cfg=json.loads((output/'signature.json').read_text())['protocol']
    checks=[];raw=[];readback=0
    for n in [3,5]:
        manifest=json.loads((output/f'bin{n}/manifest.json').read_text())
        for sid in ['590','n0','20w']:
            record=next(r for r in manifest['sessions'] if r['session']==sid)
            assert digest(output/f'bin{n}/{sid}.npz')==record['prepared_sha256']
            if n==3:
                sha=digest(record['source']);assert sha==record['source_sha256']
                raw.append(dict(session=sid,sha256=sha,unchanged=True))
            d=dict(np.load(output/f'bin{n}/{sid}.npz'))
            old=dict(np.load(reference/f'bin{n}/{sid}.npz'))
            iq,oq=np.where(np.isclose(d['q'][:,None],old['q'][None,:],rtol=0,atol=1e-10))
            ie,oe=np.where(np.isclose(d['energy'][:,None],old['energy'][None,:],rtol=0,atol=1e-6))
            assert len(ie)>0,'Cannot claim parity with an empty energy intersection'
            for i,j in zip(iq,oq):
                assert abs(d['width'][i]-old['width'])<1e-12
                np.testing.assert_array_equal(d['frames'][i][:,ie],old['frames'][j][:,oe])
            checks.append(dict(session=sid,channels=n,common_q_bins=len(iq),
                energy_channels=len(ie),exact=True if len(iq) else None,
                status='verified' if len(iq) else 'no_common_q_centre'))
            rr=[r for r in rows if r['session']==sid and int(r['high_channels'])==n]
            assert len(rr)==len(d['q'])
            for i,r in enumerate(rr):
                assert abs(float(r['q_Ainv'])-d['q'][i])<1e-12
                result=extract(d['energy'],d['frames'][i].astype(float).sum(0),
                    cfg['smooth_fwhm_meV'],cfg['peak_window_meV'],cfg['prominence_fraction'])
                expected=result['peak_meV'] if abs(d['q'][i])>1e-10 else np.nan
                np.testing.assert_allclose(float(r['peak_meV']),expected,equal_nan=True,atol=1e-10)
                readback+=1
    result=dict(previous_fixed_bin_parity=checks,raw_sources=raw,apex_rows_verified=readback,
                verification_code_sha256=digest(Path(__file__)))
    dump(output/'additional_verification.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    a=p.parse_args();verify(a.output,a.reference)
