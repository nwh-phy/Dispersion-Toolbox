"""Guard against invented dispersion and hidden q/energy averaging."""
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'case_studies/bisb2026/exploration/composition'))
from prepare import adaptive_q_bins
from ridge_extract import extract, block_weights


def test_adaptive_bins_partition_native_pixels_without_seven_channel_bins():
    for dq in [.0005,.00025]:
        for n in [3,5]:
            q,w=adaptive_q_bins(dq,.006,.020,n)
            np.testing.assert_allclose(q,-q[::-1])
            np.testing.assert_allclose(q[:-1]+w[:-1]/2,q[1:]-w[1:]/2)
            np.testing.assert_allclose(w[abs(q)<=.006+1e-12],dq)
            np.testing.assert_allclose(w[abs(q)>.006+1e-12],n*dq)
            assert q[-1]+w[-1]/2 <= .020+dq/2+1e-12


def test_independent_apices_recover_rising_flat_and_falling_peaks():
    e=np.arange(166.,2200.,16.)
    for truth in [[450,700,1100,1400],[900]*4,[1400,1100,700,450]]:
        recovered=[]
        for mu in truth:
            y=5+100*np.exp(-.5*((e-mu)/130)**2)
            recovered.append(extract(e,y)['peak_meV'])
        np.testing.assert_allclose(recovered,truth,atol=1)


def test_monotone_tail_has_no_invented_interior_peak():
    e=np.arange(166.,2200.,16.)
    for y in [e**-3,e]:
        result=extract(e,y)
        assert result['edge'] and np.isnan(result['peak_meV'])


def test_multiple_peaks_remain_candidates_and_dominance_can_switch():
    e=np.arange(166.,2200.,16.)
    low=np.exp(-.5*((e-650)/90)**2); high=np.exp(-.5*((e-1450)/110)**2)
    a=extract(e,2*low+high);b=extract(e,low+2*high)
    assert len(a['candidates'])==len(b['candidates'])==2
    assert abs(a['peak_meV']-650)<1 and abs(b['peak_meV']-1450)<1


def test_bootstrap_preserves_frame_count_and_is_reproducible():
    w=block_weights(300,100,10,123)
    np.testing.assert_equal(w.sum(1),np.full(100,300))
    np.testing.assert_array_equal(w,block_weights(300,100,10,123))
    assert np.any(w>1) and np.any(w==0)
