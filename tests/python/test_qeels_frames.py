import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src" / "python"))
from qeels_frames import (energy_axis, momentum_axis, integer_alignment,
                          apply_integer_alignment, invalid_channels,
                          noise_from_frames, circular_block_indices, dl_apex, elastic_origin)
from qeels_rebin import fractional_bin_weights


def test_explicit_axes_units_and_signed_q():
    e = energy_axis({"spatial_calibrations": [{"units": "eV", "scale": .004, "offset": -.4}]}, 110)
    assert e[100] == 0
    np.testing.assert_allclose(np.diff(e), 4)
    np.testing.assert_allclose(momentum_axis(5, 3, .0005), [-.001, -.0005, 0, .0005, .001])
    with pytest.raises(ValueError):
        energy_axis({}, 100)


def test_alignment_recovers_injected_shift_and_conserves_retained_counts():
    e = np.arange(-200, 804, 4.)
    profile = 1000 * np.exp(-.5*(e/12)**2) + 30*np.exp(-.5*((e-400)/50)**2)
    offsets = np.array([-4, -2, 0, 1, 3])
    frames = np.stack([np.roll(profile, k) for k in offsets])
    a = integer_alignment(frames, e)
    np.testing.assert_array_equal(a['offsets'], offsets)
    out = apply_integer_alignment(frames, a)
    np.testing.assert_allclose(out, np.tile(profile[a['support']], (5, 1)), atol=1e-20)
    np.testing.assert_allclose(np.sum(out),
        sum(np.sum(frames[t, a['support']+k]) for t, k in enumerate(offsets)),
        rtol=5*np.finfo(float).eps, atol=0)
    assert np.argmax(np.sum(out, axis=0)) == np.argmin(np.abs(a['energy_meV']))


def test_invalid_frame_excluded_and_no_circular_edges():
    e = np.arange(-200, 204, 4.)
    profile = 1000*np.exp(-.5*(e/8)**2)
    frames = np.stack([profile, np.roll(profile, 5), np.ones_like(profile)])
    a = integer_alignment(frames, e)
    assert a['valid'].tolist() == [True, True, False]
    out = apply_integer_alignment(frames, a)
    assert out.shape == (2, len(e)-5)


def test_measured_elastic_origin_precedes_frame_alignment():
    e = np.arange(-400, 804, 4.)
    profile = 1000*np.exp(-.5*((e+200)/12)**2)
    frames = np.stack([np.roll(profile,k) for k in [-3,0,2]])
    calibrated, origin = elastic_origin(frames, e, (-350,100))
    assert abs(origin+200)<=4
    a=integer_alignment(frames,calibrated)
    assert np.all(a['valid'])
    np.testing.assert_array_equal(a['offsets'],[-3,0,2])
    assert np.argmax(apply_integer_alignment(frames,a).sum(0)) == np.argmin(abs(a['energy_meV']))


def test_sentinel_channel_cannot_become_zero_reference():
    raw = np.ones((4, 5, 30)); raw[0, 4, 0] = 1.3e10; raw[3, 3, 2] = np.nan
    np.testing.assert_array_equal(invalid_channels(raw, 2**32-1), [False,False,False,True,True])


def test_noise_gain_recovered_without_assuming_one_electron_counts():
    rng = np.random.default_rng(8)
    expected = np.linspace(20, 200, 101)
    frames = .82*rng.poisson(expected, size=(4000,101))
    noise = noise_from_frames(frames, np.ones(4000))
    assert abs(np.median(noise['g'])-.82) < .04
    assert abs(np.median(noise['rho_time'])) < .03


def test_block_draws_preserve_adjacency_and_length():
    ix = circular_block_indices(300, 10, np.random.default_rng(1))
    assert len(ix) == 300 and np.min(ix) >= 0 and np.max(ix) < 300
    np.testing.assert_array_equal(np.diff(ix.reshape(-1,10),axis=1)%300, np.ones((30,9),int))


def test_dl_apex_is_not_oscillator_parameter():
    e = np.linspace(1, 1500, 150000)
    e0, g = 1000., 600.
    y = e*g/((e**2-e0**2)**2+e**2*g*g)
    assert abs(e[np.argmax(y)] - dl_apex(e0,g)) < .02
    assert dl_apex(e0,g) < e0


def test_fractional_q_bin_preserves_constant_counts_and_partition():
    q=np.arange(-10,11)*.0005
    w=fractional_bin_weights(q,.001,.0015,.00017)
    assert abs(w.sum()-3)<1e-12
    assert np.all((w>=0)&(w<=1+1e-12))
    a=fractional_bin_weights(q,.00025,.001,.00017)
    b=fractional_bin_weights(q,.00125,.001,.00017)
    union=fractional_bin_weights(q,.00075,.002,.00017)
    np.testing.assert_allclose(a+b,union,atol=1e-12)
    with pytest.raises(ValueError):fractional_bin_weights(q,.005,.002,.0005)
    with pytest.raises(ValueError):fractional_bin_weights(q,np.nan,.001)
    with pytest.raises(ValueError):fractional_bin_weights(np.array([0.]),0.,.001)


def test_fractional_q_bin_recovers_injected_subpixel_translation():
    q=np.arange(-30,31)*.0005
    reference=1000+1e5*q  # exactly linear cell-average density
    offset=.00017
    shifted=1000+1e5*(q-offset)
    baseline=fractional_bin_weights(q,.002,.0015)@reference
    corrected=fractional_bin_weights(q,.002,.0015,offset)@shifted
    np.testing.assert_allclose(corrected,baseline,rtol=1e-12)
