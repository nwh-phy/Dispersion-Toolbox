"""Independent observed-count apex extraction; no q-tracking prior.

The main apex is the highest significant local maximum in the prescribed B1
window. If the window maximum is at an edge, no main peak is assigned. All
other local maxima remain candidates, not inferred oscillators or shoulders.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


def smooth_counts(energy, counts, fwhm):
    step = np.median(np.diff(energy))
    if not np.allclose(np.diff(energy), step):
        raise ValueError('Uniform energy sampling is required')
    return gaussian_filter1d(np.asarray(counts, float), fwhm/2.354820045/step,
                             axis=-1, mode='nearest')


def peaks_from_smoothed(energy, smooth, window=(300, 1800), prominence=.02):
    ix = np.flatnonzero((energy >= window[0]) & (energy <= window[1]))
    e, y = energy[ix], np.asarray(smooth)[ix]
    if not np.all(np.isfinite(y)) or y.max() <= 0:
        return dict(peak_meV=np.nan, edge=True, candidates=[])
    step = np.median(np.diff(e))
    peaks, properties = find_peaks(y, prominence=prominence*y.max(),
                                   distance=max(1, int(round(64/step))))
    candidates = []
    for j, p in enumerate(peaks):
        # Three-point local parabolic interpolation. It does not impose a
        # full spectral line shape or couple adjacent q bins.
        curvature = y[p-1]-2*y[p]+y[p+1]
        offset = .5*(y[p-1]-y[p+1])/curvature if curvature < 0 else 0.
        position = e[p]+np.clip(offset, -.5, .5)*step
        candidates.append(dict(energy_meV=float(position), height=float(y[p]),
            prominence_fraction=float(properties['prominences'][j]/y.max())))
    edge = int(y.argmax()) in (0, len(y)-1)
    chosen = max(candidates, key=lambda c: c['height']) if candidates and not edge else None
    return dict(peak_meV=chosen['energy_meV'] if chosen else np.nan,
                edge=edge, candidates=candidates)


def extract(energy, counts, fwhm=100, window=(300, 1800), prominence=.02):
    return peaks_from_smoothed(energy, smooth_counts(energy, counts, fwhm), window, prominence)


def block_weights(n_frames, repeats, block, seed):
    """Circular moving-block bootstrap; one draw is shared across all q bins."""
    rng = np.random.default_rng(seed)
    starts = rng.integers(n_frames, size=(repeats, int(np.ceil(n_frames/block))))
    indices = ((starts[..., None]+np.arange(block)) % n_frames).reshape(repeats, -1)[:, :n_frames]
    weights = np.zeros((repeats, n_frames), dtype=float)
    for i, ix in enumerate(indices):
        weights[i] = np.bincount(ix, minlength=n_frames)
    return weights


def bootstrap_peaks(energy, frames, cfg, seed):
    weights = block_weights(frames.shape[1], cfg['bootstrap_repeats'],
                            cfg['bootstrap_block_frames'], seed)
    y = np.einsum('bt,qte->bqe', weights, frames.astype(float), optimize=True)
    smooth = smooth_counts(energy, y, cfg['smooth_fwhm_meV'])
    result = np.full(smooth.shape[:2], np.nan)
    for b in range(len(smooth)):
        for q in range(smooth.shape[1]):
            result[b, q] = peaks_from_smoothed(energy, smooth[b,q], cfg['peak_window_meV'],
                cfg['prominence_fraction'])['peak_meV']
    return result
