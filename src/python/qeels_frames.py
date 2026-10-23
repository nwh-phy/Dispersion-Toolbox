"""Count-preserving processing of explicitly calibrated q-EELS frame sequences.

No sample, filename, instrument or peak-model assumptions are embedded here.
Arrays use [frame, momentum, energy]. Counts are never clipped or denoised.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import median_filter


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def energy_axis(metadata: dict, n_energy: int) -> np.ndarray:
    axes = metadata.get("spatial_calibrations", [])
    # A documented axis order is required; do not guess from array lengths.
    if not axes or axes[-1].get("units") not in {"eV", "meV"}:
        raise ValueError("Last raw dimension must have explicit eV or meV calibration")
    axis = axes[-1]
    factor = 1000.0 if axis["units"] == "eV" else 1.0
    out = (float(axis["offset"]) + np.arange(n_energy) * float(axis["scale"])) * factor
    if not np.all(np.isfinite(out)) or not np.all(np.diff(out) > 0):
        raise ValueError("Energy axis must be finite and increasing")
    return out


def momentum_axis(n_q: int, zero_column: int, dq: float) -> np.ndarray:
    if not 1 <= zero_column <= n_q or dq <= 0 or not np.isfinite(dq):
        raise ValueError("Invalid explicit q calibration")
    return (np.arange(n_q) - (zero_column - 1)) * dq


def elastic_origin(reference: np.ndarray, energy: np.ndarray,
                   search_window: tuple[float, float]) -> tuple[np.ndarray, float]:
    """Calibrate the common energy origin from a declared raw elastic region.

    The sequence-summed reference locates the elastic feature. Subsequent
    frame alignment still measures each frame independently in a narrow window.
    """
    ref = np.asarray(reference, dtype=float)
    idx = np.flatnonzero((energy >= search_window[0]) & (energy <= search_window[1]))
    if len(idx) < 3 or ref.ndim != 2 or ref.shape[1] != len(energy):
        raise ValueError('Invalid elastic-origin search region')
    profile = ref[:, idx].sum(axis=0)
    k = int(np.argmax(profile))
    if not np.all(np.isfinite(profile)) or k in (0, len(idx)-1) or profile[k] <= 3*max(np.median(profile),np.finfo(float).eps):
        raise ValueError('No unambiguous interior elastic reference in the declared region')
    origin = float(energy[idx[k]])
    return energy-origin, origin


def invalid_channels(raw: np.ndarray, threshold: float) -> np.ndarray:
    bad = np.zeros(raw.shape[1], dtype=bool)
    for frame in raw:
        bad |= np.any(~np.isfinite(frame) | (frame >= threshold), axis=-1)
    return bad


def integer_alignment(reference: np.ndarray, energy: np.ndarray,
                      window=(-100.0, 100.0)) -> dict:
    """Reference is [frame, E]; output indices preserve a common valid support.

    Matches qe_zlp_integer_align: local maxima must be interior and >3 times
    the window median. Offsets are measured positions minus median position.
    """
    ref = np.asarray(reference, dtype=float)
    energy = np.asarray(energy, dtype=float)
    window_idx = np.flatnonzero((energy >= window[0]) & (energy <= window[1]))
    if len(window_idx) < 3 or ref.ndim != 2 or ref.shape[1] != len(energy):
        raise ValueError("Invalid reference axis/window")
    local = ref[:, window_idx]
    maximum = np.argmax(local, axis=1)
    heights = local[np.arange(len(local)), maximum]
    valid = (np.all(np.isfinite(local), axis=1) & (maximum > 0)
             & (maximum < len(window_idx) - 1)
             & (heights > 3 * np.maximum(np.median(local, axis=1), np.finfo(float).eps)))
    if not np.any(valid):
        raise ValueError("No valid elastic reference frames")
    peaks = window_idx[maximum]
    # MATLAB round for nonnegative pixel indices (not NumPy bankers rounding).
    anchor = int(np.floor(np.median(peaks[valid]) + 0.5))
    offsets = peaks - anchor
    lo = max(0, -int(np.min(offsets[valid])))
    hi = min(len(energy), len(energy) - int(np.max(offsets[valid])))
    support = np.arange(lo, hi)
    return dict(valid=valid, anchor=anchor, offsets=offsets, support=support,
                energy_meV=energy[support] - energy[anchor], peak_pixels=peaks)


def apply_integer_alignment(frames: np.ndarray, alignment: dict) -> np.ndarray:
    """Return valid frames [T, ..., E] with no wrap, interpolation or padding."""
    keep = np.flatnonzero(alignment["valid"])
    return np.stack([np.take(frames[t], alignment["support"] + alignment["offsets"][t],
                             axis=-1) for t in keep])


def noise_from_frames(frames: np.ndarray, reference: np.ndarray,
                      smooth_bins=21) -> dict:
    """Conditional detector-noise estimate from common-reference-corrected frames.

    Difference variance removes slow drift; raw and corrected variance are
    retained as diagnostics. This is not a full covariance or bootstrap CI.
    """
    frames = np.asarray(frames, dtype=float)
    ref = np.asarray(reference, dtype=float)
    if frames.ndim != 2 or len(ref) != len(frames) or len(ref) < 3 or np.any(ref <= 0):
        raise ValueError("Need at least three frames and positive beam reference")
    corrected = frames * (np.mean(ref) / ref[:, None])
    mu = np.mean(frames, axis=0)
    diffvar = np.var(np.diff(corrected, axis=0), axis=0, ddof=1) / 2
    good = np.isfinite(mu) & (mu > 0)
    ratio = np.full_like(mu, np.nan)
    ratio[good] = diffvar[good] / mu[good]
    finite = np.isfinite(ratio) & (ratio > 0)
    if np.count_nonzero(finite) < 3:
        raise ValueError("Insufficient positive data to estimate noise")
    interp = np.interp(np.arange(len(mu)), np.flatnonzero(finite), ratio[finite])
    g = median_filter(interp, size=smooth_bins, mode="nearest")
    variance_sum = g * np.maximum(np.sum(frames, axis=0), 1)
    centered = corrected - np.mean(corrected, axis=0)
    time_num = np.mean(centered[:-1] * centered[1:], axis=0)
    time_den = np.mean(centered**2, axis=0)
    rho_time = np.divide(time_num, time_den, out=np.full_like(mu, np.nan), where=time_den > 0)
    return dict(sigma=np.sqrt(variance_sum), g=g, raw_ratio=np.var(frames, axis=0, ddof=1)
                / np.maximum(mu, 1e-12), corrected_ratio=np.var(corrected, axis=0, ddof=1)
                / np.maximum(mu, 1e-12), rho_time=rho_time)


def circular_block_indices(n_frames: int, block_length: int,
                           rng: np.random.Generator) -> np.ndarray:
    if not 1 <= block_length <= n_frames:
        raise ValueError("Invalid block length")
    starts = rng.integers(n_frames, size=int(np.ceil(n_frames / block_length)))
    return ((starts[:, None] + np.arange(block_length)) % n_frames).ravel()[:n_frames]


def dl_apex(e0: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    c = 2 * np.asarray(e0)**2 - np.asarray(gamma)**2
    return np.sqrt((c + np.sqrt(c**2 + 12 * np.asarray(e0)**4)) / 6)
