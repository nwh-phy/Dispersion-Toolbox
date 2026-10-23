"""Conservative subpixel bin transport for calibrated momentum axes."""
import numpy as np


def fractional_bin_weights(q_axis: np.ndarray, center: float, width: float,
                           offset: float = 0.) -> np.ndarray:
    """Conservative overlap of a translated interval with detector pixel cells.

    Counts are treated as uniform inside each pixel. Weights are overlap/dq,
    not unit-normalised; their sum is the number of channels in the bin.
    Adjacent rebinned intervals share input pixels and become correlated.
    """
    q=np.asarray(q_axis,dtype=float)
    if q.ndim!=1 or len(q)<2 or not np.all(np.isfinite(q)) or not np.all(np.isfinite([center,width,offset])):
        raise ValueError('Finite one-dimensional q axis and finite bin parameters required')
    dq=float(np.median(np.diff(q)))
    if dq<=0 or not np.allclose(np.diff(q),dq) or width<=0:
        raise ValueError('Positive width and a uniform increasing q axis required')
    left=center+offset-width/2;right=center+offset+width/2
    if left<q[0]-dq/2-1e-12 or right>q[-1]+dq/2+1e-12:
        raise ValueError('Translated bin extends outside acquired support')
    return np.maximum(0,np.minimum(q+dq/2,right)-np.maximum(q-dq/2,left))/dq
