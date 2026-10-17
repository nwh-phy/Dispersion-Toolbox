"""Three-way comparison of the 2d_h5 adaptive-binning results: 590 (1L), n0 (1L repeat), 20w (2L).
Usage: python3 b1_layer_compare.py <run_root>/stage5/   (numbers of reports/layer_1L_2L_20261017.md)
"""
import csv, sys, math
S = sys.argv[1]
sets = [('590 1L', 'e0_stage2_2d_h5_590_PL2_10w.csv'), ('n0 1L', 'e0_stage2_2d_h5_n0_PL2_10w_repeat.csv'),
        ('20w 2L', 'e0_stage2_2d_h5_no_PL2_20w_2film.csv')]
D = {}
for nm, f in sets:
    D[nm] = {round(float(r['q']), 5): r for r in csv.DictReader(open(S + f))}
    st = {}
    for r in D[nm].values(): st[r['status']] = st.get(r['status'], 0) + 1
    print(nm, 'status:', st)
def band(d, lo, hi, key, side=None):
    v = [float(r[key]) for q, r in d.items() if lo <= abs(q) < hi and r['status'] == 'converged' and (side is None or (q > 0) == (side > 0))]
    if not v: return float('nan'), float('nan'), 0
    m = sum(v) / len(v); sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else float('nan')
    return m, sd / math.sqrt(len(v)), len(v)
bands = [(0.0021, 0.0031), (0.0031, 0.0051), (0.0051, 0.0071), (0.0071, 0.0091), (0.0091, 0.0125)]
for key in ['apex_low', 'apex_high', 'apex_measured_LF', 'E0_low', 'E0_high']:
    print(f'\n{key}: mean +- sem (n) per band; then n0-590 and 20w-590')
    for lo, hi in bands:
        r = [band(D[nm], lo, hi, key) for nm, _ in sets]
        cells = ' | '.join(f'{m:6.0f}+-{s:3.0f} ({n:2d})' for m, s, n in r)
        print(f'  {lo:.4f}-{hi:.4f}: {cells} || n0-590 {r[1][0]-r[0][0]:+5.0f}   20w-590 {r[2][0]-r[0][0]:+5.0f}   20w-n0 {r[2][0]-r[1][0]:+5.0f}')
print('\n+q / -q per band (apex_low, apex_high):')
for nm, _ in sets:
    print(' ', nm)
    for lo, hi in bands:
        print(f'    {lo:.4f}-{hi:.4f}: low {band(D[nm],lo,hi,"apex_low",1)[0]:6.0f} / {band(D[nm],lo,hi,"apex_low",-1)[0]:6.0f}   high {band(D[nm],lo,hi,"apex_high",1)[0]:6.0f} / {band(D[nm],lo,hi,"apex_high",-1)[0]:6.0f}')
