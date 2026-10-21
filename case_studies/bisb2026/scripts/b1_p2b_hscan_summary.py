import csv, os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
R = os.path.join(ROOT, 'paper_results', 'b1_adaptq', '20261009T054142Z', 'p2b') + os.sep
runs = [(0.0005, 'fit_590_33db389c'), (0.001, 'fit_590_f594eda0'), (0.002, 'fit_590_679b7905'), (0.004, 'fit_590_91493ef7'), (0.0066, 'fit_590_79ee195c')]
print('h        variant        chi2     dchi2(vs best)  chi2/dof  holdout chi2/pt  holdout mean  neg.tails  lobo pred')
rows = []
for h, d in runs:
    lobo = list(csv.DictReader(open(R + d + '/p2b_lobo.csv')))[0]
    for r in csv.DictReader(open(R + d + '/p2b_fits.csv')):
        rows.append((h, r['variant'], float(r['chi2']), float(r['chi2dof']), float(r['holdout_chi2_per_point']), float(r['holdout_mean_residual_sigma']), r['negative_tails'], lobo['chi2n_predicted'] if r['variant'] == 'main' else ''))
best = min(x[2] for x in rows)
for h, v, c, cd, ho, hm, nt, lp in rows:
    print(f'{h:.4f}  {v:<13} {c:9.0f}  {c-best:+9.0f}      {cd:.3f}     {ho:8.1f}        {hm:+6.2f}      {nt:>4}     {lp[:5]}')
print('\nmain-variant curves (meV): apex and E0 of both branches at Q')
Qs = [0.001, 0.002, 0.003, 0.004, 0.005, 0.0065, 0.008]
for meas in ['apex_b1', 'E0_b1', 'apex_b2', 'E0_b2', 'G_b2']:
    print(f'  {meas}: ' + '  '.join(f'Q={q}' for q in Qs))
    for h, d in runs:
        C = list(csv.DictReader(open(R + d + '/p2b_curves_main.csv')))
        vals = []
        for q in Qs:
            k = min(range(len(C)), key=lambda i: abs(float(C[i]['Q']) - q))
            vals.append(f"{float(C[k][meas]):6.0f}+-{float(C[k][meas + '_err']):<3.0f}")
        print(f'    h={h:.4f}: ' + ' '.join(vals))
