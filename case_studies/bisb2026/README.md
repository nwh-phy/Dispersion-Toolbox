# BiSb 2026 Case Study

This directory contains the BiSb-specific analysis runners and workflow tests
used for the 2026 q-EELS case study. It is not the core toolbox API.

Core reusable code stays in `src/`. Scripts here are examples layered on top of
the toolbox and may require ignored local inputs such as `20260120 BiSb/` and
ignored outputs under `paper_results/`.

## Analysis Convention

- The primary peak-position evidence uses the no-BG, area-normalized Fano
  extraction route.
- Background-subtracted analyses are diagnostic or robustness checks.
- Labels such as `10w` and `20w` refer to defocus settings, not laser/electron
  beam power.
- `590 PL2 10w` and `n0 PL2 10w repeat` are treated as 1film repeats.
- `no PL2 20w 2film` is treated as the 2film data set.

## Commands

From the repository root:

```matlab
startup; runtests('tests')
runtests('case_studies/bisb2026/tests')
```

For individual case-study exports:

```matlab
addpath('case_studies/bisb2026/scripts')
output = run_b1_physical_fit_analysis();
output = run_b1_physical_fit_enhancements();
output = run_b3_scatter_only_export();
output = run_q015_linewidth_by_session_export();
```

Result folders remain under the ignored `paper_results/` tree so existing paths
and generated artifacts are not broken by this cleanup.
