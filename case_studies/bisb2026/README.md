# BiSb 2026 Case Study

This directory contains the BiSb-specific analysis runners and workflow tests
used for the 2026 q-EELS case study. It is not the core toolbox API.

Core reusable code stays in `src/`. Scripts here are examples layered on top of
the toolbox and may require ignored local inputs such as `20260120 BiSb/` and
ignored outputs under `paper_results/`.

## Thesis-era Analysis Convention

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

Current analysis entry points and configurations are indexed in
[`scripts/README.md`](scripts/README.md). The earlier export runners have moved
to `scripts/legacy_202605/` and are frozen; the former commands that added only
`scripts/` no longer locate those runners. Choose the current configuration and
run directory from the index rather than rerunning a historical export by default.

The GUI/API numerical acceptance tests are in `tests/test_gui_joint_fit_realdata.m`.
They compare six representative channels from the three local eq3D datasets,
including both signs of q. They require local data and report an explicit skip
when it is unavailable. Agreement checks the interface, not the physical validity
of the single-channel model.

Result folders remain under the ignored `paper_results/` tree so existing paths
and generated artifacts are not broken by this cleanup.
