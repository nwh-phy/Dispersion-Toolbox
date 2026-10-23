# q-EELS processing exploration, 2026-10-23

Completed findings: [RESULTS.md](RESULTS.md). The strongest processing result
is relative q-drift correction validated against temporal holdout spectra.
The model/uncertainty limitations and proposed GUI adoption are explicit.

B1 composition follow-up: [the research judgment](composition/REPORT.md)
and [the reproducible Python comparison](composition/README.md). Its full-Q
measurement model tests two spectral contributions against a single broad
response; spectral components are not automatically distinct plasmons.

This isolated experiment starts at research commit `3d66386`. It processes
the existing C20 MoS2/Bi-rich BiSb/MoS2 sequences; it does not assume pure Bi,
a crystal structure, independently verified layer counts, carrier density,
or an absolute dielectric function. The local thesis is an older research
baseline, not the submitted thesis or a current momentum calibration.

## Questions and fixed scope

1. Can raw frames be processed reproducibly without historical MAT outputs?
2. How much do energy drift, count binning and measured noise affect the
   broad low-loss band and its two-component description?
3. Are statistical uncertainty and sensitivity to the measurement model
   distinguishable? Do results repeat across signed q and the three sessions?

The first milestone is a tested raw-frame pipeline, with count conservation
on the retained common support, explicit axes, masks and input hashes.
The next milestone uses the existing MATLAB joint fitter for matched-window
one/two-component comparisons, sensitivity fits and temporal validation.
Independent synthesis/injection tests quantify recovery and false splitting.
Frame/block resampling is conditional on the chosen model; changing aperture
or peak shape is reported separately as model sensitivity.

## Reproduction

Use Python >=3.11 and MATLAB with Optimization Toolbox. In the repository root:

```sh
python3 -m venv .venv
.venv/bin/python3 -m pip install -r case_studies/bisb2026/exploration/requirements.txt
.venv/bin/python3 -m pytest tests/python -q
.venv/bin/python3 case_studies/bisb2026/exploration/run.py prepare \
  --data-root '/path/to/20260120 BiSb' --output '/path/to/new/output'
```

For the complete pipeline (preparation, checks, all variants, synthesis,
bootstrap, plots and final raw-input integrity audit), run:

```sh
.venv/bin/python3 case_studies/bisb2026/exploration/run_all.py \
  --data-root '/path/to/20260120 BiSb' --output '/path/to/new/output' \
  --matlab '/Applications/MATLAB_R2026a.app/bin/matlab' --jobs 3
```

Repeat with `--resume` after interruption. A changed preparation config,
raw input, fitter or cached synthetic model is rejected rather than silently
mixed into the same run. Per-fit checkpoints keep numerical failures and
boundary states. Python does frame processing, minimal physics and reporting;
MATLAB supplies the existing joint optimizer. This is not a Python-only port
of that optimizer.

Individual MATLAB stages can also be called after `startup` and adding this
directory to the path:

```matlab
bisb_explore_fit(output_dir, 'check');
bisb_explore_fit(output_dir, 'fit', '590'); % also 'n0', '20w'
bisb_explore_fit(output_dir, 'bootstrap', '590');
bisb_explore_synthetic(fullfile(output_dir, 'synthetic'));
```

`summarize.py --output ...` rebuilds tables and plots from saved fits.
`verify_run.py --output ...` re-hashes the original files and checks fit
coverage/signatures. [METHODS.md](METHODS.md) gives equations, assumptions,
energy definitions, verified literature and the route to band/DFT response.
The figures linked by RESULTS.md are local reproduction outputs; raw arrays,
fit checkpoints and generated figures are deliberately not Git assets.

The data root is an explicit read-only input. Outputs must be outside that
root. The input JSON gives energy dispersion and exposure; the checked
dataset configuration gives q calibration and its provenance. Filenames are
never used to infer physical units. Raw frames are memory-mapped, and only
selected spectra are materialised. Generated data, fits and figures remain
local under `paper_results/` and are not published by the run.

The metadata energy origin is checked against the measured elastic reference
in the declared raw search region before narrow-window frame alignment.
In these sequences the elastic peak is about 160–204 meV below the metadata
zero; blindly using a +/-100 meV window around that zero rejects normal frames.
The empirical origin and every frame shift are saved separately from the
metadata. The 30 kV beam setting is also checked against acquisition metadata.

`requirements-lock.txt` records the exact packages used in this run. The
original raw-file SHA-256 hashes, configuration hash and source hashes are
recorded in the output manifest.

The Python frame alignment is the same noncircular integer shift used by
`qe_zlp_integer_align.m`: one elastic-reference shift per frame, common to
all q. It does not flatten the q-dependent elastic line. The first milestone
also compares the independent Python and MATLAB implementations.

## Decision rules

- Fixed sampling and variants are specified before examining fit results.
- Keep failed, boundary and collapsed fits in the output tables.
- Distinguish observed count-band features, fitted component E0, component
  loss-function apex, and physical mode assignments.
- No lower chi-square alone is accepted as proof of a second physical mode.
- No forced plus/minus-q symmetry, layer dependence, or gapless dispersion.
- Reuse `qe_zlp_joint_fit` and `qe_kinematic_prefactor` without changing their
  scientific defaults globally; exploratory choices are explicit arguments.
- Commit code at milestones; do not merge the research and GUI branches.

Known prior evidence motivating the experiment: the repository's S0
simulation found model bias despite small residuals; the h-scan found
aperture-sensitive E0. Neither recorded forward fit is treated as ground
truth about the sample.
