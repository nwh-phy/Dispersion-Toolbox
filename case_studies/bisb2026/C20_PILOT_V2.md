# C20 component pilot (P0–P3)

From the repository root in MATLAB:

```matlab
addpath('case_studies/bisb2026/scripts');
out = run_b1_component_pilot_v2();
```

The entry creates a unique timestamp/config-hash/HEAD directory under
`paper_results/b1_components_v2`. An existing explicit output directory is
rejected. Native NPY/JSON sources are required and read without using or writing
source caches. Historical MAT/CSV and the 590 history are inventoried separately;
there is no silent processed-data fallback.

The registry supplies session identities and dq. The pilot uses ±0.015 Å⁻¹ and
q_skip=0.0005 Å⁻¹, and keeps signs separate. It saves N=1/3/5 bins for all sessions
and N=6/10 for the 20w physical-width comparison. Only 590 is fitted at three
preselected low/mid/high |q| locations; this is not frozen cross-session validation.

The raw adapter excludes an entire native q column containing any nonfinite or
uint32 saturation-sentinel sample before calibration. Excluded output columns are
NaN. The 2026-09-12 input audit found that an unmasked sentinel in column 461
otherwise corrupts both zero definitions. L1 retains acquisition-side detector
corrections and sums frames without per-q alignment, normalization, denoising,
presubtraction or deconvolution. Frame ZLP diagnostics are saved.

`qe_prepare_count_bins` retains the existing extractor's mean/member field
conventions and can be selected there with `binning_policy='fixed_nonoverlapping'`.
Its default legacy adaptive policy is unchanged. Missing data exclude whole q
columns; groups never cross a sign, native-column gap or invalid channel. Measurement
variance is unknown unless explicitly supplied. Member scatter is a distinct field.

`qe_compare_component_models` reuses `peak_models` and `measure_peak_fwhm` with a
bounded multistart path that never deletes a weak component or invokes an unbounded
fallback. `lorentz` retains its Drude–Lorentz meaning; `lorentz_symmetric` uses an
area-normalized symmetric Lorentzian. Both explicit n=1/2 models fit a joint power-law
background in the same window. Unweighted LS is descriptive, not a calibrated
likelihood. Internal solver scaling does not normalize or overwrite the input spectra.

The main window is 300–1800 meV, with named 300–2000/2100 sensitivities. Each model
uses 24 starts; all candidates, boundaries, solver messages and curve arrays are saved.
No Fano default, trend constraint, jump repair, downstream physical fit or index/thesis
update is executed. Scientific labels remain unresolved pending P4/P5. Native DL
Gamma, symmetric FWHM and finite-window floor-based width diagnostics are distinct.

Tests:

```matlab
startup;
r = runtests('tests/test_b1_component_pilot_v2.m');
assert(all([r.Passed]));
```

After a completed run, saved arrays can be independently checked with
`audit_b1_component_pilot_output_v2(out)`. This also reruns one fixed-start pair.
It does not establish statistical identifiability, a false-positive rate, or an
intrinsic lifetime. The broader pre-existing test failure in
`QeAssignPeakBranchesByWindowsTest/testOverlappingWindowsAssignEachPeakOnlyOnce`
is outside this pilot's local energy-order labeling path.
