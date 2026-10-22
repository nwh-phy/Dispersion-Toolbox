# GUI workflow acceptance — 2026-10-22

## Scope and branch alignment

Work is on `codex/gui-workflow-acceptance`, based on research branch
`feat/zlp-joint-fit` at `f2105ef`. The existing GUI feature commit `ee65ba5`
was applied as `977e517`. The two branches' `qe_zlp_joint_fit.m` and
`qe_kinematic_prefactor.m` were identical before integration, so the GUI uses
the existing fitting core. The research scripts, calibration fixes, and legacy
organization are retained. No main-branch merge or remote push was performed.

This stage covers the existing single-spectrum workflow. It does not add batch
joint fitting or reinterpret BiSb results.

## Reproduced failure and fixes

The new `testFailedRefitCannotAcceptThePreviousResult` failed on the existing GUI:
after a successful ZLP joint fit, an unsupported model caused the next fit to fail,
but the old pending result, enabled Accept button, and old overlay remained.
Calling Accept still added the previous point. The reproduction log is
`/tmp/dispersion_gui_repro.log` (temporary local evidence).

The GUI now clears the pending fit and overlay before a fit attempt, on selected
channel/view refresh, and when fit parameters change. Failed attempts cannot
accept the previous result. Loading a different points file also clears dispersion
curves calculated from the previous points.

History now records the normalization method and auto/manual stack-offset mode.
Restoration applies dq to the dataset and q axes before resolving the saved signed
q selection. An old file without the normalization method uses the historical
`Area` default; the missing original choice cannot be recovered from that file.
Failed restoration is reported rather than silently labelled successful.

## Behavioral coverage

The existing GUI hooks now also expose path-based save/load/export operations.
They execute the same callbacks as the buttons and bypass only file selection.
An empty hook path follows the cancellation path.

`tests/test_gui_joint_fit_smoke.m` covers:

- The existing synthetic single-peak fit, acceptance, old method, and unsupported
  peak-model checks.
- Failed refits, channel changes, and fit-setting changes clearing pending results.
- Saving, closing, and reopening with normalization method, dq, selected signed q,
  stack mode, and joint-fit settings restored.
- Legacy snapshots without the newer fields and a broken snapshot reporting failure.
- Accepting a fit into an explicit branch, saving points and history, reopening,
  then recovering signed coordinates, branch matrices, correction log, and target.
- Exporting from Comparison view while retaining physical raw data and axes,
  with preprocessed values equal to a direct preprocessing call.
- Cancelled file operations preserving curated points, correction log, and history.

`case_studies/bisb2026/tests/test_gui_joint_fit_realdata.m` compares GUI and direct
function results on six representative channels:

| Local session | q values (Å⁻¹) |
|---|---|
| 590 10w | −0.002, +0.005 |
| n0 10w repeat | −0.005, +0.002 |
| 20w 2film | −0.002, +0.010 |

Both paths use the same raw counts, two Drude-Lorentz peaks, initial E0 values,
300–1800 meV signal window, 30 kV kinematic factor, dq, and Poisson weights.
The tests compare E0, width, loss-function apex, predicted curve, validity flags,
and reduced chi-square. Local-data tests explicitly skip when input files are absent.

## Reproduction

Run from the project root in MATLAB R2026a:

```matlab
startup
core = runtests('tests');
bisb = runtests('case_studies/bisb2026/tests');
assert(all([core.Passed]) && all([bisb.Passed]));
```

For the bounded GUI acceptance alone:

```matlab
startup
runtests('tests/test_gui_joint_fit_smoke.m')
runtests('case_studies/bisb2026/tests/test_gui_joint_fit_realdata.m')
```

## Validation results

The focused GUI run passed 10/10 tests before the additional malformed-history
case and explicit legacy-normalization check were added. The standalone real-data
run passed 3/3 cases (all six channels, no skips). The final suites below include
the additional checks.

The pre-fix suite sizes at `977e517` were independently discovered as 166 core
tests and 75 case-study tests. Both baseline component runs completed without a
failure or incomplete-test summary. The baseline shell command then exited with
an array-dimension error while concatenating its two result arrays; it did not save
the planned aggregate MAT file. This is a runner-recording issue, not a passing
aggregate result. Evidence: `/tmp/dispersion_gui_baseline.log` and
`/tmp/dispersion_gui_discovery.log`.

The final runner saved the two result arrays separately and asserted every test
passed. It exited successfully:

| Suite | Before fixes (discovered count) | Final passed | Failed | Incomplete/skipped |
|---|---:|---:|---:|---:|
| `tests/` | 166 | 174/174 | 0 | 0 |
| `case_studies/bisb2026/tests/` | 75 | 78/78 | 0 | 0 |

The final run adds eight GUI behavior cases and three real-data cases. Its 252
tests all passed, including the old-method and historical analysis regression
checks. Structured results are in `/tmp/dispersion_gui_final.mat` (`core`, `bisb`);
the run log is `/tmp/dispersion_gui_final.log`. These temporary local artifacts
are not committed and may be removed by the operating system.

In the real-data comparison, the upper peak of two low-q channels reached the
signal-window boundary and was marked invalid by both interfaces. Numerical
agreement includes this rejection; it does not promote that peak to accepted
scientific evidence.

## Limits

These tests run the actual MATLAB GUI and callback logic, but bypass native file
dialogs. They do not verify manual dialog interaction, visual layout across screen
sizes, image/PDF export appearance, or every preprocessing combination.
Input round trips here use eq3D datasets; reconstruction of raw 4D crop,
alignment, and binning stages is not covered.

History and points remain two files. To recover a curated workflow, load the input
data, load points, then load history. Pending single-spectrum fits are not persisted.
No unified session format is introduced in this stage.

GUI/API numerical agreement does not establish physical identifiability. A
single-channel joint fit is not P2b's shared multi-channel acceptance model, and a
loss-function apex is not the observed-spectrum peak or model E0. Physical claims
and the thesis submission remain outside this engineering acceptance.
