# GUI session persistence acceptance — 2026-10-22

## Scope

This stage builds on `68fb373` on `codex/gui-workflow-acceptance`.
It adds complete projected q–E session persistence and accepted single-spectrum
fit provenance. It does not add batch ZLP joint fitting or reinterpret BiSb data.

## Saved state and compatibility

The History panel now has **Save Session** and **Load Session**. A version-1 MAT
file contains the `session` variable with the following fields:

- `format`, `version`, `saved_at`: explicit format and timestamp.
- `sources`: source path, byte count, SHA-256 captured when data was imported;
  an NPY companion JSON is included when present.
- `data`: the dataset, projected eq3D/physical counts, calibrated axes and stage.
- `data_hashes`: fingerprints of both counts and axes, including dq/q-zero.
- `snapshot`, `history`: the current live controls and earlier operation history.
- `analysis`: signed points, branch matrices, correction log, dispersion and
  automatic fit results, and the accepted-fit audit archive.
- `pending_fit`, `display_fit`: numeric fit output and plotted data, preserving
  an unaccepted fit or an accepted overlay without re-running optimization.

Comparison views are rebuilt from the saved physical data and settings. Graphics
handles, preview caches and solver closures are excluded. Numerical prefactor
values and the kinematic configuration remain in the joint-fit record; its
anonymous solver prefactor closure is not persisted. MATLAB value tables used by
branch assignment are preserved and checked recursively.

A session can be opened into an empty browser or passed directly as the first
argument to `interactive_qe_browser`. Old history and curated point files remain
supported separately. Old points without fit provenance keep an empty archive.
Loading a new dataset resets the old analysis state so it cannot be silently
saved as belonging to the new source.

## Integrity and failure behavior

All structural checks and source verification happen before changing the GUI.
Embedded data/calibration fingerprints must match. A source that exists at its
recorded path must still have the recorded SHA-256; a same-size modification is
also detected. An unavailable source is allowed because the projected counts are
embedded, and the browser explicitly reports offline recovery.

Controls must restore exactly. Any restoration/rebuild/rendering error rolls
back the old data, controls and displayed/pending fit. A saved pending fit must
match the selected signed q, calibration, data fingerprint, fit settings and
peak-validity acceptance mask. Accepted overlays restore with Accept disabled.

Saves validate first, write a temporary file in the target directory, then replace
the target. An invalid save leaves an existing file intact; the source data file
cannot be the save target. Cancelling a dialog does not modify live analysis.

## Fit provenance

Accept stores a numeric fit record alongside its source identity, settings at
fit time, signed q/channel, calibration, signal/fit window, input kind, timestamp,
accepted peak indices and energies. Model E0, loss-function apex, quality flags
and full predicted/component curves remain separately available in the result.
For ZLP joint the fit window includes the gain side and ZLP; its signal window is
recorded separately. Preprocessed background-plus-peaks fits identify that input
kind instead of claiming raw counts.

The archive records historical acceptance. Deleting, moving, replacing or
reassigning a point later does not make an old record a new validation of the
current point. Loading an old point file does not invent missing provenance.
`Save Curated` also includes these records in `accepted_fits`.

## Export and laptop display fixes

A visual export check reproduced blank PNG/PDF pages when an entire UIAxes was
copied into a hidden ordinary figure. `qe_export_browser_panel` now recreates
ordinary axes and copies graphics children, including hidden peak markers, while
rebuilding labels/legend and preserving scales, limits and colormap. A fixed
paper canvas preserves the selected aspect ratio without requiring newer
`exportgraphics` sizing options than the documented R2024b baseline. PDF failures
are no longer silently counted as successful panel exports.

Exported long titles wrap and reserve measured headroom. Single-fit Auto Y also
reserves room for E0/apex annotations. The control grid is scrollable so its fixed
columns remain reachable when the window is narrower than the original layout.

## Validation

MATLAB R2026a full regression passed and the runner exited successfully:

| Suite | Previous stage | This stage | Failed | Incomplete/skipped |
|---|---:|---:|---:|---:|
| `tests/` | 177 | 191/191 | 0 | 0 |
| `case_studies/bisb2026/tests/` | 78 | 78/78 | 0 | 0 |

All 269 cases passed. This adds eleven GUI behavior cases and three I/O cases.
Temporary evidence: `/tmp/dispersion_session_after_export.log` and
`/tmp/dispersion_session_after_export.mat` (`core`, `bisb`). The added-case run and
separate automatic-fit table check are in `/tmp/dispersion_session_new_cases.log`
and `/tmp/dispersion_session_auto.log`. Temporary files can be removed by the OS.

Reproduction from the project root:

```matlab
startup
core = runtests('tests');
bisb = runtests('case_studies/bisb2026/tests');
assert(all([core.Passed]) && all([bisb.Passed]));
```

The added session behavior tests exercise a nondefault dq and negative q selection,
branch correction and targeted dispersion fits, Comparison normalization,
closing/reopening into an empty browser, pending and accepted overlays,
background-plus-peaks fits, real automatic-fit detail tables, offline sources,
changed sources, corrupted counts, incompatible versions, unsupported controls,
wrong pending channels/masks, rollback after UI failure, cancelled operations,
failed-save preservation and old point-file compatibility.

Independent I/O tests check the standard SHA-256 digest for all byte values
0–255, same-size source changes, unavailable sources, and count/axis changes.
The existing real-data GUI/API regression still checks six representative
channels without reclassifying invalid peaks as scientific evidence.

Native and rendered checks used a synthetic dataset, not new research analysis:

- In interactive MATLAB, the native **Save Session** dialog was opened, cancelled,
  opened again and used to write a 1.6 MB session. A separate process reopened it
  and verified a nonempty pending fit. File: `/tmp/dispersion_session_native/synth 10w/qe_session.mat`.
- The native **Load Session** chooser was opened. Its final file selection, and
  the Export Panels selection/folder dialogs, were not completed as a single
  computer-controlled workflow. Path-based load callbacks and the panel export
  backend are exercised independently; this is not claimed as full native-dialog
  acceptance.
- The actual restored UI was inspected for fit curves and annotation headroom.
  Narrow-window behavior beyond scroll access is not comprehensively validated.
- The shared panel export backend generated single-spectrum and Comparison PNG
  and PDF files. PDF pages were rendered with Poppler for visual inspection.
  The final single-spectrum and Comparison PDFs were checked after Poppler
  rendering: curves, hidden peak markers, labels, legend, wrapped titles and
  margins were present. Temporary artifacts and harnesses are in
  `/tmp/dispersion_session_native/`.
- A formal GUI export case checks that the restored spectrum PNG is nonblank,
  retains the requested 4:3 ratio, writes a PDF and leaves the pending fit intact.

Only R2026a was executed. R2024b compatibility was retained by avoiding newer
export sizing arguments, but an installed R2024b runtime was not available for a
separate run.

## Limits

A session embeds the current projected q–E data, not the original raw 4D tensors
or a reconstruction of their crop/alignment/binning intermediate states. Missing
raw sources therefore do not prevent reviewing/refitting the saved projection,
but do prevent re-running the original raw-data pipeline from this file alone.

The current format is version 1. Unknown versions are rejected; there is no
cross-version migration promise. A matching source and GUI/API numerical
agreement do not establish physical identifiability or replace research
acceptance criteria. Scientific thesis/results interpretation is outside scope.
