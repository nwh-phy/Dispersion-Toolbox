# C20/B1 v4 actual execution and review
This is a new run. The v3 complete labels are not accepted as evidence: v3 lacked new fitted arrays, component mapping, and applied A1. No old run was overwritten.
## Inputs and identity
590 native NPY/JSON selected by session_id and extension, hashes matched parent manifest before use and raw after use. Native sequence-q-E is 300 x 512 x 1028. Invalid q mask independently reproduced parent column 461; sentinel mechanism is not inferred from its value. Native dtype/header and detector mask maxima are recorded.
Camera metadata contains is_sequence, nimages=300, exposure=10 s and frame_time=10 s. Fixed position, per-frame timestamps and sample stability are not verified. Sequence blocks are not presumed independent.
## Actual boundary audit
108 parent fits / 2592 candidates / 16848 parameter rows / 702 selected rows. Every candidate has its own raw-to-energy-order map and native-unit conversion. All selected native curves and objectives reconstruct. Legacy boundary flags: zero mismatches.
98 selected fits have at least one boundary: background r lower in 96, background B lower in 52, r upper in 1, E0 upper in 2. Categories overlap; no selected width or peak-amplitude bound hits. Background boundary does not itself disqualify peak energies, and absence of peak bounds is not identifiability.
## Real fitting and saved data
Legacy reconstruction: 12 models, 18 components, 288 candidates; selected boundary count 10. Mode 2 has independent model streams plus 12 asymmetric weak-component starts for n2. Mode 3 adds a nonnegative constant to BOTH hypotheses. Total: 36 models, 54 component records, 1008 optimized candidates, 18 separate H0 feasible witnesses.
Across solver/background comparisons, largest absolute normalized-objective change from legacy is 2.2732e-09. This does not establish decomposition uniqueness. Sparse profiles have 27 points, 0 with no converged optimizer start.
All comparisons use the same three N=3 centers, 300-1800 meV, unchanged width bounds and no physical prior. Width columns distinguish native Lorentzian FWHM/DL Gamma from zero-baseline finite-window bracketed widths. Areas refer to Wref=[300,1800] and have ordinate*meV units, not electron counts. Native A is not comparable between DL and symmetric models.
## Real frames and A1
A0 native-frame sum versus parent L1: maximum error 0. Fixed elastic window [-100,100] meV and reference |q|<=0.001 A^-1 produced 300 valid integer offsets. Correction is minus observed offset. Original member spectra, corrected members, support, all shifts, 18 actual block sums and per-frame means are in sequence_block_spectra.mat.
Reference-band cumulative zero-baseline ZLP FWHM A0/A1: 25.8157 / 23.013 meV. These are sample elastic reference diagnostics, not independent response calibration or intrinsic-width correction. Shift lag1=0.96383, linearly detrended=0.95362: descriptive quantized correlations, not effective independent frame counts.
## Minimal P5
Four named hypothetical scenarios were fixed before trials: static single, q-mixed single, overlapping broad double, background mismatch. Each has 10 smoke draws then 100 different-seed pilot draws. Detector noise was NOT fitted: independent Gaussian noise at 1 percent max ordinate is only an engineering assumption. Simulation search uses 4 independent starts plus 4 n2 extras, explicitly less than the real-data search budget. No experimental false-positive rate or CI is claimed.
The predeclared engineering split flag requires relative objective gain >0.1, reference-window component fraction >0.05 and center separation >4 meV. This is NOT a calibrated detection rule. Failures and boundary/collapse are retained. Wilson intervals describe only the hypothetical trial proportion. Generated H0 means never include the observed dual-structure mean or uncentered experimental residual.
| Pilot scenario | Trials | Engineering split flags | Failures | n2 boundary |
|---|---:|---:|---:|---:|
| static_single | 100 | 0 | 0 | 24 |
| q_mixed_single | 100 | 100 | 0 | 0 |
| overlapping_double | 100 | 100 | 0 | 0 |
| background_mismatch | 100 | 94 | 0 | 80 |
## Verification and limits
New contracts (10), simulation-generator tests (2), legacy pilot tests (7): 19 passed, 0 failed. Separate persisted-array tests recompute 36 fits and 1008 candidates and every actual block. The archive is additionally extracted and checked by the packaging entry; its result is saved alongside the ZIP.
No full repository test suite was run. The historical overlapping branch-window test is outside this local sorting chain and remains a known unaddressed regression. P5 experiment-calibrated inference, frame bootstrap, dense profile coverage, mirrored-q fits, q/time single-mode mixture exclusion, and physical attribution remain unperformed.
Minimum additional information: confirm whether all 300 frames observe the same location and whether beam/sample changes are expected; same-condition elastic/dark references are needed for intrinsic response/noise claims. Existing blocks show ordered shape/intensity changes, so do not blindly iid-resample them.
Visual spot checks: actual paired-fit legends and block/A0-A1 plots were rendered and inspected; arrays are the authoritative payload.
## Recompute without raw
From the project root: addpath('case_studies/bisb2026/scripts'); result=c20_v4_verify_packet('<extracted packet folder>');
A full new run is run_b1_component_v4(); it is not required to review this packet. Source snapshots identify core, P5 and finalization stages separately.