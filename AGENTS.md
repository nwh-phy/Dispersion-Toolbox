# Project instructions

This file is the shared development policy for Codex and Claude Code. `CLAUDE.md` imports it; keep current rules here rather than maintaining a second copy. Daily usage is in `README.md`, detailed GUI help in `docs/toolbox/gui_workflow.md`, and dated software evidence in `reports/README.md`.

# Research progress and results alignment

## Start from the graduation thesis

The user explicitly requested on 2026-09-03 that discussions of this project's research progress, analysis results, interpretation, next steps, and publication positioning begin with their graduation thesis in the adjacent project.

- Thesis project: `C:/Users/HP/Desktop/vibecoding/本科毕业论文` on the Windows machine. The Mac has a clone at `~/Documents/NJUThesis` (GitHub `nwh-phy/NJUThesis`) whose last commit is 2026-05-15, older than the submission, and it has no submission package; use it only as an entry point and check passages against the submission.
- Main source: `thesis.tex`
- Research overview: `frontmatter/abstract-zh.tex`
- Results and evidence: `chapters/03-实验结果与分析.tex`
- Completed work and open questions: `chapters/04-总结与展望.tex`
- Methods and provenance details: `chapters/02-实验方法.tex` and `appendices/appendix-a.tex`
- Submitted PDF: `终稿提交包_221240048聂文辉_20260605/221240048聂文辉终稿.pdf`
- Archived submission sources: `终稿提交包_221240048聂文辉_20260605/221240048聂文辉终稿_LaTeX源码/`

Paths below the thesis project are relative to that project. For the completed thesis baseline, prefer the 2026-06-05 submitted PDF and its matching archived sources unless the user identifies a newer revision. The listed chapter paths also apply inside the archived source directory.

The root `thesis.pdf` predates the submitted PDF. A read-only comparison found that the root abstracts, methods, and conclusions match the submission after normalizing line endings, but other files differ, including a B3 interpretation/citation passage in the results chapter. Root sources are a convenient entry point; verify the relevant passage against the submission before treating it as submitted wording or evidence.

## Reading and verification order

1. Read the relevant thesis sections before assessing what has already been completed or proposing missing analyses. Preserve the thesis's distinction between established observations, candidate explanations, and unresolved questions.
2. Incorporate the user's subsequent corrections, verified sample information, and revised plans. The thesis is the starting baseline, not a reason to revert newer explicit instructions.
3. Follow the thesis's source chain into figures, tables, analysis outputs, scripts, and raw-data records when a specific result requires verification. Older `paper_results` reports (indexed in `paper_results/00_INDEX.md`), `archive/old_notes/PROJECT_MEMORY.md`, `archive/old_notes/WORKSPACE_INDEX.md`, and other model-generated notes are supporting historical records, not the default summary of current scientific progress.
4. If the thesis and upstream evidence disagree, identify the discrepancy and verify it explicitly. Do not assume that either a document's completion status or the model that generated an analysis determines scientific correctness.
5. This alignment instruction does not authorize raw-data processing, replacement figures, thesis edits, or publication. Keep those actions within the user's requested scope and follow the thesis project's own instructions when working there.

## Language: MATLAB only

The user explicitly required on 2026-10-24 that new data processing, fitting, simulations, figures and validation use the existing MATLAB toolbox. Do not add or run Python analysis code, and reuse the existing MATLAB functions before writing a new implementation.

- The Codex Python exploration of 2026-10-23/24 was removed on 2026-10-27 at the user's request. It survives only in the local tag `archive/codex-exploration-20261024` (`git worktree add <dir> archive/codex-exploration-20261024` to look at it) and is not evidence for the MATLAB pipeline. Its MATLAB low-branch analysis was kept: `b1_low_branch_matlab*` in `case_studies/bisb2026/scripts/`.
- The figure on page 4 of the 2026-10-25 weekly report, `01_low_component_count_dispersion.png`, came from that Python run. Its outputs are kept locally in `paper_results/b1_low_branch_python_20261024/`.
- Run a batch in one MATLAB process with a cell-array task loop and a per-task `try/catch`, not one MATLAB launch per parameter combination. Keep at most 2–3 MATLAB processes at once; repeated macOS launches leave large `com.mathworks.matlab.code_sign_clone` caches.

## Code organization

Agreed with the user on 2026-10-21, after the clean-up of the BiSb scripts (commits `db55a4d`, `a1fcebd`, `1628013`). The script index is `case_studies/bisb2026/scripts/README.md`.

1. A new round of analysis adds a config function and calls the existing code. Do not copy the previous round's script and edit the copy; put what changes into the config or into a new function.
2. When a helper would be copied into a second file, move it into its own function file instead: `src/` if it is generic q-EELS processing, otherwise the top level of `case_studies/bisb2026/scripts/`.
3. Pass the variants of a run as arguments, as a modified copy of the config struct (see `b1_p2b_slot_config` and `b1_p2b_fit(run_root, cfg, ocfg)`). Do not override a function with a same-name file placed earlier on the path, and do not pass settings through environment variables.
4. Commit at the end of every analysis stage; do not leave scripts untracked across sessions. Find the project root with `bisb_find_project_root`, never an absolute `/Users/...` path. The data and `paper_results/` exist only in the main checkout, so analysis scripts run there; a worktree must not read uncommitted files of another checkout.
5. When a line of analysis is finished, move its scripts into a `legacy_*` folder (`legacy_202605/`, `legacy_c20/`) and stop editing them.
6. `src/` holds the GUI toolbox and generic functions only; BiSb-specific code lives in `case_studies/bisb2026/`.
7. Before committing a change to shared code, run `tests/` and `case_studies/bisb2026/tests` and compare with the previous result. The latest full run is 2026-11-09, after the dataset-registry refactor: 214 core and 85 BiSb tests passed, with no failures or skips, on macOS MATLAB R2026a in the main checkout. The last GUI integration baseline (2026-10-23: 205 core, 78 BiSb) is in `reports/gui_integration_20261023.md`. Do not apply these counts to older commits or data-free worktrees.
8. The momentum calibration was corrected by a factor of 10 on 2026-05-14 (dq 0.0005 1/A per pixel at 10w, 0.00025 at 20w). A q threshold in 1/A written before that date selects 10x more channels than intended; divide it by 10 when reviving older code.
9. Every acquisition is registered once in `case_studies/bisb2026/scripts/bisb_datasets.m` (folder, file, dq, q = 0 channel, slot acceptance h, SHA-256, sample, layers). Code takes the dataset id; do not write data folders, dq or q = 0 into scripts. For new data: put the raw files in a git-ignored data folder at the project root, add the registry entry, run `bisb_quicklook('<id>')` for the standard first figures.
10. Search before writing: check `src/` and the script index for an existing function, and extend it with an option or a config instead of copying its body. Raw sequences go through `qe_align_sequence_frames` + `qe_transport_frames` (cached map: `bisb_frame_map`); spectra and figures through `qe_pm_spectra`, `qe_zlp_center`, `qe_powerlaw_tail`, `qe_peak_apex`, `qe_plot_heatmap`. One-off scripts that only make a figure live next to that figure in `paper_results/`, not in the repository.

## GUI work

Updated with user approval on 2026-10-23 after integrating `gui-joint-fit`.

1. The main checkout stays on `feat/zlp-joint-fit` and is the daily entry point for both GUI and research. The accepted GUI changes are also on `main`; research-only scripts must not be merged into `main` without approval.
2. Future generic GUI work starts from `main` on a short-lived branch/worktree. After verification and user authorization, merge into `main`, then integrate into the research branch. Retire completed worktrees after checking for uncommitted and ignored files. Do not recreate the retired `.claude/worktrees/gui-joint-fit` as a permanent GUI entry point.
3. Shared fitting functions must remain consistent across the stable and research lines. Prefer integrating reviewed commits over maintaining manual copies.
4. Run both `tests/` and `case_studies/bisb2026/tests/` after integration. Worktrees without data may skip real-data tests; run those in the main checkout or through a temporary external script, and keep generated results out of git.
5. A new GUI control must stay optional when older history and session files are restored: a missing field takes the control's default.

## GUI implementation and result records

- `src/ui/qe_browser_ui.m` declares controls; browser state and callbacks live in `src/interactive_qe_browser.m`. Preserve this separation.
- In overlapping energy windows, seed-tagged peaks keep their branch; untagged peaks go to the later matching branch. This is the existing tested rule.
- ZLP-joint fits accept model E0; BG + peaks accepts the peak apex. Preserve `energy_kind` and existing session provenance. GUI/API agreement does not establish scientific model validity.
- Keep dated reports immutable as historical evidence; update their navigation index instead of rewriting old test counts as current.
- Local result classification and recovery are indexed by `paper_results/00_INDEX.md`. An old run may still be an input. Check consumers before moving it; preserve hashes, archives, and a restoration route. Raw data and thesis evidence are not disposable build artifacts.
