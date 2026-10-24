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

## Code organization

Language requirement, explicitly confirmed by the user on 2026-10-24:
new data processing, fitting, simulations, figures and validation use the
existing MATLAB toolbox. Do not add or run Python analysis code. Reuse the
existing MATLAB functions before introducing a new implementation. The
Python exploration is a frozen historical record; its test results do not
validate the MATLAB pipeline or GUI. Preserve its recorded outputs unless
the user separately authorizes cleanup.

Run a batch in one MATLAB process with a cell-array task loop and per-task
`try/catch`; do not launch MATLAB for each parameter combination. If parallel
processes are needed, keep at most 2–3 running at once. Repeated macOS launches
can leave large `com.mathworks.matlab.code_sign_clone` caches. The user has
authorized cleaning that specific cache only after all MATLAB processes exit.

Agreed with the user on 2026-10-21, after the clean-up of the BiSb scripts (commits `db55a4d`, `a1fcebd`, `1628013`). The script index is `case_studies/bisb2026/scripts/README.md`.

1. A new round of analysis adds a config function and calls the existing code. Do not copy the previous round's script and edit the copy; put what changes into the config or into a new function.
2. When a helper would be copied into a second file, move it into its own function file instead: `src/` if it is generic q-EELS processing, otherwise the top level of `case_studies/bisb2026/scripts/`.
3. Pass the variants of a run as arguments, as a modified copy of the config struct (see `b1_p2b_slot_config` and `b1_p2b_fit(run_root, cfg, ocfg)`). Do not override a function with a same-name file placed earlier on the path, and do not pass settings through environment variables.
4. Commit at the end of every analysis stage; do not leave scripts untracked across sessions. Find the project root with `bisb_find_project_root`, never an absolute `/Users/...` path. The data and `paper_results/` exist only in the main checkout, so analysis scripts run there; a worktree must not read uncommitted files of another checkout.
5. When a line of analysis is finished, move its scripts into a `legacy_*` folder (`legacy_202605/`, `legacy_c20/`) and stop editing them.
6. `src/` holds the GUI toolbox and generic functions only; BiSb-specific code lives in `case_studies/bisb2026/`.
7. Before committing a change to shared code, run `tests/` and `case_studies/bisb2026/tests` and compare with the previous result. Both suites passed in full on 2026-10-21 (159 and 75 tests).
8. The momentum calibration was corrected by a factor of 10 on 2026-05-14 (dq 0.0005 1/A per pixel at 10w, 0.00025 at 20w). A q threshold in 1/A written before that date selects 10x more channels than intended; divide it by 10 when reviving older code.

## GUI work

Agreed with the user on 2026-10-22.

1. GUI work (`src/interactive_qe_browser.m`, `src/ui/`, the browser files in `src/io/`, and their tests) goes on the branch `gui-joint-fit`, which starts from `main` and is checked out in the worktree `.claude/worktrees/gui-joint-fit`. Do not put GUI work on `feat/zlp-joint-fit` or start another GUI branch from it: the research scripts of `feat/zlp-joint-fit` must not enter `main`.
2. The main checkout stays on `feat/zlp-joint-fit`; do not switch it to another branch. Do not merge `gui-joint-fit`, `main` and `feat/zlp-joint-fit` into each other without the user's approval.
3. `src/fitting/qe_zlp_joint_fit.m` and `src/fitting/qe_kinematic_prefactor.m` exist on both branches. A change to either on one branch is copied to the other.
4. The worktree has no data, so the real-data GUI tests in `case_studies/bisb2026/tests/test_gui_joint_fit_realdata.m` skip there. Check real data with a temporary script outside the repository that reads the main checkout's data by absolute path; its results stay out of git.
5. A new GUI control must stay optional when older history and session files are restored: a missing field takes the control's default.
