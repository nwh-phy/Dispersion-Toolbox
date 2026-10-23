# GUI integration — 2026-10-23

## Git integration

With user approval, local main was fast-forwarded from bd2d7da to 01c9d1b.
The GUI branch was merged into feat/zlp-joint-fit in 00ed806 without conflicts.
Research-specific commits were not merged into main.
Commit 026bba0 updates the daily startup instructions and branch workflow.

## Validation

MATLAB R2026a, main checkout after integration:

| Suite | Passed | Failed | Incomplete |
|---|---:|---:|---:|
| tests | 205 | 0 | 0 |
| case_studies/bisb2026/tests | 78 | 0 | 0 |

The real-data GUI/API checks cover three datasets and six q channels.
The saved acceptance session loaded successfully (MERGED_SESSION_OK).
Background-model tests emitted invalid-start warnings and completed successfully.
This verifies software behavior, not the physical interpretation of every fitted peak.

## Cleanup and recovery

The completed GUI worktree was clean, backed up, and removed. Its branch is retained.
The main checkout remains on feat/zlp-joint-fit; the exploration worktree is unchanged.
Local backup branches preserve main at bd2d7da and research at 6694e56.
Full logs, MATLAB results and the worktree archive remain outside Git under
../Dispersion-Toolbox-GUI-acceptance/20261023/merge/.
No raw datasets or full session artifacts are committed.
