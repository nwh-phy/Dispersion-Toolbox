# C20（2026-09，冻结）

Codex 在 Windows 上写的 C20 v2–v6，以及 2026-10 的 v7 镜像检验。不再修改，也不要重跑：v5、v6 按固定路径读取上一版的运行目录并校验 ZIP 哈希（v5 的 `next_action.md` 写明不重跑）；v2、v4、v5 的入口把自己和依赖文件的旧路径（`case_studies/bisb2026/scripts/<名字>.m`、`src/…`）写进了溯源清单。它们的运行目录仍是后续分析的输入，见本机的 `paper_results/00_INDEX.md`。

| 版本 | 入口 | 配套 |
|---|---|---|
| v2 | `run_b1_component_pilot_v2.m` | `b1_component_config_v2.m`、`audit_b1_component_pilot_output_v2.m`、`write_b1_component_pilot_report_v2.m` |
| v3 | `run_b1_p4p5_diagnostics_v3.m` | — |
| v4 | `run_b1_component_v4.m` | `c20_v4_*.m`（11 个：读写、帧、导出、作图、目标函数剖面、模拟、P5、补充诊断、包校验、打包、收尾） |
| v5 | `run_b1_data_constrained_pilot_v5.m` | `c20_v5_*.m`（6 个） |
| v6 | `run_b1_physics_preview_v6.m` | `c20_v6_physics_figures.m`、`c20_v6_finish.m`、`c20_v6_literature.ps1` |
| v7 | `run_b1_mirror_q_v7.m` | ±q 镜像检验，`b1_binned/` 里全 q 阶段 0 的参照；用到 v4 的读写和作图函数 |
| 共用 | `c20_packet_readback.m` | 在临时目录解压 ZIP、核对清单和哈希，验证后清理 |

只给这批脚本用的函数（2026-11-10 从 `src/` 移来）：`qe_area_scopes.m`、`qe_component_mapping.m`、`qe_component_prediction.m`、`qe_fit_member_models.m`、`qe_member_prediction.m`、`qe_spectral_centroids.m`，以及三组数据的清单 `thesis_sessions.m`（路径和 dq 与 `bisb_datasets.m` 一致，由 `test_bisb_datasets.m` 核对）。

测试在 `case_studies/bisb2026/tests/`：`C20*.m`、`test_b1_component_pilot_v2.m`。

## v2 试跑说明（2026-09 原文，原 `C20_PILOT_V2.md`）

下面的路径和命令是当时的写法：脚本那时在 `scripts/` 顶层，测试在 `tests/`。

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
