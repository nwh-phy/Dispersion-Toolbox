# 论文期脚本（2026-03 至 05，冻结）

毕业论文的数值和图出自这里，不再修改。7 个 Lorentz 追踪脚本和 `run_b1_double_peak_binning_analysis.m` 靠 `addpath(script_dir)` 调用同目录的函数，所以这批整体放在一个文件夹里。结果目录见本机的 `paper_results/00_INDEX.md`"论文期证据"一节。

## 当时的口径

- 峰位的主要证据是不扣背景、面积归一化的 Fano 提取；扣背景的分析只作诊断和稳健性检验。
- `10w`、`20w` 是离焦量，不是束流功率。
- `590 PL2 10w` 和 `n0 PL2 10w repeat` 当作一层的两次重复，`no PL2 20w 2film` 是两层。

## 脚本

| 用途 | 脚本 |
|---|---|
| GUI 操作历史重现（主分析） | `run_590_gui_history_area_analysis.m`：590、n0、20w 三组，面积归一化；读取 `op_history_260506.mat` 和 `paper_results/legacy_branch_points_260521/` |
| B1 Lorentz 峰追踪 | `run_b1_lorentz_tracking_optimization.m`、`_until300`、`_upper_stability`、`_upper_rapidrise_plateau`、`_ridge_guided_optimization`、`_ridge_guided_failure_retry`、`_balanced_highqbin7_sg71`；`run_b1_lorentz_peak_evidence_audit.m` |
| B1 双峰提取 | `run_b1_double_peak_binning_analysis.m`、`run_b1_double_peak_waterfall_extraction_overlay.m` |
| 物理拟合 | `run_b1_physical_fit_analysis.m`、`run_b1_physical_fit_enhancements.m`、`run_epsilon_bg_sensitivity.m` |
| 20w | `run_20w_b1_highq_audit.m`、`run_20w_candidate_overlay_export.m`、`run_20w_lowq_gap_diagnostic.m`、`run_20w_lowq_sqrt_zero_diagnostic.m` |
| 低 q | `run_lowq_gap_comparison.m`、`run_lowq_powerlaw_deviation_export.m` |
| 线宽 | `run_current_q_linewidth_export.m`、`run_q015_linewidth_by_session_export.m`、`run_wideq_linewidth_analysis.m`、`run_lorentz_fano_comparison.m` |
| 出图 | `run_b1_current_waterfall_three_dataset_export.m`、`run_b1_qe_heatmap_three_dataset_export.m`、`run_b1_waterfall_three_dataset_comparison.m`、`run_b3_scatter_only_export.m` |
| 峰位稳健性 | `run_peak_position_robustness_matrix.m` |
| 背景扣除比较 | `compare_bg_dual_window_590.m`、`compare_bg_dual_window_all.m` |
| Do et al. 式批处理 | `run_bosman_pipeline.m`（输出到 `output/bosman_pipeline/`） |

四个 Lorentz 追踪和证据审计脚本运行时会往 `case_studies/bisb2026/RESULTS_INDEX.md` 追加一段；该文件 2026-11-10 已删，脚本在它不存在时会新建。

## 只给这批脚本用的函数

- 2026-10-21 从 `src/` 移来：`b1_double_peak_binning_extract.m` 及其子函数 `b1_double_peak_candidate_path_select.m`、`b1_double_peak_ridge_reference.m`，`b1_double_peak_repair_tracking_points.m`、`b1_double_peak_upper_stability_metrics.m`、`b1_peak_evidence_audit_classify.m`、`qe_recommend_lowq_background_configs.m`、`qe_summarize_lowq_background_eval.m`，以及早期的 10w/20w 对比脚本 `analyze_dispersion.m`。
- 2026-11-10 从 `src/` 移来：`qe_select_low_energy_branch.m`、`qe_preprocess_preset.m`、`thesis_config.m`（原 `src/thesis/`，项目根目录改由 `bisb_find_project_root` 找）。`src/thesis/` 里其余的论文流水线同日删除：它没有保留下来的输出，论文数值来自上表的 GUI 操作历史重现。

测试在 `case_studies/bisb2026/tests/`。
