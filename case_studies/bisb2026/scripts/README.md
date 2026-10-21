# BiSb case-study 脚本索引

在项目根目录的 MATLAB 里运行：`startup; addpath(genpath(fullfile(pwd,'case_studies','bisb2026')))`，再按函数名调用。脚本用 `bisb_find_project_root`（留在本目录顶层，三批共用）找根目录、自己加 `src/` 路径，结果写到 `paper_results/<名字>/<时间戳>/`（各结果文件夹的说明见本地的 `paper_results/00_INDEX.md`）。

目录（2026-10-21 整理）：

| 位置 | 内容 | 状态 |
|---|---|---|
| 本目录顶层 | 2026-10 起的 B1 分析：全 q、adaptive-q、低 q 的 P1/P1b/P2a/P2b、简化流程检验 | 在用 |
| `legacy_c20/` | C20 v2–v6（2026-09）和 v7 镜像检验 | 冻结 |
| `legacy_202605/` | 论文期脚本（2026-03 至 05），以及只给它们用的 9 个函数（2026-10-21 从 `src/` 移来，见下） | 冻结 |

冻结的两批不再修改。v2、v4、v5 的入口脚本把自己和依赖文件的旧路径（`case_studies/bisb2026/scripts/<名字>.m`、`src/b1_*.m`）写进了溯源清单，移动后不能原样重跑；要重跑先改清单里的路径。

## 当前：B1 低 q（2026-10-09 起）

| 批次 | 脚本 | 作用 |
|---|---|---|
| adaptive-q | `b1_adaptq_config.m`、`_build`、`_look`、`_fit_worker`、`_aggregate`、`_dispersion`；图 `_e0_figs`、`_fig_twopeak`、`_fig_saturation`、`_fig_kinematic`；诊断 `_lowq_diag`、`_raw_lowq` | 计数允许处用原生通道（N=1）的自适应分箱，阶段 1–6 同全 q |
| P1 | `b1_p1_config.m`、`b1_p1_lowq.m`、`b1_p1_report.m` | 低 q 复核第一步：单个原生通道的拟合 |
| P1b | `b1_p1b_config.m`、`_diag`、`_noise`、`_lowq`、`_report` | 不同 kernel 的低 q 拟合，帧数据的噪声模型 σ² = g(E)·counts |
| P2a | `b1_p2a_inputs(_config)`、`_response(_config)`、`_build_ops`、`_operator(_config)`、`_synth(_config)`、`_escale` | 590 的输入检查、slice 响应、块算子、注入回收 |
| P2b | `b1_p2b_config.m`、`b1_p2b_fit.m`、`b1_p2b_dispersion.m` | 20 个原生通道共用一个 L(Q,E) 的前向拟合，以及两支色散 |
| P2b 的 h | `b1_p2b_slot_config.m`、`b1_p2b_h66_run.m`、`b1_p2b_h66_conv.m`、`b1_p2b_hscan_run.m`、`b1_p2b_hscan_summary.py` | 狭缝接收 h = 0.0005–0.0066 的重跑与扫描。`b1_p2b_slot_config('h66')` 或 `('hscan', h)` 在 `b1_p2b_config`、`b1_p2a_operator_config` 的基础上改出两份配置，传给 `b1_p2b_fit(run_root, cfg, ocfg)` |
| 简化流程 | `b1_simple_config.m`、`b1_simple_root.m`、`b1_simple_s0.m` | SI 式简单模型在多小的 q 还可靠 |
| 单项检验 | `b1_kinematic_ridge.m`、`b1_layer_fig.m`（配 `b1_layer_compare.py`）、`b1_whole_band_check.m` | 列归一化 q-E 图的峰脊、1L/2L 对比图、不分两支的整个 B1 带 |
| 共用 | `b1_sha256.m` | 配置和输入的 SHA-256，用来给输出文件夹命名 |

报告在 `case_studies/bisb2026/reports/`。

## 当前：B1 全 q 双分量（2026-10）

拟合核心在 `src/fitting/qe_zlp_joint_fit.m` 和 `src/fitting/qe_kinematic_prefactor.m`。运行记录和启动脚本在 `paper_results/b1_fullq/20261007T171525Z/stage*/launch.sh`。

| 阶段 | 脚本 | 作用 |
|---|---|---|
| 配置 | `b1_fullq_config.m` | 共享设置（q 网格、能窗、模型、运动学因子、重采样次数） |
| 0 | `b1_fullq_stage0.m`、`b1_fullq_hscan.m`、`b1_fullq_stage0c.m` | 在 v7 的 6 个点上确定运动学因子形式、垂直接收范围 h、低 q 退化 |
| 1 | `b1_fullq_build.m`、`b1_fullq_check_build.m` | 原始帧 → 逐帧 ZLP 对齐的 N=3 分箱和噪声模型；核对帧求和谱等于 A1 全 q 图 |
| 2、4 | `b1_fullq_fit_worker.m`、`b1_fullq_run_part.m` | 每个分箱的全部拟合配置和帧重采样；把 78 个分箱分给 3 个进程 |
| 3 | `b1_fullq_synth_worker.m` | 双支合成谱验证 |
| 5、6 | `b1_fullq_aggregate.m` | 汇总、色散拟合、图和数值表 |
| 图 | `b1_fullq_fig_twopeak.m` | 双峰证明图：单谱分解和单峰/双峰残差 |

同一条线上的试跑和核对：

- `run_zlp_joint_fit_pilot.m`：ZLP 联合拟合和旧的幂律窗口拟合对比
- `run_zlp_aux_phonon_test.m`：50–150 meV 辅助峰是否含 MoS₂ 声子
- `run_zlp_ncomp_test.m`：低 q 处 ZLP 分量数
- `legacy_c20/run_b1_mirror_q_v7.m`：±q 镜像检验，b1_fullq 阶段 0 的参照（用到 C20 v4 的读写和作图函数，所以和 C20 放在一起）

## C20（2026-09，Codex 在 Windows 上写的；`legacy_c20/`，冻结）

v5、v6 按固定路径读取上一版的运行目录并校验 ZIP 哈希，不要重跑（v5 的 `next_action.md` 写明不重跑）。

| 版本 | 入口 | 配套 |
|---|---|---|
| v2 | `run_b1_component_pilot_v2.m` | `b1_component_config_v2.m`、`audit_b1_component_pilot_output_v2.m`、`write_b1_component_pilot_report_v2.m` |
| v3 | `run_b1_p4p5_diagnostics_v3.m` | — |
| v4 | `run_b1_component_v4.m` | `c20_v4_*.m`（11 个：读写、帧、导出、作图、目标函数剖面、模拟、P5、补充诊断、包校验、打包、收尾） |
| v5 | `run_b1_data_constrained_pilot_v5.m` | `c20_v5_*.m`（6 个） |
| v6 | `run_b1_physics_preview_v6.m` | `c20_v6_physics_figures.m`、`c20_v6_finish.m`、`c20_v6_literature.ps1` |
| 共用 | `c20_packet_readback.m` | 在临时目录解压 ZIP、核对清单和哈希，验证后清理 |

测试：`tests/C20*.m`、`tests/test_b1_component_pilot_v2.m`。

## 论文期（2026-03 至 05；`legacy_202605/`，冻结）

对应的测试在 `case_studies/bisb2026/tests/`。7 个 Lorentz 追踪脚本和 `run_b1_double_peak_binning_analysis.m` 靠 `addpath(script_dir)` 调用同目录的函数，所以这批整体放在一个文件夹里。

2026-10-21 从 `src/` 移来的函数（使用者只有本批脚本、C20 v2 和测试）：`b1_double_peak_binning_extract.m` 及其子函数 `b1_double_peak_candidate_path_select.m`、`b1_double_peak_ridge_reference.m`，`b1_double_peak_repair_tracking_points.m`、`b1_double_peak_upper_stability_metrics.m`、`b1_peak_evidence_audit_classify.m`、`qe_recommend_lowq_background_configs.m`、`qe_summarize_lowq_background_eval.m`，以及早期的 10w/20w 对比脚本 `analyze_dispersion.m`。

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
| 共用 | `bisb_find_project_root.m` |
