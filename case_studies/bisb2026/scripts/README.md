# BiSb 分析脚本索引

在项目根目录的 MATLAB 里先运行：

```matlab
startup
addpath(genpath(fullfile(pwd, 'case_studies', 'bisb2026')))
```

然后按函数名调用。脚本用 `bisb_find_project_root` 找根目录，结果写到 `paper_results/<名字>/…`（各结果文件夹的说明见本机的 `paper_results/00_INDEX.md`）。数据和 `paper_results/` 不随 Git 分发，只在主目录里有。

| 位置 | 内容 | 状态 |
|---|---|---|
| 本目录 | 数据集登记表、新数据的首批图、三批脚本共用的小函数 | 在用 |
| [`b1_lowq/`](#2-在用的分析b1_lowq) | 原生通道的 B1 低 q 分析：P1、P1b、P2a、P2b、简化流程、低能支、单项检验 | 在用 |
| [`b1_binned/`](#3-已做完仍被调用b1_binned) | 分箱的全 q 和 adaptive-q 双分量拟合（2026-10 上旬） | 已做完；`b1_lowq/` 读它的输出、调它的函数，不再往里加新分析 |
| [`legacy_c20/`](legacy_c20/README.md) | C20 v2–v6（2026-09）和 v7 镜像检验 | 冻结 |
| [`legacy_202605/`](legacy_202605/README.md) | 论文期脚本（2026-03 至 05） | 冻结 |

冻结的两批不再修改，文件清单和注意事项在各自的 README 里。报告在 `case_studies/bisb2026/reports/`，测试在 `case_studies/bisb2026/tests/`。

## 1. 新数据来了

原始数据放进项目根目录的数据文件夹（加进 `.gitignore`），在 `bisb_datasets.m` 里加一条，再运行 `bisb_quicklook('<id>')`。

| 本目录的函数 | 作用 |
|---|---|
| `bisb_datasets.m` | 唯一的数据集登记表：每次采集的文件夹、`.npy`、dq、q = 0 列、狭缝接收 h、SHA-256、样品和层数。分析代码按 id 取，脚本里不写路径、dq 或 q = 0。现有 id：`590_PL2_10w`（`590`）、`n0_PL2_10w_repeat`（`n0`）、`no_PL2_20w_2film`（`20w`）、`Bi1107_10w`、`Bi1107_20w` |
| `bisb_frame_map.m` | 某个数据集逐帧对齐、按束斑漂移分箱后的 q–E 计数图（另有按时间分块的求和），缓存在 `paper_results/frame_maps/` |
| `bisb_quicklook.m` | 新数据的标准首批图：两张列归一化热图（\|q\| ≤ 0.015、≤ 0.005）、堆积图、q = 0 单谱（含扣 ZLP 尾巴后的峰位）和 `summary.txt`，输出到 `paper_results/quicklook/` |
| `bisb_find_project_root.m` | 找项目根目录（含 `startup.m` 和 `src/` 的最近一层） |
| `b1_simple_root.m` | 返回根目录，并把 `src/`、本目录、`b1_lowq/`、`b1_binned/` 加进路径 |
| `b1_sha256.m` | 配置和输入的 SHA-256，用来给输出文件夹命名 |

需要别的图时，先用 `src/` 里这些通用函数，不要另写对齐、取谱或画热图的代码：

| `src/` 里的函数 | 作用 |
|---|---|
| `qe_align_sequence_frames` | 读原始序列，逐帧整数对齐 ZLP，测束斑 q 漂移 |
| `qe_transport_frames` | 把对齐后的帧按 q 分箱，跟随每帧的漂移 |
| `qe_pm_spectra` | 某些 \|q\| 处 +q、−q 两侧的平均谱 |
| `qe_zlp_center` | ZLP 中心和半高宽（P1b/P2a 的三点对数抛物线规则） |
| `qe_powerlaw_tail` | 在一个能量窗里拟合 ZLP 幂律尾巴；窗口要避开约 130 meV 的积碳声子 |
| `qe_peak_apex` | 不靠模型的峰顶位置（平滑后找最显著的局部极大，再在原始数据上用抛物线细化） |
| `qe_plot_heatmap` | 列面积归一化的 q–E 热图，按分位数定色标，导出 PNG 和 PDF（线性或 log） |

## 2. 在用的分析（`b1_lowq/`）

拟合核心在 `src/fitting/qe_zlp_joint_fit.m`、`qe_kinematic_prefactor.m`、`qe_lowq_joint_forward.m`。P1 到 simple 的输出都在 `paper_results/b1_adaptq/20261009T054142Z/` 下的同名子目录里。一轮新的分析加一个配置函数、调用已有代码，不复制上一轮的脚本。

| 分析线 | 入口 | 配置 | 输出 / 报告 |
|---|---|---|---|
| P1：单个原生通道的拟合 | `b1_p1_lowq`、`b1_p1_report` | `b1_p1_config` | `p1/` |
| P1b：不同 kernel 的低 q 拟合，帧数据的噪声模型 σ² = g(E)·counts | `b1_p1b_lowq`、`b1_p1b_diag`、`b1_p1b_noise`、`b1_p1b_report` | `b1_p1b_config` | `p1b/` |
| P2a：590 的输入检查、slice 响应、块算子、注入回收 | `b1_p2a_inputs`、`b1_p2a_response`、`b1_p2a_operator`、`b1_p2a_synth`；共用 `b1_p2a_build_ops`、`b1_p2a_escale` | 各自的 `_config` | `p2a/` |
| P2b：20 个原生通道共用一个 L(Q,E) 的前向拟合，两支色散 | `b1_p2b_fit(run_root, cfg, ocfg)`、`b1_p2b_dispersion` | `b1_p2b_config` | `p2b/fit_590_<配置哈希>/` |
| P2b 出图 | `b1_p2b_lossmap_fig`（调用 `b1_p2b_intrinsic_L`） | — | 实测计数、正向模型、本征 L(Q,E)、谱形残差四联 q–E 图，一条单谱扣除运动学因子前后的对比 |
| P2b 的狭缝接收 h | `b1_p2b_h66_run`、`b1_p2b_h66_conv`、`b1_p2b_hscan_run`、`b1_p2b_hscan_summary` | `b1_p2b_slot_config('h66')` 或 `('hscan', h)`，在 `b1_p2b_config`、`b1_p2a_operator_config` 的基础上改出两份配置 | `p2b/`；`reports/p2b_hscan_20261017.md` |
| 简化流程：SI 式简单模型在多小的 q 还可靠 | `b1_simple_s0` | `b1_simple_config` | `simple/`；`reports/simple_pipeline_S0_20261017.md` |
| 低能支：直接读原始 NPY，比较一个和两个 B1 响应 | `b1_low_branch_matlab(stage, …)`，`stage` 为 `prepare`、`pilot`、`fit`、`bootstrap`、`report`、`verify` | `b1_low_branch_matlab_config` | `paper_results/b1_low_branch_matlab_20261024/`；`reports/b1_low_branch_matlab_20261024.md` |
| 一层 / 两层对比 | `b1_layer_fig`（逐点图）、`b1_layer_compare`（分 q 段的平均表） | — | `stage5/figures/`；`reports/layer_1L_2L_20261017.md` |
| 单项检验 | `b1_kinematic_ridge`（列归一化 q–E 图的峰脊）、`b1_whole_band_check`（不分两支的整个 B1 带）、`b1_single_divide_fig`（单个通道直接除以运动学因子，不拟合） | — | — |

## 3. 已做完、仍被调用（`b1_binned/`）

两条线的运行记录、启动脚本和报告在 `paper_results/b1_fullq/20261007T171525Z/` 和 `paper_results/b1_adaptq/20261009T054142Z/`（各自的 `REPORT.md`）。

| 线 | 阶段 | 脚本 |
|---|---|---|
| 全 q 双分量（N = 3 分箱） | 配置 | `b1_fullq_config` |
| | 0：在 v7 的 6 个点上确定运动学因子形式、垂直接收范围 h、低 q 退化 | `b1_fullq_stage0`、`b1_fullq_hscan`、`b1_fullq_stage0c` |
| | 1：原始帧 → 逐帧 ZLP 对齐的分箱和噪声模型；核对帧求和谱 | `b1_fullq_build`、`b1_fullq_check_build` |
| | 2、4：每个分箱的全部拟合配置和帧重采样，78 个分箱分给 3 个进程 | `b1_fullq_fit_worker`、`b1_fullq_run_part` |
| | 3：双支合成谱验证 | `b1_fullq_synth_worker` |
| | 5、6：汇总、色散拟合、图和数值表；双峰证明图 | `b1_fullq_aggregate`、`b1_fullq_fig_twopeak` |
| adaptive-q（计数允许处用原生通道的自适应分箱） | 阶段 1–6 同全 q | `b1_adaptq_config`、`_build`、`_look`、`_fit_worker`、`_aggregate`、`_dispersion` |
| | 图 | `b1_adaptq_e0_figs`、`_fig_twopeak`、`_fig_saturation`、`_fig_kinematic` |
| | 诊断 | `b1_adaptq_lowq_diag`、`b1_adaptq_raw_lowq` |
| 试跑和核对 | ZLP 联合拟合对比旧的幂律窗口拟合；50–150 meV 辅助峰是否含 MoS₂ 声子；低 q 处 ZLP 分量数 | `run_zlp_joint_fit_pilot`、`run_zlp_aux_phonon_test`、`run_zlp_ncomp_test` |

## 写新代码时

- 新脚本叫 `b1_<分析线>_<动作>.m`，放进 `b1_lowq/`；三批都用的放本目录，叫 `bisb_*.m`。版本号不进文件名，写进配置（`cfg.version`）。
- 新测试叫 `test_<被测函数>.m`，放 `case_studies/bisb2026/tests/`。项目根目录的 `tests/` 只测 `src/`。
- 旧文件不改名：旧报告和 `paper_results/` 里的溯源记录按名字引用它们。
- 一条分析线做完后整体移进 `legacy_*`，并在那里写一份 README；前提是在用的脚本不再调用它。
