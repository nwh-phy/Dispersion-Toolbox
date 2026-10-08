# 源文件与独立计算证据索引

仅指本轮上传的v4 ZIP。源码、数值CSV均未用网络版本替换。行号从文件第一行计；长行不按分号拆分。evidence/内为对应小文件原样副本，SHA256如下。大MAT在用户原交付中，不重复打入本计划。

用户ZIP SHA256：`77dfe0c21292212469a0063406cbe279972f4e2817226e5466d5659a9fe9388c`

## 源证据

### S01 — `run_report.md`
副本：`evidence/S01/run_report.md`；32行；SHA256 `fda9d1198c35a6945f7df862d72d5dc31d4f144e9df745da6708be6a0c9ab942`。

### S02 — `provenance/source_snapshot/case_studies/bisb2026/scripts/c20_v4_export_fits.m`
副本：`evidence/S02/c20_v4_export_fits.m`；100行；SHA256 `cf45935ffadb82698a8df0d263806c613816d37a327248d771812ac402e89841`。

### S03 — `590_PL2_10w/component_parameters.csv`
副本：`evidence/S03/component_parameters.csv`；19行；SHA256 `a2760a5df6a2f8c20bee4409b46a441dfb9f33174ea96b0f77c568a124e11085`。

### S04 — `provenance/source_snapshot/src/fitting/qe_compare_component_models.m`
副本：`evidence/S04/qe_compare_component_models.m`；122行；SHA256 `b444605453a5f6bf6b2f932456fc0c858bf2ae5e1a0839581b89752826d6d433`。

### S05 — `provenance/source_snapshot/src/qe_zlp_integer_align.m`
副本：`evidence/S05/qe_zlp_integer_align.m`；21行；SHA256 `4f202697691c20307b57e61c3296d29c7371305d44fbaf0931710985abd2ea84`。

### S06 — `provenance/source_snapshot/case_studies/bisb2026/scripts/c20_v4_frames.m`
副本：`evidence/S06/c20_v4_frames.m`；84行；SHA256 `f2a0090d97185f3b7e0d235e8a4fb493062ca20749daf380257ed49350c30b02`。

### S07 — `validation/simulation_summary.csv`
副本：`evidence/S07/simulation_summary.csv`；9行；SHA256 `2541d3718809d47dc3a5ae764ea562fa98d7118eca23ec0a069850a82c32c24d`。

### S08 — `provenance/source_snapshot/c20_v4_simulate.m`
副本：`evidence/S08/c20_v4_simulate.m`；25行；SHA256 `cd6e97097f99ceace97e3f817137866b68fa62a12a1dd32674da6444b132d18b`。

### S09 — `provenance/source_snapshot/c20_v4_p5_355a5f09.m`
副本：`evidence/S09/c20_v4_p5_355a5f09.m`；67行；SHA256 `355a5f094231c1e801b682e3fefcf8d5b4725073007f82454288b046dda413cd`。

### S10 — `provenance/source_snapshot/c20_v4_profiles.m`
副本：`evidence/S10/c20_v4_profiles.m`；49行；SHA256 `6786854a0cd5a43c6e252a2910b4f7b051c021845aea0e2b2cc24262a0fb8a5b`。

### S11 — `590_PL2_10w/solver_background_differences.csv`
副本：`evidence/S11/solver_background_differences.csv`；37行；SHA256 `5b55f82b7e611a7b4e69f4d0e5d9265969e4aaee4b541b954d982327ede37b57`。

### S12 — `provenance/source_snapshot/case_studies/bisb2026/scripts/run_b1_component_v4.m`
副本：`evidence/S12/run_b1_component_v4.m`；99行；SHA256 `381939aa31bf790ff20f8e86053e3e7d366ba09b3ed854d0d0e0a2db84832149`。

### S13 — `stage_status.json`
副本：`evidence/S13/stage_status.json`；23行；SHA256 `2408562d42b88e5d77d93aaee4a72514f394d27d96ef131a6e77b3d9613ec664`。

### S14 — `590_PL2_10w/frame_qc.csv`
副本：`evidence/S14/frame_qc.csv`；301行；SHA256 `0ff8673f8312b9378bea82a0b30a6c0ed5987d458197bc16d0abd9227bedd136`。

### S15 — `590_PL2_10w/centered_bins.csv`
副本：`evidence/S15/centered_bins.csv`；10行；SHA256 `3ba84a05e43c5cee7ff65473ee13ce2d30d065d25877839ace81539e2556e393`。

### S16 — `tests/results.csv`
副本：`evidence/S16/results.csv`；18行；SHA256 `e6fd25a37b43ee4e3839994012ae3d8ec59d9914fb2ead875f3f55d785c23954`。

### S17 — `tests/simulation_tests.csv`
副本：`evidence/S17/simulation_tests.csv`；3行；SHA256 `e38db25a84fe53d22a687c6d7931c1c6584f0100673f133c70d5ebac985a7b81`。

### S18 — `590_PL2_10w/block_spectral_diagnostics.csv`
副本：`evidence/S18/block_spectral_diagnostics.csv`；19行；SHA256 `ba35c840037b3c4a101652dab2e0cb1173da49969369e531f42e426e6a246084`。

### S19 — `590_PL2_10w/A0_A1_diagnostics.json`
副本：`evidence/S19/A0_A1_diagnostics.json`；18行；SHA256 `dd998589d34863f3caf3243387c24f702b10aa1442c1e48138251823d1165481`。

### S20 — `provenance/source_snapshot/src/fitting/qe_component_mapping.m`
副本：`evidence/S20/qe_component_mapping.m`；30行；SHA256 `038d4cd77a9e91c084f4ed834a8d04e8d44460477261ea20297f7205eee384be`。

## 独立计算

A1：`independent_checks/parent_mapping_checks.csv`、`selected_parent_boundary_hits.csv`。从每个候选raw值独立重建排序逆映射、原判据上下界和原生缩放。
A2：`real_fit_checks.csv`、`real_candidate_checks.csv`、`witness_checks.csv`。36个真实谱模型、1008优化候选、18零幅度见证；重新计算公式、背景、参数排序、预测、残差、目标。
A3：`bin_checks.csv`、`block_checks.csv`、`frame_derived_diagnostics.csv`、`q_peak_channel_by_block.csv`。从15成员300帧数组独立累计、应用保存偏移、重建9bin和18块。
A4：`parent_area_scope_audit.csv`。从父模型原生参数分别积分实际fit窗和统一300–1800窗，显示字段误标的108行；本次没有重拟合。
A5：`simulation_checks.csv`、`hypothetical_double_recovery_errors.csv`。440个模拟数据集5280候选目标回算；规则、误报/恢复的工程场景统计。
A6：`profile_checks.csv`。27个点81次约束拟合的目标与约束回算；保留1次非成功试次，80次正退出约束通过。

所有上述A表相对路径均位于 `independent_checks/`。运行脚本：`check_delivery.py`；环境只需numpy/pandas/scipy。完整运行输出 `check_execution.txt`，总结 `independent_checks/CHECK_SUMMARY.json`。

## 可复核边界

本次独立计算没有调用MATLAB优化器。没有重新拟合实验谱；没有包外NPY/JSON，所以不能独立再次认证完整采集身份、各帧同一位置、独立响应或样品稳定性。父候选映射全部可查，父全部原始观测不完整，因此未称父2592目标全部独立重算。包内19项MATLAB测试通过是交付实录；本次单独执行了Python回读，二者不混称。
