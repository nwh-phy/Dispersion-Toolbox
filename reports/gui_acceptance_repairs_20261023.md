# 原生 GUI 验收补修 — 2026-10-23

本轮针对 `565cbd3` 原生验收发现的两项问题，工作分支为 `gui-joint-fit`。

## 修复

1. 初次载入先把数据的默认参考 q 区间同步到控件，再计算 Comparison。原来计算使用
   0.010–0.020 Å⁻¹，控件随后显示 0.001–0.002 Å⁻¹；保存后重新计算便出现图像变化。
   MAT 和 raw 导入入口都调整了顺序。修改 dq 后也会根据同步后的参考区间重算对照图。
2. 导出 figure、axes 和图例明确使用白底，标题、刻度、坐标标签和图例文字使用深色。
   不依赖桌面主题或 root graphics defaults；只修改导出副本，保留源图、曲线数据、颜色、
   colormap、隐藏峰标记及固定比例画布。`Theme` 属性使用存在性检查兼容旧版接口。

## 新增回归测试

`tests/test_gui_acceptance_regressions.m` 包含三个行为用例：

- 初次导入后，不按 Build Views，参考区间、Comparison 结构和实际图像 CData 在保存重开前后一致。
- 修改 dq 后重复上述检查。
- 强制深色 root 默认样式，检查导出 PNG 的白底、深色标题和坐标标签、4:3 比例及 PDF 文件；
  核对源图的颜色、图例和曲线数据未改变。

三个用例在修复前全部失败，原因与原生验收一致；修复后全部通过。

本地完整证据保存在仓库外 `../Dispersion-Toolbox-GUI-acceptance/20261023/repair/`
（路径相对于主工作区）。其中 `new_tests_before_fix.log` 和 `new_tests_after_fix.log`
分别保存失败复现及通过结果。真实数据、完整会话和导出文件不进入 git。

## 完整回归

MATLAB R2026a，修复前后运行 `tests/` 和 `case_studies/bisb2026/tests/`：

| 测试集 | 修复前 | 修复后 | 结果变化 |
|---|---|---|---|
| 核心 | 189 通过，0 失败，0 跳过 | 192 通过，0 失败，0 跳过 | 新增 3 个用例通过 |
| BiSb | 67 通过，3 失败，8 跳过 | 67 通过，3 失败，8 跳过 | 78 项逐项不变 |

三个既有失败为源码文本预期不匹配：

- `Run590GuiHistoryAreaOptionsTest/defaultsQRangeToQ015WhenUnspecified`
- `test_epsilon_bg_sensitivity_workflow/testScriptUsesCurrentThreeAreaReports`
- `test_gui_history_highq_refinement_profile/testTwentyWProfileTargetsOnlyB1LargeQ`

八个跳过用例需要工作区内的真实数据，而 GUI worktree 不含数据。本轮没有修改这些研究脚本
或降低断言。前后日志、MAT/TSV 结果和逐项比较
分别保存在 `before_tests.log`、`after_tests.log`、`{before,after}_results.mat`、
`{before,after}_{core,bisb}.tsv`、`suite_comparison.json`。

## 真实数据与原生复验

环境：macOS、MATLAB R2026a、深色界面。输入为主工作区的
`20260120 BiSb/590 PL2 10w 0.004 10sx300/eq3D.mat`。

以下点击、输入、文件选择、保存、载入和导出由 Computer Use 操作原生 GUI 完成；
MATLAB 启动、文件数值核对及 PDF 渲染由终端辅助。

| 操作 | 结果 | 证据 |
|---|---|---|
| Load Data，设置 q = ±0.01 Å⁻¹、E = 300–1800 meV，不按 Build Views | 初次 Comparison 与当前 0.001–0.002 Å⁻¹ 参考区间一致 | `01_fresh_saved.jpg` |
| Save Session，完整退出 MATLAB，新启动空白窗口，再 Load Session | 对照图中心带与两侧结构保持一致 | `02_restarted_empty.jpg`、`03_fresh_restored.jpg` |
| 恢复上一轮完整会话 | 63 个分支点、2 条修正、2 个色散拟合及待接受单谱拟合保留 | `04_full_session_restored.jpg` |
| Publication、4:3、All panels，经原生目录选择器导出 | 4 PNG + 4 PDF，GUI 确认 4 个面板导出成功 | `05_export_completed.jpg` |
| 全部 PNG 和全部 PDF 渲染页逐张检查 | 白底、深色标题/坐标/图例，曲线和峰标记可见，标题未越界 | `20261023_203018_*_4x3_publication.*`、`pdf_render/` |

四张 PNG 均为 2500×1875；四个 PDF 均为单页 600×450 pt。文件清单及尺寸见
`export_files.json`。Poppler 将所有 PDF 渲染成长边 1600 px 的图像后再做视觉检查。

仓库外的 `verify_repair_realdata.m` 另外完成 **20/20 项数值检查**：

- 使用测试钩子分别检查真实 590 数据初次导入、修改 dq 两种情况的参考区间、完整
  Comparison 结构、图像 CData 和 Physical 数据，8 项通过。这部分是自动化补充检查。
- 比较两个原生快照 `fresh_before_restart.mat` / `fresh_after_restart.mat`，以及旧完整
  快照 `../590_before_restart.mat` / 新快照 `590_after_fixed_export.mat`；每组比较
  `data`、`data_hashes`、`snapshot`、`analysis`、`pending_fit`、`display_fit`，12 项通过。
  操作历史和保存时间允许变化。结果为 `repair_numeric_checks.json` 及同名日志。

这 20 项是数值断言，不计入正式测试集的 270 个用例。

## 边界

两项指定缺陷在上述环境与 590 数据上通过复验。核心拟合函数、研究脚本和输入数据未修改；
主工作区保持 `feat/zlp-joint-fit`，未合并、推送或开 PR。

R2024b、Windows、其他真实数据格式及 raw 4D 中间态恢复未做原生复验。本轮沿用完整会话中
既有的拟合结果，没有重新做物理模型评估，也没有把 15 个待复查候选点标为科研验收完成。
上一轮记录的长文本截断、Correct Peak 按钮样式和 Fit Curated 后标记展示差异仍属于另行处理项。

原始 eq3D 输入补修前后的 SHA-256：
`4fcfc98369ab772e08d7eef39d056cfa2b58c1935d249b8fd683d97e555a0e9d`。
