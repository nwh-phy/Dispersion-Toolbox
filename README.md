# q-EELS Plasmon Dispersion Toolbox

用于动量分辨低损 EELS（q-EELS）数据导入、交互检查、峰拟合和色散提取的 MATLAB 工具箱。

## 启动

在项目主目录启动 MATLAB，运行：

```matlab
cd('path/to/Dispersion-Toolbox')
startup
interactive_qe_browser
```

打开数据或恢复完整会话：

```matlab
interactive_qe_browser('path/to/eq3D.mat')
interactive_qe_browser('path/to/qe_session.mat')
```

代码只有 `main` 一条线，工具箱和 BiSb 分析都在上面；测试过的稳定点用 `stable-YYYYMMDD` 标签标记。

要求 MATLAB R2024b+，以及 Optimization、Curve Fitting、Signal Processing、Image Processing Toolbox。
近期验证环境为 macOS / MATLAB R2026a。

## 日常流程

1. **Load Data**，核对 dq、q 范围和能量窗。
2. 选代表性 q 通道，**Fit Spectrum** 检查后再接受拟合。
3. **Auto Fit Branches** 提取候选分支，逐点复查并修正。
4. **Save Session** 保存完整工作状态，**Load Session** 继续。
5. **Export Panels → Publication** 导出 PNG 和 PDF。

`ZLP joint` 使用原始计数，接受到色散中的能量为模型 E0；`BG + peaks` 保留峰顶定义。
两种能量不要混用。详细设置、运动学因子、手动修正和脚本接口见
[GUI 使用与接口参考](docs/toolbox/gui_workflow.md)。

## 去哪里找

| 需求 | 入口 |
|---|---|
| GUI 操作、模型选项、数据格式 | [GUI 使用参考](docs/toolbox/gui_workflow.md) |
| BiSb 分析脚本与调用方法 | [脚本索引](case_studies/bisb2026/scripts/README.md) |
| 本机分析结果、依赖和归档 | `paper_results/00_INDEX.md`（本地文件，不随 Git 分发） |
| 软件验收和整理记录 | [报告索引](reports/README.md) |
| 开发约定、论文依据、分支规则 | [AGENTS.md](AGENTS.md) |

通用代码在 `src/`，第三方库在 `lib/`，BiSb 专用代码在 `case_studies/bisb2026/`。
实验数据 `20260120 BiSb/`、分析结果 `paper_results/` 和多数文献材料仅保存在本地。
`archive/` 是本地的历史记录（不随 Git 分发），日常操作从上表进入。

## 验证

```matlab
startup
runtests('tests')
runtests('case_studies/bisb2026/tests')
```

2026-11-10 整理后（提交 `576b650`，标签 `stable-20261110`），macOS / MATLAB R2026a，主目录、本地数据：

- BiSb **128/128** 通过。
- 核心里不开 GUI 窗口的 23 个文件 **103/103** 通过。
- 打开 GUI 的 8 个文件 **60/60** 通过；这一项是在提交 `a8220ec`、无数据的工作树里跑的，GUI 源码和这 8 个测试文件与 `576b650` 逐字相同。

无失败或跳过。两套测试上一次一口气全量跑是 2026-11-09（那时核心 214、BiSb 85，之后有 34 个 C20 测试从核心移到 BiSb）。
更早的 GUI 合并基线见[合并验收记录](reports/gui_integration_20261023.md)。

## 使用范围

本项目供组内学术研究使用。研究进展和结论按 AGENTS.md 的论文与后续证据链核对；
软件测试通过只验证相应的软件行为。
