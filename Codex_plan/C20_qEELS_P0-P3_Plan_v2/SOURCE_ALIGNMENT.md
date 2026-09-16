# 最新源码对齐记录

核查日期：2026-09-11。固定提交：`18af79b66cd27365fd00bcf3de26658ddd7e7d53`。

本次检查了关键入口/配置、读取、分箱、拟合、证据分类及相关测试定义；不是全仓库逐文件审计，未运行 MATLAB/实验数据。旧 `bbf24965` 审阅不能直接当成当前事实。以下链接为本次读取的固定版本；行号指源码行，不是聊天工具包装行号。

## 引用登记

- **S01 — 版本。** [提交记录](https://github.com/nwh-phy/Dispersion-Toolbox/commit/18af79b66cd27365fd00bcf3de26658ddd7e7d53)。GitHub compare 返回相对 `bbf24965b6c1ea177541bf382e0b28fd7cdafdf7` ahead 24 commits。本次未把 compare 的零行统计解释为文件无变化。
- **S02 — 论文与成果对齐。** [AGENTS.md](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/AGENTS.md)。正式提交版优先、用户后续更正优先于旧文字，根目录论文与提交版需区分。
- **S03 — 当前工具箱主线。** [README.md](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/README.md)。不预扣背景、curated工作流、Save Pts与完整live snapshot保存、通用/case测试拆分。README不是所有case脚本默认值的统一配置。
- **S04 — 当前 dq。** [infer_qe_dq_Ainv.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/io/infer_qe_dq_Ainv.m)。10w=0.0005、20w=0.00025 Å⁻¹/pixel。该helper实际按路径含10w/20w匹配，故新adapter仍须验证session身份，不能用于任意同名材料。
- **S05 — 三组session。** [thesis_sessions.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/thesis/thesis_sessions.m)。三组数据名、相对路径、角色与dq。
- **S06 — dq测试。** [test_qe_dq_calibration.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/tests/test_qe_dq_calibration.m)。验证两套步长和显式override优先；这是软件断言，不是仪器独立校准实验。
- **S07 — 输入/缓存。** [load_qe_dataset.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/io/load_qe_dataset.m#L23-L205)。显式MAT按类型读取；raw缓存验证import_provenance和crop；folder仍有processed优先，legacy auto-crop缓存仍可能被接受。证明最新代码已有改动，不应重复声称旧版所有缓存问题仍原样存在。
- **S08 — GUI-history三组复现。** [run_590_gui_history_area_analysis.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/case_studies/bisb2026/scripts/run_590_gui_history_area_analysis.m#L1-L280)。默认q范围±0.015；强制Area；590 history模板；20w高q阈值0.010和窗口1000–1700；analysis_results.mat中output结构。仅检查了这里的关键加载和配置，不将全部后续绘图与物理拟合视为本轮已验证。
- **S09 — B1双峰分箱主调用链。** [wrapper](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/case_studies/bisb2026/scripts/run_b1_double_peak_binning_analysis.m#L1-L490)；[core](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/b1_double_peak_binning_extract.m)。重点读取core 1–300、500–910、932–1435行附近：四种模式；均值谱；forced two的success检查；窗口各自单峰；旧q默认；noise/adaptive分箱；成员与预处理元数据。未从函数标题“Mandatory”推断两个物理分量已经成立。
- **S10 — 最新通用fitter。** [fit_loss_function.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/fitting/fit_loss_function.m#L21-L490)。模型参数数目动态、可配弱峰阈值、bootstrap、apex及质量字段；仍存在单初值、无边界fallback、删峰前后输出契约风险。没有实际运行确认其影响。
- **S11 — 模型定义。** [peak_models.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/fitting/peak_models.m#L23-L120)。现有fano；lorentz仍是Drude–Lorentz，不能改名偷换。`measure_peak_fwhm.m`及测试已由树/compare确认存在，执行前仍应读其实现并测边界/多峰行为。
- **S12 — 点级证据分类。** [b1_peak_evidence_audit_classify.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/b1_peak_evidence_audit_classify.m)。local support、component SSE、constraint release、robustness、宽度等经验规则；`robust_stable` 对非finite robust_delta可放行，故新结果须另存缺失测试标记。
- **S13 — 已有结果状态与审计入口。** [RESULTS_INDEX.md](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/case_studies/bisb2026/RESULTS_INDEX.md)；[run_b1_lorentz_peak_evidence_audit.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/case_studies/bisb2026/scripts/run_b1_lorentz_peak_evidence_audit.m#L1-L185)。500meV及相关tracking参考作废；正负q平均图作废；历史v14/v15/审计统计。审计入口固定目录并会更新索引，不适合未经修改就写新run。
- **S14 — 已有分箱测试。** [test_b1_double_peak_binning_workflow.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/case_studies/bisb2026/tests/test_b1_double_peak_binning_workflow.m#L1-L220)。真实数值运算合成测试包含noise bin、低q保护、metadata、highqforce；部分合成q仍为0.03–0.15等。需要补充当前物理dq及±0.015范围测试，不删除旧测试。
- **S15 — 另一套thesis配置。** [thesis_config.m](https://github.com/nwh-phy/Dispersion-Toolbox/blob/18af79b66cd27365fd00bcf3de26658ddd7e7d53/src/thesis/thesis_config.m)。它的pre_subtracted=true、low窗口500–2100等与新的B1双分量主线并非同一用途；不能只因文件名含thesis就替代提交版论文的实际证据链。

## 执行中必须核验的区别

### 1. q标定已有明确更新，但数值使用点未完全统一

读取端与登记的0.0005/0.00025作为当前项目口径；而B1默认q_skip0.005与no-bin0.05仍存在。若直接把history生成的±0.015谱交给旧double wrapper，可能既跳过大量低q谱，又没有任何可分箱点。这是从条件表达式推得的风险，不是本轮运行结果。

### 2. “已经有双峰代码”不等于“联合双分量已经可靠提取”

independent/propagated模式会调用同窗两峰；windowed模式分别拟合上下窗口单峰。后者在谱重/宽度解释方面不能直接与前者等价。首次新主线使用无趋势约束的联合同窗模型，旧candidate/path作为辅助。

### 3. 数值成功与科学可辨识不同

旧double success主要检查返回两峰、能量finite、间距未塌缩；其“成功”不包含完整双分量对单分量证据、仪器展宽、系统误差或恢复验证。继续保留数值status，另添科学验证字段。

### 4. 不再重复旧审阅中已经过时的结论

最新版已有Fano、论文session、保存更多拟合/curated状态、显式MAT与crop缓存规则改进、B1分箱和evidence审计；不要再将这些列为全新模块。对未核查的新实现只写待测，不把旧版问题未经验证地继承成当前bug。

### 5. 文档记载与真实数组始终分开

没有随仓库读取到的原始实验数据、paper_results MAT/CSV和正式论文，不会仅凭路径或索引被标成已验证。plan允许缺数据先完成合成测试，不允许制造实验结果。
