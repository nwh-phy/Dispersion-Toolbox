# 本轮审阅证据索引与行号

v3前缀指 `paper_results/b1_components_v3/20260912T035735540Z/`。

这些是用户粘贴片段的转录或上一轮附件的注明来源提取，不是从用户电脑直接复制的原始运行快照。行号以每段文件内容首行=1计；保留原长行、不按分号拆行。转录SHA只标识本包转录文本，不能证明父run实际执行身份。CSV表格已复算；原始谱、实际MAT、参数级CSV和PNG没有在本轮挂载输入中提供。本包不包含或声称重建这些缺失工件。

| 简记 | 文件 | 行数 | 来源说明 |
|---|---|---:|---|
| R | [case_studies/bisb2026/scripts/run_b1_p4p5_diagnostics_v3.m](evidence/run_b1_p4p5_diagnostics_v3.m.numbered.txt) | 35 | 本轮用户粘贴代码；Markdown转义还原，长行不拆分 |
| F | [src/fitting/qe_compare_component_models.m](evidence/qe_compare_component_models.m.numbered.txt) | 85 | 本轮粘贴与上一轮附件相同的拟合器文本；从上一轮挂载附件提取后对照本轮 |
| Q | [src/qe_prepare_count_bins.m](evidence/qe_prepare_count_bins.m.numbered.txt) | 64 | 本轮粘贴与上一轮附件相同的分箱器文本；不要与R中内联同中心代码混淆 |
| S | [v3/audit/boundary_summary.csv](evidence/boundary_summary.csv.numbered.txt) | 109 | 本轮贴出的父运行模型级汇总，不是boundary_by_parameter.csv |
| B | [v3/590_PL2_10w/centered_bins.csv](evidence/centered_bins.csv.numbered.txt) | 10 | 本轮九个同中心bin的成员表 |
| G | [v3/590_PL2_10w/frame_block_summary.csv](evidence/frame_block_summary.csv.numbered.txt) | 7 | 本轮六个分块的标量汇总 |
| H | [v3/590_PL2_10w/alignment_shifts.csv](evidence/alignment_shifts.csv.numbered.txt) | 301 | 本轮300个观测偏移；不是已应用的校正 |
| C | [v3/590_PL2_10w/component_parameters.csv](evidence/component_parameters.csv.numbered.txt) | 13 | 本轮12项模型级汇总，尽管名称如此但不含分量参数 |
| V | [v3/tests/p4p5_test_status.csv](evidence/p4p5_test_status.csv.numbered.txt) | 8 | 本轮提交的七项旧pilot测试通过记录；本审阅未执行MATLAB |
| J | [v3/stage_status.json](evidence/stage_status.json.numbered.txt) | 6 | 本轮状态声明，不代表经本审阅验收 |
| U | [v3/run_report.md](evidence/run_report.md.numbered.txt) | 5 | 本轮文字执行声明 |
| T | [tests/test_b1_component_pilot_v2.m](evidence/test_b1_component_pilot_v2.m.numbered.txt) | 113 | 上一轮用户附件《粘贴的 markdown (1)。md》；只作已有测试定义，不冒充本轮源码身份 |

## 实际已做与未做

已做：解析108+12模型汇总、逐对比较Q、九bin成员几何重算、300偏移索引/中位数/极值检查、偏移与六块表内部对账、总信号首末变化、序列lag1样本相关、静态代码审查。

未做：MATLAB执行、输入/源码SHA与用户本地比对、原谱重拟合、实际参数级候选重算、观看未提供PNG、A1变换或统计校准。2592候选/16848参数行等为结构对账预期，不是独立确认数量。

## 外部方法核对（不作为实验证据）

MathWorks官方 `lsqnonlin` 文档用于核对输出和终止条件；Protassov等2002原始方法论文用于说明新增谱分量检验的边界问题。

- https://www.mathworks.com/help/optim/ug/lsqnonlin.html
- https://arxiv.org/abs/astro-ph/0201547

## 复算

在本包目录执行 `python recompute_review_tables.py`。需要numpy与pandas。该程序仅做所贴CSV的代数核对，输出 `recomputed/`；它不是MATLAB分析入口，不会访问用户原始EELS或修改仓库。
