# C20 q-EELS：接续Codex P0–P3结果的P4/P5计划

把本文件夹放在现有 `Dispersion-Toolbox/codex_plan/` 下，Codex继续打开整个仓库，不需要新建一个与仓库隔离的空工程，也不需要复制原始数据。

默认下一轮目标：**已有边界解审计 → 同中心q分箱 → 590帧诊断 → 三个N3位置重拟合**。不是重做P0–P3，也不是立即把全q两条色散画出来。

## 文件

- `C20_P4P5_ExecPlan_v3.md`：完整设计、阶段任务、测试与交付。
- `P0P3_REVIEW_AND_HANDOFF.md`：本轮源码/报告审阅、证据与待验证项。
- `CODEX_START_PROMPT.txt`：直接发送给本地Codex的执行提示词。
- `next_stage_policy.template.yaml`：设计参数模板，不是已接入MATLAB的配置。
- `EVIDENCE_MAP.json`：用户附件hash及其中源码/报告片段位置。

父run是 `paper_results/b1_components_v2/20260912T021813081Z_31fcc543_18af79b6/`，保持只读。新结果写 `paper_results/b1_components_v3/<run_id>/`。

本包只包含计划与审阅，不包含新的实验拟合结果或可以直接替代现有程序的MATLAB实现。代码名称中的新增v3入口需由Codex先实现、测试后再调用。不要自动重跑父run的审计器覆写旧输出。

把 `CODEX_START_PROMPT.txt` 的内容发给当前Codex对话即可。下一轮返回应包括小型 `review_packet.zip` 和关键图/CSV，避免只返回成功率与文字总结。
