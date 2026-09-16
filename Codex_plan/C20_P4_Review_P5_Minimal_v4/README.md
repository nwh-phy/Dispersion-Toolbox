# C20 P4审阅与下一轮最小验证包 v4

先读 `REVIEW_P4_ACTUAL.md`；本地Codex执行用 `NEXT_EXECUTION_PLAN.md` 和 `CODEX_START_PROMPT.txt`。

`evidence/`只是用于本次审阅的用户片段转录，**不要用其中.m覆盖项目源代码，也不要把它作为新的运行程序执行**。`SOURCE_EVIDENCE_INDEX.md`说明行号和来源；`recomputed/`是本次实际Python表格复算，不是新EELS结果。`recompute_review_tables.py`可复现这些代数核对。

本次审阅确认了九个名义同中心bin和CSV汇总数字；未获得实际参数级CSV、新完整拟合MAT或原始EELS，未运行MATLAB。下一轮应修复持久化、排序映射和阶段验收后进行小型P5，不重复全q或全窗口扫描。

将本包解压到现有工程 `codex_plan/C20_P4_Review_P5_Minimal_v4/`；Codex工作区仍为整个工程。让Codex先读取实际本地AGENTS/源码和两个父run，再按提示词执行。
