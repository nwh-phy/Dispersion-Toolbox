# C20 v4审阅与v5数据约束计划

先读 `REVIEW_v4.md`，再让Codex执行 `CODEX_START_PROMPT.txt` 与 `NEXT_EXECUTION_PLAN_v5.md`。

本包包含本次独立MAT/CSV回读脚本及已运行结果、精确源码证据索引、小文件快照。没有新实验拟合结果，也没有复制巨大raw或源交付中的全部MAT。

放入现有项目 `codex_plan/C20_v4_Review_P5_Plan_v5/`，Codex工作区仍打开整个项目根目录。所有父结果只读，新输出独立到 `paper_results/b1_components_v5/`。

独立复核重跑：
```bash
python check_delivery.py /path/to/extracted_v4_packet /path/to/new_audit_output
```
复核输出必须位于输入包外；此脚本不调用MATLAB、不拟合新谱、不更改输入。
