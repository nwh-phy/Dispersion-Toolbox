# C20 q-EELS Codex 执行包 v2

基于 `nwh-phy/Dispersion-Toolbox@18af79b66cd27365fd00bcf3de26658ddd7e7d53`；编制于2026-09-11。

- `C20_qEELS_ExecPlan_v2.md`：完整分阶段任务、原始输入、源码对应、分箱与模型规范、验证和验收。
- `SOURCE_ALIGNMENT.md`：本次实际读取的固定提交来源及与旧审阅的区别。
- `input_manifest.template.yaml`：由本地Codex自动填充的输入台账，不是现有MATLAB运行参数。
- `pilot_spec.json`：首轮实现目标的机器可读规范，不是现有仓库已经支持的配置接口。
- `CODEX_START_PROMPT.txt`：直接交给本地Codex的启动提示词。

将此目录放在授权工程内的独立plans/目录，向Codex发送启动提示词，并提供实际数据根目录。现有库和原始数据不需要在聊天端重新上传。入口名 `run_b1_component_pilot_v2` 等是待实现接口，不能在实现前直接运行。

本包不含原始实验数据、已运行的拟合结果或已修改的源码。没有写入GitHub。本版替代前一版通用执行计划；旧文件继续保留用于追溯。
