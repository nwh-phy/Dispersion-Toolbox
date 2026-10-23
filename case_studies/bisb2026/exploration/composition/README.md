# B1 宽谱带构成探索

先读 [REPORT.md](REPORT.md) 的判断，再看 [METHODS.md](METHODS.md) 的假说、单位和限制。[LITERATURE.md](LITERATURE.md) 记录原始论文与 DOI。此模块属于 BiSb 个案研究，尚未作为通用 GUI 默认算法。

本轮复用上一轮已验证的能量/q 漂移轨迹，从原始 NPY 重新读取更宽能区。源数据只读，所有计数、预测、模拟和图像均放在 gitignored 的 `paper_results/`。Python 依赖沿用上级目录的 `requirements-lock.txt`。

## 重现正式结果

在项目根目录、已建立的 Python 环境中运行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python3 \
  case_studies/bisb2026/exploration/composition/run.py \
  --previous paper_results/qeels_exploration_20261023 \
  --output paper_results/b1_composition_20261023 --resume
```

这是本轮已实际执行的入口；用户无需再手动运行。新目录运行时去掉 `--resume`。默认最多两个 CPU 任务，支持中断后按数据/代码哈希续跑；代码或数据改变后会拒绝复用旧拟合。上一轮输出的构建方法见 [上级 README](../README.md)。原始数据位置从上一轮 manifest 读取，不写死家目录。

1. `prepare.py`：核验原始哈希，复用漂移轨迹，重建 8 q × 300 帧 × 200 能量通道。
2. `compare.py`：统一的 14 条件 × 4 模型 × 3 数据集，配置来源为 `protocol_v3.toml`。
3. `validate.py`：4 种已知真值、独立精细积分网格、噪声模拟与实际解收敛检查。
4. `summarize.py`：配对时间块重采样、参数和面积表、留出预测图、无色散峰脊反例。
5. `verify.py`：原始文件未变、旧处理数值一致、结果覆盖、源码签名、独立网格恢复与收敛。

`models.py` 是唯一模型实现；`io_utils.py` 只负责哈希和 JSON。保留的 `protocol.toml` / `protocol_v2.toml` 记录早期探索与模型漏洞，正式结果以 v3 为准；旧结果在本地输出的 `fits/` 和 `fits_v2/`。不要混用不同版本的残差排序。

正式产物在输出目录：`fits_v3/`、`synthetic/`、`summary/`、`convergence.json`、`verification.json`。`summary/model_scores.csv` 包含全部状态，不静默丢弃边界拟合。`summary/paired_time_bootstrap.csv` 的区间只条件于已有训练拟合；它不是材料参数的总置信区间。

科学回归检查：

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python3 -m pytest tests/python -q
```

没有修改 MATLAB 公共拟合器、GUI、研究主分支或原始数据。要移植进 GUI，先做通用接口与独立数据验收，不能将本个案的 0.3–1.8 eV 范围和邻近能带约束直接套用所有 q-EELS。
