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

真实 q–E 图、留出实测谱与分量、观测峰脊与条件裸峰顶的图片见
[REAL_DATA_FIGURES.md](REAL_DATA_FIGURES.md)。独立重现入口：

    .venv/bin/python3 case_studies/bisb2026/exploration/composition/real_data_figures.py

它从原始帧重建图片并核对既有正 q 分箱，复用 v3 拟合；不再优化参数，不以理论示意谱代替真实数据。

## 后续：低、高能成分公平放开色散

最新判断见 [DISPERSION_REPORT.md](DISPERSION_REPORT.md)，公式、参数语义和统计边界见
[DISPERSION_METHODS.md](DISPERSION_METHODS.md)。这轮保留冻结的 `models.py`，
`dispersion_models.py` 调用同一个测量算子和背景函数，实现有序的 A/B/C/D 双成分假设与单振子 S。
`dispersion_protocol.toml` 管理主比较和敏感性条件；`dispersion_block10.toml` 只改变时间块检查。

以下是已经实际运行的复现入口，用户无需再手动执行。在根目录的现有 `.venv` 中：

```bash
for sid in 590 n0 20w; do
  for stage in main sensitivity bootstrap; do
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python3 \
      case_studies/bisb2026/exploration/composition/dispersion_run.py \
      --source paper_results/b1_composition_20261023 \
      --output paper_results/b1_ordered_dispersion_20261023 \
      --session "$sid" --stage "$stage"
  done
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python3 \
    case_studies/bisb2026/exploration/composition/dispersion_block_check.py \
    --source paper_results/b1_composition_20261023 \
    --reference paper_results/b1_ordered_dispersion_20261023 \
    --output paper_results/b1_ordered_dispersion_block10_20261023 --session "$sid"
done
.venv/bin/python3 case_studies/bisb2026/exploration/composition/dispersion_summary.py \
  --source paper_results/b1_composition_20261023 --output paper_results/b1_ordered_dispersion_20261023
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python3 \
  case_studies/bisb2026/exploration/composition/dispersion_verify.py \
  --source paper_results/b1_composition_20261023 --output paper_results/b1_ordered_dispersion_20261023 \
  --block-output paper_results/b1_ordered_dispersion_block10_20261023
```

每个已完成拟合都会保存独立结果；数据或计算源码签名不一致时拒绝复用缓存。
长时间块检查通过相对符号链接引用主拟合，并另外核验其哈希。
摘要脚本仅读取结果，不调整曲线以得到预期趋势。它把振子 E0、模型损失峰顶与实际计数分开绘制。

## 热图亮脊与运动学因子的复核

[KINEMATIC_AUDIT.md](KINEMATIC_AUDIT.md) 对照用户提供的 590 三联图，
区分整条总谱亮脊的 1.05 eV 上升与 n0 局部低能 E₀ 的 283 meV 变化。
`kinematic_audit.py` 复现旧图、检查实际狭缝下的除法诊断，并在冻结原始计数上
比较 K=1 和只冻结 K 能量依赖的拟合；`kinematic_audit_verify.py` 检查独立精细
积分与多种峰位定义。旧代码、旧结果保持冻结。

## 高 q 的 3/5 通道合并

[QBIN_REPORT.md](QBIN_REPORT.md) 比较三个数据集在 ±.010、±.015 Å⁻¹ 的
3/5 原始 q 通道合并，能量 bin 保持 16 meV。`qbin_protocol.toml` 指定范围，
`qbin_run.py` 复用准备流程并做独立单谱的单/双成分比较，`qbin_summary.py`
检查复现与积分并绘图。5 通道提高计数信噪比约 25%–28%，部分分峰仍对窗口敏感。
