# B1 两种纵向响应：推导和 CPU 模拟

阅读 [推导](DERIVATION.md) 和 [结果](RESULTS.md)。本轮比较两个空间分离的 Drude 通道，与一个 Drude 通道加动态电子极化的自洽响应。没有新实验，也没有把经验分量直接宣布为两个已证实 plasmon。

代码复用现有 Python 环境及 q-EELS 接受积分算子，不修改 MATLAB 或 GUI。

不依赖本地实验文件的演示，在仓库根目录运行：

    .venv/bin/python3 case_studies/bisb2026/exploration/two_mode_theory/run.py \
      --standalone --output paper_results/two_mode_demo

使用已完成 B1 分析中的能量尺度、几何及分量比例：

    .venv/bin/python3 case_studies/bisb2026/exploration/two_mode_theory/run.py \
      --composition paper_results/b1_composition_20261023 \
      --output paper_results/b1_two_mode_theory_20261023

验证：

    .venv/bin/python3 -m pytest tests/python -q

尚无虚拟环境时，依赖清单位于上一级 requirements.txt；用 python3 创建本地 venv，不全局安装。

输出包括三组 PNG/PDF 图、CSV 参数扫描/色散/权重、summary.json 和输入/代码 SHA-256 清单。数值参数在 [config.toml](config.toml) 中，模型方程只有一份实现 [models.py](models.py)。

独立演示采用配置中冻结的参考 E₀ 和几何尺度，且不画实验分量比例。完整模式读取已冻结的经验拟合和分箱文件，不读取原始采集文件，不改动其它分支。两种模式均完成 180 组参数扫描。

扫描明确列出 Δ 的九个位置、r₀ 的五个值及 γ 的四个值，全部结果保留。示例 Δ=1.05 eV、r₀=500 Å、γ=0.65 eV 没有因扫描结果而调整；它没有精准匹配实测低能比例。扫描包络是不同假设的范围，不是统计置信区间。

参考点两条频率由输入能量解析校准，其它动量上的谱形和权重是条件预测。本轮未将这些物理模型直接拟合原始计数来进行模型优选。
