# 2D-Metal-Plasmon Toolbox

MATLAB toolbox for analyzing **plasmon dispersion** in atomically thin quasi-2D metals from momentum-resolved EELS (q-EELS) data.

Implements the quasi-2D plasmon model from [da Jornada et al., *Nat. Commun.* **11**, 1013 (2020)](https://doi.org/10.1038/s41467-020-14826-8):

$$E(q) = \sqrt{\frac{A \cdot |q|}{(1+\varepsilon_s)/2 + \rho_0 \cdot |q|}}$$

where $\rho_0$ is the interband screening length and $A$ is the Drude weight parameter.

---

## Quick Start

```matlab
% Launch with preprocessed eq3D data
interactive_qe_browser("path/to/eq3D.mat");

% Launch with previously saved operation history
interactive_qe_browser("path/to/eq3D.mat", "path/to/op_history.mat");
```

**输入格式**：仅支持 `eq3D.mat` 文件（包含 `eq3D` 三维矩阵和 `dq` 动量分辨率）。

---

## Toolbox Architecture

```
┌─────────────────────────────────────────────────────────┐
│                interactive_qe_browser.m                 │
│           (GUI: 可视化 · 拟合 · 导出 · 历史)             │
├──────────┬────────────────┬──────────────────────────────┤
│ Data I/O │   Analysis     │   Fitting & Modeling         │
├──────────┼────────────────┼──────────────────────────────┤
│load_qe_  │build_          │fit_loss_function             │
│dataset   │comparison_qe   │  (Lorentz + power-law BG)    │
│          │                │fit_quasi2d_plasmon            │
│make_qe_  │                │  (色散关系拟合 → ρ₀, A)       │
│struct    │                │                              │
└──────────┴────────────────┴──────────────────────────────┘
```

---

## Core Functions

### Data Loading

| Function | Description |
|----------|-------------|
| `load_qe_dataset` | 加载 `eq3D.mat` 并构建标准化 q-E 数据结构，自动计算 q 轴（Å⁻¹）|
| `make_qe_struct` | 从 intensity 矩阵、energy 和 dq 构建标准 q-E struct |

### On-Axis Signal Separation

| Function | Description |
|----------|-------------|
| `build_comparison_qe` | 使用对称 off-axis 参考带重建 on-axis 信号，分离等离激元激发与 ZLP 尾巴。支持 ZLP 面积/峰值归一化 |

### Spectral Fitting

| Function | Description |
|----------|-------------|
| `fit_loss_function` | 单谱 Drude-Lorentz 拟合：$S(E) = B_0 E^{-\alpha} + \sum_i \text{Lorentz}_i$。支持手动指定初始峰位 |
| `fit_quasi2d_plasmon` | 拟合 quasi-2D 等离激元色散关系，提取 $\rho_0$（屏蔽长度）和 $A$（Drude 权重），支持正负 q 同时拟合 |

---

## Interactive q-E Browser

**`interactive_qe_browser.m`** 是 toolbox 的核心交互界面，提供四面板实时可视化和完整的交互式分析流程。

### 界面布局

| Panel | 位置 | 功能 |
|-------|------|------|
| **q-E Heatmap** | 左上 | 动量-能量色散图，点击选择 q channel，叠加已拟合峰位 |
| **Single Spectrum** | 右上 | 当前 q channel 的能谱，支持 Lorentz 拟合叠加显示 |
| **Comparison Map** | 左下 | off-axis 归一化信号图，去除 on-axis ZLP 干扰 |
| **Dispersion** | 右下 | 色散关系散点图 + quasi-2D 拟合曲线 |

### 控件功能

#### Row 1: 数据加载与流程控制
| 控件 | 功能 |
|------|------|
| **Load eq3D** | 加载 `eq3D.mat` 文件 |
| **Build Views** | 构建归一化 off-axis 比较图 |
| **↩ Undo / Reset To** | 撤销操作 / 重置到指定阶段 |
| **View Mode** | Physical（物理 q 轴）/ Normalized（归一化）显示 |

#### Row 2: 动量轴设置
| 控件 | 功能 |
|------|------|
| **dq (1/A)** | 显示当前动量分辨率 |
| **Ovr dq** | 手动覆盖 dq 值 |
| **q start / q end / q step** | 设定 q 轴显示与分析范围 |
| **Exp q equal.** | 指数 q 等化（非均匀在均匀 q 网格上重采样）|

#### Row 3: 能量轴与 Y 轴
| 控件 | 功能 |
|------|------|
| **E min / E max** | 能量范围（meV）|
| **Auto Y / Y min / Y max** | 自动/手动 Y 轴范围 |
| **Y scale** | linear / log 选择 |

#### Row 4: 显示与平滑
| 控件 | 功能 |
|------|------|
| **Trace** | 显示模式：`display`（原始）/ `background-subtracted` / `curvature` |
| **Offset** | waterfall 偏移量 |
| **Base λ** | ALS 基线平滑参数 |
| **Smooth** | 平滑方法：`gaussian` / `off` / `sgolay` |
| **Width** | 平滑窗口宽度 |
| **SG ord / SG len** | Savitzky-Golay 阶数和帧长 |

#### Row 5: 归一化与背景减除
| 控件 | 功能 |
|------|------|
| **ZLP min / max** | ZLP 归一化能量窗口 |
| **Ref\|q\|min / max** | 对称参考带的 q 范围 |
| **Area Norm** | 面积归一化开关 |
| **NormE lo / hi** | 归一化能量积分范围 |
| **BG Sub** | 背景减除开关 |
| **BG Method** | 背景模型：`Power` / `ExpPoly3` / `Pearson` |

#### Row 6: 降噪与反卷积
| 控件 | 功能 |
|------|------|
| **Denoise** | 降噪开关 |
| **Method** | `Wiener2D` / `SavGol` |
| **Noise σ** | 噪声水平（0 = 自动估计）|
| **Deconv** | Lucy-Richardson 反卷积开关 |
| **LR iter** | 反卷积迭代次数 |

#### Row 7: 色散数据管理与导出
| 控件 | 功能 |
|------|------|
| **Pick Peaks** | 在 q-E 热图上交互式点选峰位 |
| **Undo Pt / Clear Pts** | 撤销/清除峰位点 |
| **Save Pts / Load Pts** | 保存/加载色散数据点 |
| **Split 3** | 将散点按能量自动分为多支（branch） |
| **Fit Model** | 对色散数据进行 quasi-2D plasmon 拟合 → 提取 $\rho_0$, $E_{\text{flat}}$ |
| **📷 Export** | 导出面板图像（PNG 300dpi + PDF 矢量） |
| **Auto Fit ω(q)** | 自动批量 Lorentz 拟合所有 q channel |
| **Export Ratio** | 导出图像宽高比：`1:1` / `4:3` / `3:2` / `16:9` / `2:1` |

#### Row 8: Lorentz 拟合参数
| 控件 | 功能 |
|------|------|
| **Prominence** | `findpeaks` 最小突出度 |
| **Smooth W** | 拟合前平滑窗口 |
| **Max Peaks** | 最大峰数 |
| **Guesses** | 手动指定初始峰位（逗号分隔 meV 值，如 `800,2500,3200`） |
| **Fit Spectrum** | 对当前单谱执行 Lorentz + power-law 拟合 |
| **Accept Fit** | 接受拟合结果，将峰位加入色散图 |
| **Show Γ** | 打开独立窗口显示 Γ(q) 和 Γ(E) 分析图 |

### 侧边栏：操作历史 (Operation History)
- 右侧面板记录所有操作步骤
- **💾 Save / 📂 Load**: 保存/读取完整操作历史到 `.mat` 文件
- **Clear History**: 清空历史记录
- 点击历史条目可回溯到任意操作步骤

---

## 工作流程

### 1. 数据加载
```matlab
interactive_qe_browser("eq3D.mat");
% 自动检测 dq 并显示 q-E 热图
```

### 2. 预处理设置
```
调整 q 范围 (q start/end)、能量范围 (E min/max)
→ 可选：勾选 Denoise / BG Sub / Deconv / Area Norm
→ 可选：调整 Smooth 参数
```

### 3. On-Axis 信号分离
```
设置 Ref|q|min 和 Ref|q|max → 点击 Build Views
→ Comparison Map 显示去除 on-axis ZLP 后的纯等离激元信号
```

### 4. 交互式 Lorentz 拟合（逐谱）
```
1. 点击 q-E 热图选择一个 q channel
2. 在 Guesses 输入估计峰位，如 "800, 2500, 3200"
3. 调整 Prominence / Smooth W / Max Peaks
4. 点击 Fit Spectrum → 检查拟合叠加曲线
   - 红色实线: 总拟合
   - 灰色点线: power-law background
   - 彩色虚线: 各 Lorentz 峰分量
5. 满意后点击 Accept Fit → 峰位添加到色散图
6. 切换下一个 q channel，重复 2-5
```

### 5. 自动批量拟合
```
在 Guesses 中填入全局峰位猜测 → 点击 Auto Fit ω(q)
→ 自动对所有显示范围内的 q channel 进行 Lorentz 拟合
→ 结果按 branch 分类显示在色散图上
```

### 6. 色散模型拟合
```
点击 Fit Model → 对当前色散数据进行 quasi-2D plasmon 拟合
→ 提取 ρ₀ (screening length)、E_flat (无色散极限能量)
→ 拟合曲线叠加显示在 Dispersion 面板
```

### 7. 线宽分析
```
完成 Auto Fit 后 → 点击 Show Γ
→ 独立窗口显示各 branch 的 Γ(q) 和 Γ(E)
→ 用于分析等离激元阻尼机制
```

### 8. 图像导出
```
选择导出宽高比 → 点击 📷 Export → 选择面板 → 选择保存路径
→ 同时输出 PNG (300 dpi) + PDF (矢量格式)
→ 图表内容按所选宽高比自动拉伸
```

---

## Fitting Model Details

### Single-Spectrum Loss Function

$$S(E) = B_0 \cdot E^{-\alpha} + \sum_{i=1}^{N} \frac{A_i \cdot E \cdot \Gamma_i}{(E^2 - \omega_{p,i}^2)^2 + E^2 \cdot \Gamma_i^2}$$

- **Power-law background** ($B_0 \cdot E^{-\alpha}$): 显式建模 ZLP 尾巴，避免 Lorentz 峰扭曲
- **Drude-Lorentz peaks**: 每个峰由 ($\omega_p$, $\Gamma$, $A$) 参数化

### Quasi-2D Plasmon Dispersion

$$E(q) = \sqrt{\frac{A \cdot |q|}{(1+\varepsilon_s)/2 + \rho_0 \cdot |q|}}$$

**物理含义**:
- 小 $|q|$: $E \propto \sqrt{|q|}$ — 类 2D plasmon 行为
- 大 $|q|$: $E \to E_{\text{flat}} = \sqrt{A/\rho_0}$ — 无色散极限
- $\rho_0$: 带间屏蔽长度，反映非局域介电屏蔽效应强度
- $q_c = (1+\varepsilon_s)/(2\rho_0)$: 色散→无色散转折动量

---

## Dependencies

- **MATLAB R2024b+** (需要 `uigridlayout`, `uifigure`, `exportgraphics` 支持)
- **Optimization Toolbox** (`lsqcurvefit`)
- **Signal Processing Toolbox** (`findpeaks`, `smoothdata`)
- [Nion-EELS 3D Toolbox](../Nion-EELS数据处理toolbox/) (可选，用于 BG Sub 中的 ALS 背景减除)

---

## File Structure

```
2D-metal-plasmon/
├── interactive_qe_browser.m     # 交互式 GUI 主入口
├── load_qe_dataset.m            # eq3D 数据加载 → 标准化 q-E struct
├── make_qe_struct.m             # q-E 数据结构构造函数
├── build_comparison_qe.m        # On-axis 信号分离
├── fit_loss_function.m          # 单谱 Lorentz + power-law 拟合
├── fit_quasi2d_plasmon.m        # Quasi-2D 色散拟合 (ρ₀, A)
└── references/                  # 参考文献 markdown
```

---

## References

> da Jornada, F. H., Xian, L., Rubio, A., & Louie, S. G. (2020). Universal slow plasmons and giant field enhancement in atomically thin quasi-two-dimensional metals. *Nature Communications*, **11**, 1013.

> Do, H. T. B., Zhao, M., Li, P., et al. (2025). Slow and highly confined plasmons observed in atomically thin TaS₂. *Nature Communications*, **16**, 5801.
