# q-EELS Plasmon Dispersion Toolbox

面向课题组的 MATLAB 工具箱，用于从 **Nion 原始 4D-EELS 数据**（或已预处理的 eq3D 数据）中提取等离激元色散关系。适用于多种二维金属薄膜材料体系。

支持完整的数据处理链：原始数据导入 → 帧累加 → q 裁切 → ZLP 对齐 → 能量标定 → 交互式分析 → 色散拟合 → 出图。

---

## 谁需要这个工具箱

- 做 q-EELS（动量分辨 EELS）实验，需要从 q-E 图中提取等离激元色散关系
- 数据来自 Nion 电镜，格式为 `.npy + .json` 或 Nion 4D Toolbox 处理后的 `EELS4D` `.mat` 文件
- 已有课题组 Nion 工具箱处理过的 `eq3D.mat` 数据

---

## 快速开始

### 首次使用：初始化路径

```matlab
cd('path/to/2D-metal-plasmon')   % cd 到项目根目录
startup                           % 自动将 src/ 加入 MATLAB 路径
```

> **注意**：处理 `.npy` 原始数据需要预先将 [Nion EELS 4D Toolbox](https://github.com/Gao-Group) 加入 MATLAB 路径（包含 `readNPYtoEELS`、`Align4DEELS` 等函数）。如果仅使用已处理的 `eq3D.mat` 文件，则不需要该 Toolbox。

### 方式一：从原始 Nion 数据开始（推荐）

```matlab
% 直接传入 .npy 文件路径
interactive_qe_browser("path/to/Sequence EELS Image 273.npy")

% 或传入 Nion 4D Toolbox 预对齐的 .mat（含 EELS4D 对象）
interactive_qe_browser("path/to/Sequence EELS Image 273_Align.mat")
```

首次加载时，工具箱会自动完成以下预处理：
1. **帧累加**：将 N 帧序列（如 300 帧）求和
2. **q 裁切**：自动检测有效 q 范围（或弹出对话框手动设置）
3. **ZLP 对齐**：移植自 Nion 4D Toolbox 的 `Align4DEELS` 算法
4. **能量标定**：从 `.json` 元数据或 `EELS4D.ene` 构建物理能量轴
5. **缓存**：处理结果自动保存为 `eq3D_processed.mat`，下次加载秒开

### 方式二：从已有 eq3D 数据开始

```matlab
% 传入 eq3D.mat 文件
interactive_qe_browser("path/to/eq3D.mat")

% 传入文件夹——自动查找 eq3D_processed.mat 或 eq3D.mat
interactive_qe_browser("path/to/data_folder/")
```

### 方式三：恢复之前的分析会话

```matlab
% 加载数据 + 恢复操作历史
interactive_qe_browser("path/to/eq3D.mat", "path/to/op_history.mat")
```

---

## 与 Nion 工具箱的关系

本工具箱**不替代** Nion EELS 工具箱，而是接在它的下游：

```
  Nion 采集                  Nion 4D Toolbox（可选）           本工具箱
┌──────────┐    .npy+.json   ┌───────────────────┐  _Align.mat  ┌─────────────────────┐
│ Nion     │ ─────────────→  │ readNPYtoEELS     │ ──────────→  │                     │
│ Swift    │                 │ Align4DEELS       │              │  load_raw_session    │
│          │                 │ Denoise4DEELS     │              │  (帧累加+对齐+标定)  │
└──────────┘                 │ Squeeze4DEELS     │              │         ↓            │
                             └───────────────────┘              │  interactive_qe_     │
                                                                │  browser             │
  如果你已经有 eq3D.mat：                                        │  (交互式分析)       │
     dealTBG 或其他脚本生成 ─────────────────────────────────→  │         ↓            │
                                                                │  色散拟合 + 出图     │
                                                                └─────────────────────┘
```

**关键区别**：
- Nion 工具箱：专注于 4D 数据的对齐、去噪、q 畸变校正、降维（4D → 3D）
- 本工具箱：专注于 2D q-E 谱的**物理分析**——峰拟合、色散提取、模型拟合

**你可以跳过 Nion 工具箱**：直接将 `.npy` 传给本工具箱，内置的 `load_raw_session` 已经移植了核心的 ZLP 对齐算法。

---

## 完整工作流

### Step 1：加载数据

点击 GUI 左上角的 **「Load Data」** 按钮，支持三种数据格式：

| 格式 | 文件类型 | 说明 |
|------|----------|------|
| **原始 Nion** | `.npy` + `.json` | 从 Nion Swift 直接导出的序列数据 |
| **Nion 4D** | `_Align.mat` | 经 `Align4DEELS` 对齐后的 `EELS4D` 对象 |
| **eq3D** | `eq3D.mat` | 已预处理的 2D (E × q) 矩阵 |

加载时自动识别格式，首次处理的结果会被缓存（`eq3D_processed.mat`），后续加载瞬间完成。

### Step 2：预处理

在 GUI 控件区可按需启用：

- **Despike**：宇宙射线点缺陷移除（中值滤波 + 5σ 阈值检测）
- **Denoise**：Wiener 2D 去噪 / Savitzky-Golay 滤波 / **SVD/PCA 截断降噪**（自动拐点检测或手动设定保留分量数）
- **Deconv**：Lucy-Richardson 迭代退卷积 / **Fourier-Log 多次散射消除**（Egerton §4.2）
- **BG Sub**：背景扣除（详见下方）
- **Build Views**：用对称的离轴参考带重建在轴信号

#### 背景扣除（BG Sub）

基于 Fung et al. (2020) 的诊断透明框架和 Ruishi Qi (Nion EELS Toolbox) 的鲁棒拟合策略，实现了多模型、双窗口、迭代重拟合的背景扣除引擎。

**6 种背景模型**：

| 模型 | 函数形式 | 适用场景 | MATLAB 实现 |
|------|----------|----------|-------------|
| **Power** (默认) | A·E^B | 标准 EELS 背景（最常用） | `fit(e, s, 'power1')` |
| **Power2** | a·E^b + c | 带常数偏移的幂律 | `fittype('a*x^b + c')` |
| **ExpPoly3** | exp(p₃·E³+p₂·E²+p₁·E+p₀) | 快速衰减的背景 | `fit(e_norm, log(s), 'poly3')` |
| **Pearson** | exp(p₂·(logE)²+p₁·logE+p₀) | Log-log 二次多项式 | `fit(le_norm, log(s), 'poly2')` |
| **PearsonVII** | I·w^(2m)/(w²+(2^(1/m)-1)·(2E-2t₀)²)^m | ZLP 物理尾部模型 | `lsqcurvefit` + 非对称惩罚 |
| **Exp2** | a·e^(bE) + c·e^(dE) | 双指数衰减背景 | `fit(e, s, 'exp2')` |

> **推荐**：对于大多数 EELS 数据，默认 **Power** 模型即可得到良好结果（R² > 0.95）。当 ZLP 尾部形状需要更精确建模时，使用 **PearsonVII**。

**双窗口拟合**（Dual Window）：

在前景信号的两侧分别设定拟合窗口，用两端的背景特征锚定拟合曲线：

```
窗口 1 (bg_win_lo)：[50, 300] meV   ← ZLP 和信号之间的前景区
窗口 2 (bg_win_hi)：[3000, 3500] meV ← 信号之后的高能区
```

勾选 GUI 中的 **Dual** 复选框启用。适合信号两侧都有可靠背景参考的情况。

**迭代重拟合**（Iterative Refit）：

勾选 **Iter** 复选框启用。逻辑为：
1. 首次拟合后，检查信号区域是否有 `spec − bg < 0`（过扣除）
2. 如果存在负值 → 增大这些区域的拟合权重 → 重新拟合
3. 防止背景曲线过高导致有效信号变成负值

**拟合诊断**（自动显示在 FitInfoLabel）：

| 指标 | 含义 | 良好标准 |
|------|------|----------|
| **R²** | 决定系数（拟合优度） | > 0.9 |
| **RMSE** | 均方根残差 | 越小越好 |
| **h** | Fung h-参数 = (I_B + var(I_B))/I_B | — |
| **SNR** | 信噪比 = I_k / √(I_k + h·I_B) | > 3 |

**可视化**：在单谱模式（SingleAxes）勾选 BG Sub 后，自动叠加三线图：
- 🔵 蓝线 — 原始谱（扣除前）
- 🔴 红色虚线 — 背景拟合曲线
- 🟢 绿线 — 扣除后信号
- 浅红色填充 — 拟合窗口区域

**编程接口**（用于批处理脚本）：

```matlab
opts.do_bg_sub    = true;
opts.bg_method    = 'Power';           % 模型名称
opts.bg_win_lo    = [50, 300];         % 窗口 1 [meV]
opts.bg_win_hi    = [3000, 3500];      % 窗口 2 [meV]（仅双窗口）
opts.bg_iterative = true;             % 启用迭代重拟合

% 单输出（向后兼容）
qe_out = qe_preprocess(qe_in, opts);

% 双输出（获取诊断）
[qe_out, bg_diag] = qe_preprocess(qe_in, opts);
% bg_diag(qi).rsquare, .rmse, .h_param, .snr, .bg_curve, .residuals
```

### Step 3：峰拟合

1. **手动拟合**：在 q-E 热图上点击选择 q 通道 → 在 Guesses 栏输入峰位（如 `800,2500`）→ Fit Spectrum
2. **自动拟合**：点击 **Auto Fit ω(q)** 批量拟合所有 q 通道

拟合模型为 Drude-Lorentz 损失函数（默认，另支持 Gaussian、Voigt、Damped HO 峰型）：

```
S(E) = B₀·E^(-α) + Σᵢ [ Aᵢ·E·Γᵢ / ((E² - ωₚᵢ²)² + E²·Γᵢ²) ]
```

其中每个 Lorentz 峰由三个物理参数描述：峰位 `ωₚ`、线宽 `Γ`、幅度 `A`。

### Step 4：色散拟合

拟合好各 q 通道的峰位后：

1. **Split** 分支（将等离激元峰与带间跃迁分开）
2. **Fit Dispersion** → 选择色散模型
3. **Γ & Amp** → 打开线宽分析面板（6 子图）：
   - Γ(q)、Γ(E) — 线宽随动量/能量的变化
   - Amplitude(q)、Amplitude(E) — 模型振幅
   - Raw A(q) — 原始峰高 log-log 图 + 幂律拟合
   - Q-factor(q) — 等离激元品质因子 ωₚ/Γ

内置色散模型包括：

| 模型 | 公式 | 适用场景 |
|------|------|----------|
| Quasi-2D Plasmon | E = sqrt(A·abs(q) / (ε_bg + ρ₀·abs(q))) | 准二维金属薄膜（ε_bg 默认=1，即 suspended） |
| Acoustic Linear | E = v_s·abs(q) | 声学声子/线性色散 |
| Optical Constant | E = ω₀ | 光学声子平带（无色散） |
| Optical Quadratic | E = ω₀ − β·q² | 光学声子抛物色散 |

### Step 5：出图

点击 **Export** 按钮，支持选择面板和宽高比，输出 300 dpi PNG + 矢量 PDF。

---

## 文件结构

```
├── startup.m                    # 路径初始化（cd 到项目根目录后运行）
│
├── src/                         # 所有 MATLAB 源码
│   ├── interactive_qe_browser.m #   主入口：状态管理、回调、绘图（~2850 行）
│   ├── qe_preprocess.m          #   预处理管线：归一化 → 去噪 → 背景扣除 → 反卷积
│   ├── qe_auto_fit.m            #   批量多峰 Drude-Lorentz 拟合
│   ├── qe_plot_helpers.m        #   绘图工具类（分支配色、散点/误差棒、色散叠加）
│   ├── qe_gamma_dashboard.m     #   线宽 Γ(q)/Γ(E)、振幅、品质因子分析面板
│   ├── build_comparison_qe.m    #   在轴信号分离（对称离轴参考）
│   ├── analyze_dispersion.m     #   自动化批处理脚本
│   │
│   ├── io/                      #   数据导入模块
│   │   ├── load_qe_dataset.m    #     数据加载路由（自动识别格式）
│   │   ├── load_raw_session.m   #     统一原始数据导入（.npy + 4D .mat）
│   │   ├── load_npy_session.m   #     NPY 导入流水线
│   │   ├── read_npy.m           #     NPY 文件读取器（自包含）
│   │   └── make_qe_struct.m     #     q-E 数据结构构造器
│   │
│   ├── fitting/                 #   拟合 & 色散分析
│   │   ├── fit_loss_function.m  #     单谱 Drude-Lorentz 拟合
│   │   ├── peak_models.m        #     峰型模型（Lorentz/Gaussian/Voigt/Damped HO）
│   │   ├── propagate_seed_peaks.m #   峰位跟踪/传播
│   │   ├── fit_quasi2d_plasmon.m #    准二维等离激元色散拟合
│   │   ├── fit_dispersion_generic.m # 通用色散拟合接口
│   │   └── dispersion_models.m  #     色散关系模型库（quasi-2D/linear/power-law）
│   │
│   └── ui/                      #   GUI UI 布局
│       └── qe_browser_ui.m      #     声明式 UI 布局（纯控件，不访问状态）
│
├── figures/                     # 导出图片 & 分析结果
├── docs/                        # 个人文档（开题报告、周报等）
├── references/                  # 参考文献（Markdown 格式）
├── 20260120 Bi/                 # 实验数据
└── graphene/                    # 石墨烯数据
```

---

## 数据文件说明

数据文件（`.mat`, `.npy`, `.json` 等）不包含在 Git 仓库中。需要在本地放置：

```
├── your_sample/
│   ├── Sequence EELS Image XXX.npy      # 原始 Nion 数据
│   ├── Sequence EELS Image XXX.json     # 元数据（能量标定等）
│   ├── eq3D.mat                         # 已有的预处理数据（可选）
│   └── eq3D_processed.mat              # 本工具箱自动生成的缓存
```

---

## dq（动量分辨率）标定

q 轴的物理标定依赖 Nion 的离焦设置：

| 离焦 | dq (Å⁻¹/pixel) | 路径关键字 |
|------|-----------------|-----------|
| 20w (200 μm) | 0.0025 | 路径含 `20w` |
| 10w (100 μm) | 0.005 | 路径含 `10w` |

如果路径中没有离焦信息，GUI 会弹出对话框让你手动输入。也可以在 GUI 控件区的 **Ovr dq** 栏直接覆盖。

---

## 环境要求

- **MATLAB R2024b+**
- Optimization Toolbox（拟合必需）
- **Curve Fitting Toolbox**（背景扣除必需）
- Signal Processing Toolbox（去噪、滤波）
- Image Processing Toolbox（Lucy-Richardson 反卷积；有手动 fallback）
- **不依赖** Nion EELS 工具箱（原始数据导入已内置）

---

## License

This project is provided for academic and research use within the group.
