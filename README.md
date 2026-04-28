# q-EELS Plasmon Dispersion Toolbox

用于 **动量分辨低损 EELS (q-EELS)** 数据处理与色散分析的 MATLAB 工具包。从原始数据出发，通过 GUI 完成去噪、背景扣除、峰拟合和色散提取。

---

## 安装

```matlab
cd('path/to/Dispersion-Toolbox')
startup   % 自动加载 src/ 和 lib/ 路径
```

如需外部 Nion 工具箱，复制 `toolbox_config.m.template` → `toolbox_config.m` 并设置路径。

**环境要求**：MATLAB R2024b+，需要 Optimization Toolbox、Curve Fitting Toolbox、Signal Processing Toolbox、Image Processing Toolbox。

---

## 支持的数据格式

| 输入 | 说明 |
|---|---|
| `eq3D.mat` | 已预处理的 2D q-E 矩阵（变量 `a3` + `e`） |
| `.npy` | Nion Swift 原始数据（自动查找同名 `.json` 元数据） |
| 3D / 4D `.mat` | 兼容的 EELS 对象或数值数组 |
| 文件夹 | 自动查找 `eq3D_processed.mat` → `eq3D.mat` → 最大 `.npy` |

---

## GUI 工作流

启动 GUI：

```matlab
interactive_qe_browser("path/to/data")
```

也可以恢复之前的操作历史：

```matlab
interactive_qe_browser("path/to/data", "path/to/op_history.mat")
```

以下是典型的数据处理流程：

### Step 1：加载数据

- 点击 **Load Data** 选择数据文件，或直接在启动时传入路径
- 检查 **dq (1/A)** 是否正确，如需修改可在 **Ovr dq** 中覆盖
- 加载后自动显示 q-E map

### Step 2：调整显示范围

- **q start / q end**：设置感兴趣的 q 范围
- **E min / E max**：设置能量窗口（同时用于拟合范围）
- **Trace**：切换显示模式（原始 / 背景扣除后 / 曲率）

### Step 3：预处理

按需开启以下选项（顺序：归一化 → 去噪 → 背景扣除 → 反卷积）：

| 控件 | 作用 | 建议 |
|---|---|---|
| **Normalize** + `ZLP Peak` | 用 ZLP 峰高归一化，校正束流差异 | 一般选 ZLP Peak |
| **Denoise** + `Wiener2D` | 2D 自适应去噪 | 大多数情况足够 |
| **QE BG** + `Auto` | 自动选择最佳模型扣除 ZLP 尾巴 | 推荐先用 Auto |
| **BG W1** | 背景拟合窗口（meV） | 应在真实峰之前 |
| **Dual** + **BG W2** | 双窗口锚定（峰前后各取一段） | 有明确的高能侧无峰区时使用 |
| **Deconv** | Lucy-Richardson 反卷积 | 可选，注意不要过度 |

> **关键**：背景窗口不要覆盖真实的 loss 峰！先在不同 q 通道检查扣除效果。

### Step 4：拟合色散

有两种方式提取 ω(q)：

**方式 A — 自动拟合**（推荐首选）：

1. 设置 **Prominence**、**Max Peaks**（每个 q 通道最多几个峰）
2. 选择 **Peak Model**（一般用 Lorentz）
3. 点击 **Auto Fit ω(q)**
4. 检查散点是否合理

**方式 B — Seed propagation**（峰位漂移大时更稳定）：

1. 在 **Guesses** 文本框输入峰位估计（如 `800,2500`），或点击 **Pick Guesses** 在谱上点选
2. 点击 **Auto Fit ω(q)**
3. 算法会从 seed 点逐通道追踪峰位

**方式 C — 手动选点**：

1. 点击 **Pick Peaks**，在 q-E map 上逐点选取
2. 用 **Split 3** 分离不同 branch

### Step 5：色散拟合

1. 选择 **Disp Model**（如 Quasi-2D Plasmon）
2. 点击 **Fit Dispersion**
3. 拟合曲线和参数会显示在色散图中

### Step 6：诊断检查

- **Γ & Amp**：查看线宽 Γ(q)、振幅 A(q) 和品质因子的 dashboard
- **I_kin Map**：生成动能校正 loss map

### Step 7：导出

| 控件 | 导出内容 |
|---|---|
| **Export** (📷) | 当前显示的图像（q-E map / 单谱 / 色散图） |
| **Export Data** (📦) | 预处理后的数据 `.mat`（含 axes、预处理参数、背景诊断） |
| **Export Disp.** | 色散拟合结果 |
| **Save Pts** | 峰位点数据 |
| History 面板 **Save** | 操作历史 `op_history.mat`，可在下次加载时恢复 |

---

## 推荐流程

对新数据，建议按以下顺序检查：

1. 确认数据加载正确，检查能量轴和 q 轴
2. 在不同 q 通道浏览原始谱，了解数据质量
3. 开启 Normalize (ZLP Peak) + Denoise (Wiener2D)
4. 开启 QE BG (Auto)，检查几个代表性 q 通道的背景扣除效果
5. 调整背景窗口，避免过扣除（切换到 `background-subtracted` trace 查看）
6. 先做单谱拟合（**Fit Spectrum**），确认峰型合理
7. 再做全局自动拟合（**Auto Fit ω(q)**）
8. 检查 branch assignment 后拟合色散
9. 保存操作历史，确保可复现

---

## 背景扣除模型

| 模型 | 简述 |
|---|---|
| `Power` | 幂律背景（最简单） |
| `Power2` | 幂律 + 常数项 |
| `ExpPoly3` | log-intensity 多项式 |
| `Pearson` | log-log 二次经验模型 |
| `PearsonVII` | Pearson-VII ZLP 尾模型 |
| `Exp2` | 双指数 |
| `DualWindow` | 双窗口 side-band 锚定 |
| `Auto` | 自动对候选模型打分，选最优 |

---

## 峰型与色散模型

**峰型模型**（`Peak Model`）：

| 模型 | 适用场景 |
|---|---|
| Lorentz (Drude) | 金属等离激元（默认推荐） |
| Gaussian | 分辨率限制的窄峰 |
| Pseudo-Voigt | Lorentz + Gaussian 混合 |
| Damped HO | 有限温度声子 |

**色散模型**（`Disp Model`）：

| 模型 | 物理含义 |
|---|---|
| Quasi-2D Plasmon | √(Aq/(ε+ρ₀q))，2D 等离激元色散 |
| Acoustic Linear | E = v_s·q，声学支 |
| Optical Constant | E = ω₀，光学常数支 |
| Optical Quadratic | E = ω₀ − βq²，负色散光学支 |

---

## 脚本接口（可复现）

GUI 适合交互检查；可复现流程建议用脚本：

```matlab
startup
dataset = load_qe_dataset("path/to/data", dq_Ainv);
qe = dataset.qe;

% 预处理
opts = struct();
opts.do_normalize = true;  opts.norm_method = 'ZLP Peak';
opts.norm_min = -50;       opts.norm_max = 50;
opts.do_denoise = true;    opts.denoise_method = 'Wiener2D';  opts.denoise_sigma = 0;
opts.do_bg_sub = true;     opts.bg_method = 'Auto';
opts.bg_win_lo = [50, 300]; opts.bg_win_hi = [];
opts.do_deconv = false;
[qe_pp, bg_diag] = qe_preprocess(qe, opts);

% 拟合
fit_opts = struct();
fit_opts.E_min = 300;      fit_opts.E_max = 3000;
fit_opts.q_start = min(qe_pp.q_Ainv);
fit_opts.q_end = max(qe_pp.q_Ainv);
fit_opts.peak_model = 'lorentz';
fit_opts.max_peaks = 2;
fit_opts.energy_mask = qe_pp.energy_meV >= fit_opts.E_min & qe_pp.energy_meV <= fit_opts.E_max;
fit_opts.energy_axis = qe_pp.energy_meV(fit_opts.energy_mask);
results = qe_auto_fit(qe_pp, qe_pp, fit_opts);

% 色散拟合
fit_result = fit_dispersion_generic( ...
    results.branches{1}(:,1), results.branches{1}(:,2), ...
    'model', 'quasi2d_plasmon');
```

参数应由具体数据决定，以上仅为示例。

---

## 测试

```matlab
cd('path/to/Dispersion-Toolbox')
startup
runtests('tests')
```

---

## License

This project is provided for academic and research use within the group.
