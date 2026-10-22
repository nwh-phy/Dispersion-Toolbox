# q-EELS Plasmon Dispersion Toolbox

用于 **动量分辨低损 EELS (q-EELS)** 数据处理与色散分析的 MATLAB 工具包。从原始数据出发，通过 GUI 完成去噪、可选背景扣除、峰拟合和色散提取。

> **当前 BiSb 主线（2026-05-05）**：主色散结论优先使用不扣背景的 q-E 数据提取峰位。背景扣除保留为可选的稳健性检查和诊断工具，不作为证明 quasi-2D plasmon dispersion 的必要前提。

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

按需开启以下选项（顺序：归一化 → 去噪 → 可选背景扣除 → 反卷积）：

| 控件 | 作用 | 建议 |
|---|---|---|
| **Normalize** + `ZLP Peak` | 用 ZLP 峰高归一化，校正束流差异 | 一般选 ZLP Peak |
| **Denoise** + `Wiener2D` | 2D 自适应去噪 | 大多数情况足够 |
| **QE BG** + `Auto` | 自动选择最佳模型扣除 ZLP 尾巴 | BiSb 主线默认关闭，仅作稳健性检查 |
| **BG W1** | 背景拟合窗口（meV） | 应在真实峰之前 |
| **Dual** + **BG W2** | 双窗口锚定（峰前后各取一段） | 有明确的高能侧无峰区时使用 |
| **Deconv** | Lucy-Richardson 反卷积 | 可选，注意不要过度 |

> **关键**：BiSb 主线先不扣背景提取峰位；如果为了稳健性检查而开启背景扣除，背景窗口不要覆盖真实的 loss 峰。

### Step 4：拟合色散

有四种方式提取或修正 ω(q)：

**方式 A — 自动拟合**（推荐首选）：

1. 设置 **Prominence**、**Max Peaks**（每个 q 通道最多几个峰）
2. 选择 **Peak Model**（一般用 Drude-Lorentz）
3. 点击 **Auto Fit Branches**
4. 检查散点和 branch assignment 是否合理

**方式 B — Seed propagation**（峰位漂移大时更稳定）：

1. 在 **Guesses** 文本框输入峰位估计（如 `800,2500`），或点击 **Pick Guesses** 在谱上点选
2. 点击 **Auto Fit Branches**
3. 算法会从 seed 点逐通道追踪峰位

**方式 C — 自动拟合 + 局部人工矫正**（论文级 curated workflow 推荐）：

1. 先点击 **Auto Fit Branches**，把自动结果作为 branch 底稿
2. 在 q-E map 中选中需要复查的 q 通道，并在单谱窗口确认正确峰位
3. 用 correction target 下拉框选择 **Auto**、**Branch N** 或 **New Branch N+1**
   - **Auto**：自动选择当前 q 附近、能量最接近的 branch 点
   - **Branch N**：强制把本次矫正写入指定 branch
   - **New Branch N+1**：把本次矫正写入新 branch
4. 点击 **Correct Auto**，再在单谱窗口点击正确峰位；被矫正点会以黄色菱形标出
5. 如误点，点击 **Undo Corr** 撤销最近一次 correction
6. 也可先用 **Fit Spectrum** 拟合当前单谱，再用 **Accept Fit** 把拟合峰位写回所选 target branch

**单谱拟合方法**（`Fit Method`，最后一行控件）：

- **BG + peaks**（默认，旧方法）：在预处理后的谱上拟合幂律背景加若干峰（`fit_loss_function`）。
- **ZLP joint**：直接拟合所选通道的原始计数，ZLP 和损失峰一起拟合，不扣背景，不受归一化、去噪、背景扣除设置影响（`qe_fit_spectrum_joint` → `qe_zlp_joint_fit`）。拟合窗为 max(E min, 250) 到 E max；峰数取 **Max Peaks**（恰好这么多个，不是上限）；峰形只能选 **Drude-Lorentz** 或 **Sym. Lorentz**。**Kinematic factor**（默认开）给损失峰乘上该通道的 2D 运动学因子（`qe_kinematic_prefactor`，按 **Beam kV**），小 q 处不开会让 E0 系统性偏低（30 kV、q = 0.002 Å⁻¹、Γ = 400 meV 的合成谱上偏约 70 meV）；**Poisson weights** 按 σ = √counts 加权。信息栏给出 E0、峰顶（损失函数最大值，阻尼 DL 的峰顶低于 E0）和 χ²_red，比较 1 峰和 2 峰时看 χ²_red。叠加图里灰线是实际拟合的原始计数（按当前显示缩放）。这种拟合还没有误差棒。

**方式 D — 手动选点**（旧式 fallback）：

1. 点击 **Pick Peaks**，在 q-E map 上逐点选取
2. 用 **Split 3** 分离不同 branch

### Step 5：色散拟合

1. 选择 **Disp Model**（如 Quasi-2D Plasmon）
2. 选择要拟合的 branch，点击 **Fit Curated**
3. 拟合曲线和参数会显示在色散图中

### Step 6：诊断检查

- **Γ & Amp**：查看线宽 Γ(q)、振幅 A(q) 和品质因子的 dashboard
- **I_kin Map**：生成动能校正 loss map

### Step 7：导出

| 控件 | 导出内容 |
|---|---|
| **Export Panels** | 当前图像（q-E map / 单谱 / 色散图）；Publication 模式同时写 PNG 和 PDF，保留所选比例 |
| **Export Data** (📦) | 预处理后的数据 `.mat`（含 axes、预处理参数、背景诊断） |
| **Export Disp.** | 色散拟合结果 |
| **Save Session / Load Session** | 一份 `.mat` 保存/恢复当前 q–E 数据、标定和控件、分支点/修正、自动/色散拟合、单谱拟合记录及当前曲线 |
| **Save Curated / Load Curated**（原 Save Pts / Load Pts） | 单独交换峰位点：signed `q_Ainv`、`energy_meV`、`branches`、`corrections`、GUI `snapshot` 和新增的 `accepted_fits` 来源记录 |
| History 面板 **Save** | 操作历史 `op_history.mat`，可在下次加载时恢复 |

**日常保存与重开优先使用 Save Session / Load Session。** 会话文件内嵌当前投影后的物理 q–E 计数和坐标，不必重新选择数据、点文件和历史文件。也可以直接运行 `interactive_qe_browser('qe_session.mat')` 打开会话。尚未接受的单谱拟合会恢复曲线和 Accept 状态；已接受且仍在显示的曲线也会恢复，Accept 保持禁用。

会话格式有明确版本号，保存时用同目录临时文件完成写入后再替换目标。加载会检查数据/标定指纹和记录的来源文件 SHA-256；可访问的来源若有变化，会拒绝加载并保留当前工作。来源暂不可访问时，可使用内嵌 q–E 数据恢复，并在界面明确提示。保存的是 q–E 工作阶段，不包含原始 4D 张量及其裁剪/对齐/分箱的完整中间态。

每次接受单谱拟合还会保留完整数值结果、signed q、dq、拟合时的设置、实际信号窗/拟合窗、输入类型、来源指纹、有效性标记和所接受能量的定义（损失函数峰顶或模型 E0）。`accepted_fits` 是历史接受记录；之后删除、重分配或修正点，不会让旧记录变成当前点的重新验证。旧点文件没有这些信息时保持为空，不补造来源。

旧的点文件与 History 文件仍可单独使用：先加载数据，再 **Load Curated**，最后 History **Load**。History 只保存参数，保存时会追加当前 live snapshot；不含归一化方法的旧历史按原默认 `Area` 恢复。

会话格式、测试结果及原生对话框验证边界见 [会话验收报告](reports/gui_session_acceptance_20261022.md)。

重新拟合、切换通道或修改拟合/预处理设置会使待接受的单谱结果失效。拟合失败后，**Accept Fit** 不会继续接受上一次结果。
调整导出模式或导出比例会保留拟合曲线和待接受结果。历史文件缺少必需字段时会在修改界面前拒绝加载；恢复过程若遇到无效参数等错误，会回退数据和控件，并清除待接受拟合，避免使用与当前标定不一致的旧结果。

---

## 推荐流程

对新数据，建议按以下顺序检查：

1. 确认数据加载正确，检查能量轴和 q 轴
2. 在不同 q 通道浏览原始谱，了解数据质量
3. 开启 Normalize (ZLP Peak) + Denoise (Wiener2D)，保持 **QE BG** 关闭
4. 先做单谱拟合（**Fit Spectrum**），确认不扣背景时峰位仍清楚
5. 再做全局自动拟合（**Auto Fit Branches**，目前使用 BG + peaks 路径；**ZLP joint** 用于单谱 **Fit Spectrum**）
6. 逐 q 检查 branch assignment；对少数明显错误点，用 **Correct Auto** + target branch 下拉框做局部矫正，必要时用 **Undo Corr** 回退误操作
7. 用 **Save Session** 保存整个工作状态；如需单独交换分支点，再用 **Save Curated**
8. 可选：开启 QE BG 作为稳健性检查，确认主 branch 峰位没有被背景模型显著改变

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
| Drude-Lorentz | 金属等离激元（默认推荐）；A·E·Γ/((E²−E0²)²+E²Γ²)，代码里叫 `lorentz` |
| Sym. Lorentz | 以 E0 为中心的对称 Lorentz 峰（`lorentz_symmetric`） |
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
opts.do_bg_sub = false;    opts.bg_method = 'Auto';  % BG is optional robustness check
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

GUI 的单通道 ZLP joint 与 BiSb 的 P2b 多通道前向模型是不同方法。前者的 `apex_energy_meV` 表示损失函数峰顶，`omega_p` 表示模型 E0；它们也不等同于观测谱上的峰顶。数值验收只证明 GUI 和函数接口计算一致，不证明低 q 的模型参数具有物理唯一性。

---

## 测试

```matlab
cd('path/to/Dispersion-Toolbox')
startup
runtests('tests')
```

---

## Case Studies

Reusable toolbox code lives in `src/`; core tests are run with:

```matlab
startup; runtests('tests')
```

BiSb-specific analysis examples live in `case_studies/bisb2026/`. These scripts
may require ignored local data and outputs such as `20260120 BiSb/` and
`paper_results/`. Run their workflow tests separately:

```matlab
runtests('case_studies/bisb2026/tests')
```

---

## License

This project is provided for academic and research use within the group.
