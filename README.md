# q-EELS Plasmon Dispersion Toolbox

这是一个用于 **q-EELS / 动量分辨低损 EELS** 数据处理与色散分析的 MATLAB 工具包。README 只描述当前仓库源码中已经实现的通用能力；具体样品路径、实验窗口、特定分析参数和本地批处理脚本不作为通用说明的一部分。

工具包的核心目标是：把不同来源的 q-EELS 数据统一成标准 `qeData` 结构，然后通过 GUI 或脚本完成预处理、峰拟合、色散拟合和结果导出。

```text
supported input data
        │
        ▼
load_qe_dataset / load_raw_session
        │
        ▼
qeData: intensity(E,q), energy_meV, q_Ainv, metadata
        │
        ├── interactive_qe_browser       # GUI workflow
        │
        └── script API                   # batch / reproducible workflow
                │
                ▼
        qe_preprocess
                │
                ▼
        fit_loss_function / qe_auto_fit
                │
                ▼
        dispersion fitting / diagnostics / export
```

---

## 1. 功能概览

当前源码提供以下能力：

### 数据读取与标准化

- 读取已有二维 q-E 数据文件 `eq3D.mat`。
- 读取 raw Nion `.npy` 数据，并可使用同名 `.json` 元数据建立能量轴。
- 尝试读取兼容的 3D / 4D `.mat` 数据，包括可识别的 EELS 对象或数值数组。
- 将所有支持格式转换为统一的 `qeData` struct。
- 对 raw 数据生成 `eq3D_processed.mat` 缓存，加快后续加载。

### 预处理

由 `qe_preprocess.m` 实现：

- 可选 cosmic-ray / spike removal。
- 强度归一化：`Area` 或 `ZLP Peak`。
- 去噪：`Wiener2D`、`SavGol`、`BM3D`。
- 准弹性 / ZLP-tail 背景扣除。
- 可选 Lucy-Richardson 反卷积。
- 返回每个 q 通道的背景扣除诊断信息。

### 拟合与色散分析

- 单谱多峰拟合：`fit_loss_function.m`。
- 自动跨 q 通道拟合：`qe_auto_fit.m`。
- Seed peak propagation：`propagate_seed_peaks.m`。
- 峰型模型：Drude-Lorentz、Gaussian、Pseudo-Voigt、Damped Harmonic Oscillator。
- 色散模型：Quasi-2D Plasmon、Acoustic Linear、Optical Constant、Optical Quadratic。
- 线宽、振幅和品质因子诊断：`qe_gamma_dashboard.m`。
- q 相关 loss-map 诊断：`qe_loss_map.m`。

### GUI 与导出

`interactive_qe_browser.m` 提供交互式工作流：

- 数据加载与 q-E map 浏览。
- q/E 显示范围、trace mode、smoothing、normalization、background、deconvolution 等参数调节。
- 手动 peak picking、单谱拟合、自动拟合、branch reassignment、色散拟合。
- 操作历史保存/恢复。
- 数据和图像导出。

---

## 2. 安装与初始化

在 MATLAB 中进入仓库根目录并运行：

```matlab
cd('path/to/Dispersion-Toolbox')
startup
```

`startup.m` 会把 `src/` 和 `lib/` 加入 MATLAB path。

如果需要外部 Nion / EELS 工具箱路径，可复制模板：

```text
toolbox_config.m.template -> toolbox_config.m
```

并在本地文件中设置：

```matlab
cfg.nion_toolbox_root = 'path/to/external/toolbox';
```

`toolbox_config.m` 是本地配置文件，不应提交到 Git。

---

## 3. 快速开始

### 3.1 使用 GUI

```matlab
cd('path/to/Dispersion-Toolbox')
startup
interactive_qe_browser("path/to/data_folder_or_file")
```

支持直接传入：

```matlab
interactive_qe_browser("path/to/eq3D.mat")
interactive_qe_browser("path/to/raw_data.npy")
interactive_qe_browser("path/to/compatible_3d_or_4d_data.mat")
```

恢复之前保存的 GUI 操作历史：

```matlab
interactive_qe_browser("path/to/eq3D.mat", "path/to/op_history.mat")
```

### 3.2 使用脚本接口读取数据

```matlab
dataset = load_qe_dataset("path/to/data_folder");
qe = dataset.qe;

size(qe.intensity)   % [n_energy, n_q]
qe.energy_meV        % energy axis, meV
qe.q_Ainv            % momentum axis, Å^-1
```

若无法自动确定动量步长，显式传入 `dq_Ainv`：

```matlab
dataset = load_qe_dataset("path/to/eq3D.mat", dq_Ainv);
```

---

## 4. 支持的数据格式

主入口是：

```matlab
dataset = load_qe_dataset(path, dqOverride)
```

| 输入 | 源码行为 |
|---|---|
| 文件夹 | 优先查找 `eq3D_processed.mat`，其次 `eq3D.mat`，最后选择最大的 `.npy` 文件 |
| `eq3D.mat` | 读取变量 `a3` 和 `e`，可选读取 `q` |
| `.npy` | 使用内置 `read_npy.m` 读取 raw array，并查找同名 `.json` 元数据 |
| 3D / 4D `.mat` | 尝试识别 EELS 对象或大尺寸数值数组，并走 raw import pipeline |

### `eq3D.mat` 最小要求

```matlab
a3  % 2D intensity matrix, [energy × q]
e   % energy axis, usually in meV
q   % optional q-channel vector
```

### raw import pipeline

`load_raw_session.m` 对 raw 数据执行：

1. 读取数据。
2. 对 frame / spatial 维度求和，得到二维 `[E × q]` 矩阵。
3. 自动或手动裁切 q 维度。
4. 对不同 q 通道做 ZLP alignment。
5. 从元数据或数据本身建立能量轴。
6. 建立动量轴所需的 `dq_Ainv`。
7. 保存 `eq3D_processed.mat` 缓存。

---

## 5. 标准数据结构 `qeData`

`make_qe_struct.m` 生成的 `qeData` 是后续 GUI、预处理和拟合的共同数据结构。

| 字段 | 含义 |
|---|---|
| `intensity` | 2D intensity matrix, `[n_energy × n_q]` |
| `energy_meV` | 能量轴，单位 meV |
| `q_channel` | q 通道编号 |
| `q_zero_index` | q=0 对应的通道索引 |
| `q_Ainv` | 动量轴，单位 Å^-1 |
| `q_abs_Ainv` | `abs(q_Ainv)` |
| `dq_Ainv` | 动量步长，Å^-1 / pixel |
| `source_path` | 数据来源路径 |
| `source_kind` | 数据来源类型 |
| `label` | 显示标签 |
| `view_kind`, `stage_name`, `note` | GUI / processing metadata |

只要数据被转换为 `qeData`，后续处理就不再依赖原始文件格式。

---

## 6. GUI 工作流

启动：

```matlab
interactive_qe_browser("path/to/data_folder_or_file")
```

GUI 布局在 `src/ui/qe_browser_ui.m`，状态和回调在 `src/interactive_qe_browser.m`。

### 6.1 数据加载与显示

主要控件：

| 控件 | 功能 |
|---|---|
| `Load Data` | 加载文件或文件夹 |
| `Build Views` | 构建 comparison / normalized view |
| `↩ Undo` | 撤销到上一个操作历史 |
| `Reset To` | 重置到指定 stage（raw / eq3d / crop / align / bin） |
| `Physical / Normalized` | 切换显示视图 |
| `Ovr dq` | 覆盖动量步长 dq（Å⁻¹/pixel） |
| `q start`, `q end`, `q step` | q 范围与采样 |
| `E min`, `E max` | 能量显示和拟合范围 |
| `Trace` | `display`、`background-subtracted`、`curvature` |
| `Offset` | 堆叠谱偏移量 |
| `Smoothing` | `gaussian`、`off`、`sgolay`；可选 `Width`、`SG ord`、`SG len` 子参数 |
| `Auto Y`, `Y min`, `Y max`, `Y scale` | 单谱 Y 轴范围和比例设置 |

### 6.2 预处理控件

| 控件 | 对应功能 |
|---|---|
| `Normalize` | 启用强度归一化 |
| `Area / ZLP Peak` | 选择归一化方式 |
| `Denoise` | 启用去噪 |
| `Wiener2D / SavGol / BM3D` | 选择去噪方法 |
| `Noise σ` | 去噪噪声参数（0 = 自动估计） |
| `QE BG` | 启用准弹性背景扣除 |
| `BG Method` | 选择背景模型或 `Auto` |
| `BG W1` (lo / hi) | 第一背景拟合窗口上下界 |
| `Dual` + `BG W2` (lo / hi) | 启用第二背景拟合窗口及上下界 |
| `Iter` | 启用迭代重拟合 |
| `Deconv` | 启用 Lucy-Richardson 反卷积 |
| `LR iter` | Lucy-Richardson 迭代次数 |

### 6.3 Peak picking 与拟合

| 控件 | 功能 |
|---|---|
| `Pick Peaks` | 手动选取色散点 |
| `Undo Pt` / `Clear Pts` | 撤销或清除已选取的手动点 |
| `Save Pts / Load Pts` | 保存或加载手动点 |
| `Fit Model` | 对手动点拟合色散模型 |
| `Fit Spectrum` | 对当前 q 通道做单谱拟合 |
| `Accept Fit` | 接受当前单谱拟合结果 |
| `Auto Fit ω(q)` | 跨 q 通道自动拟合峰位 |
| `Guesses` | 文本框，手动输入 seed 峰位（如 `800,2500,3200`） |
| `Pick Guesses` | 在谱上交互式选择 seed 峰位 |
| `Max Shift` | Seed propagation 相邻 q 通道最大偏移量 |
| `Reassign Pts` | 重新分配点到分支 |
| `Split 3` | 分离 peak branches |

### 6.4 色散与诊断

| 控件 | 功能 |
|---|---|
| `Peak Model` | Lorentz、Gaussian、Voigt、Damped HO |
| `Disp Model` | Quasi-2D Plasmon、Acoustic Linear、Optical Constant、Optical Quadratic |
| `Fit Dispersion` | 拟合选定色散模型 |
| `Γ & Amp` | 显示 linewidth / amplitude / quality-factor dashboard |
| `I_kin Map` | 显示 q-dependent loss-map 诊断 |
| `Export Disp.` | 导出色散图 |
| `Export Data` | 导出预处理后的 q-E 数据（含 axes、预处理参数、背景诊断） |

### 6.5 保存和导出

| 输出 | 控件 | 内容 |
|---|---|---|
| 操作历史 | History panel `Save / Load` | `op_history.mat`，保存参数和操作快照 |
| 预处理数据 | `Export Data` | intensity、axes、原始数据、预处理参数、背景诊断 |
| 图像 | `Export` | q-E map、comparison map、single spectrum、dispersion 等 |
| 手动点 | `Save Pts` | 手动点或已接受拟合点 |

---

## 7. 脚本接口示例

GUI 适合交互检查；可复现实验流程建议用脚本调用核心函数。

### 7.1 读取和预处理

```matlab
cd('path/to/Dispersion-Toolbox')
startup

dataset = load_qe_dataset("path/to/data_folder", dq_Ainv);
qe = dataset.qe;

opts = struct();
opts.do_despike = false;
opts.do_normalize = true;
opts.norm_method = 'ZLP Peak';
opts.norm_min = -50;
opts.norm_max = 50;
opts.do_denoise = true;
opts.denoise_method = 'Wiener2D';
opts.denoise_sigma = 0;
opts.do_bg_sub = true;
opts.bg_method = 'Auto';
opts.bg_candidate_methods = {'Power','ExpPoly3','Pearson','PearsonVII','Power2','Exp2'};
opts.bg_win_lo = [50, 300];
opts.bg_win_hi = [];
opts.bg_iterative = false;
opts.do_deconv = false;

[qe_pp, bg_diag] = qe_preprocess(qe, opts);
```

### 7.2 自动峰拟合

```matlab
fit_opts = struct();
fit_opts.E_min = 300;
fit_opts.E_max = 3000;
fit_opts.q_start = min(qe_pp.q_Ainv);
fit_opts.q_end = max(qe_pp.q_Ainv);
fit_opts.prominence = 0.001;
fit_opts.smooth_width = 25;
fit_opts.max_peaks = 2;
fit_opts.peak_model = 'lorentz';
fit_opts.pre_subtracted = opts.do_bg_sub;
fit_opts.guesses = [];
fit_opts.seed_idx = [];
fit_opts.max_shift = 80;
fit_opts.R2_threshold = 0.3;
fit_opts.verbose = true;
fit_opts.energy_mask = qe_pp.energy_meV >= fit_opts.E_min & qe_pp.energy_meV <= fit_opts.E_max;
fit_opts.energy_axis = qe_pp.energy_meV(fit_opts.energy_mask);

raw_opts = opts;
raw_opts.do_normalize = false;
qe_raw = qe_preprocess(qe, raw_opts);

results = qe_auto_fit(qe_pp, qe_raw, fit_opts);
branches = results.branches;
```

`qe_auto_fit.m` 返回：

```text
results.all_peaks    % 经 R² threshold 过滤后的峰矩阵
results.branches     % cell array of per-branch subsets
results.fit_details  % per-channel fit result structs
results.n_success    % 成功拟合的 q 通道数
results.used_seed    % 是否使用了 seed propagation
```

`results.all_peaks` 的列定义取决于工作模式：

**Blind mode**（无 seed guesses）— Nx12 矩阵：

```text
[q, omega_p, gamma, R2, amplitude,
 E_ci_lo, E_ci_hi, G_ci_lo, G_ci_hi, A_ci_lo, A_ci_hi,
 raw_peak_height]
```

**Seed mode**（有 seed guesses）— 列结构由 `propagate_seed_peaks` 决定，第 6 列为 `branch_id`，分支提取使用前 5 列 `[q, E, Γ, R², A]`。

### 7.3 色散拟合和诊断图

```matlab
fit_result = fit_dispersion_generic( ...
    branches{1}(:,1), branches{1}(:,2), ...
    'model', 'quasi2d_plasmon');

fig1 = qe_gamma_dashboard(results.all_peaks, branches);
```

如果要使用 q-dependent loss-map / intensity diagnostics，可参考：

```matlab
phys = qe_physics_extract(branches, struct('branch_index', 1));
fig2 = qe_loss_map(qe_pp, phys, struct('mode', 'experimental'));
```

示例数值仅用于展示 API。实际能量窗口、q 窗口、背景窗口和拟合阈值应由当前数据决定。

---

## 8. 背景扣除

`qe_preprocess.m` 的背景扣除是准弹性 / ZLP-tail 的 **subtraction** 步骤，不是 normalization。支持模型来自源码中的 `local_background_model_candidates` 和相关 dispatch 逻辑：

| Method | 说明 |
|---|---|
| `Power` | 幂律背景 |
| `Power2` | 幂律 + 常数项 |
| `ExpPoly3` | log-intensity polynomial model |
| `Pearson` | log-log quadratic empirical model |
| `PearsonVII` | Pearson-VII-style tail model |
| `Exp2` | 双指数背景 |
| `DualWindow` | 双窗口 side-band anchoring |
| `Auto` | 对候选模型打分并为每个 q 通道选择模型 |

常用选项：

```matlab
opts.do_bg_sub = true;
opts.bg_method = 'Auto';
opts.bg_candidate_methods = {'Power','Exp2','ExpPoly3','Pearson','PearsonVII'};
opts.bg_win_lo = [E1, E2];
opts.bg_win_hi = [];        % or [E3, E4]
opts.bg_iterative = false;
```

可返回的 `bg_diag` 字段包括：

| 字段 | 含义 |
|---|---|
| `selected_method` | 实际使用的背景模型 |
| `rsquare`, `rmse` | 背景窗口拟合质量 |
| `bg_curve` | 拟合背景曲线 |
| `residuals` | 拟合残差 |
| `h_param`, `snr` | 背景/信号诊断量 |
| `candidate_methods`, `candidate_scores`, `candidate_details` | `Auto` 模型选择信息 |
| `linear_rmse` | 线性残差诊断 |
| `neg_fraction`, `neg_area_fraction`, `neg_peak_fraction` | 过扣除诊断 |
| `bg_fraction` | 扣除背景相对强度 |

推荐做法：先在多个 q 通道上检查背景曲线与扣除后谱，再批量拟合峰位。背景窗口不应覆盖待分析的真实峰。

---

## 9. 动量标定

动量轴在 `make_qe_struct.m` 中按下式构造：

```matlab
q_Ainv = ((1:n_q) - q_zero_index) * dq_Ainv;
```

其中：

- `q_zero_index` 由 ZLP 附近强度自动估计，或由输入 metadata 指定。
- `dq_Ainv` 是每个 q pixel 对应的 Å^-1 步长。

如果自动推断不适用于当前数据，应显式传入：

```matlab
dataset = load_qe_dataset("path/to/data", dq_Ainv);
```

或在 GUI 中用 **Ovr dq** 覆盖。

---

## 10. 源码结构

```text
├── README.md
├── startup.m
├── toolbox_config.m.template
│
├── src/
│   ├── interactive_qe_browser.m
│   ├── qe_preprocess.m
│   ├── qe_auto_fit.m
│   ├── qe_physics_extract.m
│   ├── qe_loss_map.m
│   ├── qe_gamma_dashboard.m
│   ├── qe_prepare_stacked_spectra.m
│   ├── build_comparison_qe.m
│   ├── qe_recommend_lowq_background_configs.m
│   ├── qe_summarize_lowq_background_eval.m
│   ├── analyze_dispersion.m             # 多 session 自动化色散分析脚本
│   ├── qe_plot_helpers.m                # 可复用绘图工具类（branch scatter / fit overlay / map display）
│   │
│   ├── io/
│   │   ├── load_qe_dataset.m
│   │   ├── load_raw_session.m
│   │   ├── load_npy_session.m
│   │   ├── read_npy.m
│   │   ├── make_qe_struct.m
│   │   └── resolve_q_crop_bounds.m
│   │
│   ├── fitting/
│   │   ├── fit_loss_function.m
│   │   ├── fit_quasi2d_plasmon.m
│   │   ├── fit_dispersion_generic.m
│   │   ├── dispersion_models.m
│   │   ├── peak_models.m
│   │   ├── propagate_seed_peaks.m
│   │   └── measure_peak_height.m
│   │
│   └── ui/
│       └── qe_browser_ui.m
│
├── scripts/      # example / project-specific batch scripts; review before reuse
├── tests/        # MATLAB unit tests
├── docs/         # 开题报告、周报等项目文档
├── figures/      # 分析输出图像和结果
├── references/   # 参考文献资料
└── lib/BM3D/     # bundled BM3D implementation
```

`src/` 是通用工具包主体。`scripts/` 中的文件更适合作为批处理示例或项目内实验脚本，复用前应检查路径、窗口和参数。

---

## 11. 环境要求

建议环境：

- MATLAB R2024b 或更新版本。
- Optimization Toolbox：非线性拟合。
- Curve Fitting Toolbox：背景模型拟合。
- Signal Processing Toolbox：平滑、滤波和去噪。
- Image Processing Toolbox：Lucy-Richardson 反卷积相关功能。

`lib/BM3D/` 中包含 BM3D 去噪实现及相关 mex 文件。

---

## 12. 测试

在 MATLAB 中运行：

```matlab
cd('path/to/Dispersion-Toolbox')
startup
runtests('tests')
```

当前测试文件覆盖：

- `measure_peak_height`
- pre-subtracted auto fit
- background subtraction
- physics extraction
- stacked spectra preparation
- background configuration recommendation
- low-q background summary utilities
- real-data background selection smoke test（需要对应本地数据时才适用）

---

## 13. 使用建议

对新数据，建议按以下顺序检查：

1. 确认 raw / processed 数据能正确加载为 `qeData`。
2. 检查能量轴、q 轴和 ZLP alignment。
3. 先在 GUI 中检查代表性 q 通道。
4. 调整 normalization、denoising 和 background 参数。
5. 检查 `bg_diag`，避免明显过扣除。
6. 先做单谱拟合，再做跨 q 自动拟合。
7. 检查 branch assignment 后再拟合色散。
8. 保存 `op_history.mat` 或脚本参数，保证结果可复现。

---

## License

This project is provided for academic and research use within the group.
