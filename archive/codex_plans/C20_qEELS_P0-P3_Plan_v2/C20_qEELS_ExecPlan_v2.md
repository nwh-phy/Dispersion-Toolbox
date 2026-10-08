# C20 q-EELS：基于最新版 Dispersion-Toolbox 的双分量分析执行计划

版本：2.0｜编制日期：2026-09-11  
源码基线：`nwh-phy/Dispersion-Toolbox@18af79b66cd27365fd00bcf3de26658ddd7e7d53`  
该提交时间：2026-09-11 05:53:08 UTC；提交信息：`Sync q calibration, branch tracking, and analysis exports`。  
对比基线：`bbf24965b6c1ea177541bf382e0b28fd7cdafdf7`；新基线领先 24 个提交。[S01]

**文件性质：交给本地 Codex 的执行规范，不是已经运行的新分析结果。** 本轮已在线检查版本、关键源码、测试定义和结果索引；没有读取用户本地实验数组，没有运行 MATLAB，没有修改或推送 GitHub 仓库。结果索引中的统计仅为仓库记载，不视为本轮重算或验证。

本文件取代上一份 v1 计划作为新一轮执行入口。v1 和旧审阅保留作历史记录，不再并行执行相互冲突的默认配置。客户端按用户要求选用相应模型与推理强度；本计划不编造模型 ID、CLI 参数或 API 支持情况。

## 0. 目标与不可改变的科学口径

先在原 B1 能区提取两个重叠的正谱学分量 P1/P2，报告各自的峰位、宽度、有限能窗积分及不确定度，再研究其色散和可能的激发归属。P1/P2 不是旧 B1/B2 的改名，也不预先等同于两种 plasmon。用户观察到的是两个正峰合并，而不是凹陷型反共振；Fano 不作为新主线的预设机制。

C20 是本科毕业设计样品；用户最新明确更正是 Bi-rich，Bi:Sb 约 70:30 原子比。不得采用更早文档中的反向比例，也不得把旧目录中的 `Bi` 自动理解为纯 Bi。研究对象默认是完整的 MoS₂–metal–MoS₂ 封装体系。文件名 `1film/2film` 不自动构成受控单层/双层证据；不得默认启用厚度因子 1/2 的物理约束。

**成功标准不是所有 q 上都有两条光滑线，而是分清哪些分箱确实支持两个分量，以及每个参数能被确定到什么程度。** 可以得到峰位可辨识而线宽不可辨识的结果；也可以得到部分区间未分辨的结果。不能通过压窄线宽上限、强制分离、插值补点或筛掉不顺眼的点制造双分支。

### 0.1 当前仓库已经有的内容，不能再当作从零开发任务

- 论文与数据登记入口：`AGENTS.md`、`src/thesis/thesis_sessions.m`、`thesis_config.m`。
- 统一 dq 推断：`src/io/infer_qe_dq_Ainv.m`；配有标定回归测试。
- GUI history/curated points 与三组数据复现：`run_590_gui_history_area_analysis.m`。
- B1 分箱和双分量候选提取：`src/b1_double_peak_binning_extract.m`。
- 独立、传播、分窗口、ridge-guided candidate-path 等已有追踪模式。
- 点级诊断：候选、路径、失败、repair、分箱来源及 `b1_peak_evidence_audit_classify.m`。
- Fano、峰顶位置、峰顶区间、峰质量诊断和 FWHM 辅助模块。
- 通用 `tests/` 与 `case_studies/bisb2026/tests/` 两层测试。[S02–S12]

执行策略是**复用已有计算与诊断，补足输入一致性、分箱统计、联合分解与验证**，不全面重写或转语言，不开发新 GUI，不再用“不断加 tracking variant 直到 score≤300”作为科研目标。

### 0.2 论文、历史索引与权限

先读取仓库 `AGENTS.md`，再读取本地正式毕业论文相关段落，按“已完成观察—候选解释—尚待验证”登记。正式基线优先 2026-06-05 提交 PDF 及匹配归档源码；根目录 `thesis.pdf` 不是提交版。论文是起点，但用户后续更正优先。不能把旧 AI 报告当作当前成果的权威摘要。[S02]

原始数据只读；不覆盖论文、旧输出、当前 curated 入口或废弃目录。读取历史脚本不代表获准执行其中所有写文件、更新索引或物理拟合动作。所有新输出必须写入独立 run 目录。不得 reset --hard、clean、force push、自动上传实验数据、修改全局模型设置或擅自升级整个 MATLAB 环境。

## 1. 最新源码对应的重点决策

| 主题 | 已核对的最新版行为 | 本轮执行决定 |
|---|---|---|
| dq | `infer_qe_dq_Ainv` 与 `thesis_sessions` 均为 10w=0.0005、20w=0.00025 Å⁻¹/pixel；测试有相同断言 | 以此作为当前项目配置。核对加载后轴与旧 MAT/CSV，不再把旧 0.005/0.0025 当成同等有效默认，也不声称本轮独立重做了仪器标定 |
| q 阈值迁移 | 双峰 wrapper/core 仍有范围 ±0.15、q_skip=0.005、low_q_no_bin=0.05 等；history 主入口已用 ±0.015 | 把阈值迁移和保存结果轴检查列为 P1。不能对所有数字机械除 10 |
| 当前主线 | README 明确不预先扣背景提取主峰，背景扣除用于诊断 | 新主线不预扣；拟合中仍建模背景，不把“不扣背景”误解为“背景恒为零” |
| 旧双峰输入 | wrapper 读 `analysis_results.mat/output.qe_pp`、`branch1_points.csv`，优先 `wideq030` 目录 | 历史复现继续读原输入；新分析显式选定最小处理输入，禁止自动借用不明标定的旧 wide-q 文件 |
| 旧 raw 命名 | wrapper 的 raw 重建只关掉归一化，可能仍去噪/扣背景；失败还可回退 qe_pp | 新主线缺 raw 就标明缺失，不能静默把 qe_pp 冒充原始计数 |
| 相邻 q 组合 | core 按组求均值、保存成员、两侧分别处理；有噪声阈值和低 q 保护 | 保留这些功能；补充显式固定 N 分箱、sum/mean/方差/边界/有效成员。不要重复造一个不兼容的分箱流程 |
| 双分量模式 | `independent_double_peak` 是同窗两分量拟合；`windowed_branch_tracking` 对上下窗口各做一次单峰拟合 | 后者仅作候选追踪/初值诊断，不能当作完整同谱联合分解或直接解释两分量面积 |
| 峰形 | `lorentz` 仍为 Drude–Lorentz，最新版已有 `fano` | 不重定义旧名字。保留 DL 对照，新增明确命名的标准对称 Lorentzian 作为首轮有效线形 |
| 默认行为 | 双峰 core 默认 fano；wrapper `runFits=true` 可接物理拟合 | 所有新入口显式指定峰形，提取阶段禁物理拟合，不调用 wrapper 默认全流程 |
| 历史结果 | 500 meV 下界和其派生参考、±q 平均图已标废弃；已有 v14/v15 和证据审计 | 不恢复废弃参考；主能窗下界从 300 meV 开始，正负 q 分开。历史评分及 evidence 标签只是诊断 |

来源与关键代码位置见同包 `SOURCE_ALIGNMENT.md`，所有源码链接固定到本次 SHA。[S03–S15]

## 2. 原始输入及已有中间结果：按用途区分，先自动发现

### 2.1 三组数据的准确登记

下列是最新版 `thesis_sessions.m` 和 history 入口的登记值；**路径是相对当前工程根目录的线索，不代表这些未跟踪文件已经随 GitHub 下载**。[S05, S08]

| session_id | 当前目录线索（`20260120 BiSb/` 下） | 当前 dq / Å⁻¹·pixel⁻¹ | 本轮用途 |
|---|---|---:|---|
| `590_PL2_10w` | `590 PL2 10w 0.004 10sx300/` | 0.0005 | 方法开发与首轮试点 |
| `n0_PL2_10w_repeat` | `n0 pl2 10w 0.004 10s x300/` | 0.0005 | 新配置冻结后的重复一致性检查 |
| `no_PL2_20w_2film` | `no pl2 20w 0.004 10sx300 2film/` | 0.00025 | 更细采样与区域比较 |

重复数据过去已进入多轮共同调参，不能宣称它从未被观察过或构成严格的盲测。新一轮冻结设置后再跑它，称为冻结方案下的一致性验证。

### 2.2 必需与可选输入

**A. 新分量提取的最小可启动包。** 至少一组完整数值谱、能量轴及 signed-q 标定来源、当前仓库和可以识别处理层级的记录。优先 `.npy` + 同名 `.json` 或原始多帧 MAT；也接受经核查的 `eq3D.mat`（`a3/e/q`），但必须说明它已否经过累计、裁剪、对齐、去噪、归一化等。不是强制所有原始帧齐备才能开始；只得到二维谱时先做二维有效谱参数，限制误差与寿命解释。

**B. 必须同时查找的采集与标定信息。** 实际 shape、dtype、轴顺序和单位；E 标定/零点；q 中心、步长、原始索引、方向和裁剪位置；帧数、曝光、束流、坏点/饱和、暗场/增益、探测器是否计数或 ADU。无法获得的信息保持 unknown，不从目录名自动填成已验证采集参数。json 中的单位必须解析，而不是只凭字符串含 eV 就把 meV 再乘 1000。

**C. 旧结果复现输入。** 对应的 `eq3D.mat`、`eq3D_processed.mat` 与 `import_provenance`；`op_history_260506.mat`、其他 `op_history*.mat`；Save Pts 的 curated points、corrections/snapshot；GUI 导出的 `*_preprocessed*.mat`。优先读取文件内部字段，不按文件名、mtime 或目录标签判定版本。

**D. 最新双峰脚本实际需要的中间输入。** 当前 wrapper 依赖旧 `paper_results/<session_tag>/analysis_results.mat` 中的 `output.qe_pp`、`output.dataset.qe`、`output.preprocess_opts`、`output.snap`，以及同目录 `branch1_points.csv`。最新 history 脚本将 590 的 `op_history_260506.mat` 作为另两组的模板，且对 20w 有高-q 局部再拟合与可能的点替换记录；这种跨组模板使用是已存在规则，不要偷偷换成每组目录下另一个过时 history。[S08]

需要发现的旧目录：

```text
paper_results/590_gui_history_area_260506/
paper_results/n0_PL2_10w_gui_history_area_260506/
paper_results/no_PL2_20w_2film_gui_history_area_260506_highq_refined/
各自的 *_wideq030/ 目录（仅作为需检查的已有输入，不自动优先）
paper_results/00_CURRENT_B1_DOUBLE_PEAK_260509/
paper_results/b1_lorentz_tracking_peak_evidence_audit_260511/
```

读取双峰旧结果时继续追踪：`b1_double_peak_binning_results.mat`、lower/upper/combined points、binning_map、noise_profile、fit_failures、repair/exclusion、candidate/path-selection CSV。若缺失它们，新提取仍可从 A 启动，但不能声称已复现历史双峰。

**E. 仪器响应与物理解释输入。** 同采集条件的真空/参考 ZLP、能量分辨、动量接受核；区域对应图、MoS₂ 封装层数/取向/厚度资料、HAADF/NBED/EDS 等。缺少仪器响应时先报观测/有效 FWHM，不反卷积出本征寿命；缺少结构对照不阻止分峰，但限制归属。

**F. 正式论文。** 旧本地根目录线索：`C:/Users/HP/Desktop/vibecoding/本科毕业论文/`；提交 PDF 相对路径：`终稿提交包_221240048聂文辉_20260605/221240048聂文辉终稿.pdf`，以及匹配的归档源文件。重点读方法、结果、总结和附录。路径迁移到 Mac/其他盘符后由用户授权根目录解析，不全盘扫描。

### 2.3 台账和缺失处理

本地执行第一步自动生成 `input_manifest.resolved.yaml`；随包模板只是台账，不是现有 MATLAB 已支持的配置文件。不要要求用户提前填完它。对实际采用的文件保存 SHA256、源文件/实际读取文件、轴与单位、处理层级、字段位置、数据与输出映射；metadata 不足时登记影响和可继续阶段。

数据层使用：L0=原生帧；L1=累计/已校准物理 q–E；L2=归一化/去噪等历史处理谱；L3=由参考谱重建或仅供显示的数组。L1 不是自动等于独立 Poisson 计数；L2/L3 不得重新命名为 L0。

## 3. 新一轮实施路径与文件职责

保留 MATLAB 为主，不为统一语言移植。`src/` 放可复用运算，`case_studies/bisb2026/scripts/` 放本项目配置与薄入口。先检查是否已有等效 helper；满足要求就复用，不重复建立平行框架。

拟新增的最小接口（**这些名称是实施目标，当前基线尚不存在，不得假装已可调用**）：

| 拟新增/扩展接口 | 职责及约束 |
|---|---|
| `b1_component_config_v2.m` | session/输入模式、单位、能窗、分箱、拟合、验证、输出根统一配置；源于 registry，不再复制会话常数 |
| `run_b1_component_pilot_v2.m` | headless 入口；支持 validate_inputs_only / pilot / full 和独立 output_root；不自动调用物理模型、旧索引更新或复制当前结果 |
| `qe_prepare_count_bins.m` 或已有 helper 的等效扩展 | signed-q 连续有效通道的固定分箱、sum/mean/方差、边界/成员/掩码；可作为现有 extractor 的显式 extraction-units 输入 |
| `qe_compare_component_models.m` 或 fitter 的向后兼容扩展 | 同窗 n=1/n=2、多初值、背景与峰联合、噪声/响应、统一结果表；不夹杂绘图 |
| `qe_validate_component_models.m` | 单峰误拆、双峰恢复、配置扰动、bootstrap/profile 与分类输出 |

现有 `b1_double_peak_binning_extract` 的成员来源表、失败表、fit_details、seed/ridge 候选及 overlay 思路继续利用。必要时仅抽出当前任务涉及的小函数，先写行为保持测试，再重构。每个提交只改一类东西：输入/标定、分箱、模型、验证或导出；不要一次同时改所有规则。

## 4. P0：版本、论文和已有流程冻结

1. 记录 `git status --short`、HEAD、branch、diff、实际 MATLAB 版本和 toolboxes；不覆盖用户未提交修改。如果本地 HEAD 已超过本计划基线，先列出相关差异，再决定配置是否仍适用；不要强制切回旧 SHA。
2. 阅读 `AGENTS.md`、正式论文、case README/RESULTS_INDEX 和关键脚本。将现有工作登记为“已有功能/已有历史输出/此次尚未执行”，不重复报告成新增成果。
3. 选择一个 590 代表输出建立历史只读基准，记录数据、快照、峰形、q 标定、人工修改和输出对应关系。用已有输出先建立基准，不为了逐字节重现全部旧报告阻塞新分峰。
4. 运行并保存通用测试；case 测试先区分纯函数、静态源码断言、需本地数据/会写结果的测试。无数据或无 MATLAB 的测试标明 skipped/not-run，不把它们算 passed。先读测试，避免测试自动覆盖旧 `paper_results`。
5. `which run_thesis_pipeline -all` 等检查同名入口解析；只使用明确路径。本轮不默认执行论文全流程。

已有安全检查入口（在已确认可写临时输出的 MATLAB 环境中）：

```matlab
startup
which load_qe_dataset -all
which b1_double_peak_binning_extract -all
r_core = runtests('tests');
% case 测试按需分组，在确认本地数据依赖与写出行为后执行。
```

P0 输出：`environment.txt`、`repo_state.json`、`baseline_inventory.md`、`test_baseline/`、初版输入台账。

## 5. P1：统一 q 口径与实际输入，防止新版读取配旧版阈值

### 5.1 当前配置与旧结果迁移

当前项目采用 10w=0.0005、20w=0.00025 Å⁻¹/pixel。程序一致性测试不能代替独立仪器标定证据；但没有相反的实际证据时，不能继续把旧 0.005/0.0025 作为新配置候选反复阻塞分析。[S04, S05, S06]

对每个实际加载的 L1/L2 数组检查：signed q 是否有序；`median(diff(q))` 与登记值是否一致；q_zero、source_channel、裁剪范围是否可追溯；旧 CSV/MAT 的 q 是否指向同一批原始列。读 `analysis_results.mat` 不会自动经过新推断函数，必须单独核查。

若旧数组 q 相差 10 倍，只能在证明同一原始通道映射后派生一份新轴，不原位覆写；必须同时审查旧 q-dependent 筛选、分箱、平滑和参考带。**旧参数点简单除以 10，不等于分析已经按新口径重做。**

### 5.2 阈值迁移清单

逐项导出当前值、来源、对应原始列、影响、新值及理由：qRange、q_skip、lowQNoBin、highQForceBin、fitDenoiseQStart/End、各 trackingHighQ、trend anchor/smallQ/plateau、reference q bands、max_q_gap、plot/export q limits、evidence classifier 的 small-q 切换点。

已定位风险：双峰 wrapper/core 的 `q_skip=0.005` 比 thesis/history 的 0.0005 大十倍；若分析 |q|≤0.015，而 `low_q_no_bin=0.05` 原样保留，则所有候选点都受低-q 保护，不发生真正分箱。[S07, S09]

新试点显式采用 q_range=[−0.015,0.015] 与 q_skip=0.0005，裁到实际可用范围；这是当前项目分析基线，不等于仪器分辨率。若真实中心束掩码需要更宽排除区，按实测掩码另立具名敏感性版本，不能按拟合是否成功选掩码。固定 N 分箱测试关闭旧 adaptive 低-q 保护，但仍遵守中心束/无效通道掩码。自适应分箱在固定方案通过后单独验证。

**不要全仓库搜索替换除 10。** ±q、像素数、能量阈值、dimensionless penalty、纯合成测试坐标不能混改。先生成 `q_rule_migration.csv`。旧历史数值保存在 legacy 配置；新主线只读统一配置。

### 5.3 输入与缓存

显式 `eq3D.mat` 当前已不会被邻近缓存无条件替代；显式 raw 裁剪也增加了 provenance 校验。保留并运行 `test_load_qe_dataset_cache`，不要按旧审阅重复“修复”已变更行为。[S07]

仍需针对本轮增加显式 raw/no-cache、只读 source/独立 cache 输出策略。避免读取 raw 后自动在原目录重写 `eq3D_processed.mat`。缓存身份至少纳入文件内容/元数据、轴解释、裁剪、对齐算法版本/参数。未知旧缓存只可当已处理输入，不能冒充完整可复现 raw 导入。

P1 验收：三个 registry、实际加载轴及新 runtime 配置一致；在真实 dq 的合成网格上证明中心掩码和 N>1 分箱确实产生预期成员；不确定数据被标为 unresolved/legacy，不混入主线。

## 6. P2：建立两套互不混淆的数据路径

### 6.1 历史复现路径

保留实际 history 的 Area、去噪、background、Fano 或 DL 选择、20w refinement 和 curated edits。只有明确的新 run 目录才执行复现；不能调用会写原结果或替换索引的默认脚本。GUI 展示平滑、AsLS trace 和参考重建视图只作诊断，不能自动变成定量输入。

`run_590_gui_history_area_analysis` 会强制 Area 归一化，并读取 590 快照给其他 session 作模板；`src/thesis/thesis_config` 又是另一套含背景扣除、低支 500 meV 筛选的 baseline。二者不是可互换的“论文最终分析”。先与提交版证据链对照，当前分峰不能自动调用后一配置。[S08, S15]

### 6.2 新定量路径

从核对后的 L0/L1 出发，保留完整能量范围供校准/ZLP诊断，再选拟合窗口。首轮不逐-q归一化、不做 Wiener/BM3D/PCA/NMF 跨 q 去噪、不做 Lucy–Richardson 反卷积、不预先扣背景、不用 reference-reconstructed 中心谱。

允许必要的增益/暗场、坏点掩码、裁剪和经验证的能量对齐；每步记账。不能移动各 q 的非弹性峰令它们重合。旧 raw importer 的帧累计及对齐顺序要检查；有帧时检查时间漂移和束照演化，无帧时不要伪造重复维。

原生单位沿用 MATLAB 代码的 **meV、Å⁻¹**。v1 提过内部切换 eV，本版取消这一非必要改动；对外展示 eV 时单独换算并带字段单位。输入 json 或 MAT 单位未知不能静默猜。

计数/ADU/已校正浮点数据的噪声模型分别记录。raw 缺失时，先在 L2 上复现/探索；文件名明确 `processed_domain_exploratory`，不输出绝对谱重或本征线宽结论。

## 7. P3：相邻 q 分箱与首轮三个代表位置

### 7.1 分箱规则

先实现并对照 N=1、3、5 的固定非重叠分箱；N=7及更大值以后按需要测试，不直接继承历史 bin7/bin11 或 sg71/sg111 的最优标签。分组前按 signed q、中心束掩码、坏点/缺测间隙和采集区域拆段；不得跨零点或跨缺测原生列。尾部不足 N 保留 partial bin 并报告实际大小，不补零。

对均匀采样的等权组同时保存：

\[
C_B(E)=\sum_{i\in B} C_i(E),\quad
\bar C_B(E)=C_B(E)/N,\quad
V_B(E)=\sum_i V_i(E)
\]

最后一式只在相互独立时使用；有相关项则 V_B=∑V_i+2∑Cov_ij，或传播 Σ_out=TΣ_inTᵀ/采用端到端帧重采样。均值的方差为 V_B/N²。sum 和 mean 的预测值与权重同步缩放后应给出等价峰形估计；幅度因子按定义变换，不能混用。

当前 core 的 mean 不是天然错误；要补的是统计定义和输入层标记。禁止把逐谱 Area/ZLP Peak 归一化后均值当成原始计数叠加；禁止将 N 个探测器 q 像素合并记成 N 倍曝光时间。

每个 bin 保存 bin_id、所有原生成员、signed q centers/edges、q_left/right/center、总宽度、每个能量像素有效成员、sum/mean、方差/协方差来源、掩码、处理层、partial 标志、归一化和插值记录。已有 `source_q_count/source_q_index/source_q_Ainv` 保留为兼容字段。

`omitnan` 不能让各能量点偷偷来自不同数量的谱：缺失点用明确有效成员/权重计入前向模型，或采用共同有效能量掩码。q_center 为几何或已知采集权重中心，不用待拟合峰高倒推出模式相关的“实测 q”。不重新调用自动找 q=0 来改变标定。

### 7.2 不同采样组的可比性

10w 的 N=3、5 对应像素边界总宽度 0.0015、0.0025 Å⁻¹；20w 分别用 N=6、10才匹配。首末中心距离是 (N−1)dq，与总宽度 Ndq 区分。20w 同样可以跑 N=1/3/5 诊断，但跨组比较必须注明不同宽度。

这些宽度不是仪器动量分辨率。横向条表示接受区间，q 零点/刻度不确定度另列；不得把接受宽度除 √N 当成统计 q 误差。替代分箱起点移动一个原生通道仅作敏感性检查；不同 N 和重叠替代网格来自同一原数据，不当独立重复共同计入似然。

### 7.3 首轮实际交付

按原始谱质量而非双峰拟合成功率，选择低/中/高 |q| 三个代表位置，保存其原始成员和分箱谱。随后扩大至 12–20 个覆盖双肩清楚/严重重叠、正负 q 和低计数的代表 bin。

先完成无物理约束的同能窗单/双分量小规模试拟合。若数据不足以拟合，仍应交付正确的输入、分箱和明确原因，不通过放宽所有规则把失败改成成功。

## 8. P4：同谱联合提取与可辨识性检验

### 8.1 模型与能窗

主模型为“背景 + 两个正的标准对称 Lorentzian 分量”，同时保留一个分量及当前 Drude–Lorentz 的 n=1/n=2 对照。标准 Lorentzian 单位面积模板：

\[
L(E;E_j,\Gamma_j)=\frac{1}{\pi}\frac{\Gamma_j/2}{(E-E_j)^2+(\Gamma_j/2)^2}.
\]

采用新模型名 `lorentz_symmetric`（实施时注册并测试），绝不改写旧 `lorentz` 的公式/参数含义。DL 的 E0/Gamma/A 与标准 Lorentzian 的峰中心/FWHM/面积不是一一同名对应。[S11]

使用 **300–1800 meV** 作为首轮窗口，300–2000、300–2100 meV 作具名敏感性版本；上界扩大需考虑旧其他能区的尾部。不得将已废弃 500 meV 下界或其派生追踪参考用于当前主结果。若实测能量覆盖不够则报 actual window；若低分量贴近300边缘，标为 edge-limited，必要时在可信背景模型下增加更低下界的独立探索，不截掉它再宣称分量消失。

### 8.2 不预扣背景，但要联合拟合背景

新谱先保持未扣背景状态，在同一统计模型中拟合背景与峰。可从现有幂律+峰开始；增加低阶缓变基线/实测 ZLP 尾只依据数据范围与物理条件。n=1 与 n=2 必须使用相同能窗、背景自由度、噪声与响应。

不要设置 `pre_subtracted=true` 来伪装成“常数背景模式”；扩展明确的 `baseline_mode`/输入层字段。历史安全 cap、Auto 背景评分、负值修正只进入历史复现/单项敏感性版本，不作为新统计模型里的隐含规则。背景扣除后为负的观测不能裁成零。

### 8.3 仪器与 q 平均前向模型

若有可信能量响应，模型用响应卷积后与观测比较；能量采样较粗时对像素能量边界积分。对宽 q bin，模型必须与成员/权重用同一求和算子。

先在 bin 层提取有效参数；检查 bin 内峰位置变化是否可忽略。必要时局部用 E_j(q)=E_j(q_B)+v_j(q−q_B) 等低阶展开计算所有成员期望，再求和；不强制 √q 或 rapid-rise/plateau。原始单通道的 q 接受已纳入响应时，不再次重复卷积其像素孔径。

无可靠响应时先报有效/观测宽度。不能把 Lorentzian 宽度按高斯方差直接相减，也不能把叠加展宽、未分离分量和内禀阻尼混为一谈。

### 8.4 拟合实现

- 显式 n=1/n=2，而非只设置 max_peaks=2。复用当前双初值入口时检查实际拟合峰数，并保留 collapse/failure；双分量失败不自动变成成功单分量。
- 单分量是独立对照模型，不是失败回退；H0/H1结果并列记录。
- 首轮主线关闭弱峰硬删除，`min_peak_amplitude_fraction=0`，但必须再检验零幅度、退化与参数可辨识性；这不是宣布所有弱分量都可信。
- 多初值覆盖峰位、间距、宽度、面积比及背景。建议每个双分量代表 bin 24–40 组起点，资源不足可先12组并明确标记；加入旧单峰中心±split只是其中一族，不能是唯一初值。
- 拟合前验证 lb≤p0≤ub、finite、尺度与条件数；初值放在合法内部。记录 exitflag、边界命中、目标函数及所有局部解。不静默切换到无边界 fminsearch；失败原因要保留。
- 允许能量排序作标签定义，但不将数值最小间距当分辨率。无法区分时输出 collapse/unresolved；不人为硬分开。
- 将总曲线、所有分量、背景、残差、目标函数建立同一契约。当前通用 fitter 仍有删峰前评分与删峰后曲线不一致的风险，先加回归测试，再修复或以关闭删除的独立路径规避并记录。[S10]
- `windowed_branch_tracking` 的上下两次单峰可作初值，不作为同谱双分量联合面积结果；最终必须回到完整共同窗口联合重拟合。[S09]

### 8.5 参数与分支追踪

每个分量输出模型原生参数和曲线派生量：E0、峰顶 E_apex、观测/模型 FWHM、峰高、固定有限能窗积分、总积分参数（若有）、区间、面积比及相关性。raw 局部峰高含其他分量尾部，不能直接等于独立分量强度。沿用新增 apex/FWHM helper，但数值网格加密应做收敛检查。[S10–S11]

先独立分解，再做 q 连续性辅助。只在同一侧、合理 gap 内一对一分支匹配；两标签不能指向同一个候选。P1/P2 默认局部低/高能标签，不自动等于穿越处的绝热本征模身份。

传播/ridge/path、manual anchor、quality retry、repair 均保留来源。趋势 penalty 关闭为 baseline；需要开启时提供约束释放对照。不得默认上支必须高-q平台，也不得为了 score≤300 选取路径。已有点级 evidence 分类可继续作诊断接口，但其经验阈值不是测量置信概率；缺失 robustness 等测试不能解释为测试通过。[S12–S13]

## 9. P5：验证，而不是只把曲线画顺

### 9.1 合成谱恢复与误报

用实测采样、计数/噪声、背景和响应建立两类真值：单分量、双分量。至少包含：单峰随 q 移动后分箱产生偏斜、时间漂移/多域混合、背景尾模型不匹配、弱第二峰、两峰距离变化、宽度差异、非均匀强度、有限分辨率和真实 q 缺口。逐一标明哪些是与样品有关的可检验场景，不能把模拟中的结构当成样品事实。

使用完全相同的分析和选择链处理合成谱；若使用候选路径/约束/筛选，也要包含在误报测试里。推荐 pilot 每类约100次，冻结方案后关键情形≥500次；报告实际重复次数与 Monte Carlo 误差，不能把经验通过率写成普适保证。

比较单/双模型时，不直接使用普通 F-test/标准似然比 χ² 阈值来证明加峰；第二峰面积为零时其位置/宽度无法定义。AIC/BIC只作辅助，主证据是校准后的误报率、恢复偏差、可辨识区间、原谱重建及重复一致性。

### 9.2 不确定度与模型敏感性

优先从原生重复帧重采样整个校准—分箱—拟合流程；有时间关联就按时间块。只有二维谱时用已说明噪声模型的参数 bootstrap/profile-likelihood；不要把经滤波、归一化的残差当作独立同分布样本随意打散。区间必须注明只包含统计误差还是也传播了背景、响应、标定和分箱选择。

做背景方案、E上界、N与bin起点、初值族、标准Lorentzian/DL、原谱/历史去噪、约束有无、±q分别等敏感性分析。N=1低计数未能独立分开，不自动否定N=3/5结果；但N更大时须排查单峰色散混合误拆。若模型明显不适用，不以有限宽度界把所有残差塞进第二峰。

### 9.3 验收标签

保留已有 `data_supported/tracking_assisted/suspicious` 字段作为 historical_diagnostic_class，另建新结果的参数级结论：

- `resolved_at_bin_scale`：该接受区间支持双分量，指定参数通过恢复/扰动检验。
- `assisted_or_model_dependent`：依赖连续性、特定峰形、背景或强平滑；单独标记。
- `unresolved_or_invalid`：分量退化、边界、输入不明、伪双峰风险高或关键检验缺失。

每个参数还分 energy/width/area_identifiable，不用一个总标签替代。分类阈值先在合成谱试点中固定，并记录误报/偏差；未分辨是合法科研结果，不为了验收强求 resolved 数量。

旧 RESULTS_INDEX 记载 v14 score=300、v15=326，以及 v15 证据审计252点中77支持、24辅助、151可疑；这些属于不同历史版本/经验诊断，不作为本轮目标、对照真值或新统计结果。[S13]

## 10. P6：冻结配置，扩展三组，再讨论准粒子

先冻结 590 试点方案，再处理 repeat 与20w；匹配物理 bin 宽度，保存各组自己的响应/曝光/区域属性。正负 q 不平均；比较同侧对应区间，再讨论对称性。不得把1film/2film强制等同厚度比2，不自动复用物理拟合中的 `epsilon_bg=4.5` 或其他材料参数。

只有完成分量提取和验证后，才做可选物理模型：按测得的能量、色散、宽度、谱重、方向与区域依赖提出候选。检索原始论文并建立逐项对照表：完整MoS₂–metal–MoS₂结构优先，Bi/Sb/BiSb相近体系作为次级参考；局域纳米结构共振不能仅凭能量相近等同传播q分支。没有找到完全相同结构时明确写相似程度，不伪造对应文献。

新一轮不默认调用 `run_b1_physical_fit_analysis/enhancements`。原 wrapper 的 `runFits` 默认为true，而其含义是后接物理拟合；提取试点必须显式false或由无该副作用的新入口执行。Fano只允许后续对旧单峰流程作同输入回代，检验历史有效宽度，不能据拟合函数名称判定准粒子。[S09]

## 11. 最小测试与交付契约

### 11.1 必须新增或扩展的测试

1. 最新实际dq测试网格：±0.015、10w/20w步长、qskip0.0005；N=3确实发生，低-q区域不因旧0.005/0.05阈值全部丢失/保护。
2. 输入身份：指定eq3D不被替代；显式raw/no-cache只读；相同mtime不同内容/metadata不误用；已有L2不会伪装L0。
3. 分箱：sum守恒、mean和方差尺度、一对一成员、signed侧分开、partial bin、缺口/NaN、q中心/边界与接受宽度、20w matched width。
4. 同谱模型：n=1/n=2均显式；sum/mean域变换拟合一致；初值和边界合法；线形名称/单位不混；无效响应不填默认。
5. 输出一致：背景+分量=总预测、观测−预测=残差、目标函数重算一致；弱峰筛选前后保留不同字段；CI fallback带标志。
6. 一对一匹配、无跨gap和±q混合、分窗口单峰只作初始化、失效点不自动插值、manual correction与原点可回溯。
7. 单真峰分箱不误拆的负对照、双真峰恢复、宽度/面积参数退化、约束释放和模型依赖。
8. 再执行得到一致数值；新输出不覆盖旧结果/索引/论文；失败/跳过不计入通过数。

已有测试先保持，不删除不通过测试以凑绿色。旧测试中500meV或大q的合成场景可继续做算法回归，不等于允许新实验主线回到废弃窗口；新增真实标定尺度测试补足覆盖。

### 11.2 每次run需要的输出

```text
paper_results/b1_components_v2/<run_id>/
  config_resolved.mat + config_resolved.json
  input_manifest.resolved.yaml
  repo_state.json + environment.txt
  q_rule_migration.csv + data_lineage.json
  tests/ + logs/
  <session_id>/
    bins.csv + binned_spectra.mat
    single_component_fits.csv
    double_component_candidates.csv
    component_parameters.csv
    model_comparison.csv
    fit_failures.csv
    parameter_identifiability.csv
    historical_rule_sensitivity.csv
    fit_details.mat
    figures/ (raw/member/bin、两分量/背景/总谱/残差、signed-q色散与参数)
  validation/ (null_false_split、recovery、uncertainty、sensitivity)
  run_report.md + DECISIONS.md
```

run_id包含UTC时间、配置hash与源码短SHA；同run输出已存在时不覆盖。图中标明数据层、bin成员数/宽度、峰形、能窗、背景与分类；PNG用于查看，实际参数与数组必须同时提供。无法跑到后期时只输出已执行部分，不制造空结果伪装完成。

## 12. 第一轮给 Codex 的明确任务边界

**本轮优先完成 P0–P3：实际输入/轴审计、必要的配置修复、带测试的分箱、三个代表位置的试拟合。随后在条件允许时继续代表性bin验证。** 不要耗尽预算在全仓库重构、GUI美化、物理拟合或反复重写计划。

若本地没有实际数据：照样完成源码映射、标定与分箱合成测试、入口/配置实现；报告缺失路径，不输出实验峰位。如果只有处理谱：先跑明确标记的历史复现/探索，同时指出要升级计数域推断缺什么。缺ZLP/响应时仍允许提取观测参数，不编造本征线宽。

每阶段更新下面的记录。仅将实际执行成功标为完成；预计耗时不是完成证据。

| 阶段 | 初始状态 | 最小验收产物 |
|---|---|---|
| P0 版本/论文/输入/基线测试 | 未执行 | repo与输入台账、基线来源、测试实录 |
| P1 q/阈值/缓存一致性 | 未执行 | migration表、真实dq网格测试、明确输入层 |
| P2 最小处理与历史对照 | 未执行 | 一组可追溯L1/或清楚标注L2及处理链 |
| P3 分箱与三个代表位置 | 未执行 | N=1/3/5成员与噪声、模型/背景/残差 |
| P4 联合分解与区间 | 未执行 | 代表bin多初值候选与参数级不确定度 |
| P5 误报/恢复/稳健性 | 未执行 | 合成验证和冻结分类规则 |
| P6 三组验证/物理讨论 | 未执行 | signed-q比较、文献候选和待检验问题 |

收尾报告必须写：改了哪些文件；用了哪些实际输入和dq；测试通过/失败/跳过；真正跑了哪些bin和模型；哪些参数能报告；仍未解决的输入/科学问题；下一条可以在本地直接执行的命令。不以“代码已写好”替代实际试跑，不以“任务复杂”自动结束在新计划上。

## 13. 本计划依据

[S01–S15] 为固定提交源码/文档及测试定义；详见 `SOURCE_ALIGNMENT.md`。以上数学误差传播来自明确的线性变换与噪声假设，新增模块、统计检验和阶段门槛是本计划提出的实施规范，不声称仓库已实现或验证。输入台账中的未知项不得用历史AI文本补齐。
