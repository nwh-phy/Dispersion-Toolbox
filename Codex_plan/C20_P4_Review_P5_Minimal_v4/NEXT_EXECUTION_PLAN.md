# C20/B1：P4实际交付补全与最小P5执行规范 v4

性质：待本地Codex实施和运行的计划，不是已完成结果。本计划接续而不覆盖v2/v3。审阅依据是用户粘贴的实际v3入口/拟合器/CSV；不把最新远端代码自动当成父run的执行源码。详情见REVIEW_P4_ACTUAL.md与SOURCE_EVIDENCE_INDEX.md。

## 0. 本轮只解决四件事

A. 读父MAT完成真实的参数级边界与raw→ordered映射；B. 把12个N3拟合的全部数组保存出来，补数值嵌套/失败路径与少量背景对照；C. 证明帧与ZLP校正的实际含义，输出真实的分块B1谱；D. 在这条固定小流程上开始单模式误拆与双分量恢复的最小P5，而非全q扫描。

不要为了形式上补齐上一长计划而做GUI、全仓库重构、全窗口/N/峰形的笛卡尔积或下游物理拟合。此次成功不要求两个分量必须分开，而要求可说明“哪些参数由数据约束、哪些由边界/背景决定”。Fano、√q、预设准粒子与厚度比不进入此次主线。

科学口径延续：C20是用户更正的Bi-rich约Bi:Sb=70:30；完整MoS2–metal–MoS2；P1/P2是旧B1内部有效分量，不是旧B1/B2标签。文件夹名不是材料成分、帧语义、层数或剂量的证明。

## 1. 输入、只读约束与实际身份

两个父目录保持只读：

```text
P0-P3: paper_results/b1_components_v2/20260912T021813081Z_31fcc543_18af79b6/
P4:    paper_results/b1_components_v3/20260912T035735540Z/
新run: paper_results/b1_components_v4/<UTC>_<configHash8>_<codeHash8>/
```

输入按用途获取：

| 用途 | 优先实际输入 | 缺失行为 |
|---|---|---|
| 边界审计 | v2/590_PL2_10w/fit_details.mat、all_multistart_candidates.csv、model_comparison.csv、config/source_snapshot | 缺MAT可按完整p/lb/ub/scale等字段恢复，但字段不足不得造结果 |
| 同中心重拟合 | v3/590_PL2_10w/centered_binned_spectra.mat，v2/590_PL2_10w/L1_minimal.mat | 重新从已验证L1构造9个bin，不读旧GUI处理谱代替 |
| 旧12项对账 | v3/component_parameters.csv、centered_bins.csv、run_report与stage_status | 该CSV没有真正分量参数，只作模型级对照 |
| 帧语义和A1 | v2 manifest中的590 NPY和同名JSON、raw_quality、invalid_native_q、相关采集记录 | 从session_id与文件角色查找，禁止sessions(1).files(1) |
| 已有新MAT | 若本地有v3 fit_details或候选另存，查其生成代码、时间和hash | 不假定存在；不存在就只重跑12模型来恢复可审计输出 |
| 仪器参考 | 同条件ZLP/暗场/增益/束流/接受核 | missing_optional，限制本征解释，不阻塞有效谱诊断 |

读取父输入哈希并重算实际使用的raw/JSON/L1/MAT身份；父manifest含不存在的可选文件时保留missing。不要通过修改父manifest“修好”匹配。历史op_history等不得成为这一新分析的强制依赖。

查询`git -C <root>`状态、staged与unstaged差异；记录实际`which -all`解析路径。至少快照入口、fitter、registry/peak_models、分箱/帧helper、审计/写出helper及实际测试。保存执行源码集合hash，不只HEAD。用户未提交源码保留；不reset/clean/push。不在原始数据目录写缓存。

本轮代码复用现有MATLAB模块，优先扩展现有v3入口或抽出少数短helper。先实现后才提供可执行MATLAB命令，不把本计划中的建议接口名当已存在命令。

## 2. A：先审已有候选，不先重跑2592次

### A1. 完整性与字段口径

先从MAT实际统计54个key、108任务和各自起点。24起点配置下预期2592候选，5/8参数下预期16848参数行、702选中参数行；任何差异列出原因，不填充假行。每项成功/失败、参数向量长度、上下界、finite和selected索引都独立检查。

CSV统一使用`bin_N_requested`、`source_q_count`、`n_components`，不能把f.n_components又命名为N。记录窗口实际上下界、原key、signed q、成员、scale、ampunit、parent文件hash。

### A2. 对每个候选重建映射

对p(3:end) reshape为n×3；按第一列E排序获得order，建立逆映射raw_to_order(order)=1:n。原参数idx≥3对应raw_component=ceil((idx-2)/3)。输出原槽位与最终energy-order标签两列。

映射必须在每个候选上独立计算，不复用selected_start的order。候选交换对称、相等E0和A≈0分别标记；能量排序只是标签约定，不能据此证明两个模态可区分。

对选中解按排序后参数重建曲线，与父f.parameters/f.components/f.prediction逐点核对，既检查映射也检查单位。背景参数ordered_component=0/NA，不挂到P1。

### A3. 边界与单位

原始p/lb/ub/distance保留，另存native值及单位：E/width乘1000回meV，peakA乘ampunit*scale，背景幅度乘scale，r不变。距离按同一线性换算；Inf上界保持Inf，不赋予虚假命中。

复算原fitter的上下界容差与原布尔标志，比较mismatch；新增工程容差需具名，不能改变旧标准后声称“触边率下降”。失败/NaN候选boundary_status=not_assessable，不等于无触边。

分类至少包括背景B下界、r下/上界、分量A下界、E0下/上界、width下/上界、equal-center/零幅度退化。用背景边界/峰边界正交字段，不把背景B=0自动判成峰位失败。A接近0时区分仅数值系数阈值与在参考窗中的实际分量贡献；后者未统计校准不能命名检出阈值。

### A4. 产物与门槛

输出真正的boundary_by_parameter_native.csv、boundary_type_summary.csv、candidate_solution_families.csv、parent_contract_checks.csv。解族包含全部候选的目标、exitflag、finite、是否选中、排序参数；无原J时scaled-J诊断标未运行，可对少数selected重新算差分J，不造旧运行最优性。

必须自动核对按parameter映射重组出的selected-any-boundary与父CSV的98，且能说清每类具体数量。所有统计类别可重叠，合计不需要等于98。若无法重构，先输出不一致表，不把这个阶段写complete。

## 3. B：恢复可复核的12个拟合，再做最小数值/背景检查

### B1. 保存比再次提高拟合次数更优先

先查v3是否有独立保存的新拟合MAT；没有则从9个bin内的三个N3重新运行原12任务，保持300–1800、两峰形、n1/n2、旧24起点和旧seed+ti对账。保存完整对象，不要求重新分析其他N/q/window。这是新run中的reconstruction，不声称找回了原来未保存的瞬时对象。

每个key独立保存checkpoint，即使后续图像/CSV失败也不丢失已完成拟合。元数据记录真实数据hash、mask、sum/mean尺度、原生E、模型函数hash与实际选项。

`model_comparison.csv`应有12任务记录；`component_parameters.csv`应有18分量记录（失败任务仍用明确NA行而非伪造值）；`all_candidates.csv`在旧24起点流程应有288候选；额外witness和新增起点另用candidate_type、candidate_id计数。`fit_details.mat`包含全部curves/parameters/candidates。

### B2. H0嵌入和求解诊断

同一数据/scale/window/model/bg/bounds下，把H0的p原样放入H1第一分量，增加合法center/width和精确A2=0，直接计算曲线和目标。要求Q_witness≈Q_H0。

分别保存Q1_best_optimized、Q2_best_optimized、Q2_witness、Q2_best_feasible和selection_type。若只因witness使Q2≤Q1，不把这个恒真的构造当作H1优化成功；保留optimized nesting violation标志。witness不标exitflag>0或“两分量发现”。

保留旧起点表作为legacy对照，新增模型专用随机流/持久起点表；增加或删除H0起点不应改变H1既有起点。为n2追加小规模非对称宽度、弱面积和split起点（例如4个受控加8个独立随机），只有诊断需要才扩大。不得无限加起点凑低SSE。

保存output.firstorderopt/iterations/funcCount、exitflag、原p、边界、lambda（适用时），有限但未成功且目标更低的候选一并展示。KKT/边界最优性需结合算法，不用全空间梯度非零直接判错。cond(J)同时记录列尺度；背景零幅度导致的结构秩亏与峰参数相关性分开。

不强制实现variable projection。若仍有明显求解差异，再用固定非线性参数下的线性系数剖面作小规模交叉检查。其active-set和数值导数须测试。

### B3. 背景对照不能继续停留在配置名称

先完成相同旧B0模型的求解对照，再增加B1=power_law_plus_nonnegative_constant。H0/H1共同增加C，不负背景抵消巨峰。r=0时B/C重复列、B=0时r未定义需标注。

以三个N3位置为上限，每模型/背景均保存，不同时扩窗口、缩width上界或固定峰间距。新背景额外最多12个任务，且在代码接口与输出中可核对baseline_mode真正生效。触边类型若指出其他更相关机制，登记理由再调整，不按哪套分得漂亮选背景。

当前width上界5内部单位=5000meV保留。仅在实际命中且有诊断价值时用10000meV作具名敏感性；若参数继续追随上界，应判window-unconstrained而不是压回更窄界。

### B4. 共同物理口径与失败输出

保存native_E0/native_width/native_A，以及零基线曲线峰顶和参考窗[300,1800]meV积分。半高交点未被观测窗包围时width_bracketed=false。finite-floor FWHM保留legacy字段但不用于寿命。DL-Gamma不是标准Lorentzian FWHM，native_A不可跨模型比较。

内部1000/ampunit*scale是数值反变换；不能把它误修成没有缩放。当前point-sampled预测不额外乘4。后续像素平均作为独立测试，对观测纵轴、模板A和dE同时检查。

访问selected_start前检查success与有效索引；失败写字段NA和原因，不abort整个run、不把失败从汇总删除。图使用显式句柄obs/bg/P1/P2/total，正确数目legend，附真实q成员/窗口与处理级别。

## 4. C：同中心与帧检查的最小补全

### C1. 分箱保持成功的部分

保留现有三个中心231、251、261及其N1/3/5成员。扩展统一helper，保存sum/mean、source_q_count、native members、q_left/right/width、mask、partial与processing_level。缺方差就NaN；成员散布不能充当measurement variance。

必须测试的是新实际调用路径，不是仅调用旧qe_prepare_count_bins一次。新helper接受显式q_axis与source_channel，不依赖当前位置索引=原生列假设。对这次均匀网格检查目标中心、步长、中心掩码和有效性，不重估q零点。

保留镜像N3作为具名待执行任务，先可构造原谱对照，不为了完成清单马上增加所有模型。不同N共享数据，不拼为独立样本。N1/5的拟合敏感性排在模型输出与帧检查之后。

### C2. 源码语义和原始数据一致性

从manifest按590 session_id解析raw与JSON；核对真实shape/dtype/axes/is_sequence、时间或空间语义、曝光/时序。`is_sequence=true`若不足以证明同位置等条件，继续读采集记录；无法核实集中问用户一次，不从目录名10sx300代填。

继承父invalid_native_q及实际阈值掩码；计算零损统计时采用显式valid_q成员列表，不把未验证的硬编码461视为通用规则。可以在求和中排除列，但输出invalid mask，不把零填数组重新当观测。读取一次raw并保存小型派生帧×能量×代表q块，避免每个fit重读大文件。

重算A0累计并与父L1在共同轴/掩码上对账；不一致先停真实数据推断。没有新证据不修改父461掩码。

### C3. A1真正实现而不是只改文件名

用固定的ZLP能区（初始可取父全局零点附近[-100,100]meV，执行前检查其支持）和有足够弹性信号的预先声明q参考带估每帧共同偏移。选择原则来自ZLP有效性，不来自B1拟合改善。输出window/reference、raw ZLP位置、有效性与offset；无效帧不得自动置0。

先整数非循环平移：measured_offset=observed-reference，correction=-measured_offset。用已知±位移、零位移、缺ZLP和端点样本测试符号与共同支持。对共同保留支持比较计数，不要求裁剪后的和等于含已丢边缘的全原数组。A0/A1拟合比较使用相同能量支持。

A1有效性判断用ZLP宽度/峰位重复和合成真值；不能用B1是否更分开决定接受。亚像素插值、静态δ(q)和前向漂移响应只在必要时追加，明确其协方差与未知区间。

### C4. 真正的6块B1谱与统计目标

先把1–50…251–300称为sequence blocks，帧语义核实后才称time blocks。对三个目标N3至少保存18条块谱及对应成员/每块有效帧数，同时保存同一完整序列的A0/A1。

对时间重复且曝光相同，输出block_sum和per_frame_mean两种表示以避免50帧块与300帧累计的幅度混淆。若有真实曝光记录，可按剂量/曝光定义比较；不能用逐块Area归一化掩盖谱形演化。对每块拟合仅在需要时做；先看原谱及共同模板的幅度/形状残差。

给出ZLP/强度随序列、去趋势前后相关性与B1残差变化；一阶相关大不能直接用AR(1)公式套有效样本数。六块不是天然独立；块长50不是已选择的bootstrap长度。若非平稳，先定义分段或时间平均目标。

## 5. D：最小P5，不让修复无限拖延科学验证

### D0. 两级前进门槛

工程级模拟可在帧语义尚未完备时运行，但只能写algorithmic_smoke，不报告实验误报率。真实数据统计标签升级需要A/B输出可复核、raw与帧语义可解释、具体统计目标/噪声依据已声明。

边界仍存在不禁止模拟检测极限；禁止的是把带未验证噪声、边界或原始映射不明的结果称resolved。没有MATLAB或raw的环境不能由文字计划冒充运行。

### D1. 先冻结主分析，少开组合

主分析保持三个预先给定N3位置、300–1800窗、标准Lorentzian有效模型；DL为具名敏感性。不把“哪个模型/N/窗口最容易显著”选择后忽略搜索过程。若经诊断改变主背景，写冻结记录和选择理由；开发模拟种子与正式验证种子独立。

### D2. 最少情形与次数

先各10例端到端smoke，运行通过再扩大：

1. 单静态峰+所选背景：验证双分量优化/判据的基线误拆。
2. 每个原生q单模式，有数据约束的局部色散/强度变化，并按真实N3成员合并；时间漂移按A0/A1规则处理。参数范围未知则做具名场景，不说它是真实系统。
3. 两个真分量：覆盖实际边界审计附近的间距、弱分量比和宽度比，不只复用容易的750/1250例子。
4. 依据A/B/C定位结果再加最相关的一种竞争情形：背景不足或时间状态变化。不能所有问题都跳过，只检验完美H0。

主情形smoke通过后，每个关键情形至少100独立生成的模拟试点；正式性能再用独立seed扩至500或更多。每个情形保存truth、生成器、noise provenance、actual_trials、fail/collapse/边界率和Monte Carlo区间。10或100次零误报不代表0风险。

检验对象应明确为“某固定分析规则是否把本来一个局部模式误报为两个有效分量”。即使拒绝静态单峰，也不能直接排除q/时间混合或认定两种准粒子共存。

### D3. 噪声与重采样的底线

完整帧向量保持E×q相关性。帧稳定时用适当时间块；非平稳时分段/显式状态模型。先定义统计目标是累计谱、每帧平均谱还是某段谱。

令u_t(E)=mean_{i in bin}X_t(E,q_i)，实际Y(E)=sum_t u_t(E)。独立同分布条件下Var(Y)=T Var(u_t)；时间相关时加跨帧协方差。不能用Var(u_t)/T当累计谱方差。300帧估376维能量协方差秩≤299，不能直接逆；先用未加权拟合配块重采样或受验证简化权重。

参数重复性bootstrap可重采样实际帧/块；H0模拟必须从单模式均值出发。经验噪声用围绕已定义段内/经验均值中心化的整帧残差块，检查异方差；不能把含系统第二峰的H0拟合残差未经处理加回H0。校准固定的结果标conditional_on_fixed_calibration，考虑校准误差时重估对应校准。

ΔQ=Q1−Q2只在同一观测/窗口/权重下比较；普通F/χ²阈值不适用于直接宣称第二峰成立。H0模拟标定阈值或尾概率需要声明nuisance取值范围。多个N/window不是独立实验，不能汇总成108次独立验证。

### D4. 最小参数可辨识性检查

三个固定N3位置在主背景/模型下，先做ΔE、参考窗分量比、一个关键width的少量profile-objective点（固定该参数后重优化其余）。无噪声校准不叫95%profile-likelihood。允许零幅度、合并和边界结果保留；profile仅被人工界截住就报告未约束。

真正得到统计区间时保存所有bootstrap失败/collapse/零幅度结果；不得只对“两个峰都成功”的幸存样本算窄CI。energy/width/area分别给状态。等中心不同width并不自动等于单分量；H0按A=0等实际定义嵌套。

## 6. 本轮必须运行的新路径测试

至少覆盖以下场景，测试状态绑定实际执行代码hash：

- raw_component_1高能且width触上界、raw_component_2低能：上界应归到最终P2；反向、单峰、等中心、NaN/失败、背景触边独立测试。
- 已验证L1的同中心N1/3/5正确成员和q边界；缺列/invalid_q/qskip/不均匀q轴拒绝或清楚标无效；sum/mean和方差尺度正确；测试调用新实际helper。
- H0精确嵌入H1曲线/目标一致；witness和优化成功标签不同；n1起点数变化不扰动n2既有起点。
- 全部失败或某bin无效时仍保存失败结果，不能索引NaN selected_start或无条件写complete。
- 12任务、18分量、常规288候选（按实际配置推导）的输出回读；总曲线/残差/SSE与native参数曲线重算一致。
- B0/B1配置确实改变设计矩阵，H0/H1自由度公平；sum/mean输入与nativeA缩放一致；符号/单位不能误乘4。
- A1正负/零位移、invalid ZLP、共同支持/非循环边缘；无效列不产生伪ZLP；帧×q维解释不靠shape大小猜。
- 图例句柄P1/P2/total数量与实际曲线一致；status从工件+测试判定，不从函数未报错直接complete。
- H0生成不把实际双结构均值加回；算法ic smoke和实际noise-calibrated结果标记不同。

已有七个pilot测试仍执行并保留；它们不替代这些新测试。宽域branch-window旧失败若不在调用链，记录隔离理由，不删测试。不要写全仓库全绿，除非真的全量执行通过。

## 7. 交付、状态和本轮停止条件

```text
new_run/
  provenance/    # 两个父run、实际输入/源码hash、which、配置、环境
  audit/
    candidate_inventory.csv
    boundary_by_parameter_native.csv
    boundary_type_summary.csv
    component_order_checks.csv
    parent_contract_checks.csv
  590_PL2_10w/
    centered_bins.csv + centered_binned_spectra.mat
    frame_semantics.json + frame_qc.csv + alignment_shifts.csv
    sequence_block_spectra.mat + A0_A1_diagnostics.csv
    fit_details.mat
    model_comparison.csv
    component_parameters.csv
    all_candidates.csv
    nested_objective_checks.csv
    parameter_identifiability.csv
    figures/
  validation/   # 仅实际执行的smoke/null/recovery/profile
  tests/        # 新旧测试明细与执行命令
  stage_status.json
  run_report.md
  review_packet.zip
```

小MAT纳入review_packet：至少三个原谱/mean、所有12拟合数组、order与全部288候选（通常不大）、18条真实块谱、ZLP位移、q成员。父2592候选可以给压缩CSV，不要只交汇总。巨型raw不打包。

README必须列出包内文件、大小、hash、每阶段已做/未做；打包后实际解压到临时目录验读，重算12拟合曲线/参数契约。不能说“大MAT本地保存”却实际没有save语句。

分开标记：
- input_validity
- nominal_bin_geometry
- sequence_semantics
- A1_status = estimated_only / tested_applied / not_possible
- numerical_status
- model_adequacy
- energy/width/area_identifiability
- origin_status
- P5_status

本轮最低交付A、B1/B2、C1/C2和实际小型smoke。metadata可验证且ZLP可用时继续C3/C4及B3，并开始D的100例关键场景。遇真正缺失字段仅限制相应结论；不能只再次写计划结束，也不能不做就标complete。

停止扩规模的条件是明确的具体障碍：参数表/输入无法回读、嵌套见证不等价、A1无有效参考、统计帧语义未知、拟合在极宽分量与背景之间无约束。对每项给出一个最小补证/补测，而不是笼统“数据不足”。新P5工程模拟可继续估计方法边界，真实统计声明不越级。
