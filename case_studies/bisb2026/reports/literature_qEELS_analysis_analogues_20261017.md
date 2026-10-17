# 和我们 B1 分析做法相近的文献

整理日期：2026-10-17。

目的：找出和我们做法相近的工作，用来定位分析里的盲点：
- 运动学整形与假色散；
- 动量分辨、接收范围与色散形状；
- 层数系列；
- 二维金属层 plasmon。

运动学因子本身的理论依据见 `docs/references/kinematic_factor/README.md`（9 篇，本表不重复，只在相关处引用）。

"核实程度"一栏：
- **摘要**：读过 arXiv 摘要；
- **Zotero**：已在用户的 Zotero 库里；
- **书目**：只核实了题目、期刊、卷页，没读摘要或正文。

## 1. 运动学"整形"会造出假的峰位或色散

对应：峰脊检验（`b1_kinematic_ridge.m`，不色散模型乘运动学因子即可复现原始热图上的"√q"）、正向模型。

| 文献 | 核实 | 做了什么 | 和我们的关系 |
|---|---|---|---|
| A. Guandalini et al., *Phys. Rev. B* **111**, L041401 (2025)，"Direct observation of the vanishing electron energy loss spectroscopy cross section in graphene" | Zotero（QT9ZARZJ）；PDF 已在 `kinematic_factor/` | 二维截面 ∝ q²/(q²+q_z²)²·Re σ(ω)/ω，整形函数 A = q⁴/(q²+q_z²)²。30 keV 时小 q 下强烈整形；以往"最低 q 的峰最强"来自 Δq > q_z 时混进的大 q 电子 | 峰脊检验就是这一点在我们数据上的实验演示 |
| J. Hong, A. Guandalini, W. Wu, H. Sun, F. Wu, S. Chen, C. Ma, K. Suenaga, T. Pichler, F. Mauri, arXiv:2510.08397 (2025)，"Ultraviolet optical conductivity, exciton fine-structure and dispersion of freestanding monolayer h-BN" | 摘要 | 单层到少层 hBN 的高动量分辨 q-EELS；用光学电导近似（OCA）从 q-EELS 反推 q = 0 的 σ(ω)；小 q 下第一亮激子带有两个分量 A、A′，BSE 没有预言 | 思路和我们相同：不直接读峰位，而是用模型换算回本征量。可与师兄的 hBN 文章对照（同一材料，处理方式不同） |
| A. Guandalini, A. Ferretti, D. Varsano, P. Barone, F. Mauri, arXiv:2609.29789 (2026)，"Unified description of exciton, phonon and plasmon dispersions in 2D materials from an optical conductivity approximation" | 摘要 | OCA：只用 q = 0 的光学电导就能算出有限 q 的纵向响应（q-EELS）；二维材料中激子、声子、plasmon 的线性色散主要来自二维宏观库仑相互作用的形式，而不是能带色散；复现了 hBN 的 LO 声子、亮激子和石墨烯 π plasmon | 可直接借用的正向框架：有了夹层的 σ(ω)，就能预测 B1 两支的 q 依赖 |
| H. T. B. Do et al., *Nat. Commun.* **16**, 5801 (2025) | 已有 | 运动学因子加完整接受范围积分 | 我们的基线 |

## 2. 测得的色散形状取决于动量分辨和分析方法

对应：slot 接收范围 h 的问题（`p2b_hscan_20261017.md`）、毕设热图上的"√q"。

| 文献 | 核实 | 内容 | 和我们的关系 |
|---|---|---|---|
| M. K. Kinyanjui et al., *EPL* **97**, 57005 (2012)，DOI 10.1209/0295-5075/97/57005 | 书目 | 石墨烯 π plasmon 测得线性色散，归因于局域场效应 | 争论的一方 |
| S. C. Liou, C.-S. Shie, C. H. Chen, R. Breitwieser, W. W. Pai, G. Y. Guo, M.-W. Chu, *Phys. Rev. B* **91**, 045418 (2015)，arXiv:1412.8538 | Zotero（EJ96S87X）；摘要 | q 分辨提高到约 0.001 Å⁻¹ 后，π plasmon 在小 q 呈 √q，推翻此前的线性色散报告；并区分了 π plasmon 和 Dirac 二维 plasmon | 小 q 下的色散形状由动量分辨决定，和我们 h 的处境相同 |
| F. J. Nelson, J.-C. Idrobo, J. D. Fite, Z. L. Mišković, S. J. Pennycook, S. T. Pantelides, J. U. Lee, A. C. Diebold, *Nano Lett.* **14**, 3827 (2014)，"Electronic excitations in graphene in the 1–50 eV range: the π and π+σ peaks are not plasmons" | 书目（作者取自师兄 SI 的 ref. [7]） | 认为单层石墨烯 EELS 中的 π、π+σ 峰是带间跃迁，不是 plasmon | "峰 = plasmon"不能默认；Novko 等后来有反驳 |
| V. U. Nazarov, *New J. Phys.* **17**, 073018 (2015)，arXiv:1506.05220 | 已有；摘要 | 准二维晶体不能直接用 −Im 1/ε 解释 EELS；同一能区里，色散的集体激发（plasmon）和不色散的单粒子带间跃迁同时存在、互相影响 | 和我们"一支平、一支色散"的格局直接对应 |
| arXiv:2302.06367，"Excitonic effects in energy loss spectra of freestanding graphene" | 书目 | 自由石墨烯损失谱中激子效应显著改变能隙和 π plasmon 的位置、色散和线形 | 同上 |

## 3. 层数、厚度系列

对应：一层与两层的对比（`layer_1L_2L_20261017.md`）。

| 文献 | 核实 | 内容 | 和我们的关系 |
|---|---|---|---|
| J. Hong, M. K. Svendsen, M. Koshino, T. Pichler, H. Xu, K. Suenaga, K. S. Thygesen, arXiv:2203.05088 (2022)，"Momentum-dependent oscillator strength crossover of excitons and plasmons in two-dimensional PtSe₂" | 摘要；Zotero 没有 | q-EELS 加第一性原理，从单层到体相：有限 q 下由两支随 q 往高能走的带间 plasmon 主导；单层缺少长程屏蔽，长波 plasmon 形成不了，小 q 激发全是激子性的，损失谱与光学谱重合 | 层数系列加"带间 plasmon / 激子"区分的范例，和我们的问题最接近 |
| H. T. B. Do et al. 2025（TaS₂ 单层、双层）；何沛一等（hBN 1L–25L，师兄文章及 SI） | 已有 | 以层数为主要变量 | 做法范例 |
| "Thickness dependence of surface plasmon damping and dispersion in ultrathin Ag films", *Phys. Rev. B* **72**, 205405 (2005) | 书目 | HREELS 加 STM：膜越厚色散越强，阻尼边随厚度变化，归因于垂直方向的量子阱态 | 提示一层、两层除了比能量，还要比高能支色散的斜率 |
| "Plasmonics in atomically-thin crystalline silver films", arXiv:1901.07739 | 书目 | 原子级薄的晶态 Ag 膜的 plasmon | 金属超薄膜的对照 |
| A. Rodríguez Echarri et al., *PRR* **2**, 023096 (2020) | 已有 | 超薄 Ag、Au 膜的 EELS 理论；Ag 少数原子层时局域经典模型仍适用 | 理论对照 |

## 4. 二维金属层 plasmon：经典测量与综述

| 文献 | 核实 | 内容 | 和我们的关系 |
|---|---|---|---|
| T. Nagao, T. Hildebrandt, M. Henzler, S. Hasegawa, *Phys. Rev. Lett.* **86**, 5747 (2001)，"Dispersion and damping of a two-dimensional plasmon in a metallic surface-state band" | 书目 | Si(111)-√3×√3-Ag 表面态能带中二维 plasmon 的色散与阻尼（HREELS） | √q 二维 plasmon 的经典实验 |
| H. Pfnür, C. Tegenkamp, L. Vattuone, "Plasmons in one and two dimensions", arXiv:1701.05049 | 摘要 | 综述：二维、三维电子体系的耦合，以及面内各向异性（台阶阵列），都会把 √q 色散"拉直"；准一维体系的色散形状取决于电荷分布和链间耦合 | α-Bi(110) 面内各向异性；夹层中 BiSb 与 MoS₂ 的耦合 |
| H. Pfnür, T. Langer, J. Baringhaus, C. Tegenkamp, arXiv:1101.2768，"Monopole and multipole plasmons in a two-dimensional system" | 摘要片段 | 单层石墨烯（ELS-LEED）中观察到不止一支色散的损失峰 | 多支二维 plasmon 的先例 |
| C. Elgvin et al. (2026)，"Advances in momentum-resolved EELS of phonons, excitons and plasmons in 2D materials and their heterostructures" | Zotero（364IDGPI） | 二维材料及异质结 q-EELS 方法综述 | 方法总览 |
| M. McCauley et al., arXiv:2410.23983 (2024)，"Control of 2D plasmons in the topological insulator Bi₂Se₃ with highly crystalline C₆₀ overlayers" | 摘要 | Bi₂Se₃/C₆₀ 异质结：界面 π 态杂化，二维 plasmon 被淬灭；用 q-EELS 测色散，DFT 用界面电荷偶极解释 | 封装界面会改变二维 plasmon，和"界面激发"候选有关 |

## 5. 对我们解释的启示

- **一个判据**：Nelson 2014、Nazarov 2015、Hong 2022、Guandalini 2025/2026 指向同一点：二维体系小 q 下 q-EELS 主要测的是光学电导 Re σ(ω)。不随 q 变的峰对应带间（单粒子或激子性）跃迁，随 q 往高能走的峰才是 plasmon。
- **套到 B1 上**（候选解释，和当前"两支都是 plasmon"的汇报口径不同，由用户决定是否采用）：
  - 低能支平坦、有隙、对层数不敏感，和 Re σ 中的一个带间吸收峰一致；
  - 高能支随 q 往高能走（峰脊检验和拟合都支持），符合 plasmon 的特征。
- **可做的检验**：
  1. 红外或近红外吸收、椭偏测夹层的光学电导，看 0.6–0.7 eV 有没有吸收峰；
  2. 用测得或算得的 σ(ω) 按 OCA 预测两支的 q 依赖，和 P2b 及高 q 结果比较；
  3. 一层和两层比较高能支色散的斜率，不只比能量（参考 Ag 膜的厚度依赖）；
  4. 确定 slot 相对 Bi 晶轴的方向，考虑面内各向异性对色散形状的影响（Pfnür 综述）。

## 6. 还没做的

- PDF 没有下载。建议先补进 Zotero 的是：Hong 2022（PtSe₂）、Hong 2025（hBN OCA）、Guandalini 2026（OCA）、Nelson 2014。
- 标"书目"的条目只核实了题目和卷页，作者和结论以原文为准。
