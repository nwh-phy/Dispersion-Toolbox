# B1 双响应的两个最小理论模型

本推导回答：在不新增实验、不指定 BiSb 晶相和载流子参数的条件下，什么机制可以产生两种纵向电子响应，并让它们在 q-EELS 中可见？

这是长波、电静近似下的条件模型，不是该异质结构的完整 RPA/DFT 结果。旧论文的宽 B1 plasmon-like candidate 是历史起点；这里使用修正后的 q 标定及后续原始帧分析。Mac 上的论文源码不是已核验的最终提交稿。

## 1. 从电子响应到集体模

定义 E=ℏω，q 为面内动量模长。外场 U_ext 引起密度变化 δn=Π₀U_tot；诱导电荷又产生库仑势 Vδn。因此

\[
\delta n=\Pi_0(U_{\rm ext}+V\delta n),\qquad
\chi=(1-\Pi_0V)^{-1}\Pi_0.
\]

集体模的候选位置由介电矩阵的零点确定：

\[
\det(1-V\Pi_0)=0.
\]

这一步把电子运动与它自身产生的电场联系起来。仅给损失谱相加两个 Lorentz 峰，没有自动施加这个自洽条件。

二维抛物带自由载流子的长波、高相速度近似为

\[
\Pi_i(q,E)=\frac{\hbar^2 n_iq^2/m_i}{E(E+i\gamma_i)},\qquad
V_{ij}(q)=\frac{e^2}{2\epsilon_0\kappa q}e^{-q|z_i-z_j|}.
\]

γ 的单位是能量；这个近似没有电子–空穴连续谱和 Landau damping。定义

\[
\alpha_i=\frac{\hbar^2e^2 n_i}{2\epsilon_0\kappa m_i},
\qquad [\alpha_i]={\rm eV^2\,\mathring A},\qquad A_i(q)=\alpha_iq .
\]

所有 α 都是待约束的电荷刚度，不分别指定 n、κ 或 m。数值代码检查了 SI 与 eV/Å 单位转换。

## 2. 两个空间分离的 Drude 导电片

令两个通道相距 d，忽略隧穿，采用均匀、无损背景。在零阻尼极限，

\[
(E^2-A_1)(E^2-A_2)-A_1A_2e^{-2qd}=0,
\]

\[
E_\pm^2=\frac{A_1+A_2\pm
\sqrt{(A_1-A_2)^2+4A_1A_2e^{-2qd}}}{2}.
\]

相同导电片给出同相与反相模式：

\[
E_\pm^2=\alpha q(1\pm e^{-qd}),\qquad
E_+\sim\sqrt{2\alpha q},\quad E_-\sim q\sqrt{\alpha d}\quad(qd\ll1).
\]

这与双层电子气的长波 RPA 结果一致。[Nakamura & Tagawa, JPSJ 58, 646 (1989), DOI:10.1143/JPSJ.58.646](https://www.jstage.jst.go.jp/article/jpsj1946/58/2/58_2_646/_article/-char/en)

### 能量比的条件约束

记 R=E₊²/E₋²。对任意正 A₁、A₂，相同权重使能量比最小，因为

\[
\frac{(A_1-A_2)^2+4A_1A_2e^{-2qd}}{(A_1+A_2)^2}
=1-\frac{4A_1A_2}{(A_1+A_2)^2}(1-e^{-2qd})\ge e^{-2qd}.
\]

若把两个经验 E₀ 当成无阻尼模式尺度，则

\[
d\ge-\frac1q\ln\frac{R-1}{R+1}.
\]

代入 q=0.005 Å⁻¹、E₋=0.83545 eV、E₊=1.41736 eV，得到 d≥145 Å，即约 14.5 nm。相同导电片达到下限；改变两片的 Drude 权重不能压小这个下限。

**这是简单模型的条件约束，不是样品厚度测量。** 输入 E₀ 来自阻尼经验拟合，尚未被证明为介电本征频率；有限接受范围须在直接拟合原始谱时重新处理。非均匀介电环境、隧穿、有限厚度三维电荷分布等不受这个特定模型约束。

### 两个本征模不保证都能被电子束看见

法向穿透的电子在不同高度产生相位不同的外势：

\[
\phi_{\rm ext}(q,z,E)\propto
\frac{e^{iEz/\hbar v}}{q^2+q_E^2},\qquad q_E=E/(\hbar v).
\]

法向 EELS 需要这一空间相位投影，见 [Rodríguez Echarri et al., PRResearch 2, 023096 (2020), Eq. 13–15](https://www.nanophotonics.org/assets/pdfs/paper344.pdf)。

归一化电荷坐标的动力学矩阵和束流驱动向量为

\[
C=\begin{pmatrix}A_1&\sqrt{A_1A_2}e^{-qd}\\
\sqrt{A_1A_2}e^{-qd}&A_2\end{pmatrix},\qquad
f=(\sqrt{A_1},\,\sqrt{A_2}e^{iEd/\hbar v})^T .
\]

投影损失为

\[
\mathcal L=-\operatorname{Im}
f^\dagger[E^2I-C+iE\,{\rm diag}(\gamma_1,\gamma_2)]^{-1}f .
\]

相同导电片的模权重含 1±cos(Ed/ℏv)。当 Ed/ℏv 很小时，反相运动的净驱动几乎抵消。30 kV 时 ℏv≈648 eV Å。在上述 14.5 nm 示例中，低能模相位约 0.187 rad。计入两模不同的频率和运动学权重，弱阻尼、点 q 极限的低能模积分占比仅约 1.86%。

所以计算必须包含电子束耦合。弱阻尼的 1.86% 不等同于宽峰拟合分账，也不排除非对称通道、复杂屏蔽或其它几何的双模。

同一平面内两个纯 Drude 通道是重要极限：d=0 时两项极化率相加，只给出一个有净电荷权重的有限频率模式，E²=(α₁+α₂)q。多能带体系的额外声学模需要更完整的非局域响应，不能仅凭“两种载流子”自动计数为两个 plasmon。

## 3. 载流子与动态电子极化共同形成两支纵向模

真实电子系统还有有限能量电子跃迁。用一个正强度 Lorentz 极化近似其中一组跃迁，并把其它较远跃迁记为静态背景：

\[
\epsilon(q,E)=b(q)
-\frac{a(q)}{E(E+i\gamma_D)}
-\frac{c(q)}{E^2-\Delta^2+i\gamma_L E},
\qquad b=1+r_0q,\quad a=\alpha_Dq,\quad c=\alpha_Lq.
\]

Δ 是未加本模型自洽库仑反馈时的电子共振能量；未指定为某条 BiSb 或 MoS₂ 跃迁。r₀ 是其它极化通道的有效屏蔽长度，不能与显式 Lorentz 通道重复计数。r₀、Δ、α_D、α_L 均未知。

静态非局域屏蔽使二维 plasmon 色散变平的依据见 [da Jornada et al., Nature Communications 11, 1013 (2020), DOI:10.1038/s41467-020-14826-8](https://pmc.ncbi.nlm.nih.gov/articles/PMC7035343/)。其材料、q 尺度和参数不能直接移植到 BiSb；这里额外保留动态极化通道，是本轮提出的最小模型。

零阻尼下，置 A=a/b、B=c/b、x=E²，ε=0 给出

\[
x^2-(\Delta^2+A+B)x+A\Delta^2=0,
\]

\[
\boxed{E_\pm^2=\frac{\Delta^2+A+B\pm
\sqrt{(\Delta^2+A+B)^2-4A\Delta^2}}{2}}.
\]

等价的耦合坐标本征值问题为

\[
H=\begin{pmatrix}A&\sqrt{AB}\\
\sqrt{AB}&\Delta^2+B\end{pmatrix}.
\]

耦合 √(AB) 由自洽响应固定，不是额外任意调整的参数。A、B>0 时，

\[
E_-<\Delta<E_+,\quad E_-^2+E_+^2=\Delta^2+A+B,\quad E_-^2E_+^2=A\Delta^2.
\]

裸电子共振处于两支纵向模之间。小 q 下低支从零出发，高支靠近 Δ；强屏蔽使 A、B 在较大 q 趋于常数，两支都可趋平。

### 频率和权重一起计算

当 γ_D=γ_L=γ 时，记 u_j 为 H 的正交本征矢、g=(√A,√B)，有精确分解

\[
\frac1\epsilon=\frac1b+\sum_{j=\pm}
\frac{S_j}{E^2-E_j^2+i\gamma E},
\qquad S_j=\frac{|g^Tu_j|^2}{b}\ge0,
\]

\[
-\operatorname{Im}\epsilon^{-1}
=\sum_j S_j\frac{\gamma E}{(E^2-E_j^2)^2+\gamma^2E^2}.
\]

两个损失振子的振幅由同一套参数决定，不能独立逐 q 调整。正残数保证被动响应；阻尼不同时直接计算复介电函数，不沿用这个简单正峰分解。

数学模型的能量加权和规则为

\[
\int_0^\infty dE\,E[-\operatorname{Im}\epsilon^{-1}]
=\frac\pi2\,\frac{(\alpha_D+\alpha_L)q}{(1+r_0q)^2}.
\]

这是频率无关背景 b 下的有效模型和规则，不是完整材料的总电子数和规则。真实高频下须恢复远处跃迁的频率依赖。

### 什么是输入，什么是预测

给定参考 q 的两条目标模式尺度和处于其间的 Δ，可以反解

\[
A_{\rm ref}=\frac{E_-^2E_+^2}{\Delta^2},\qquad
B_{\rm ref}=E_-^2+E_+^2-\Delta^2-A_{\rm ref}.
\]

本轮用此式匹配经验 E₀ 的平均值，再计算其它 q 的能量、残数、计数谱和分量比例。**参考点两条能量属于校准，不能作为模型正确的新增证据。** 仅靠两条能量不能确定 Δ 和 r₀。

把所有极化引起的库仑项乘 λ，即 (a,c,r₀)→λ(a,c,r₀)，连续追踪两支模。λ→0 时低支回到零，高支回到 Δ。高支有电子跃迁起源，也受到库仑重整化；“带间辅助的纵向集体响应”比“第二个独立自由载流子 plasmon”更明确。完整材料中的集体性仍须能带极化和电荷振型来检验。

## 4. 与 q-EELS 几何连接

相同电静近似下，每单位面内动量面积的谱形使用

\[
I_i(E)\propto R_E*
\int d^2Q\,W_i(\mathbf Q)
\frac{Q}{(Q^2+q_E^2)^2}\mathcal L(Q,E).
\]

没有额外插入径向 q 雅可比；W 是探测器动量面积上的接受权重，环积分与像素积分不能混用。共面混合模型采用 \(\mathcal L=-\operatorname{Im}(1/\epsilon)\)，双导电片采用空间相位投影。

数值复用上一轮完整 Q 接收积分及能量卷积算子。计数谱逐 q 归一化，只比较谱形；低能比例为两个计算出的 B1 模在 0.5–1.8 eV 的积分比例，没有独立分量振幅。算子的每 q 数值归一化不改变这个比例，不能据图比较绝对计数。

示例用 590 的记录几何：slot 半宽 0.0066 Å⁻¹、q 分箱宽 0.0015 Å⁻¹、探针 σ=0.000414 Å⁻¹、Gaussian 能量核 FWHM=23.013 meV。没有新增真实仪器核、解卷积或多重散射分析。

## 5. 阻尼和可见峰数

E_j 为无阻尼本征能量。相同 γ 下的复极点为

\[
E_{\rm pole}=-i\gamma/2\pm\sqrt{E_j^2-\gamma^2/4}.
\]

损失峰顶取决于线形，计数峰顶进一步受到运动学和孔径影响，三者分别报告。两支复响应可以合并成一条宽计数峰，局部极大值数不能直接计数集体模。

模型省略真实电子–空穴连续谱、q/能量依赖阻尼、各向异性、SOC 能带、隧穿、迟滞及封装中的空间非局域响应。极小 q 的电静曲线仅作数学极限检查，不是靠近光线区域的定量预测。
