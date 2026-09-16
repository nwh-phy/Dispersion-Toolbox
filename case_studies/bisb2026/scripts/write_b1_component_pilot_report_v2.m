function write_b1_component_pilot_report_v2(out)
% Summarize only completed artifacts; refuse to replace an existing report.
out=char(out); report=fullfile(out,'run_report.md'); assert(~isfile(report),'Report exists');
cfg=jsondecode(fileread(fullfile(out,'config_resolved.json')));
m=jsondecode(fileread(fullfile(out,'input_manifest.resolved.yaml')));
dest=fullfile(out,'590_PL2_10w'); c=readtable(fullfile(dest,'model_comparison.csv'));
p=readtable(fullfile(dest,'component_parameters.csv')); sel=readtable(fullfile(dest,'representative_selection.csv'));
id=p(:,{'key','n','component'});
id.energy_identifiability=repmat("not_assessed_P4_P5",height(id),1);
id.width_identifiability=id.energy_identifiability; id.area_identifiability=id.energy_identifiability;
id.interval_status=repmat("not_computed",height(id),1);
writetable(id,fullfile(dest,'parameter_identifiability.csv'));
fid=fopen(report,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# C20 P0–P3 实际执行记录\n\n');
fprintf(fid,'本轮完成输入与 q 规则核对、原生帧最小累计、固定分箱和三个位置的同窗试拟合。未进入 P4/P5 的可辨识性验证或下游物理拟合。\n\n');
fprintf(fid,'## 实际输入\n\n| session | dq / Å⁻¹ | 原生 q=0 列 | 能量范围 / meV | 排除原生列 | 分箱数 |\n|---|---:|---:|---|---|---:|\n');
for i=1:numel(m.sessions)
 s=m.sessions(i); a=s.actual;
 fprintf(fid,'| %s | %.7g | %d | %s | %s | %d |\n',s.session_id,a.dq_Ainv,a.q_zero_native_channel,mat2str(a.energy_range_meV),mat2str(a.invalid_native_q),s.bin_count);
end
fprintf(fid,'\n各原生数组均为 300×512×1028，dE=4 meV。输入 SHA256、字段、历史轴与 CSV 栅格检查见 input_manifest.resolved.yaml（JSON 兼容 YAML 1.2）。\n\n');
fprintf(fid,'在首次未掩码试跑中，探测器饱和值错误主导 ZLP 校准。该 run 已中止并标记 INVALID_RUN.md。本轮先按显式阈值排除任一帧含 uint32 饱和值/非有限值的整列 q，再累计和定标；排除列在输出中保存为 NaN。原始文件未改写。\n\n');
fprintf(fid,'L1 保留采集端已有校正；新处理只有帧累计、轴转换和零点定义，无逐-q归一化、去噪、预扣背景、反卷积或逐-q能量对齐。帧 ZLP 漂移范围：\n\n');
for i=1:numel(m.sessions)
 s=m.sessions(i); q=jsondecode(fileread(fullfile(out,s.session_id,'raw_quality.json')));
 z=q.frame.zlp_energy_pixel;
 fprintf(fid,'- %s：%g–%g 像素，跨度 %g meV；这是漂移诊断，不是仪器分辨率。\n',s.session_id,min(z),max(z),(max(z)-min(z))*4);
end
fprintf(fid,'\n历史 qe_pp 已使用登记 dq；历史裁剪列到原生列的逐列映射尚未证明，因此未转换历史点或将旧谱用作新主输入。\n\n');
fprintf(fid,'## 分箱和实际拟合\n\n采用 ±0.015 Å⁻¹、q_skip=0.0005；按 signed-q 和原生连续有效列固定分箱。10w 的 N=1/3/5 独立保存；20w 另有 N=6/10 等物理宽度分箱。sum、mean、成员散布、有效成员数、边界及 partial 标记见各组 binned_spectra.mat 和 bins.csv。\n\n');
fprintf(fid,'探测器校正/插值与帧漂移使独立 Poisson 假设未经验证。measurement variance 保留 NaN；member_scatter_variance 仅表示成员谱之间的散布，不能当作测量方差或置信区间。\n\n');
fprintf(fid,'| N | 位置 | 实际 q / Å⁻¹ | bin id |\n|---:|---:|---:|---:|\n');
for i=1:height(sel), fprintf(fid,'| %d | %d | %+.5f | %d |\n',sel.N(i),sel.representative(i),sel.actual_q(i),sel.bin_id(i)); end
fprintf(fid,'\n位置在拟合前按低/中/高 |q| 目标指定，选择同侧最近的完整有效 bin，未按双峰拟合成功率挑选。\n\n');
fprintf(fid,'实际执行 %d 个单/双模型任务（%d 数值收敛、%d 无收敛解，%d 个选中解命中边界），每个模型 %d 个初值。\n\n',height(c),nnz(c.success),nnz(~c.success),nnz(c.boundary),cfg.n_starts);
fprintf(fid,'模型为 lorentz_symmetric 与保留原公式的 lorentz（Drude–Lorentz），均联合拟合相同自由度的幂律背景。主窗 300–1800 meV；300–2000、300–2100 是具名敏感性。目标函数是未加权最小二乘描述量，不是经标定的似然；不以 SSE 改善证明双分量。\n\n');
fprintf(fid,'保存了 %d 行分量参数。native_width 对标准 Lorentzian 是模型 FWHM，对 DL 是 Gamma；有限窗局部 floor 的 FWHM 另列为诊断，不能混同。面积仅是给定模型在有限窗内的积分。\n\n',height(p));
fprintf(fid,'所有科学标签维持 unresolved_or_invalid。没有本轮统计置信区间、恢复校准后的误报率、参数可辨识区间或本征寿命。多初值、边界、Jacobian 条件数和残差用于发现退化，不替代这些验证。\n\n');
fprintf(fid,'## 验证记录与限制\n\n测试实录见 tests/；首轮另有 baseline_inventory.md。逐谱输出契约检查了背景+分量=总预测、观测−总预测=残差、重新计算 SSE 与优化目标一致。source_integrity.txt 记录输入内容哈希复核。\n\n');
fprintf(fid,'未进行帧对齐敏感性、时间相关噪声估计、P4/P5 bootstrap/profile/误报校准、冻结方案跨组拟合或物理归属。中心排除宽度目前仅为明确配置，未完成独立仪器接受核和更宽中心束掩码的验证。论文及旧结果/索引未覆盖。\n\n');
fprintf(fid,'下一次重跑（项目根目录 MATLAB，自动生成新目录）：\n\n```matlab\naddpath(''case_studies/bisb2026/scripts''); out = run_b1_component_pilot_v2();\n```\n');
end
