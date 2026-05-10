function out = run_b1_lorentz_tracking_optimization(options)
%RUN_B1_LORENTZ_TRACKING_OPTIMIZATION Optimize B1 Lorentz double-peak tracking.
%
% This diagnostic keeps existing Fano / single-peak B1 outputs untouched. It
% runs Lorentz-only B1 double-peak extraction variants, scores continuity
% metrics, and records the recommended tracking version.

arguments
    options.baseDateTag {mustBeTextScalar} = "260510_lorentz_tracking_optimized"
    options.runExtraction (1,1) logical = true
    options.runOverlay (1,1) logical = true
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
addpath(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

base_tag = char(string(options.baseDateTag));
results_root = fullfile(project_root, 'paper_results');
summary_dir = fullfile(results_root, ...
    sprintf('b1_lorentz_tracking_optimization_%s', base_tag));
if ~isfolder(summary_dir)
    mkdir(summary_dir);
end

variants = local_variant_specs(base_tag);
metrics = local_empty_metrics_table();
variant_outputs = struct('id', {}, 'date_tag', {}, 'overlay_png', {}, ...
    'overlay_pdf', {}, 'score', {});

for vi = 1:numel(variants)
    variant = variants(vi);
    fprintf('\n=== B1 Lorentz tracking %s (%s) ===\n', ...
        variant.id, variant.date_tag);

    if options.runExtraction
        local_run_variant(project_root, variant, base_tag);
    end

    session_metrics = local_finalize_variant(project_root, variant);
    metrics = [metrics; session_metrics]; %#ok<AGROW>

    overlay_png = "";
    overlay_pdf = "";
    if options.runOverlay
        overlay = run_b1_double_peak_waterfall_extraction_overlay( ...
            dateTag=variant.date_tag);
        overlay_png = string(overlay.png);
        overlay_pdf = string(overlay.pdf);
    end

    variant_outputs(end + 1).id = variant.id; %#ok<AGROW>
    variant_outputs(end).date_tag = variant.date_tag;
    variant_outputs(end).overlay_png = overlay_png;
    variant_outputs(end).overlay_pdf = overlay_pdf;
    variant_outputs(end).score = local_variant_score(session_metrics);
end

metrics_path = fullfile(summary_dir, 'b1_lorentz_tracking_metrics.csv');
writetable(metrics, metrics_path);

recommendation = local_choose_recommendation(variant_outputs);
readme_path = fullfile(summary_dir, ...
    'README_lorentz_tracking_optimization.md');
local_write_readme(readme_path, base_tag, variants, metrics, ...
    recommendation);

if strlength(recommendation.overlay_png) > 0 && isfile(recommendation.overlay_png)
    copyfile(recommendation.overlay_png, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_overlay_three_dataset_comparison.png'));
end
if strlength(recommendation.overlay_pdf) > 0 && isfile(recommendation.overlay_pdf)
    copyfile(recommendation.overlay_pdf, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_overlay_three_dataset_comparison.pdf'));
end

out = struct();
out.summary_dir = summary_dir;
out.metrics_csv = metrics_path;
out.readme = readme_path;
out.recommended = recommendation;
out.variant_outputs = variant_outputs;

fprintf('\nB1 Lorentz tracking optimization complete.\n');
fprintf('  Summary: %s\n', summary_dir);
fprintf('  Recommended: %s (%s), score %.3g\n', recommendation.id, ...
    recommendation.date_tag, recommendation.score);
end


function variants = local_variant_specs(base_tag)
variants = repmat(local_empty_variant(), 1, 6);
variants(1) = local_make_variant('v1_independent', ...
    sprintf('%s_v1_independent', base_tag), ...
    'independent_double_peak', NaN);
variants(2) = local_make_variant('v2_split120', ...
    sprintf('%s_v2_split120', base_tag), ...
    'independent_double_peak', 120);
variants(3) = local_make_variant('v2_split160', ...
    sprintf('%s_v2_split160', base_tag), ...
    'independent_double_peak', 160);
variants(4) = local_make_variant('v2_split220', ...
    sprintf('%s_v2_split220', base_tag), ...
    'independent_double_peak', 220);
variants(5) = local_make_variant('v3_propagated', ...
    sprintf('%s_v3_propagated', base_tag), ...
    'propagated_double_peak', [120 160 220]);
variants(6) = local_make_variant('v4_windowed', ...
    sprintf('%s_v4_windowed', base_tag), ...
    'windowed_branch_tracking', [120 160 220]);
end


function variant = local_empty_variant()
variant = struct('id', '', 'date_tag', '', 'tracking_mode', '', ...
    'fallback_splits', NaN);
end


function variant = local_make_variant(id, tag, tracking_mode, splits)
variant = local_empty_variant();
variant.id = id;
variant.date_tag = tag;
variant.tracking_mode = tracking_mode;
variant.fallback_splits = splits;
end


function local_run_variant(project_root, variant, base_tag)
session_requests = {'590', 'n0', '20w'};
if strcmp(variant.tracking_mode, 'windowed_branch_tracking')
    ref_tag = sprintf('%s_v3_propagated', base_tag);
    for si = 1:numel(session_requests)
        ref_dir = local_session_dir(project_root, session_requests{si}, ref_tag);
        lower = readtable(fullfile(ref_dir, 'b1_double_peak_lower_points.csv'));
        upper = readtable(fullfile(ref_dir, 'b1_double_peak_upper_points.csv'));
        local_call_extraction(session_requests{si}, variant, lower, upper);
    end
else
    local_call_extraction('all', variant, table(), table());
end
end


function local_call_extraction(session_request, variant, lower_ref, upper_ref)
run_b1_double_peak_binning_analysis(session_request, ...
    runFits=false, ...
    outputDateTag=variant.date_tag, ...
    qRangeOverride_Ainv=[-0.15 0.15], ...
    peakModelOverride='lorentz', ...
    trackingMode=variant.tracking_mode, ...
    fallbackSplitCandidatesMeV=variant.fallback_splits, ...
    maxTrackingShiftMeV=180, ...
    trackingWindowHalfWidthMeV=220, ...
    trackingWindowHighQHalfWidthMeV=300, ...
    trackingWindowHighQAbsAinv=0.09, ...
    referenceLowerPoints=lower_ref, ...
    referenceUpperPoints=upper_ref, ...
    waterfallStartMeV=250, ...
    waterfallEndMeV=1600, ...
    waterfallNormMode='area', ...
    waterfallAreaNormWindowMeV=[50 3800], ...
    waterfallResidual=true, ...
    waterfallGain=2, ...
    fitDenoiseMethod='sgolay', ...
    fitDenoiseProfile='adaptive_absq', ...
    fitDenoiseLowWindow=11, ...
    fitDenoiseHighWindow=51, ...
    fitDenoiseQStartAinv=0.07, ...
    fitDenoiseQEndAinv=0.15, ...
    fitDenoiseOrder=3, ...
    highQForceBinAbsAinv=0.09, ...
    binSize=5);
end


function metrics = local_finalize_variant(project_root, variant)
session_requests = {'590', 'n0', '20w'};
metrics = local_empty_metrics_table();
for si = 1:numel(session_requests)
    session_dir = local_session_dir(project_root, session_requests{si}, ...
        variant.date_tag);
    local_copy_variant_named_csv(session_dir, variant.id);
    session_metrics = local_metrics_one_session(session_dir, variant);
    writetable(session_metrics, fullfile(session_dir, sprintf( ...
        'b1_double_peak_lorentz_tracking_%s_tracking_metrics.csv', ...
        variant.id)));
    metrics = [metrics; session_metrics]; %#ok<AGROW>
end
end


function local_copy_variant_named_csv(session_dir, variant_id)
copies = { ...
    'b1_double_peak_combined_q_points.csv', 'combined_points.csv'; ...
    'b1_double_peak_lower_points.csv', 'lower_points.csv'; ...
    'b1_double_peak_upper_points.csv', 'upper_points.csv'; ...
    'b1_double_peak_fit_failures.csv', 'failures.csv'};
for i = 1:size(copies, 1)
    src = fullfile(session_dir, copies{i, 1});
    if isfile(src)
        dst = fullfile(session_dir, sprintf( ...
            'b1_double_peak_lorentz_tracking_%s_%s', ...
            variant_id, copies{i, 2}));
        copyfile(src, dst);
    end
end
end


function metrics = local_metrics_one_session(session_dir, variant)
lower = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_lower_points.csv'));
upper = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_upper_points.csv'));
failures = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_fit_failures.csv'));

session_key = local_session_key_from_dir(session_dir);
failure_count = height(failures);
lower_row = local_branch_metric_row(variant, session_key, 'lower', lower, ...
    failure_count);
upper_row = local_branch_metric_row(variant, session_key, 'upper', upper, ...
    failure_count);
metrics = [lower_row; upper_row];
end


function row = local_branch_metric_row(variant, session_key, branch, points, ...
    failure_count)
n_points = height(points);
jump_count = 0;
edge_count = 0;
if n_points >= 1 && any(strcmp(points.Properties.VariableNames, 'energy_meV'))
    edge_count = sum(points.energy_meV < 650 | points.energy_meV > 1900, ...
        'omitnan');
end
if n_points >= 2 && any(strcmp(points.Properties.VariableNames, 'q_Ainv')) && ...
        any(strcmp(points.Properties.VariableNames, 'energy_meV'))
    [~, order] = sort(points.q_Ainv);
    energy = points.energy_meV(order);
    jump_count = sum(abs(diff(energy)) > 250, 'omitnan');
end

weight = 1;
if strcmp(session_key, 'no_PL2_20w_2film')
    weight = 2;
end
raw_score = 10 * jump_count + 3 * edge_count + 5 * failure_count;
weighted_score = weight * raw_score;
row = table({variant.id}, {variant.date_tag}, {session_key}, {branch}, ...
    n_points, jump_count, edge_count, failure_count, weight, ...
    raw_score, weighted_score, ...
    'VariableNames', {'variant_id', 'date_tag', 'session_key', ...
    'branch', 'n_points', 'large_jump_count', 'edge_count', ...
    'failure_count', 'session_weight', 'raw_score', 'weighted_score'});
end


function score = local_variant_score(metrics)
if isempty(metrics) || height(metrics) == 0
    score = Inf;
else
    score = sum(metrics.weighted_score, 'omitnan');
end
end


function recommendation = local_choose_recommendation(variant_outputs)
scores = [variant_outputs.score];
[~, idx] = min(scores);
recommendation = variant_outputs(idx);
end


function local_write_readme(path, base_tag, variants, metrics, recommendation)
fid = fopen(path, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# B1 Lorentz 双峰追踪优化记录\n\n');
fprintf(fid, '- 基准标签：`%s`\n', base_tag);
fprintf(fid, '- 强制峰型：Lorentz，不使用 Fano，不允许退回单峰。\n');
fprintf(fid, '- q 规则：signed-q，默认范围 `[-0.15, 0.15] A^-1`，不做 `+q/-q` 平均。\n');
fprintf(fid, '- 显示增强：50-3800 meV 面积归一化，250-1600 meV 显示，residual，gain=2。\n');
fprintf(fid, '- 本轮只做提取与追踪诊断，不运行 lower/upper physical fit。\n\n');

fprintf(fid, '## 版本\n\n');
for i = 1:numel(variants)
    fprintf(fid, '- `%s`：`%s`，tracking=`%s`，fallback split=`%s` meV\n', ...
        variants(i).id, variants(i).date_tag, variants(i).tracking_mode, ...
        mat2str(variants(i).fallback_splits));
end

fprintf(fid, '\n## 推荐版本\n\n');
fprintf(fid, '推荐 `%s`（`%s`），总分 %.3g。\n\n', recommendation.id, ...
    recommendation.date_tag, recommendation.score);
fprintf(fid, '评分依据是 lower/upper 分支沿 signed-q 的大跳点数量、');
fprintf(fid, '贴近边界点数量、失败点数量；20w 数据权重加倍，因为它噪声最强。\n\n');

fprintf(fid, '## Metrics 摘要\n\n');
fprintf(fid, '| version | session | branch | points | jumps | edge | failures | weighted score |\n');
fprintf(fid, '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |\n');
for i = 1:height(metrics)
    fprintf(fid, '| %s | %s | %s | %d | %d | %d | %d | %.3g |\n', ...
        local_cell_text(metrics.variant_id, i), ...
        local_cell_text(metrics.session_key, i), ...
        local_cell_text(metrics.branch, i), ...
        metrics.n_points(i), metrics.large_jump_count(i), ...
        metrics.edge_count(i), metrics.failure_count(i), ...
        metrics.weighted_score(i));
end
end


function tbl = local_read_table_if_exists(path)
if isfile(path)
    tbl = readtable(path);
else
    tbl = table();
end
end


function metrics = local_empty_metrics_table()
metrics = table('Size', [0 11], ...
    'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'variant_id', 'date_tag', 'session_key', ...
    'branch', 'n_points', 'large_jump_count', 'edge_count', ...
    'failure_count', 'session_weight', 'raw_score', 'weighted_score'});
end


function session_dir = local_session_dir(project_root, session_request, tag)
switch session_request
    case '590'
        folder = sprintf('590_gui_history_area_260506_b1_double_peak_binning_%s', tag);
    case 'n0'
        folder = sprintf('n0_PL2_10w_gui_history_area_260506_b1_double_peak_binning_%s', tag);
    case '20w'
        folder = sprintf(['no_PL2_20w_2film_gui_history_area_260506_', ...
            'highq_refined_b1_double_peak_binning_%s'], tag);
    otherwise
        error('run_b1_lorentz_tracking_optimization:UnknownSession', ...
            'Unknown session "%s".', session_request);
end
session_dir = fullfile(project_root, 'paper_results', folder);
end


function key = local_session_key_from_dir(session_dir)
[~, name] = fileparts(session_dir);
if startsWith(name, '590_')
    key = '590_PL2_10w';
elseif startsWith(name, 'n0_')
    key = 'n0_PL2_10w_repeat';
elseif startsWith(name, 'no_PL2_20w')
    key = 'no_PL2_20w_2film';
else
    key = name;
end
end


function text = local_cell_text(value, row)
if iscell(value)
    text = char(value{row});
elseif isstring(value)
    text = char(value(row));
else
    text = char(string(value(row)));
end
end
