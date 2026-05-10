function out = run_b1_lorentz_tracking_balanced_highqbin7_sg71(options)
%RUN_B1_LORENTZ_TRACKING_BALANCED_HIGHQBIN7_SG71 Stronger high-q B1 tracking.
%
% This diagnostic keeps the Lorentz double-peak model and signed-q rules. It
% first creates a propagated trend with stronger high-q binning/SG denoise,
% then runs windowed branch tracking from that trend. It does not run any
% physical fit.

arguments
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

base_tag = '260510_lorentz_tracking_balanced_highqbin7_sg71';
v5a_tag = [base_tag '_v5a_propagated'];
v5b_tag = [base_tag '_v5b_windowed'];
results_root = fullfile(project_root, 'paper_results');
summary_dir = fullfile(results_root, ...
    'b1_lorentz_tracking_optimization_260510_balanced_highqbin7_sg71');
if ~isfolder(summary_dir)
    mkdir(summary_dir);
end

if options.runExtraction
    fprintf('\n=== B1 Lorentz balanced high-q bin7 SG71 V5a propagated ===\n');
    run_b1_double_peak_binning_analysis('all', ...
        runFits=false, ...
        outputDateTag=v5a_tag, ...
        qRangeOverride_Ainv=[-0.15 0.15], ...
        peakModelOverride='lorentz', ...
        trackingMode='propagated_double_peak', ...
        fallbackSplitCandidatesMeV=[120 160 220], ...
        maxTrackingShiftMeV=180, ...
        waterfallStartMeV=250, ...
        waterfallEndMeV=1600, ...
        waterfallNormMode='area', ...
        waterfallAreaNormWindowMeV=[50 3800], ...
        waterfallResidual=true, ...
        waterfallGain=2, ...
        fitDenoiseMethod='sgolay', ...
        fitDenoiseProfile='adaptive_absq', ...
        fitDenoiseLowWindow=11, ...
        fitDenoiseHighWindow=71, ...
        fitDenoiseQStartAinv=0.07, ...
        fitDenoiseQEndAinv=0.15, ...
        fitDenoiseOrder=3, ...
        highQForceBinAbsAinv=0.08, ...
        binSize=7);

    fprintf('\n=== B1 Lorentz balanced high-q bin7 SG71 V5b windowed ===\n');
    sessions = {'590', 'n0', '20w'};
    for si = 1:numel(sessions)
        ref_dir = local_session_dir(project_root, sessions{si}, v5a_tag);
        lower = readtable(fullfile(ref_dir, 'b1_double_peak_lower_points.csv'));
        upper = readtable(fullfile(ref_dir, 'b1_double_peak_upper_points.csv'));
        run_b1_double_peak_binning_analysis(sessions{si}, ...
            runFits=false, ...
            outputDateTag=v5b_tag, ...
            qRangeOverride_Ainv=[-0.15 0.15], ...
            peakModelOverride='lorentz', ...
            trackingMode='windowed_branch_tracking', ...
            fallbackSplitCandidatesMeV=[120 160 220], ...
            maxTrackingShiftMeV=180, ...
            trackingWindowHalfWidthMeV=220, ...
            trackingWindowHighQHalfWidthMeV=300, ...
            trackingWindowHighQAbsAinv=0.09, ...
            referenceLowerPoints=lower, ...
            referenceUpperPoints=upper, ...
            waterfallStartMeV=250, ...
            waterfallEndMeV=1600, ...
            waterfallNormMode='area', ...
            waterfallAreaNormWindowMeV=[50 3800], ...
            waterfallResidual=true, ...
            waterfallGain=2, ...
            fitDenoiseMethod='sgolay', ...
            fitDenoiseProfile='adaptive_absq', ...
            fitDenoiseLowWindow=11, ...
            fitDenoiseHighWindow=71, ...
            fitDenoiseQStartAinv=0.07, ...
            fitDenoiseQEndAinv=0.15, ...
            fitDenoiseOrder=3, ...
            highQForceBinAbsAinv=0.08, ...
            binSize=7);
    end
end

metrics = local_collect_metrics(project_root, base_tag, {v5a_tag, v5b_tag}, ...
    {'v5a_propagated', 'v5b_windowed'});
metrics_path = fullfile(summary_dir, 'b1_lorentz_tracking_metrics.csv');
writetable(metrics, metrics_path);

overlay = struct('png', "", 'pdf', "");
if options.runOverlay
    run_b1_double_peak_waterfall_extraction_overlay(dateTag=v5a_tag);
    overlay = run_b1_double_peak_waterfall_extraction_overlay(dateTag=v5b_tag);
end
if strlength(string(overlay.png)) > 0 && isfile(overlay.png)
    copyfile(overlay.png, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_overlay_three_dataset_comparison.png'));
end
if strlength(string(overlay.pdf)) > 0 && isfile(overlay.pdf)
    copyfile(overlay.pdf, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_overlay_three_dataset_comparison.pdf'));
end

old_v4_score = 537;
v5b_score = local_variant_score(metrics, 'v5b_windowed');
if v5b_score < old_v4_score
    recommendation = 'v5b_windowed';
else
    recommendation = 'keep_previous_v4_windowed';
end
readme_path = fullfile(summary_dir, 'README_lorentz_tracking_optimization.md');
local_write_readme(readme_path, metrics, old_v4_score, v5b_score, ...
    recommendation);

out = struct();
out.summary_dir = summary_dir;
out.metrics_csv = metrics_path;
out.readme = readme_path;
out.v5a_tag = v5a_tag;
out.v5b_tag = v5b_tag;
out.v5b_score = v5b_score;
out.old_v4_score = old_v4_score;
out.recommendation = recommendation;

fprintf('\nB1 balanced high-q bin7 SG71 tracking complete.\n');
fprintf('  Summary: %s\n', summary_dir);
fprintf('  V5b score: %.3g; previous V4 score: %.3g; recommendation: %s\n', ...
    v5b_score, old_v4_score, recommendation);
end


function metrics = local_collect_metrics(project_root, base_tag, tags, ids)
metrics = local_empty_metrics_table();
session_requests = {'590', 'n0', '20w'};
for vi = 1:numel(tags)
    variant = struct('id', ids{vi}, 'date_tag', tags{vi});
    for si = 1:numel(session_requests)
        session_dir = local_session_dir(project_root, session_requests{si}, ...
            tags{vi});
        local_copy_variant_named_csv(session_dir, base_tag, ids{vi});
        session_metrics = local_metrics_one_session(session_dir, variant);
        writetable(session_metrics, fullfile(session_dir, sprintf( ...
            'b1_double_peak_lorentz_tracking_%s_%s_tracking_metrics.csv', ...
            base_tag, ids{vi})));
        metrics = [metrics; session_metrics]; %#ok<AGROW>
    end
end
end


function local_copy_variant_named_csv(session_dir, base_tag, variant_id)
copies = { ...
    'b1_double_peak_combined_q_points.csv', 'combined_points.csv'; ...
    'b1_double_peak_lower_points.csv', 'lower_points.csv'; ...
    'b1_double_peak_upper_points.csv', 'upper_points.csv'; ...
    'b1_double_peak_fit_failures.csv', 'failures.csv'};
for i = 1:size(copies, 1)
    src = fullfile(session_dir, copies{i, 1});
    if isfile(src)
        dst = fullfile(session_dir, sprintf( ...
            'b1_double_peak_lorentz_tracking_%s_%s_%s', ...
            base_tag, variant_id, copies{i, 2}));
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


function score = local_variant_score(metrics, variant_id)
rows = strcmp(metrics.variant_id, variant_id);
if ~any(rows)
    score = Inf;
else
    score = sum(metrics.weighted_score(rows), 'omitnan');
end
end


function local_write_readme(path, metrics, old_v4_score, v5b_score, ...
    recommendation)
fid = fopen(path, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# B1 Lorentz balanced high-q bin7 SG71 tracking\n\n');
fprintf(fid, '- 强制峰型：Lorentz；不使用 Fano；不允许退回单峰。\n');
fprintf(fid, '- q 规则：signed-q；默认范围 `[-0.15, 0.15] A^-1`；不做 `+q/-q` 平均。\n');
fprintf(fid, '- 低 q 保护：`|q| <= 0.05 A^-1` 不强制 bin。\n');
fprintf(fid, '- 高 q 强化：`|q| >= 0.08 A^-1` 强制 7 点相邻 q-bin。\n');
fprintf(fid, '- 降噪：adaptive SG，low window=11，high window=71，过渡区 `0.07-0.15 A^-1`。\n');
fprintf(fid, '- waterfall：50-3800 meV area norm，250-1600 meV display，residual，gain=2。\n');
fprintf(fid, '- 本轮只做提取诊断，不运行 physical fit。\n\n');

fprintf(fid, '## Recommendation\n\n');
fprintf(fid, '- Previous V4 score: %.3g\n', old_v4_score);
fprintf(fid, '- V5b score: %.3g\n', v5b_score);
fprintf(fid, '- Recommendation: `%s`\n\n', recommendation);

fprintf(fid, '## Metrics\n\n');
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
        error('run_b1_lorentz_tracking_balanced_highqbin7_sg71:UnknownSession', ...
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
