function out = run_b1_lorentz_tracking_upper_rapidrise_plateau(options)
%RUN_B1_LORENTZ_TRACKING_UPPER_RAPIDRISE_PLATEAU V15 upper trend retry.
%
% This retry keeps the clean 300 meV start window, mandatory Lorentz double
% peaks, signed q, no +q/-q averaging, no single-peak fallback, and no
% physical fit. It adds an upper-branch trend shape cost: small-q rapid
% rise is allowed, while isolated jumps and high-q wandering are penalized.

arguments
    options.targetScore (1,1) double = 300
    options.runExtraction (1,1) logical = true
    options.runOverlay (1,1) logical = true
    options.manualAnchorPolicy {mustBeTextScalar} = "template_only"
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
addpath(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end
if ~strcmp(char(string(options.manualAnchorPolicy)), 'template_only')
    error('run_b1_lorentz_tracking_upper_rapidrise_plateau:ManualAnchors', ...
        'manualAnchorPolicy must be template_only.');
end

summary_dir = fullfile(project_root, 'paper_results', ...
    'b1_lorentz_tracking_optimization_260510_upper_rapidrise_plateau');
if ~isfolder(summary_dir)
    mkdir(summary_dir);
end

variant = local_variant();
fprintf('\n=== B1 upper rapidrise-plateau retry %s ===\n', variant.id);
if options.runExtraction
    local_run_variant(variant);
end
[metrics, upper_diag, path_selection] = local_finalize_variant( ...
    project_root, summary_dir, variant);
score = local_variant_score(metrics, variant.id);
gate_passed = score <= options.targetScore;
overlay = struct('png', "", 'pdf', "");
if options.runOverlay
    overlay = run_b1_double_peak_waterfall_extraction_overlay( ...
        dateTag=variant.final_tag);
end
best = local_variant_output(variant, score, gate_passed, overlay);

metrics_path = fullfile(summary_dir, ...
    'b1_lorentz_tracking_upper_rapidrise_plateau_metrics.csv');
writetable(metrics, metrics_path);
upper_diag_path = fullfile(summary_dir, ...
    'b1_lorentz_tracking_upper_rapidrise_plateau_upper_jump_gamma_R2.csv');
writetable(upper_diag, upper_diag_path);
path_selection_path = fullfile(summary_dir, ...
    'upper_rapidrise_plateau_path_selection.csv');
writetable(path_selection, path_selection_path);
trend_diag_path = fullfile(summary_dir, 'upper_trend_diagnostics.csv');
writetable(local_upper_trend_diagnostics(path_selection), trend_diag_path);
variant_path = fullfile(summary_dir, ...
    'b1_lorentz_tracking_upper_rapidrise_plateau_variants.csv');
writetable(local_variant_outputs_table(best), variant_path);

if strlength(best.overlay_png) > 0 && isfile(best.overlay_png)
    copyfile(best.overlay_png, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_upper_rapidrise_plateau_overlay_three_dataset_comparison.png'));
end
if strlength(best.overlay_pdf) > 0 && isfile(best.overlay_pdf)
    copyfile(best.overlay_pdf, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_upper_rapidrise_plateau_overlay_three_dataset_comparison.pdf'));
end
trend_overlay = local_plot_trend_overlay(project_root, summary_dir, variant);

if gate_passed
    local_write_final_candidate_marker(summary_dir, best, options.targetScore);
else
    local_write_manual_anchor_template(project_root, summary_dir, best);
end
readme_path = fullfile(summary_dir, ...
    'README_lorentz_tracking_upper_rapidrise_plateau.md');
local_write_readme(readme_path, metrics, best, options.targetScore, gate_passed);
local_register_outputs(project_root, summary_dir, best, true);
local_update_results_index(project_root, summary_dir, best, gate_passed);

out = struct();
out.summary_dir = summary_dir;
out.metrics_csv = metrics_path;
out.upper_diagnostics_csv = upper_diag_path;
out.path_selection_csv = path_selection_path;
out.upper_trend_diagnostics_csv = trend_diag_path;
out.upper_trend_overlay_png = trend_overlay.png;
out.upper_trend_overlay_pdf = trend_overlay.pdf;
out.variants_csv = variant_path;
out.readme = readme_path;
out.best = best;
out.best_score = best.score;
out.metric_gate_passed = gate_passed;

fprintf('\nB1 upper rapidrise-plateau retry complete.\n');
fprintf('  Summary: %s\n', summary_dir);
fprintf('  Best: %s (%s), score %.3g, gate %d\n', ...
    best.id, best.final_tag, best.score, gate_passed);
fprintf('  no physical fit was run\n');
end


function variant = local_variant()
variant = struct();
variant.id = 'v15_upper_rapidrise_plateau';
variant.final_tag = '260510_lorentz_tracking_v15_upper_rapidrise_plateau';
variant.b1_energy_window = [300 1800];
variant.bin_default = 9;
variant.sg_default = 91;
variant.bin_20w = 11;
variant.sg_20w = 111;
variant.candidate_for_recommendation = true;
end


function local_run_variant(variant)
sessions = {'590', 'n0', '20w'};
for si = 1:numel(sessions)
    session_request = sessions{si};
    [bin_size, sg_high] = local_session_parameters(session_request, variant);
    run_b1_double_peak_binning_analysis(session_request, ...
        runFits=false, ...
        outputDateTag=variant.final_tag, ...
        qRangeOverride_Ainv=[-0.15 0.15], ...
        b1EnergyWindowOverrideMeV=[300 1800], ...
        peakModelOverride='lorentz', ...
        trackingMode='ridge_guided_candidate_path', ...
        fallbackSplitCandidatesMeV=[40 60 80 120 160 220 280 340 420], ...
        maxTrackingShiftMeV=260, ...
        trackingWindowHalfWidthMeV=320, ...
        trackingWindowHighQHalfWidthMeV=420, ...
        trackingWindowHighQAbsAinv=0.09, ...
        trackingWindowInvalidFallback='independent_double_peak', ...
        upperQualityRetry=true, ...
        upperMaxGammaOverE=1.4, ...
        upperMaxGammaMeV=1600, ...
        upperRetryWindowHalfWidthMeV=120, ...
        upperRetryHighQHalfWidthMeV=180, ...
        upperRetryHighQAbsAinv=0.09, ...
        ridgeSmoothWindow=11, ...
        candidatePathUpperMediumJumpThresholdMeV=local_upper_threshold(session_request), ...
        candidatePathLargeJumpThresholdMeV=250, ...
        candidatePathEdgeEnergyMeV=[650 1900], ...
        candidatePathUpperTrendMode='rapidrise_plateau', ...
        candidatePathUpperTrendAnchorQAbsAinv=0.005, ...
        candidatePathUpperTrendSmallQAbsAinv=0.02, ...
        candidatePathUpperTrendPlateauQAbsAinv=0.06, ...
        waterfallStartMeV=250, ...
        waterfallEndMeV=1600, ...
        waterfallNormMode='area', ...
        waterfallAreaNormWindowMeV=[50 3800], ...
        waterfallResidual=true, ...
        waterfallGain=2, ...
        fitDenoiseMethod='sgolay', ...
        fitDenoiseProfile='adaptive_absq', ...
        fitDenoiseLowWindow=11, ...
        fitDenoiseHighWindow=sg_high, ...
        fitDenoiseQStartAinv=0.07, ...
        fitDenoiseQEndAinv=0.15, ...
        fitDenoiseOrder=3, ...
        fitMinPeakAmplitudeFraction=0, ...
        highQForceBinAbsAinv=0.075, ...
        lowQNoBinAbsAinv=0.05, ...
        binSize=bin_size, ...
        enableJumpRepair=false);
end
end


function [bin_size, sg_high] = local_session_parameters(session_request, variant)
if strcmp(session_request, '20w')
    bin_size = variant.bin_20w;
    sg_high = variant.sg_20w;
else
    bin_size = variant.bin_default;
    sg_high = variant.sg_default;
end
end


function threshold = local_upper_threshold(session_request)
if strcmp(session_request, '20w')
    threshold = 120;
else
    threshold = 150;
end
end


function [metrics, upper_diag, path_selection] = local_finalize_variant(project_root, ...
    summary_dir, variant)
metrics = local_empty_metrics_table();
upper_diag = table();
path_selection = table();
sessions = {'590', 'n0', '20w'};
for si = 1:numel(sessions)
    session_dir = local_session_dir(project_root, sessions{si}, ...
        variant.final_tag);
    if ~isfolder(session_dir)
        mkdir(session_dir);
    end
    local_copy_variant_named_csv(session_dir, variant.id);
    [rows, diag_rows] = local_metrics_one_session(session_dir, variant);
    path_rows = local_path_selection_one_session(session_dir, variant);
    writetable(rows, fullfile(session_dir, sprintf( ...
        'b1_double_peak_lorentz_tracking_upper_rapidrise_plateau_%s_tracking_metrics.csv', ...
        variant.id)));
    writetable(diag_rows, fullfile(session_dir, sprintf( ...
        'b1_double_peak_lorentz_tracking_upper_rapidrise_plateau_%s_upper_jump_gamma_R2.csv', ...
        variant.id)));
    writetable(path_rows, fullfile(session_dir, sprintf( ...
        'b1_double_peak_lorentz_tracking_upper_rapidrise_plateau_%s_path_selection.csv', ...
        variant.id)));
    metrics = [metrics; rows]; %#ok<AGROW>
    upper_diag = [upper_diag; diag_rows]; %#ok<AGROW>
    path_selection = [path_selection; path_rows]; %#ok<AGROW>
end
if ~isempty(upper_diag)
    writetable(upper_diag, fullfile(summary_dir, sprintf( ...
        'b1_double_peak_lorentz_tracking_upper_rapidrise_plateau_%s_upper_jump_gamma_R2.csv', ...
        variant.id)));
end
end


function local_copy_variant_named_csv(session_dir, variant_id)
copies = { ...
    'b1_double_peak_combined_q_points.csv', 'combined_points.csv'; ...
    'b1_double_peak_lower_points.csv', 'lower_points.csv'; ...
    'b1_double_peak_upper_points.csv', 'upper_points.csv'; ...
    'b1_double_peak_fit_failures.csv', 'failures.csv'; ...
    'b1_double_peak_repair_log.csv', 'repair_log.csv'; ...
    'b1_double_peak_exclusion_points.csv', 'exclusion_points.csv'; ...
    'b1_double_peak_binning_map.csv', 'binning_map.csv'; ...
    'b1_double_peak_lorentz_candidate_points.csv', 'candidate_points.csv'; ...
    'b1_double_peak_candidate_path_selection.csv', 'path_selection.csv'};
for i = 1:size(copies, 1)
    src = fullfile(session_dir, copies{i, 1});
    if isfile(src)
        dst = fullfile(session_dir, sprintf( ...
            'b1_double_peak_lorentz_tracking_upper_rapidrise_plateau_%s_%s', ...
            variant_id, copies{i, 2}));
        copyfile(src, dst);
    end
end
end


function path_rows = local_path_selection_one_session(session_dir, variant)
path_rows = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_candidate_path_selection.csv'));
if isempty(path_rows) || height(path_rows) == 0
    return
end
session_key = local_session_key_from_dir(session_dir);
n = height(path_rows);
path_rows = addvars(path_rows, repmat({variant.id}, n, 1), ...
    repmat({variant.final_tag}, n, 1), repmat({session_key}, n, 1), ...
    'Before', 1, 'NewVariableNames', ...
    {'variant_id', 'date_tag', 'session_key'});
end


function [metrics, upper_diag] = local_metrics_one_session(session_dir, variant)
lower = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_lower_points.csv'));
upper = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_upper_points.csv'));
failures = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_fit_failures.csv'));
repair_log = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_repair_log.csv'));
exclusions = local_read_table_if_exists(fullfile(session_dir, ...
    'b1_double_peak_exclusion_points.csv'));

session_key = local_session_key_from_dir(session_dir);
failure_count = height(failures);
repair_count = height(repair_log);
exclusion_count = max(height(exclusions), failure_count);
metrics = [ ...
    local_branch_metric_row(variant, session_key, 'lower', lower, ...
    failure_count, repair_count, exclusion_count); ...
    local_branch_metric_row(variant, session_key, 'upper', upper, ...
    failure_count, repair_count, exclusion_count)];
upper_diag = local_upper_diagnostics(variant, session_key, upper);
end


function row = local_branch_metric_row(variant, session_key, branch, ...
    points, failure_count, repair_count, exclusion_count)
n_points = height(points);
jump_count = 0;
edge_count = 0;
if n_points >= 1 && ismember('energy_meV', points.Properties.VariableNames)
    edge_count = sum(points.energy_meV < 650 | points.energy_meV > 1900, ...
        'omitnan');
end
if n_points >= 2 && all(ismember({'q_Ainv', 'energy_meV'}, ...
        points.Properties.VariableNames))
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
upper = local_empty_upper_metric_values();
if strcmp(branch, 'upper')
    upper = b1_double_peak_upper_stability_metrics(points, session_key);
end
row = table({variant.id}, {variant.final_tag}, ...
    variant.candidate_for_recommendation, ...
    variant.b1_energy_window(1), variant.b1_energy_window(2), ...
    {session_key}, {branch}, ...
    n_points, jump_count, edge_count, failure_count, repair_count, ...
    exclusion_count, weight, raw_score, weighted_score, ...
    upper.upper_medium_jump_threshold_meV(1), ...
    upper.upper_large_jump_threshold_meV(1), ...
    upper.upper_medium_jump_count(1), ...
    upper.upper_large_jump_count(1), ...
    upper.upper_max_abs_jump_meV(1), ...
    upper.upper_overbroad_count(1), ...
    upper.upper_low_r2_count(1), ...
    upper.upper_stability_score(1), ...
    'VariableNames', {'variant_id', 'date_tag', ...
    'candidate_for_recommendation', 'b1_energy_min_meV', ...
    'b1_energy_max_meV', 'session_key', 'branch', ...
    'n_points', 'large_jump_count', 'edge_count', ...
    'failure_count', 'repair_count', 'exclusion_count', ...
    'session_weight', 'raw_score', 'weighted_score', ...
    'upper_medium_jump_threshold_meV', ...
    'upper_large_jump_threshold_meV', 'upper_medium_jump_count', ...
    'upper_large_jump_count', 'upper_max_abs_jump_meV', ...
    'upper_overbroad_count', 'upper_low_r2_count', ...
    'upper_stability_score'});
end


function values = local_empty_upper_metric_values()
values = table({''}, 0, NaN, 180, 0, 0, NaN, 0, 1.4, 1600, 0, 0.50, ...
    NaN, NaN, 0, 'VariableNames', {'session_key', 'n_upper_points', ...
    'upper_medium_jump_threshold_meV', ...
    'upper_large_jump_threshold_meV', 'upper_medium_jump_count', ...
    'upper_large_jump_count', 'upper_max_abs_jump_meV', ...
    'upper_overbroad_count', 'upper_max_gamma_over_E', ...
    'upper_max_gamma_meV', 'upper_low_r2_count', ...
    'upper_low_r2_threshold', 'upper_median_gamma_over_E', ...
    'upper_median_R2', 'upper_stability_score'});
end


function diag_rows = local_upper_diagnostics(variant, session_key, upper)
if isempty(upper) || height(upper) == 0
    diag_rows = table();
    return
end
[~, order] = sort(double(upper.q_Ainv));
upper = upper(order, :);
n = height(upper);
delta = [NaN; abs(diff(double(upper.energy_meV)))];
threshold = b1_double_peak_upper_stability_metrics(upper, session_key). ...
    upper_medium_jump_threshold_meV(1);
gamma_over_E = double(upper.gamma_meV) ./ max(abs(double(upper.energy_meV)), eps);
is_medium_jump = delta > threshold;
is_large_jump = delta > 180;
is_overbroad = gamma_over_E > 1.4 | double(upper.gamma_meV) > 1600;
source_mode = local_text_column(upper, 'source_mode');
source_q_Ainv = local_text_column(upper, 'source_q_Ainv');
repair_source = local_text_column(upper, 'repair_source');
repair_detail = local_text_column(upper, 'repair_detail');
diag_rows = table(repmat({variant.id}, n, 1), ...
    repmat({variant.final_tag}, n, 1), ...
    repmat(variant.candidate_for_recommendation, n, 1), ...
    repmat(variant.b1_energy_window(1), n, 1), ...
    repmat(variant.b1_energy_window(2), n, 1), ...
    repmat({session_key}, n, 1), ...
    upper.q_Ainv, upper.q_abs_Ainv, upper.energy_meV, delta, ...
    upper.gamma_meV, gamma_over_E, upper.R2, is_medium_jump, ...
    is_large_jump, is_overbroad, source_mode, ...
    upper.source_q_count, source_q_Ainv, repair_source, repair_detail, ...
    'VariableNames', {'variant_id', 'date_tag', ...
    'candidate_for_recommendation', 'b1_energy_min_meV', ...
    'b1_energy_max_meV', 'session_key', 'q_Ainv', 'q_abs_Ainv', ...
    'energy_meV', 'delta_from_previous_meV', 'gamma_meV', ...
    'gamma_over_E', 'R2', 'is_upper_medium_jump', ...
    'is_upper_large_jump_180', 'is_overbroad_upper_peak', ...
    'source_mode', 'source_q_count', 'source_q_Ainv', ...
    'repair_source', 'repair_detail'});
end


function diag = local_upper_trend_diagnostics(path_selection)
if isempty(path_selection) || height(path_selection) == 0
    diag = table();
    return
end
required = {'variant_id', 'date_tag', 'session_key', 'q_Ainv', ...
    'q_abs_Ainv', 'candidate_source', 'upper_energy_meV', ...
    'upper_delta_from_previous_meV', 'upper_trend_region', ...
    'upper_trend_mode', 'transition_cost', 'upper_trend_penalty_meV', ...
    'node_cost', 'path_cost'};
keep = intersect(required, path_selection.Properties.VariableNames, ...
    'stable');
diag = path_selection(:, keep);
if ismember('upper_delta_from_previous_meV', diag.Properties.VariableNames)
    delta = abs(double(diag.upper_delta_from_previous_meV));
    diag.is_upper_medium_jump_120 = delta > 120;
    diag.is_upper_large_jump_250 = delta > 250;
end
end


function score = local_variant_score(metrics, variant_id)
rows = strcmp(metrics.variant_id, variant_id);
if ~any(rows)
    score = Inf;
else
    score = sum(metrics.weighted_score(rows), 'omitnan');
end
end


function output = local_variant_output(variant, score, gate, overlay)
upper_score = NaN;
try
    project_root = bisb_find_project_root(fileparts(mfilename('fullpath')));
    metrics = table();
    sessions = {'590', 'n0', '20w'};
    for si = 1:numel(sessions)
        session_dir = local_session_dir(project_root, sessions{si}, ...
            variant.final_tag);
        upper = local_read_table_if_exists(fullfile(session_dir, ...
            'b1_double_peak_upper_points.csv'));
        metrics = [metrics; b1_double_peak_upper_stability_metrics( ...
            upper, local_session_key_from_dir(session_dir))]; %#ok<AGROW>
    end
    upper_score = sum(metrics.upper_stability_score, 'omitnan');
catch
    upper_score = Inf;
end
output = struct('id', variant.id, 'final_tag', variant.final_tag, ...
    'candidate_for_recommendation', variant.candidate_for_recommendation, ...
    'b1_energy_min_meV', variant.b1_energy_window(1), ...
    'b1_energy_max_meV', variant.b1_energy_window(2), ...
    'score', score, 'upper_stability_score', upper_score, ...
    'metric_gate_passed', gate, 'overlay_png', string(overlay.png), ...
    'overlay_pdf', string(overlay.pdf));
end


function tbl = local_variant_outputs_table(output)
tbl = table({output.id}, {output.final_tag}, ...
    output.candidate_for_recommendation, output.b1_energy_min_meV, ...
    output.b1_energy_max_meV, output.score, ...
    output.upper_stability_score, output.metric_gate_passed, ...
    string(output.overlay_png), string(output.overlay_pdf), ...
    'VariableNames', {'variant_id', 'date_tag', ...
    'candidate_for_recommendation', 'b1_energy_min_meV', ...
    'b1_energy_max_meV', 'score', 'upper_stability_score', ...
    'metric_gate_passed', 'overlay_png', 'overlay_pdf'});
end


function local_write_final_candidate_marker(summary_dir, best, target_score)
path = fullfile(summary_dir, ...
    'UPPER_RAPIDRISE_PLATEAU_METRIC_GATE_PASSED_NO_FIT.txt');
fid = fopen(path, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Metric gate passed. no physical fit was run.\n');
fprintf(fid, 'variant=%s\n', best.id);
fprintf(fid, 'date_tag=%s\n', best.final_tag);
fprintf(fid, 'score=%.12g\n', best.score);
fprintf(fid, 'target_score=%.12g\n', target_score);
end


function local_write_manual_anchor_template(project_root, summary_dir, best)
sessions = {'590', 'n0', '20w'};
anchor = table('Size', [0 9], ...
    'VariableTypes', {'cell', 'double', 'cell', 'double', 'double', ...
    'double', 'double', 'cell', 'cell'}, ...
    'VariableNames', {'session_request', 'q_Ainv', 'branch_label', ...
    'manual_lower_min_meV', 'manual_lower_max_meV', ...
    'manual_upper_min_meV', 'manual_upper_max_meV', 'reason', 'detail'});
for si = 1:numel(sessions)
    session_dir = local_session_dir(project_root, sessions{si}, best.final_tag);
    failures = local_read_table_if_exists(fullfile(session_dir, ...
        'b1_double_peak_fit_failures.csv'));
    if all(ismember({'q_Ainv', 'status', 'detail'}, ...
            failures.Properties.VariableNames))
        for i = 1:height(failures)
            anchor = [anchor; table(sessions(si), failures.q_Ainv(i), ...
                {'both'}, NaN, NaN, NaN, NaN, failures.status(i), ...
                failures.detail(i), ...
                'VariableNames', anchor.Properties.VariableNames)]; %#ok<AGROW>
        end
    end
end
writetable(anchor, fullfile(summary_dir, ...
    'b1_double_peak_manual_anchor_template_upper_rapidrise_plateau.csv'));
end


function paths = local_plot_trend_overlay(project_root, summary_dir, variant)
paths = struct('png', "", 'pdf', "");
sessions = {'590', 'n0', '20w'};
labels = {'590 10w defocus 1film', 'n0 10w defocus repeat 1film', ...
    '20w defocus 2film'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1650 520]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(layout, 'B1 upper rapid-rise / plateau candidate path');
for si = 1:numel(sessions)
    ax = nexttile(layout);
    session_dir = local_session_dir(project_root, sessions{si}, ...
        variant.final_tag);
    candidates = local_read_table_if_exists(fullfile(session_dir, ...
        'b1_double_peak_lorentz_candidate_points.csv'));
    selected = local_read_table_if_exists(fullfile(session_dir, ...
        'b1_double_peak_candidate_path_selection.csv'));
    hold(ax, 'on');
    local_patch_q_regions(ax, [-0.15 0.15], [300 1800]);
    if ~isempty(candidates) && height(candidates) > 0
        scatter(ax, candidates.q_Ainv, candidates.upper_energy_meV, 12, ...
            [0.65 0.65 0.65], 'filled', 'MarkerFaceAlpha', 0.20);
    end
    if ~isempty(selected) && height(selected) > 0
        [~, order] = sort(double(selected.q_Ainv));
        selected = selected(order, :);
        plot(ax, selected.q_Ainv, selected.upper_energy_meV, ...
            '-o', 'Color', [0.85 0.00 0.45], 'MarkerFaceColor', ...
            [0.85 0.00 0.45], 'LineWidth', 1.6, 'MarkerSize', 4);
        local_mark_rejected_candidates(ax, candidates, selected);
    end
    title(ax, labels{si}, 'Interpreter', 'none');
    xlabel(ax, 'signed q (A^{-1})');
    ylabel(ax, 'upper energy (meV)');
    xlim(ax, [-0.15 0.15]);
    ylim(ax, [300 1800]);
    grid(ax, 'on');
    box(ax, 'on');
end
paths.png = string(fullfile(summary_dir, ...
    'b1_upper_rapidrise_plateau_trend_overlay.png'));
paths.pdf = string(fullfile(summary_dir, ...
    'b1_upper_rapidrise_plateau_trend_overlay.pdf'));
exportgraphics(fig, char(paths.png), 'Resolution', 220);
exportgraphics(fig, char(paths.pdf), 'ContentType', 'vector');
close(fig);
end


function local_patch_q_regions(ax, x_range, y_range)
yl = y_range;
patch(ax, [-0.02 0.02 0.02 -0.02], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.95 0.97 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
patch(ax, [-0.06 -0.02 -0.02 -0.06], [yl(1) yl(1) yl(2) yl(2)], ...
    [1.00 0.97 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
patch(ax, [0.02 0.06 0.06 0.02], [yl(1) yl(1) yl(2) yl(2)], ...
    [1.00 0.97 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
patch(ax, [x_range(1) -0.06 -0.06 x_range(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.95 1.00 0.95], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.25);
patch(ax, [0.06 x_range(2) x_range(2) 0.06], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.95 1.00 0.95], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.25);
end


function local_mark_rejected_candidates(ax, candidates, selected)
if isempty(candidates) || isempty(selected) || height(candidates) == 0 || ...
        height(selected) == 0
    return
end
selected_q = double(selected.q_Ainv);
selected_e = double(selected.upper_energy_meV);
rejected_q = [];
rejected_e = [];
selected_ids = selected.selected_candidate_id;
for i = 1:height(candidates)
    if any(selected_ids == candidates.candidate_id(i))
        continue
    end
    [distance, idx] = min(abs(selected_q - candidates.q_Ainv(i)));
    if distance < 1e-9 && ...
            abs(candidates.upper_energy_meV(i) - selected_e(idx)) > 250
        rejected_q(end + 1, 1) = candidates.q_Ainv(i); %#ok<AGROW>
        rejected_e(end + 1, 1) = candidates.upper_energy_meV(i); %#ok<AGROW>
    end
end
if ~isempty(rejected_q)
    scatter(ax, rejected_q, rejected_e, 30, 'x', 'MarkerEdgeColor', ...
        [0.85 0.10 0.10], 'LineWidth', 1.1);
end
end


function local_write_readme(path, metrics, best, target_score, gate)
fid = fopen(path, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# B1 Lorentz upper rapidrise-plateau retry\n\n');
fprintf(fid, '- Model: mandatory Lorentz double peak; no single-peak fallback.\n');
fprintf(fid, '- q rule: signed q only; no +q/-q averaging.\n');
fprintf(fid, '- B1 fit window: 300-1800 meV.\n');
fprintf(fid, '- Upper trend: small-q rapid rise is allowed; isolated jumps and high-q wandering are penalized.\n');
fprintf(fid, '- Split candidates: 40, 60, 80, 120, 160, 220, 280, 340, 420 meV.\n');
fprintf(fid, '- Physical fit status: no physical fit was run.\n');
fprintf(fid, '- Target score: %.3g. Best: `%s`, tag `%s`, score %.3g, gate `%d`.\n\n', ...
    target_score, best.id, best.final_tag, best.score, gate);
fprintf(fid, '## Metrics\n\n');
fprintf(fid, '| session | branch | points | jumps>250 | edge | failures | weighted score |\n');
fprintf(fid, '| --- | --- | ---: | ---: | ---: | ---: | ---: |\n');
for i = 1:height(metrics)
    fprintf(fid, '| %s | %s | %d | %d | %d | %d | %.3g |\n', ...
        local_cell_text(metrics.session_key, i), ...
        local_cell_text(metrics.branch, i), metrics.n_points(i), ...
        metrics.large_jump_count(i), metrics.edge_count(i), ...
        metrics.failure_count(i), metrics.weighted_score(i));
end
end


function local_register_outputs(project_root, summary_dir, best, register_best)
results_root = fullfile(project_root, 'paper_results');
current_dir = fullfile(results_root, '00_CURRENT_B1_DOUBLE_PEAK_260509');
if ~isfolder(current_dir)
    mkdir(current_dir);
end
date_dir = fullfile(results_root, '00-by_date', '2026-05-10');
if ~isfolder(date_dir)
    mkdir(date_dir);
end
targets = {summary_dir};
names = {'70_lorentz_tracking_upper_rapidrise_plateau_summary'};
if register_best
    sessions = {'590', 'n0', '20w'};
    for si = 1:numel(sessions)
        targets{end + 1} = local_session_dir(project_root, sessions{si}, ...
            best.final_tag); %#ok<AGROW>
        names{end + 1} = sprintf( ...
            '%02d_%s_lorentz_tracking_upper_rapidrise_plateau_candidate_extraction', ...
            70 + si, sessions{si}); %#ok<AGROW>
    end
end
for i = 1:numel(targets)
    local_create_or_update_junction(fullfile(current_dir, names{i}), targets{i});
    local_create_junction_if_missing(fullfile(date_dir, local_file_name(targets{i})), ...
        targets{i});
    local_append_manifest(project_root, '2026-05-10', 'current_work_date', ...
        local_file_name(targets{i}), fullfile(date_dir, local_file_name(targets{i})), ...
        targets{i});
end
end


function local_update_results_index(project_root, summary_dir, best, gate)
index_path = fullfile(project_root, 'case_studies', 'bisb2026', ...
    'RESULTS_INDEX.md');
if isfile(index_path)
    src = fileread(index_path);
else
    src = "# BiSb results index" + newline;
end
summary_text = strrep(summary_dir, '\', '/');
block = sprintf(['<!-- B1_UPPER_RAPIDRISE_PLATEAU_260510_START -->\n', ...
    '## 2026-05-10 B1 Lorentz upper rapidrise-plateau retry\n\n', ...
    '- Summary: `%s`\n', ...
    '- Best tag: `%s`\n', ...
    '- Best variant: `%s`\n', ...
    '- Best score: `%.12g`\n', ...
    '- Gate <=300: `%d`\n', ...
    '- Physical fit: no physical fit was run in this retry round.\n', ...
    '- Trend rule: small-q rapid rise allowed, isolated jumps and high-q wandering penalized.\n\n', ...
    '<!-- B1_UPPER_RAPIDRISE_PLATEAU_260510_END -->'], ...
    summary_text, best.final_tag, best.id, best.score, gate);
pattern = '<!-- B1_UPPER_RAPIDRISE_PLATEAU_260510_START -->[\s\S]*?<!-- B1_UPPER_RAPIDRISE_PLATEAU_260510_END -->';
if ~isempty(regexp(src, pattern, 'once'))
    src = regexprep(src, pattern, block);
else
    src = sprintf('%s\n\n%s\n', strtrim(src), block);
end
fid = fopen(index_path, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strtrim(src));
end


function local_create_or_update_junction(link_path, target_path)
if ~isfolder(target_path)
    return
end
current_target = local_junction_target(link_path);
if strcmpi(current_target, target_path)
    return
end
if ~isempty(current_target)
    status = system(sprintf('cmd /c rmdir "%s"', link_path));
    if status ~= 0
        return
    end
elseif isfolder(link_path) || isfile(link_path)
    return
end
local_create_junction_if_missing(link_path, target_path);
end


function target = local_junction_target(link_path)
target = '';
if ~(isfolder(link_path) || isfile(link_path))
    return
end
[parent, name] = fileparts(link_path);
[status, out] = system(sprintf('cmd /c dir /AL "%s"', parent));
if status ~= 0
    return
end
expr = [regexptranslate('escape', name), '\s+\[([^\]]+)\]'];
tok = regexp(out, expr, 'tokens', 'once');
if ~isempty(tok)
    target = tok{1};
end
end


function local_create_junction_if_missing(link_path, target_path)
if isfolder(link_path) || isfile(link_path) || ~isfolder(target_path)
    return
end
[parent, ~] = fileparts(link_path);
if ~isfolder(parent)
    mkdir(parent);
end
cmd = sprintf('cmd /c mklink /J "%s" "%s"', link_path, target_path);
system(cmd);
end


function local_append_manifest(project_root, date_text, category, name, ...
    bydate_path, target_path)
manifest = fullfile(project_root, 'paper_results', '00-by_date', ...
    '_manifest.csv');
line = sprintf('"%s","%s","%s","%s","%s","%s"', date_text, category, ...
    name, char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    bydate_path, target_path);
if isfile(manifest)
    src = fileread(manifest);
    if contains(src, target_path)
        return
    end
    fid = fopen(manifest, 'a', 'n', 'UTF-8');
else
    fid = fopen(manifest, 'w', 'n', 'UTF-8');
    fprintf(fid, '"date","category","name","timestamp","by_date_path","target_path"\n');
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', line);
end


function tbl = local_read_table_if_exists(path)
if isfile(path)
    tbl = readtable(path);
else
    tbl = table();
end
end


function metrics = local_empty_metrics_table()
metrics = table('Size', [0 24], ...
    'VariableTypes', {'cell', 'cell', 'logical', 'double', 'double', ...
    'cell', 'cell', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', ...
    'double'}, ...
    'VariableNames', {'variant_id', 'date_tag', ...
    'candidate_for_recommendation', 'b1_energy_min_meV', ...
    'b1_energy_max_meV', 'session_key', 'branch', ...
    'n_points', 'large_jump_count', 'edge_count', ...
    'failure_count', 'repair_count', 'exclusion_count', ...
    'session_weight', 'raw_score', 'weighted_score', ...
    'upper_medium_jump_threshold_meV', ...
    'upper_large_jump_threshold_meV', 'upper_medium_jump_count', ...
    'upper_large_jump_count', 'upper_max_abs_jump_meV', ...
    'upper_overbroad_count', 'upper_low_r2_count', ...
    'upper_stability_score'});
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
        error('run_b1_lorentz_tracking_ridge_guided_failure_retry:UnknownSession', ...
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


function values = local_text_column(tbl, name)
if ~ismember(name, tbl.Properties.VariableNames)
    values = repmat({''}, height(tbl), 1);
    return
end
raw = tbl.(name);
if iscell(raw)
    values = cellfun(@(x) char(string(x)), raw, 'UniformOutput', false);
elseif isstring(raw)
    values = cellstr(raw);
else
    values = cellstr(string(raw));
end
values = values(:);
end


function name = local_file_name(path_value)
[~, name, ext] = fileparts(char(string(path_value)));
name = [name ext];
end
