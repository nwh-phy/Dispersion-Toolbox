function out = run_b1_lorentz_tracking_until300(options)
%RUN_B1_LORENTZ_TRACKING_UNTIL300 Iterate B1 Lorentz tracking to score gate.
%
% This controller keeps the mandatory Lorentz double-peak model and signed-q
% rule. It continues through low-cost repair and parameter variants until the
% weighted extraction score is at or below the hard target score, then runs a
% sandbox fit only for that final candidate.

arguments
    options.targetScore (1,1) double = 300
    options.runExtraction (1,1) logical = true
    options.runOverlay (1,1) logical = true
    options.runSandboxFit (1,1) logical = true
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
addpath(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

base_tag = '260510_lorentz_tracking_until300';
summary_dir = fullfile(project_root, 'paper_results', ...
    'b1_lorentz_tracking_optimization_260510_until300');
if ~isfolder(summary_dir)
    mkdir(summary_dir);
end

variants = local_variant_specs(base_tag);
all_metrics = local_empty_metrics_table();
variant_outputs = local_empty_variant_output();
best_score = Inf;
best_output = local_empty_variant_output();

for vi = 1:numel(variants)
    variant = variants(vi);
    fprintf('\n=== B1 Lorentz until300 %s ===\n', variant.id);
    if options.runExtraction
        local_run_variant(project_root, variant);
    end
    metrics = local_finalize_variant(project_root, variant);
    all_metrics = [all_metrics; metrics]; %#ok<AGROW>
    score = local_variant_score(metrics, variant.id);
    overlay = struct('png', "", 'pdf', "");
    if options.runOverlay
        overlay = run_b1_double_peak_waterfall_extraction_overlay( ...
            dateTag=variant.final_tag);
    end
    output = local_variant_output(variant, score, overlay);
    variant_outputs(end + 1) = output; %#ok<AGROW>
    if score < best_score
        best_score = score;
        best_output = output;
    end
    if best_score <= options.targetScore
        fprintf('  Score %.3g reached target %.3g; stopping variants.\n', ...
            best_score, options.targetScore);
        break
    end
end

metrics_path = fullfile(summary_dir, 'b1_lorentz_tracking_until300_metrics.csv');
writetable(all_metrics, metrics_path);
variant_path = fullfile(summary_dir, 'b1_lorentz_tracking_until300_variants.csv');
writetable(local_variant_outputs_table(variant_outputs), variant_path);

if strlength(best_output.overlay_png) > 0 && isfile(best_output.overlay_png)
    copyfile(best_output.overlay_png, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_overlay_three_dataset_comparison.png'));
end
if strlength(best_output.overlay_pdf) > 0 && isfile(best_output.overlay_pdf)
    copyfile(best_output.overlay_pdf, fullfile(summary_dir, ...
        'b1_double_peak_lorentz_tracking_overlay_three_dataset_comparison.pdf'));
end

fit_outputs = struct();
if best_score <= options.targetScore
    best_label = 'final_candidate';
    local_write_candidate_marker(summary_dir, best_output, options.targetScore);
    if options.runSandboxFit
        fit_outputs = local_run_final_candidate_fit(project_root, ...
            best_output.final_tag);
    end
else
    best_label = 'manual_anchor_needed';
    local_write_manual_anchor_template(project_root, summary_dir, best_output);
end

readme_path = fullfile(summary_dir, 'README_lorentz_tracking_until300.md');
local_write_readme(readme_path, all_metrics, variant_outputs, best_output, ...
    options.targetScore, best_label);
local_register_outputs(project_root, summary_dir, best_output);

out = struct();
out.summary_dir = summary_dir;
out.metrics_csv = metrics_path;
out.variants_csv = variant_path;
out.readme = readme_path;
out.best = best_output;
out.best_score = best_score;
out.target_score = options.targetScore;
out.fit_outputs = fit_outputs;

fprintf('\nB1 Lorentz until300 complete.\n');
fprintf('  Summary: %s\n', summary_dir);
fprintf('  Best: %s (%s), score %.3g, target %.3g, status %s\n', ...
    best_output.id, best_output.final_tag, best_score, options.targetScore, ...
    best_label);
end


function variants = local_variant_specs(base_tag)
variants = repmat(local_empty_variant(), 1, 5);
variants(1) = local_make_variant('v6_windowfallback_jumprepair', ...
    base_tag, 7, 71, 260, 340, 220, 0.08);
variants(2) = local_make_variant('v7_bin7_sg91_w300_shift240', ...
    base_tag, 7, 91, 300, 380, 240, 0.08);
variants(3) = local_make_variant('v7_bin9_sg91_w320_shift260', ...
    base_tag, 9, 91, 320, 420, 260, 0.075);
variants(4) = local_make_variant('v8_trendconstrained_bin9_sg101_w340', ...
    base_tag, 9, 101, 340, 460, 280, 0.075);
variants(5) = local_make_variant('v9_anchor_template_no_fit', ...
    base_tag, 9, 111, 360, 500, 300, 0.07);
end


function variant = local_empty_variant()
variant = struct('id', '', 'trend_tag', '', 'final_tag', '', ...
    'bin_size', NaN, 'sg_high_window', NaN, 'window_half', NaN, ...
    'highq_window_half', NaN, 'max_shift', NaN, 'highq_force_abs', NaN);
end


function variant = local_make_variant(id, base_tag, bin_size, sg_high, ...
    window_half, highq_window_half, max_shift, highq_force_abs)
variant = local_empty_variant();
variant.id = id;
variant.trend_tag = sprintf('%s_%s_trend', base_tag, id);
variant.final_tag = sprintf('%s_%s_windowed', base_tag, id);
variant.bin_size = bin_size;
variant.sg_high_window = sg_high;
variant.window_half = window_half;
variant.highq_window_half = highq_window_half;
variant.max_shift = max_shift;
variant.highq_force_abs = highq_force_abs;
end


function local_run_variant(project_root, variant)
local_call_extraction('all', variant.trend_tag, 'propagated_double_peak', ...
    variant, table(), table());
sessions = {'590', 'n0', '20w'};
for si = 1:numel(sessions)
    ref_dir = local_session_dir(project_root, sessions{si}, variant.trend_tag);
    lower = readtable(fullfile(ref_dir, 'b1_double_peak_lower_points.csv'));
    upper = readtable(fullfile(ref_dir, 'b1_double_peak_upper_points.csv'));
    local_call_extraction(sessions{si}, variant.final_tag, ...
        'windowed_branch_tracking', variant, lower, upper);
end
end


function local_call_extraction(session_request, tag, tracking_mode, variant, ...
    lower_ref, upper_ref)
run_b1_double_peak_binning_analysis(session_request, ...
    runFits=false, ...
    outputDateTag=tag, ...
    qRangeOverride_Ainv=[-0.15 0.15], ...
    b1EnergyWindowOverrideMeV=[300 2100], ...
    peakModelOverride='lorentz', ...
    trackingMode=tracking_mode, ...
    fallbackSplitCandidatesMeV=[120 160 220], ...
    maxTrackingShiftMeV=variant.max_shift, ...
    trackingWindowHalfWidthMeV=variant.window_half, ...
    trackingWindowHighQHalfWidthMeV=variant.highq_window_half, ...
    trackingWindowHighQAbsAinv=0.09, ...
    trackingWindowInvalidFallback='independent_double_peak', ...
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
    fitDenoiseHighWindow=variant.sg_high_window, ...
    fitDenoiseQStartAinv=0.07, ...
    fitDenoiseQEndAinv=0.15, ...
    fitDenoiseOrder=3, ...
    highQForceBinAbsAinv=variant.highq_force_abs, ...
    binSize=variant.bin_size, ...
    enableJumpRepair=true, ...
    largeJumpThresholdMeV=250);
end


function metrics = local_finalize_variant(project_root, variant)
metrics = local_empty_metrics_table();
sessions = {'590', 'n0', '20w'};
for si = 1:numel(sessions)
    session_dir = local_session_dir(project_root, sessions{si}, ...
        variant.final_tag);
    local_copy_variant_named_csv(session_dir, variant.id);
    rows = local_metrics_one_session(session_dir, variant);
    writetable(rows, fullfile(session_dir, sprintf( ...
        'b1_double_peak_lorentz_tracking_until300_%s_tracking_metrics.csv', ...
        variant.id)));
    metrics = [metrics; rows]; %#ok<AGROW>
end
end


function local_copy_variant_named_csv(session_dir, variant_id)
copies = { ...
    'b1_double_peak_combined_q_points.csv', 'combined_points.csv'; ...
    'b1_double_peak_lower_points.csv', 'lower_points.csv'; ...
    'b1_double_peak_upper_points.csv', 'upper_points.csv'; ...
    'b1_double_peak_fit_failures.csv', 'failures.csv'; ...
    'b1_double_peak_repair_log.csv', 'repair_log.csv'; ...
    'b1_double_peak_exclusion_points.csv', 'exclusion_points.csv'};
for i = 1:size(copies, 1)
    src = fullfile(session_dir, copies{i, 1});
    if isfile(src)
        dst = fullfile(session_dir, sprintf( ...
            'b1_double_peak_lorentz_tracking_until300_%s_%s', ...
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
end


function row = local_branch_metric_row(variant, session_key, branch, points, ...
    failure_count, repair_count, exclusion_count)
n_points = height(points);
jump_count = 0;
edge_count = 0;
if n_points >= 1 && any(strcmp(points.Properties.VariableNames, 'energy_meV'))
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
row = table({variant.id}, {variant.final_tag}, {session_key}, {branch}, ...
    n_points, jump_count, edge_count, failure_count, repair_count, ...
    exclusion_count, weight, raw_score, weighted_score, ...
    'VariableNames', {'variant_id', 'date_tag', 'session_key', ...
    'branch', 'n_points', 'large_jump_count', 'edge_count', ...
    'failure_count', 'repair_count', 'exclusion_count', ...
    'session_weight', 'raw_score', 'weighted_score'});
end


function score = local_variant_score(metrics, variant_id)
rows = strcmp(metrics.variant_id, variant_id);
if ~any(rows)
    score = Inf;
else
    score = sum(metrics.weighted_score(rows), 'omitnan');
end
end


function fit_outputs = local_run_final_candidate_fit(project_root, final_tag)
datasets = local_fit_datasets(project_root, final_tag);
fit_outputs = struct();
fit_outputs.lower_physical = run_b1_physical_fit_analysis( ...
    datasets=datasets, ...
    branchFileName='b1_double_peak_lower_points.csv', ...
    outputTag='b1_double_peak_lorentz_tracking_until300_final_candidate_fit_lower', ...
    filePrefix='b1_double_peak_lorentz_until300_lower');
fit_outputs.upper_physical = run_b1_physical_fit_analysis( ...
    datasets=datasets, ...
    branchFileName='b1_double_peak_upper_points.csv', ...
    outputTag='b1_double_peak_lorentz_tracking_until300_final_candidate_fit_upper', ...
    filePrefix='b1_double_peak_lorentz_until300_upper');
fit_outputs.lower_enhancements = run_b1_physical_fit_enhancements( ...
    datasets=datasets, ...
    branchFileName='b1_double_peak_lower_points.csv', ...
    outputTag='b1_double_peak_lorentz_tracking_until300_final_candidate_fit_lower_enhancements', ...
    filePrefix='b1_double_peak_lorentz_until300_lower');
fit_outputs.upper_enhancements = run_b1_physical_fit_enhancements( ...
    datasets=datasets, ...
    branchFileName='b1_double_peak_upper_points.csv', ...
    outputTag='b1_double_peak_lorentz_tracking_until300_final_candidate_fit_upper_enhancements', ...
    filePrefix='b1_double_peak_lorentz_until300_upper');
end


function datasets = local_fit_datasets(project_root, tag)
requests = {'590', 'n0', '20w'};
keys = {'590_PL2_10w', 'n0_PL2_10w_repeat', 'no_PL2_20w_2film'};
labels = {'590 10w defocus 1film', ...
    'n0 10w defocus repeat 1film', '20w defocus 2film'};
classes = {'1film', '1film', '2film'};
factors = [1, 1, 2];
colors = {[0.120, 0.470, 0.900], [0.160, 0.500, 0.220], ...
    [0.930, 0.280, 0.300]};
datasets = struct('session_key', {}, 'session_label', {}, ...
    'input_dir', {}, 'thickness_class', {}, 'thickness_factor', {}, ...
    'color', {}, 'marker', {});
for i = 1:numel(requests)
    datasets(end + 1) = struct( ... %#ok<AGROW>
        'session_key', keys{i}, ...
        'session_label', labels{i}, ...
        'input_dir', local_session_dir(project_root, requests{i}, tag), ...
        'thickness_class', classes{i}, ...
        'thickness_factor', factors(i), ...
        'color', colors{i}, ...
        'marker', 'o');
end
end


function local_write_candidate_marker(summary_dir, best_output, target_score)
path = fullfile(summary_dir, 'FINAL_CANDIDATE_REACHED_300_GATE.txt');
fid = fopen(path, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Final candidate reached score gate.\n');
fprintf(fid, 'variant=%s\n', best_output.id);
fprintf(fid, 'date_tag=%s\n', best_output.final_tag);
fprintf(fid, 'score=%.12g\n', best_output.score);
fprintf(fid, 'target_score=%.12g\n', target_score);
end


function local_write_manual_anchor_template(project_root, summary_dir, best_output)
sessions = {'590', 'n0', '20w'};
anchor = table('Size', [0 6], ...
    'VariableTypes', {'cell', 'double', 'cell', 'double', 'double', 'cell'}, ...
    'VariableNames', {'session_request', 'q_Ainv', 'branch_label', ...
    'manual_lower_meV', 'manual_upper_meV', 'reason'});
for si = 1:numel(sessions)
    session_dir = local_session_dir(project_root, sessions{si}, ...
        best_output.final_tag);
    failures = local_read_table_if_exists(fullfile(session_dir, ...
        'b1_double_peak_fit_failures.csv'));
    for i = 1:height(failures)
        anchor = [anchor; table(sessions(si), failures.q_Ainv(i), ...
            {'both'}, NaN, NaN, failures.status(i), ...
            'VariableNames', anchor.Properties.VariableNames)]; %#ok<AGROW>
    end
end
writetable(anchor, fullfile(summary_dir, ...
    'b1_double_peak_manual_anchor_template_until300.csv'));
end


function local_write_readme(path, metrics, variant_outputs, best_output, ...
    target_score, best_label)
fid = fopen(path, 'w', 'n', 'UTF-8');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# B1 Lorentz tracking until300\n\n');
fprintf(fid, '- Model: mandatory Lorentz double peak.\n');
fprintf(fid, '- q rule: signed-q only, no +q/-q averaging.\n');
fprintf(fid, '- Score formula: 10*jumps + 3*edge + 5*failures; 20w weight=2.\n');
fprintf(fid, '- Hard fit gate: total score <= %.3g.\n', target_score);
fprintf(fid, '- Best status: `%s`.\n', best_label);
fprintf(fid, '- Best variant: `%s`, tag `%s`, score %.3g.\n\n', ...
    best_output.id, best_output.final_tag, best_output.score);

fprintf(fid, '## Variants\n\n');
fprintf(fid, '| variant | tag | score | overlay |\n');
fprintf(fid, '| --- | --- | ---: | --- |\n');
for i = 1:numel(variant_outputs)
    fprintf(fid, '| %s | %s | %.3g | %s |\n', variant_outputs(i).id, ...
        variant_outputs(i).final_tag, variant_outputs(i).score, ...
        local_file_name(variant_outputs(i).overlay_png));
end

fprintf(fid, '\n## Metrics\n\n');
fprintf(fid, '| variant | session | branch | points | jumps | edge | failures | repairs | exclusions | weighted score |\n');
fprintf(fid, '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n');
for i = 1:height(metrics)
    fprintf(fid, '| %s | %s | %s | %d | %d | %d | %d | %d | %d | %.3g |\n', ...
        local_cell_text(metrics.variant_id, i), ...
        local_cell_text(metrics.session_key, i), ...
        local_cell_text(metrics.branch, i), metrics.n_points(i), ...
        metrics.large_jump_count(i), metrics.edge_count(i), ...
        metrics.failure_count(i), metrics.repair_count(i), ...
        metrics.exclusion_count(i), metrics.weighted_score(i));
end
end


function local_register_outputs(project_root, summary_dir, best_output)
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
names = {'31_lorentz_tracking_until300_summary'};
sessions = {'590', 'n0', '20w'};
for si = 1:numel(sessions)
    targets{end + 1} = local_session_dir(project_root, sessions{si}, ...
        best_output.final_tag); %#ok<AGROW>
    names{end + 1} = sprintf('%02d_%s_lorentz_tracking_until300_best_extraction', ...
        31 + si, sessions{si}); %#ok<AGROW>
end
fit_dirs = { ...
    'b1_double_peak_lorentz_tracking_until300_final_candidate_fit_lower', ...
    'b1_double_peak_lorentz_tracking_until300_final_candidate_fit_upper', ...
    'b1_double_peak_lorentz_tracking_until300_final_candidate_fit_lower_enhancements', ...
    'b1_double_peak_lorentz_tracking_until300_final_candidate_fit_upper_enhancements'};
for i = 1:numel(fit_dirs)
    target = fullfile(results_root, fit_dirs{i});
    if isfolder(target)
        targets{end + 1} = target; %#ok<AGROW>
        names{end + 1} = sprintf('%02d_lorentz_tracking_until300_%s', ...
            34 + i, strrep(fit_dirs{i}, ...
            'b1_double_peak_lorentz_tracking_until300_', '')); %#ok<AGROW>
    end
end
for i = 1:numel(targets)
    local_create_junction_if_missing(fullfile(current_dir, names{i}), targets{i});
    local_create_junction_if_missing(fullfile(date_dir, local_file_name(targets{i})), ...
        targets{i});
    local_append_manifest(project_root, '2026-05-10', 'current_work_date', ...
        local_file_name(targets{i}), fullfile(date_dir, local_file_name(targets{i})), ...
        targets{i});
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


function outputs = local_empty_variant_output()
outputs = struct('id', {}, 'final_tag', {}, 'score', {}, ...
    'overlay_png', {}, 'overlay_pdf', {});
end


function output = local_variant_output(variant, score, overlay)
if strlength(string(overlay.png)) == 0
    project_root = bisb_find_project_root(fileparts(mfilename('fullpath')));
    overlay_dir = fullfile(project_root, 'paper_results', ...
        sprintf(['b1_double_peak_waterfall_extraction_overlay_three_', ...
        'dataset_comparison_%s'], variant.final_tag));
    png = fullfile(overlay_dir, ...
        'b1_double_peak_waterfall_extraction_overlay_three_dataset_comparison.png');
    pdf = fullfile(overlay_dir, ...
        'b1_double_peak_waterfall_extraction_overlay_three_dataset_comparison.pdf');
    if isfile(png)
        overlay.png = string(png);
    end
    if isfile(pdf)
        overlay.pdf = string(pdf);
    end
end
output = struct('id', variant.id, 'final_tag', variant.final_tag, ...
    'score', score, 'overlay_png', string(overlay.png), ...
    'overlay_pdf', string(overlay.pdf));
end


function tbl = local_variant_outputs_table(outputs)
if isempty(outputs)
    tbl = table();
    return
end
overlay_png = string({outputs.overlay_png}).';
overlay_pdf = string({outputs.overlay_pdf}).';
tbl = table({outputs.id}.', {outputs.final_tag}.', [outputs.score].', ...
    overlay_png, overlay_pdf, ...
    'VariableNames', {'variant_id', 'date_tag', 'score', ...
    'overlay_png', 'overlay_pdf'});
end


function tbl = local_read_table_if_exists(path)
if isfile(path)
    tbl = readtable(path);
else
    tbl = table();
end
end


function metrics = local_empty_metrics_table()
metrics = table('Size', [0 13], ...
    'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double'}, ...
    'VariableNames', {'variant_id', 'date_tag', 'session_key', ...
    'branch', 'n_points', 'large_jump_count', 'edge_count', ...
    'failure_count', 'repair_count', 'exclusion_count', ...
    'session_weight', 'raw_score', 'weighted_score'});
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
        error('run_b1_lorentz_tracking_until300:UnknownSession', ...
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


function name = local_file_name(path_value)
[~, name, ext] = fileparts(char(string(path_value)));
name = [name ext];
end
