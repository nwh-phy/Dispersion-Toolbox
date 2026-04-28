function summary = qe_summarize_lowq_background_eval(qe_focus, qe_bg, bg_diag, opts)
%QE_SUMMARIZE_LOWQ_BACKGROUND_EVAL Structured low-q background-evaluation summary.
%   SUMMARY = QE_SUMMARIZE_LOWQ_BACKGROUND_EVAL(QE_FOCUS, QE_BG, BG_DIAG)
%   computes region-wise branch-preservation and oversubtraction summaries for
%   a low-q focused qe dataset after background subtraction.
%
%   Optional opts fields:
%     sub_energy_min_meV          default 50
%     center_qmax                 default 0.01
%     positive_shoulder_qrange    default [0.02 0.05]
%     negative_shoulder_qrange    default [-0.05 -0.02]
%     branch1_meV                 default [600 1700]
%     branch2_meV                 default [1700 2600]
%     branch3_meV                 default [2900 3500]

if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = local_apply_defaults(opts);

energy = qe_focus.energy_meV(:);
q_axis = qe_focus.q_Ainv(:)';
sub_mask = energy > opts.sub_energy_min_meV;
proc = qe_bg.intensity(sub_mask, :);
sub_energy = energy(sub_mask);

summary = struct();
summary.focus_channel_count = numel(q_axis);
summary.actual_q_min = local_safe_min(q_axis);
summary.actual_q_max = local_safe_max(q_axis);
summary.evaluation_ranges = struct( ...
    'sub_energy_min_meV', opts.sub_energy_min_meV, ...
    'center_qmax', opts.center_qmax, ...
    'positive_shoulder_qrange', opts.positive_shoulder_qrange, ...
    'negative_shoulder_qrange', opts.negative_shoulder_qrange);
summary.branch_windows = struct( ...
    'branch1_meV', opts.branch1_meV, ...
    'branch2_meV', opts.branch2_meV, ...
    'branch3_meV', opts.branch3_meV);

full_idx = 1:numel(q_axis);
center_idx = find(abs(q_axis) <= opts.center_qmax);
pos_idx = find(q_axis >= opts.positive_shoulder_qrange(1) & q_axis <= opts.positive_shoulder_qrange(2));
neg_idx = find(q_axis >= opts.negative_shoulder_qrange(1) & q_axis <= opts.negative_shoulder_qrange(2));

summary.regions = struct();
summary.regions.full_low_q = local_region_summary(full_idx, q_axis, proc, sub_energy, bg_diag, opts);
summary.regions.center_core = local_region_summary(center_idx, q_axis, proc, sub_energy, bg_diag, opts);
summary.regions.positive_shoulder = local_region_summary(pos_idx, q_axis, proc, sub_energy, bg_diag, opts);
summary.regions.negative_shoulder = local_region_summary(neg_idx, q_axis, proc, sub_energy, bg_diag, opts);

if isfield(qe_focus, 'q_zero_index') && isscalar(qe_focus.q_zero_index) && ...
        isfinite(qe_focus.q_zero_index)
    candidate_q0_idx = round(qe_focus.q_zero_index);
else
    candidate_q0_idx = NaN;
end
if isfinite(candidate_q0_idx) && candidate_q0_idx >= 1 && candidate_q0_idx <= numel(q_axis)
    q0_idx = candidate_q0_idx;
elseif ~isempty(q_axis)
    [~, q0_idx] = min(abs(q_axis));
else
    q0_idx = NaN;
end

if ~isfinite(q0_idx) || isempty(bg_diag) || q0_idx > numel(bg_diag)
    summary.q0 = local_empty_q0_summary();
    return;
end

q0_diag = bg_diag(q0_idx);
q0_proc = proc(:, q0_idx);
summary.q0 = struct();
summary.q0.q_index = q0_idx;
summary.q0.q_Ainv = q_axis(q0_idx);
summary.q0.selected_method = q0_diag.selected_method;
summary.q0.neg_fraction = q0_diag.neg_fraction;
summary.q0.neg_area_fraction = q0_diag.neg_area_fraction;
summary.q0.neg_peak_fraction = q0_diag.neg_peak_fraction;
summary.q0.bg_fraction = q0_diag.bg_fraction;
summary.q0.linear_rmse = q0_diag.linear_rmse;
summary.q0.branch1_peak = local_branch_peak(q0_proc, sub_energy, opts.branch1_meV);
summary.q0.branch2_peak = local_branch_peak(q0_proc, sub_energy, opts.branch2_meV);
summary.q0.branch3_peak = local_branch_peak(q0_proc, sub_energy, opts.branch3_meV);
summary.q0.processed_min_after50 = local_safe_min(q0_proc);
summary.q0.negative_points_after50 = sum(q0_proc < 0);
if isfield(q0_diag, 'candidate_methods') && ~isempty(q0_diag.candidate_methods)
    summary.q0.candidate_methods = q0_diag.candidate_methods;
    if isfield(q0_diag, 'candidate_scores')
        summary.q0.candidate_scores = q0_diag.candidate_scores;
    end
    if isfield(q0_diag, 'candidate_details')
        summary.q0.candidate_linear_rmse = local_struct_field_array(q0_diag.candidate_details, 'linear_rmse');
        summary.q0.candidate_neg_fraction = local_struct_field_array(q0_diag.candidate_details, 'neg_fraction');
        summary.q0.candidate_neg_area_fraction = local_struct_field_array(q0_diag.candidate_details, 'neg_area_fraction');
        summary.q0.candidate_neg_peak_fraction = local_struct_field_array(q0_diag.candidate_details, 'neg_peak_fraction');
        summary.q0.candidate_bg_fraction = local_struct_field_array(q0_diag.candidate_details, 'bg_fraction');
    end
end
end


function opts = local_apply_defaults(opts)
opts = local_set_default(opts, 'sub_energy_min_meV', 50);
opts = local_set_default(opts, 'center_qmax', 0.01);
opts = local_set_default(opts, 'positive_shoulder_qrange', [0.02 0.05]);
opts = local_set_default(opts, 'negative_shoulder_qrange', [-0.05 -0.02]);
opts = local_set_default(opts, 'branch1_meV', [600 1700]);
opts = local_set_default(opts, 'branch2_meV', [1700 2600]);
opts = local_set_default(opts, 'branch3_meV', [2900 3500]);
end


function opts = local_set_default(opts, field_name, default_value)
if ~isfield(opts, field_name) || isempty(opts.(field_name))
    opts.(field_name) = default_value;
end
end


function region = local_region_summary(idx, q_axis, proc, sub_energy, bg_diag, opts)
region = struct();
region.channel_count = numel(idx);
if isempty(idx)
    region.q_min = NaN;
    region.q_max = NaN;
    region.selected_methods = {};
    region.selected_method_counts = [];
    region.neg_fraction_mean = NaN;
    region.neg_area_fraction_mean = NaN;
    region.neg_peak_fraction_mean = NaN;
    region.bg_fraction_mean = NaN;
    region.linear_rmse_mean = NaN;
    region.branch1_peak_mean = NaN;
    region.branch1_peak_min = NaN;
    region.branch2_peak_mean = NaN;
    region.branch2_peak_min = NaN;
    region.branch3_peak_mean = NaN;
    region.branch3_peak_min = NaN;
    region.processed_min_after50 = NaN;
    region.negative_point_count_after50 = 0;
    return;
end

region.q_min = min(q_axis(idx));
region.q_max = max(q_axis(idx));
diag_idx = idx(idx <= numel(bg_diag));
if isempty(diag_idx)
    region.selected_methods = {};
    region.selected_method_counts = [];
    region.neg_fraction_mean = NaN;
    region.neg_area_fraction_mean = NaN;
    region.neg_peak_fraction_mean = NaN;
    region.bg_fraction_mean = NaN;
    region.linear_rmse_mean = NaN;
else
    selected = string({bg_diag(diag_idx).selected_method});
    unique_methods = unique(selected);
    counts = zeros(1, numel(unique_methods));
    for i = 1:numel(unique_methods)
        counts(i) = sum(selected == unique_methods(i));
    end
    region.selected_methods = cellstr(unique_methods);
    region.selected_method_counts = counts;

    neg_fraction = [bg_diag(diag_idx).neg_fraction];
    neg_area_fraction = [bg_diag(diag_idx).neg_area_fraction];
    neg_peak_fraction = [bg_diag(diag_idx).neg_peak_fraction];
    bg_fraction = [bg_diag(diag_idx).bg_fraction];
    linear_rmse = [bg_diag(diag_idx).linear_rmse];
    region.neg_fraction_mean = mean(neg_fraction, 'omitnan');
    region.neg_area_fraction_mean = mean(neg_area_fraction, 'omitnan');
    region.neg_peak_fraction_mean = mean(neg_peak_fraction, 'omitnan');
    region.bg_fraction_mean = mean(bg_fraction, 'omitnan');
    region.linear_rmse_mean = mean(linear_rmse, 'omitnan');
end

region_proc = proc(:, idx);
[region.branch1_peak_mean, region.branch1_peak_min] = local_region_branch_stats(region_proc, sub_energy, opts.branch1_meV);
[region.branch2_peak_mean, region.branch2_peak_min] = local_region_branch_stats(region_proc, sub_energy, opts.branch2_meV);
[region.branch3_peak_mean, region.branch3_peak_min] = local_region_branch_stats(region_proc, sub_energy, opts.branch3_meV);
region.processed_min_after50 = local_safe_min(region_proc(:));
region.negative_point_count_after50 = sum(region_proc(:) < 0);
end


function [peak_mean, peak_min] = local_region_branch_stats(region_proc, sub_energy, window_meV)
mask = sub_energy >= window_meV(1) & sub_energy <= window_meV(2);
if isempty(region_proc) || ~any(mask)
    peak_mean = NaN;
    peak_min = NaN;
    return;
end
branch_proc = region_proc(mask, :);
per_channel_peak = max(branch_proc, [], 1);
peak_mean = mean(per_channel_peak, 'omitnan');
peak_min = min(per_channel_peak, [], 'omitnan');
end


function peak_value = local_branch_peak(trace, sub_energy, window_meV)
mask = sub_energy >= window_meV(1) & sub_energy <= window_meV(2);
if isempty(trace) || ~any(mask)
    peak_value = NaN;
    return;
end
peak_value = max(trace(mask));
end


function out = local_struct_field_array(details, field_name)
out = [];
if isempty(details)
    return;
end
if ~isstruct(details) || ~isfield(details, field_name)
    return;
end
out = [details.(field_name)];
end


function value = local_safe_min(x)
if isempty(x)
    value = NaN;
else
    value = min(x);
end
end


function value = local_safe_max(x)
if isempty(x)
    value = NaN;
else
    value = max(x);
end
end


function q0 = local_empty_q0_summary()
q0 = struct( ...
    'q_index', NaN, ...
    'q_Ainv', NaN, ...
    'selected_method', '', ...
    'neg_fraction', NaN, ...
    'neg_area_fraction', NaN, ...
    'neg_peak_fraction', NaN, ...
    'bg_fraction', NaN, ...
    'linear_rmse', NaN, ...
    'branch1_peak', NaN, ...
    'branch2_peak', NaN, ...
    'branch3_peak', NaN, ...
    'processed_min_after50', NaN, ...
    'negative_points_after50', 0);
end
