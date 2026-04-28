function results = compare_bg_dual_window_all()
%COMPARE_BG_DUAL_WINDOW_ALL Compare background models across three Bi datasets.
% Evaluates single-window and Nion-inspired dual-window settings within
% |q| <= 0.15 A^-1 across all available eq3D datasets under 20260120 Bi.
%
% Usage:
%   results = compare_bg_dual_window_all();
%
% Prints JSON to stdout for checkpointing.

project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));

datasets = {
    struct('name', '590_PL2_10w', 'path', fullfile(project_root, '20260120 Bi', '590 PL2 10w 0.004 10sx300', 'eq3D.mat')), ...
    struct('name', 'n0_pl2_10w', 'path', fullfile(project_root, '20260120 Bi', 'n0 pl2 10w 0.004 10s x300', 'eq3D.mat')), ...
    struct('name', 'no_pl2_20w_2film', 'path', fullfile(project_root, '20260120 Bi', 'no pl2 20w 0.004 10sx300 2film', 'eq3D.mat')) ...
};

configs = {
    struct('name','single_auto_core4', 'method','Auto', 'candidates',{{'Power','Exp2','ExpPoly3','Pearson'}}, 'win_lo',[50 300], 'win_hi',[], 'group_qmax', NaN, 'group_ranges', zeros(0,2)), ...
    struct('name','single_auto_plus_vii', 'method','Auto', 'candidates',{{'Power','Exp2','ExpPoly3','Pearson','PearsonVII'}}, 'win_lo',[50 300], 'win_hi',[], 'group_qmax', NaN, 'group_ranges', zeros(0,2)), ...
    struct('name','dual_auto_plus_vii', 'method','Auto', 'candidates',{{'Power','Exp2','ExpPoly3','Pearson','PearsonVII'}}, 'win_lo',[50 300], 'win_hi',[2800 3400], 'group_qmax', NaN, 'group_ranges', zeros(0,2)), ...
    struct('name','dual_auto_plus_vii_group_0p01', 'method','Auto', 'candidates',{{'Power','Exp2','ExpPoly3','Pearson','PearsonVII'}}, 'win_lo',[50 300], 'win_hi',[2800 3400], 'group_qmax', 0.01, 'group_ranges', zeros(0,2)), ...
    struct('name','dual_auto_plus_vii_group_0p01_pos_shoulder', 'method','Auto', 'candidates',{{'Power','Exp2','ExpPoly3','Pearson','PearsonVII'}}, 'win_lo',[50 300], 'win_hi',[2800 3400], 'group_qmax', 0.01, 'group_ranges', [0.02 0.05]), ...
    struct('name','dual_exp2', 'method','Exp2', 'candidates',{{'Exp2'}}, 'win_lo',[50 300], 'win_hi',[2800 3400], 'group_qmax', NaN, 'group_ranges', zeros(0,2)), ...
    struct('name','dual_exppoly3', 'method','ExpPoly3', 'candidates',{{'ExpPoly3'}}, 'win_lo',[50 300], 'win_hi',[2800 3400], 'group_qmax', NaN, 'group_ranges', zeros(0,2)), ...
    struct('name','dual_pearson', 'method','Pearson', 'candidates',{{'Pearson'}}, 'win_lo',[50 300], 'win_hi',[2800 3400], 'group_qmax', NaN, 'group_ranges', zeros(0,2)), ...
    struct('name','dual_pearsonvii', 'method','PearsonVII', 'candidates',{{'PearsonVII'}}, 'win_lo',[50 300], 'win_hi',[2800 3400], 'group_qmax', NaN, 'group_ranges', zeros(0,2)) ...
};

results = cell(1, numel(datasets));
for di = 1:numel(datasets)
    results{di} = evaluate_dataset(datasets{di}, configs);
end

fprintf('%s\n', jsonencode(results));
end

function dataset_result = evaluate_dataset(dataset_cfg, configs)
dataset = load_qe_dataset(dataset_cfg.path);
qe = dataset.qe;
q_axis = qe.q_Ainv(:)';
focus_idx = find(abs(q_axis) <= 0.15);

qe_focus = qe;
qe_focus.intensity = qe.intensity(:, focus_idx);
qe_focus.q_channel = qe.q_channel(focus_idx);
qe_focus.q_Ainv = qe.q_Ainv(focus_idx);
qe_focus.q_abs_Ainv = qe.q_abs_Ainv(focus_idx);
[~, local_zero] = min(abs(qe_focus.q_Ainv));
qe_focus.q_zero_index = local_zero;

energy = qe_focus.energy_meV(:);
sub_mask = energy > 50;
branch1_mask = energy >= 600 & energy <= 1700;
branch2_mask = energy >= 1700 & energy <= 2600;
branch3_mask = energy >= 2900 & energy <= 3500;
[~, q0_idx] = min(abs(qe_focus.q_Ainv));

comparisons = cell(1, numel(configs));
for ci = 1:numel(configs)
    cfg = configs{ci};
    opts = struct();
    opts.do_despike = false;
    opts.do_normalize = false;
    opts.do_denoise = false;
    opts.do_bg_sub = true;
    opts.bg_method = cfg.method;
    opts.bg_candidate_methods = cfg.candidates;
    opts.bg_win_lo = cfg.win_lo;
    opts.bg_win_hi = cfg.win_hi;
    opts.bg_iterative = false;
    opts.do_deconv = false;
    if isfield(cfg, 'group_qmax') && isfinite(cfg.group_qmax) && cfg.group_qmax > 0
        opts.bg_auto_group_qmax = cfg.group_qmax;
    end
    if isfield(cfg, 'group_ranges') && ~isempty(cfg.group_ranges)
        opts.bg_auto_group_ranges = cfg.group_ranges;
    end

    [qe_bg, bg_diag] = qe_preprocess(qe_focus, opts);
    summary = qe_summarize_lowq_background_eval(qe_focus, qe_bg, bg_diag);

    neg_fraction = [bg_diag.neg_fraction];
    neg_area_fraction = [bg_diag.neg_area_fraction];
    neg_peak_fraction = [bg_diag.neg_peak_fraction];
    bg_fraction = [bg_diag.bg_fraction];
    linear_rmse = [bg_diag.linear_rmse];
    selected = string({bg_diag.selected_method});
    unique_methods = unique(selected);
    counts = zeros(size(unique_methods));
    for i = 1:numel(unique_methods)
        counts(i) = sum(selected == unique_methods(i));
    end

    proc = qe_bg.intensity(sub_mask, :);
    q0_proc = proc(:, q0_idx);

    item = struct();
    item.name = cfg.name;
    item.method = cfg.method;
    item.win_lo = cfg.win_lo;
    item.win_hi = cfg.win_hi;
    item.group_qmax = cfg.group_qmax;
    item.group_ranges = cfg.group_ranges;
    item.evaluation_ranges = summary.evaluation_ranges;
    item.branch_windows = summary.branch_windows;
    item.regions = summary.regions;
    item.q0_summary = summary.q0;
    item.selected_methods = cellstr(unique_methods);
    item.selected_method_counts = counts;
    item.neg_fraction_mean = mean(neg_fraction, 'omitnan');
    item.neg_fraction_median = median(neg_fraction, 'omitnan');
    item.neg_area_fraction_mean = mean(neg_area_fraction, 'omitnan');
    item.neg_area_fraction_median = median(neg_area_fraction, 'omitnan');
    item.neg_peak_fraction_mean = mean(neg_peak_fraction, 'omitnan');
    item.neg_peak_fraction_median = median(neg_peak_fraction, 'omitnan');
    item.bg_fraction_mean = mean(bg_fraction, 'omitnan');
    item.linear_rmse_mean = mean(linear_rmse, 'omitnan');
    item.channels_neg_fraction_gt_0p10 = sum(neg_fraction > 0.10);
    item.channels_neg_area_gt_0p05 = sum(neg_area_fraction > 0.05);
    item.branch1_peak_mean = mean(max(proc(branch1_mask(sub_mask), :), [], 1), 'omitnan');
    item.branch2_peak_mean = mean(max(proc(branch2_mask(sub_mask), :), [], 1), 'omitnan');
    item.branch3_peak_mean = mean(max(proc(branch3_mask(sub_mask), :), [], 1), 'omitnan');
    item.q0_selected_method = bg_diag(q0_idx).selected_method;
    item.q0_neg_fraction = bg_diag(q0_idx).neg_fraction;
    item.q0_neg_area_fraction = bg_diag(q0_idx).neg_area_fraction;
    item.q0_neg_peak_fraction = bg_diag(q0_idx).neg_peak_fraction;
    item.q0_branch1_peak = max(q0_proc(branch1_mask(sub_mask)));
    item.q0_branch2_peak = max(q0_proc(branch2_mask(sub_mask)));
    item.q0_branch3_peak = max(q0_proc(branch3_mask(sub_mask)));
    item.q0_processed_min_after50 = min(q0_proc);
    item.q0_negative_points_after50 = sum(q0_proc < 0);
    if strcmp(cfg.method, 'Auto')
        item.q0_candidate_methods = bg_diag(q0_idx).candidate_methods;
        item.q0_candidate_scores = bg_diag(q0_idx).candidate_scores;
        item.q0_candidate_linear_rmse = [bg_diag(q0_idx).candidate_details.linear_rmse];
    end
    comparisons{ci} = item;
end

dataset_result = struct();
dataset_result.dataset_name = dataset_cfg.name;
dataset_result.data_path = dataset_cfg.path;
dataset_result.focus_channel_count = numel(focus_idx);
dataset_result.actual_q_min = min(qe_focus.q_Ainv);
dataset_result.actual_q_max = max(qe_focus.q_Ainv);
dataset_result.comparisons = comparisons;
dataset_result.screening = qe_recommend_lowq_background_configs(comparisons);
end
