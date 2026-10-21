function screening = qe_recommend_lowq_background_configs(comparisons, opts)
%QE_RECOMMEND_LOWQ_BACKGROUND_CONFIGS Recommend low-q background contracts.
%   SCREENING = QE_RECOMMEND_LOWQ_BACKGROUND_CONFIGS(COMPARISONS) ranks config
%   comparisons by center/shoulder conservativeness and suggests whether the
%   grouped center / grouped shoulder contracts should be enabled.

if nargin < 2 || isempty(opts)
    opts = struct();
end
opts = local_apply_defaults(opts);
items = local_normalize_comparisons(comparisons);

names = cellfun(@(item) item.name, items, 'UniformOutput', false);
center_neg_area = cellfun(@(item) local_metric(item, {'regions','center_core','neg_area_fraction_mean'}), items);
center_neg_peak = cellfun(@(item) local_metric(item, {'regions','center_core','neg_peak_fraction_mean'}), items);
center_branch1 = cellfun(@(item) local_metric(item, {'regions','center_core','branch1_peak_mean'}), items);
pos_neg_area = cellfun(@(item) local_metric(item, {'regions','positive_shoulder','neg_area_fraction_mean'}), items);
pos_neg_peak = cellfun(@(item) local_metric(item, {'regions','positive_shoulder','neg_peak_fraction_mean'}), items);
neg_neg_area = cellfun(@(item) local_metric(item, {'regions','negative_shoulder','neg_area_fraction_mean'}), items);
neg_neg_peak = cellfun(@(item) local_metric(item, {'regions','negative_shoulder','neg_peak_fraction_mean'}), items);
q0_neg_area = cellfun(@(item) local_metric(item, {'q0_summary','neg_area_fraction'}), items);

screening = struct();
screening.available_config_names = names;
screening.rankings = struct();
screening.rankings.center_conservative_names = names(local_rank([center_neg_area(:), center_neg_peak(:), q0_neg_area(:), -center_branch1(:)]));
screening.rankings.positive_shoulder_conservative_names = names(local_rank([pos_neg_area(:), pos_neg_peak(:), center_neg_area(:), -center_branch1(:)]));
screening.rankings.negative_shoulder_conservative_names = names(local_rank([neg_neg_area(:), neg_neg_peak(:), center_neg_area(:), -center_branch1(:)]));
screening.rankings.center_branch1_names = names(local_rank([-center_branch1(:), center_neg_area(:), center_neg_peak(:)]));

finite_center = isfinite(center_neg_area);
if any(finite_center)
    best_center = min(center_neg_area(finite_center));
    safe_mask = finite_center & center_neg_area <= best_center * opts.center_safe_multiplier;
else
    safe_mask = false(size(center_neg_area));
end
safe_names = names(safe_mask);
if any(safe_mask)
    safe_order = local_rank([-center_branch1(safe_mask).', center_neg_area(safe_mask).', center_neg_peak(safe_mask).']);
    safe_names = safe_names(safe_order);
end
screening.safe_center_names = safe_names;
if ~isempty(safe_names)
    screening.recommended_center_name = safe_names{1};
elseif ~isempty(screening.rankings.center_conservative_names)
    screening.recommended_center_name = screening.rankings.center_conservative_names{1};
else
    screening.recommended_center_name = '';
end

screening.decisions = struct();
screening.decisions.center_grouping = local_group_decision(items, names, ...
    opts.center_group_baseline_name, opts.center_group_candidate_name, ...
    opts.min_center_neg_area_improvement_frac, opts.max_center_branch1_drop_frac, ...
    'center_core', 'center_core', opts.max_center_branch1_drop_frac);
screening.decisions.positive_shoulder_grouping = local_group_decision(items, names, ...
    opts.positive_shoulder_group_baseline_name, opts.positive_shoulder_group_candidate_name, ...
    opts.min_positive_shoulder_neg_area_improvement_frac, opts.max_center_neg_area_penalty_frac, ...
    'positive_shoulder', 'center_core', opts.max_positive_shoulder_branch1_drop_frac);

screening.decisions.center_grouping.enabled = screening.decisions.center_grouping.pairwise_beneficial && ...
    strcmp(screening.recommended_center_name, screening.decisions.center_grouping.candidate_name);
screening.decisions.positive_shoulder_grouping.enabled = screening.decisions.positive_shoulder_grouping.pairwise_beneficial && ...
    strcmp(screening.recommended_center_name, screening.decisions.positive_shoulder_grouping.baseline_name);

if screening.decisions.positive_shoulder_grouping.enabled
    screening.recommended_positive_shoulder_name = screening.decisions.positive_shoulder_grouping.candidate_name;
else
    screening.recommended_positive_shoulder_name = screening.recommended_center_name;
end
end


function opts = local_apply_defaults(opts)
opts = local_set_default(opts, 'center_safe_multiplier', 1.25);
opts = local_set_default(opts, 'center_group_baseline_name', 'dual_auto_plus_vii');
opts = local_set_default(opts, 'center_group_candidate_name', 'dual_auto_plus_vii_group_0p01');
opts = local_set_default(opts, 'positive_shoulder_group_baseline_name', 'dual_auto_plus_vii_group_0p01');
opts = local_set_default(opts, 'positive_shoulder_group_candidate_name', 'dual_auto_plus_vii_group_0p01_pos_shoulder');
opts = local_set_default(opts, 'min_center_neg_area_improvement_frac', 0.05);
opts = local_set_default(opts, 'max_center_branch1_drop_frac', 0.05);
opts = local_set_default(opts, 'min_positive_shoulder_neg_area_improvement_frac', 0.05);
opts = local_set_default(opts, 'max_center_neg_area_penalty_frac', 0.10);
opts = local_set_default(opts, 'max_positive_shoulder_branch1_drop_frac', 0.05);
end


function opts = local_set_default(opts, field_name, default_value)
if ~isfield(opts, field_name) || isempty(opts.(field_name))
    opts.(field_name) = default_value;
end
end


function items = local_normalize_comparisons(comparisons)
if isempty(comparisons)
    items = {};
    return;
end
if iscell(comparisons)
    items = comparisons;
elseif isstruct(comparisons)
    items = num2cell(comparisons);
else
    error('qe_recommend_lowq_background_configs:InvalidInput', ...
        'comparisons must be a cell array or struct array.');
end
end


function order = local_rank(metric_matrix)
if isempty(metric_matrix)
    order = [];
    return;
end
metric_matrix(~isfinite(metric_matrix)) = inf;
[~, order] = sortrows(metric_matrix);
end


function value = local_metric(item, path)
value = NaN;
current = item;
for i = 1:numel(path)
    key = path{i};
    if ~isstruct(current) || ~isfield(current, key)
        return;
    end
    current = current.(key);
end
if isnumeric(current) && isscalar(current)
    value = double(current);
end
end


function decision = local_group_decision(items, names, baseline_name, candidate_name, min_improve_frac, max_penalty_frac, primary_region_name, secondary_region_name, max_primary_branch1_drop_frac)
baseline_idx = find(strcmp(names, baseline_name), 1);
candidate_idx = find(strcmp(names, candidate_name), 1);

decision = struct();
decision.available = ~isempty(baseline_idx) && ~isempty(candidate_idx);
decision.enabled = false;
decision.pairwise_beneficial = false;
decision.baseline_name = baseline_name;
decision.candidate_name = candidate_name;
decision.primary_neg_area_improvement_frac = NaN;
decision.center_neg_area_improvement_frac = NaN;
decision.center_branch1_change_frac = NaN;
decision.primary_branch1_change_frac = NaN;
decision.center_neg_area_change_frac = NaN;

if ~decision.available
    return;
end

baseline_primary = local_metric(items{baseline_idx}, {'regions', primary_region_name, 'neg_area_fraction_mean'});
candidate_primary = local_metric(items{candidate_idx}, {'regions', primary_region_name, 'neg_area_fraction_mean'});
baseline_center = local_metric(items{baseline_idx}, {'regions', secondary_region_name, 'neg_area_fraction_mean'});
candidate_center = local_metric(items{candidate_idx}, {'regions', secondary_region_name, 'neg_area_fraction_mean'});
baseline_branch1 = local_metric(items{baseline_idx}, {'regions', 'center_core', 'branch1_peak_mean'});
candidate_branch1 = local_metric(items{candidate_idx}, {'regions', 'center_core', 'branch1_peak_mean'});
baseline_primary_branch1 = local_metric(items{baseline_idx}, {'regions', primary_region_name, 'branch1_peak_mean'});
candidate_primary_branch1 = local_metric(items{candidate_idx}, {'regions', primary_region_name, 'branch1_peak_mean'});

decision.primary_neg_area_improvement_frac = local_improvement_frac(baseline_primary, candidate_primary);
decision.center_neg_area_improvement_frac = local_improvement_frac(local_metric(items{baseline_idx}, {'regions', 'center_core', 'neg_area_fraction_mean'}), ...
    local_metric(items{candidate_idx}, {'regions', 'center_core', 'neg_area_fraction_mean'}));
decision.center_branch1_change_frac = local_change_frac(baseline_branch1, candidate_branch1);
decision.primary_branch1_change_frac = local_change_frac(baseline_primary_branch1, candidate_primary_branch1);
decision.center_neg_area_change_frac = local_change_frac(baseline_center, candidate_center);

if strcmp(primary_region_name, 'center_core')
    decision.pairwise_beneficial = decision.primary_neg_area_improvement_frac >= min_improve_frac && ...
        decision.center_branch1_change_frac >= -max_penalty_frac;
else
    decision.pairwise_beneficial = decision.primary_neg_area_improvement_frac >= min_improve_frac && ...
        decision.center_neg_area_change_frac <= max_penalty_frac && ...
        decision.primary_branch1_change_frac >= -max_primary_branch1_drop_frac;
end
end


function frac = local_improvement_frac(baseline_value, candidate_value)
if ~isfinite(baseline_value) || ~isfinite(candidate_value)
    frac = NaN;
    return;
end
if abs(baseline_value) < 1e-12
    if abs(candidate_value) < 1e-12
        frac = 0;
    else
        frac = NaN;
    end
    return;
end
scale = abs(baseline_value);
frac = (baseline_value - candidate_value) / scale;
end


function frac = local_change_frac(baseline_value, candidate_value)
if ~isfinite(baseline_value) || ~isfinite(candidate_value)
    frac = NaN;
    return;
end
if abs(baseline_value) < 1e-12
    if abs(candidate_value) < 1e-12
        frac = 0;
    else
        frac = NaN;
    end
    return;
end
scale = abs(baseline_value);
frac = (candidate_value - baseline_value) / scale;
end
