function result = qe_select_low_energy_branch(peaks, opts)
%QE_SELECT_LOW_ENERGY_BRANCH Select a robust low-energy q-EELS branch.
%   RESULT = QE_SELECT_LOW_ENERGY_BRANCH(PEAKS, OPTS) filters candidate peaks
%   to the configured low-energy branch, chooses the best continuous segment
%   on each signed-q side, and only enables symmetry averaging when the two
%   sides agree in the overlapping |q| range.
%
%   PEAKS follows the qe_auto_fit convention:
%   [q, E, Gamma, R2, A, E_lo, E_hi, G_lo, G_hi, A_lo, A_hi, raw_height]

if nargin < 2 || isempty(opts)
    opts = thesis_config().branch;
end
opts = local_apply_defaults(opts);

if isempty(peaks)
    peaks = zeros(0, 12);
end
if size(peaks, 2) < 12
    peaks(:, end+1:12) = NaN;
end

[low_candidates, rejected] = local_filter_low_candidates(peaks, opts);
[deduped, duplicate_rejected] = local_keep_best_per_signed_q(low_candidates, opts);
rejected = [rejected; duplicate_rejected];

neg = local_select_best_side_segment(deduped(deduped(:,1) < 0, :), opts, -1);
pos = local_select_best_side_segment(deduped(deduped(:,1) > 0, :), opts, 1);
combined = local_evaluate_combined_branch(neg, pos, opts);

if combined.trusted
    physics_branch = combined.rows;
    selected_side = 0;
    selected_fit = combined.fit;
    selected_fit_R2 = combined.fit_R2;
    symmetry_average = true;
else
    selected = local_choose_best_side(neg, pos);
    if isempty(selected.rows)
        error('qe_select_low_energy_branch:noLowBranch', ...
            'No usable low-energy branch segment remains after filtering.');
    end
    physics_branch = selected.rows;
    selected_side = selected.side;
    selected_fit = selected.fit;
    selected_fit_R2 = selected.fit_R2;
    symmetry_average = false;
end

physics_branch = local_sort_by_signed_q(physics_branch);

result = struct();
result.physics_branch = physics_branch;
result.accepted = physics_branch;
result.low_candidates = low_candidates;
result.deduped_candidates = deduped;
result.rejected = rejected;
result.symmetry_average = symmetry_average;
result.negative = neg;
result.positive = pos;
result.combined = combined;
result.quality = struct( ...
    'selected_side', selected_side, ...
    'selected_fit_R2', selected_fit_R2, ...
    'combined_fit_R2', combined.fit_R2, ...
    'symmetry_trusted', combined.trusted, ...
    'symmetry_overlap_n', combined.overlap_n, ...
    'symmetry_median_abs_delta_meV', combined.median_abs_delta_meV, ...
    'symmetry_median_rel_delta', combined.median_rel_delta, ...
    'selected_fit', selected_fit);
result.summary = struct( ...
    'candidate_n', size(low_candidates, 1), ...
    'deduped_n', size(deduped, 1), ...
    'physics_n', size(physics_branch, 1), ...
    'negative_segment_n', size(neg.rows, 1), ...
    'positive_segment_n', size(pos.rows, 1), ...
    'rejected_n', height(rejected));
result.opts = opts;
end


function opts = local_apply_defaults(opts)
if ~isfield(opts, 'q_skip_Ainv'), opts.q_skip_Ainv = 0.005; end
if ~isfield(opts, 'min_R2'), opts.min_R2 = 0.3; end
if ~isfield(opts, 'max_gamma_ratio'), opts.max_gamma_ratio = 2.0; end
if ~isfield(opts, 'score_column'), opts.score_column = 12; end
if ~isfield(opts, 'low') || ~isfield(opts.low, 'energy_window_meV')
    opts.low.energy_window_meV = [500 2100];
end
if ~isfield(opts.low, 'name'), opts.low.name = 'low'; end
if ~isfield(opts.low, 'max_delta_meV'), opts.low.max_delta_meV = 350; end
if ~isfield(opts.low, 'max_q_gap_Ainv'), opts.low.max_q_gap_Ainv = 0.015; end
if ~isfield(opts, 'selection'), opts.selection = struct(); end
opts.selection = local_set_default(opts.selection, 'min_segment_points', 4);
opts.selection = local_set_default(opts.selection, 'min_fit_points', 3);
opts.selection = local_set_default(opts.selection, 'q_min_fit_Ainv', 0.01);
opts.selection = local_set_default(opts.selection, 'epsilon_s', 1);
opts.selection = local_set_default(opts.selection, 'rho0_init', 25);
opts.selection = local_set_default(opts.selection, 'min_side_fit_R2', 0.70);
opts.selection = local_set_default(opts.selection, 'min_symmetry_overlap', 4);
opts.selection = local_set_default(opts.selection, 'q_match_tolerance_Ainv', 1e-9);
opts.selection = local_set_default(opts.selection, 'max_symmetry_abs_delta_meV', 150);
opts.selection = local_set_default(opts.selection, 'max_symmetry_rel_delta', 0.15);
opts.selection = local_set_default(opts.selection, 'max_combined_R2_drop', 0.10);
end


function s = local_set_default(s, name, value)
if ~isfield(s, name)
    s.(name) = value;
end
end


function [low_candidates, rejected] = local_filter_low_candidates(peaks, opts)
rejected = local_empty_rejected_table();
if isempty(peaks)
    low_candidates = peaks;
    return
end

q_abs = abs(peaks(:,1));
gamma_ratio = peaks(:,3) ./ max(peaks(:,2), eps);
bad_q = q_abs < opts.q_skip_Ainv;
bad_r2 = peaks(:,4) < opts.min_R2;
bad_gamma = gamma_ratio > opts.max_gamma_ratio;

quality_bad = bad_q | bad_r2 | bad_gamma;
if any(bad_q)
    rejected = [rejected; local_reject_rows(peaks(bad_q,:), 'quality', 'q_skip')];
end
if any(bad_r2 & ~bad_q)
    rejected = [rejected; local_reject_rows(peaks(bad_r2 & ~bad_q,:), 'quality', 'low_R2')];
end
if any(bad_gamma & ~bad_q & ~bad_r2)
    rejected = [rejected; local_reject_rows(peaks(bad_gamma & ~bad_q & ~bad_r2,:), 'quality', 'overdamped')];
end

usable = peaks(~quality_bad, :);
win = opts.low.energy_window_meV;
in_low = usable(:,2) >= win(1) & usable(:,2) <= win(2);
if any(~in_low)
    rejected = [rejected; local_reject_rows(usable(~in_low,:), opts.low.name, 'outside_low_window')];
end
low_candidates = usable(in_low, :);
end


function [kept, rejected] = local_keep_best_per_signed_q(candidates, opts)
rejected = local_empty_rejected_table();
if isempty(candidates)
    kept = candidates;
    return
end

q_values = unique(candidates(:,1));
kept = [];
for i = 1:numel(q_values)
    group = candidates(candidates(:,1) == q_values(i), :);
    scores = local_peak_score(group, opts);
    [~, best_idx] = max(scores);
    kept = [kept; group(best_idx, :)]; %#ok<AGROW>

    reject_idx = true(size(group, 1), 1);
    reject_idx(best_idx) = false;
    if any(reject_idx)
        rejected = [rejected; local_reject_rows(group(reject_idx,:), opts.low.name, 'duplicate_signed_q')]; %#ok<AGROW>
    end
end
kept = local_sort_by_signed_q(kept);
end


function score = local_peak_score(rows, opts)
r2 = rows(:,4);
damping = min(1, rows(:,2) ./ max(rows(:,3), 1));
amp = rows(:,5);
if isfield(opts, 'score_column') && opts.score_column <= size(rows,2)
    raw_h = rows(:, opts.score_column);
    use_raw = isfinite(raw_h) & raw_h > 0;
    amp(use_raw) = raw_h(use_raw);
end
amp_scale = max(abs(amp(isfinite(amp))), eps);
score = max(r2, 0) .* max(damping, 0) .* max(amp ./ amp_scale, eps);
score(~isfinite(score)) = 0;
end


function side = local_select_best_side_segment(rows, opts, side_value)
side = local_empty_side(side_value);
if isempty(rows)
    return
end

[~, order] = sort(abs(rows(:,1)), 'ascend');
rows = rows(order, :);
segments = local_split_continuous_segments(rows, opts);
side.n_segments = numel(segments);

best_score = -Inf;
best = local_empty_side(side_value);
for i = 1:numel(segments)
    candidate = local_score_segment(segments{i}, opts, side_value);
    if candidate.score > best_score
        best_score = candidate.score;
        best = candidate;
    end
end
side = best;
side.n_segments = numel(segments);
end


function segments = local_split_continuous_segments(rows, opts)
segments = {};
if isempty(rows)
    return
end

current = rows(1, :);
for i = 2:size(rows, 1)
    q_gap = abs(abs(rows(i,1)) - abs(rows(i-1,1)));
    energy_jump = abs(rows(i,2) - rows(i-1,2));
    if q_gap > opts.low.max_q_gap_Ainv || energy_jump > opts.low.max_delta_meV
        segments{end+1} = current; %#ok<AGROW>
        current = rows(i, :);
    else
        current = [current; rows(i, :)]; %#ok<AGROW>
    end
end
segments{end+1} = current;
end


function segment = local_score_segment(rows, opts, side_value)
segment = local_empty_side(side_value);
segment.rows = rows;
segment.n = size(rows, 1);
segment.q_span_Ainv = max(abs(rows(:,1))) - min(abs(rows(:,1)));
segment.mean_R2 = mean(rows(:,4), 'omitnan');
[fit, fit_R2, fit_n] = local_fit_segment(rows, opts);
segment.fit = fit;
segment.fit_R2 = fit_R2;
segment.fit_n = fit_n;

if fit_n < opts.selection.min_fit_points || segment.n < opts.selection.min_segment_points
    fit_term = -100;
else
    fit_term = fit_R2;
end
if ~isfinite(fit_term)
    fit_term = -100;
end

segment.score = 1000 * fit_term + 20 * fit_n + ...
    100 * segment.q_span_Ainv + max(segment.mean_R2, 0);
end


function [fit, fit_R2, fit_n] = local_fit_segment(rows, opts)
fit = struct();
fit_R2 = NaN;
q_abs = abs(rows(:,1));
mask = q_abs >= opts.selection.q_min_fit_Ainv & ...
    isfinite(q_abs) & isfinite(rows(:,2)) & q_abs > 0 & rows(:,2) > 0;
fit_n = sum(mask);
if fit_n < opts.selection.min_fit_points
    return
end

q_fit = q_abs(mask);
E_fit = rows(mask, 2);
weights = rows(mask, 4);
weights(~isfinite(weights) | weights <= 0) = 1;

try
    evalc(['fit = fit_quasi2d_plasmon(q_fit, E_fit, ', ...
        'confidence=weights, ', ...
        'epsilon_s=opts.selection.epsilon_s, ', ...
        'rho0_init=opts.selection.rho0_init);']);
    fit_R2 = fit.R_squared;
catch
    fit = struct();
    fit_R2 = NaN;
end
end


function combined = local_evaluate_combined_branch(neg, pos, opts)
combined = struct();
combined.rows = local_sort_by_signed_q([neg.rows; pos.rows]);
combined.fit = struct();
combined.fit_R2 = NaN;
combined.trusted = false;
combined.overlap_n = 0;
combined.median_abs_delta_meV = NaN;
combined.median_rel_delta = NaN;

if isempty(neg.rows) || isempty(pos.rows)
    return
end

[overlap_n, abs_delta, rel_delta] = local_symmetry_deltas(neg.rows, pos.rows, opts);
combined.overlap_n = overlap_n;
combined.median_abs_delta_meV = median(abs_delta, 'omitnan');
combined.median_rel_delta = median(rel_delta, 'omitnan');

if neg.fit_R2 < opts.selection.min_side_fit_R2 || pos.fit_R2 < opts.selection.min_side_fit_R2
    return
end

if overlap_n < opts.selection.min_symmetry_overlap || ...
        combined.median_abs_delta_meV > opts.selection.max_symmetry_abs_delta_meV || ...
        combined.median_rel_delta > opts.selection.max_symmetry_rel_delta
    return
end

[fit, fit_R2] = local_fit_segment(combined.rows, opts);
combined.fit = fit;
combined.fit_R2 = fit_R2;
min_side_R2 = min(neg.fit_R2, pos.fit_R2);
combined.trusted = isfinite(fit_R2) && ...
    fit_R2 >= min_side_R2 - opts.selection.max_combined_R2_drop;
end


function [overlap_n, abs_delta, rel_delta] = local_symmetry_deltas(neg_rows, pos_rows, opts)
tol = opts.selection.q_match_tolerance_Ainv;
neg_key = round(abs(neg_rows(:,1)) ./ tol) .* tol;
pos_key = round(abs(pos_rows(:,1)) ./ tol) .* tol;
[common_key, neg_idx, pos_idx] = intersect(neg_key, pos_key);
overlap_n = numel(common_key);
if overlap_n == 0
    abs_delta = NaN;
    rel_delta = NaN;
    return
end

E_neg = neg_rows(neg_idx, 2);
E_pos = pos_rows(pos_idx, 2);
abs_delta = abs(E_neg - E_pos);
rel_delta = abs_delta ./ max((abs(E_neg) + abs(E_pos)) ./ 2, eps);
end


function selected = local_choose_best_side(neg, pos)
if isempty(neg.rows)
    selected = pos;
    return
end
if isempty(pos.rows)
    selected = neg;
    return
end
if neg.score >= pos.score
    selected = neg;
else
    selected = pos;
end
end


function side = local_empty_side(side_value)
side = struct();
side.side = side_value;
side.rows = zeros(0, 12);
side.fit = struct();
side.fit_R2 = NaN;
side.fit_n = 0;
side.n = 0;
side.n_segments = 0;
side.q_span_Ainv = NaN;
side.mean_R2 = NaN;
side.score = -Inf;
end


function rows = local_sort_by_signed_q(rows)
if isempty(rows)
    return
end
[~, order] = sort(rows(:,1), 'ascend');
rows = rows(order, :);
end


function rejected = local_reject_rows(rows, branch, reason)
rejected = local_empty_rejected_table();
if isempty(rows)
    return
end
n = size(rows, 1);
rejected = table(rows(:,1), rows(:,2), repmat({branch}, n, 1), repmat({reason}, n, 1), ...
    'VariableNames', {'q_Ainv', 'energy_meV', 'branch', 'reason'});
end


function rejected = local_empty_rejected_table()
rejected = table(zeros(0,1), zeros(0,1), cell(0,1), cell(0,1), ...
    'VariableNames', {'q_Ainv', 'energy_meV', 'branch', 'reason'});
end
