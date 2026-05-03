function branches = assign_qe_branches(peaks, opts)
%ASSIGN_QE_BRANCHES Deterministic thesis branch assignment for q-EELS peaks.
%   BRANCHES = ASSIGN_QE_BRANCHES(PEAKS, OPTS) separates candidate peaks into
%   low/high branches using explicit energy windows, per-signed-q competition,
%   continuity constraints, and basic quality filters.  The function is meant
%   to replace ad-hoc k-means / global-gap / fixed-threshold branch splitting
%   in thesis-facing scripts.
%
%   PEAKS uses the qe_auto_fit convention:
%   [q, E, Gamma, R2, A, E_lo, E_hi, G_lo, G_hi, A_lo, A_hi, raw_height]

if nargin < 2 || isempty(opts)
    opts = thesis_config().branch;
end
if isempty(peaks)
    peaks = zeros(0, 12);
end
if size(peaks, 2) < 12
    peaks(:, end+1:12) = NaN;
end

[peaks, pre_rejected] = local_quality_filter(peaks, opts);
all_windows = [opts.low.energy_window_meV; opts.high.energy_window_meV];
outside_mask = ~local_in_any_window(peaks(:,2), all_windows);
outside_rejected = local_reject_rows(peaks(outside_mask, :), 'none', 'outside_energy_window');
usable = peaks(~outside_mask, :);

[low, low_rej] = local_assign_one_branch(usable, opts.low, opts);
[high, high_rej] = local_assign_one_branch(usable, opts.high, opts);

branches = struct();
branches.low = low;
branches.high = high;
branches.labels = {'low', 'high'};
branches.rejected = [pre_rejected; outside_rejected; low_rej; high_rej];
branches.summary = struct( ...
    'low_n', size(low.accepted, 1), ...
    'high_n', size(high.accepted, 1), ...
    'rejected_n', height(branches.rejected));
end


function [peaks_out, rejected] = local_quality_filter(peaks, opts)
rejected = local_empty_rejected_table();
if isempty(peaks)
    peaks_out = peaks;
    return
end

bad_q = abs(peaks(:,1)) < opts.q_skip_Ainv;
bad_r2 = peaks(:,4) < opts.min_R2;
bad_gamma = peaks(:,3) ./ max(peaks(:,2), eps) > opts.max_gamma_ratio;
bad = bad_q | bad_r2 | bad_gamma;

if any(bad_q)
    rejected = [rejected; local_reject_rows(peaks(bad_q,:), 'quality', 'q_skip')];
end
if any(bad_r2 & ~bad_q)
    rejected = [rejected; local_reject_rows(peaks(bad_r2 & ~bad_q,:), 'quality', 'low_R2')];
end
if any(bad_gamma & ~bad_q & ~bad_r2)
    rejected = [rejected; local_reject_rows(peaks(bad_gamma & ~bad_q & ~bad_r2,:), 'quality', 'overdamped')];
end
peaks_out = peaks(~bad, :);
end


function tf = local_in_any_window(energy, windows)
tf = false(size(energy));
for i = 1:size(windows, 1)
    tf = tf | (energy >= windows(i,1) & energy <= windows(i,2));
end
end


function [branch, rejected] = local_assign_one_branch(peaks, branch_opts, global_opts)
label = branch_opts.name;
win = branch_opts.energy_window_meV;
in_window = peaks(:,2) >= win(1) & peaks(:,2) <= win(2);
candidates = peaks(in_window, :);
rejected = local_empty_rejected_table();

if isempty(candidates)
    branch = local_empty_branch(label, win);
    return
end

[candidates, dup_rej] = local_keep_best_per_signed_q(candidates, label, global_opts);
rejected = [rejected; dup_rej];

accepted = [];
for sign_value = [-1 1]
    if sign_value < 0
        side_mask = candidates(:,1) < 0;
    else
        side_mask = candidates(:,1) > 0;
    end
    side = candidates(side_mask, :);
    if isempty(side)
        continue
    end
    [~, order] = sort(abs(side(:,1)), 'ascend');
    side = side(order, :);
    side_accept = [];
    for i = 1:size(side, 1)
        row = side(i, :);
        if isempty(side_accept)
            side_accept = row;
            continue
        end
        prev = side_accept(end, :);
        energy_jump = abs(row(2) - prev(2));
        q_gap = abs(abs(row(1)) - abs(prev(1)));
        if energy_jump > branch_opts.max_delta_meV || q_gap > branch_opts.max_q_gap_Ainv
            rejected = [rejected; local_reject_rows(row, label, 'continuity_jump')]; %#ok<AGROW>
        else
            side_accept = [side_accept; row]; %#ok<AGROW>
        end
    end
    accepted = [accepted; side_accept]; %#ok<AGROW>
end

if ~isempty(accepted)
    [~, order] = sort(accepted(:,1), 'ascend');
    accepted = accepted(order, :);
end

branch = struct();
branch.name = label;
branch.energy_window_meV = win;
branch.accepted = accepted;
branch.rejected = rejected;
branch.symmetrized = zeros(0, 12);
branch.options = branch_opts;
end


function branch = local_empty_branch(label, win)
branch = struct();
branch.name = label;
branch.energy_window_meV = win;
branch.accepted = zeros(0, 12);
branch.rejected = local_empty_rejected_table();
branch.symmetrized = zeros(0, 12);
branch.options = struct();
end


function [kept, rejected] = local_keep_best_per_signed_q(candidates, label, opts)
rejected = local_empty_rejected_table();
if isempty(candidates)
    kept = candidates;
    return
end
q_values = unique(candidates(:,1));
kept = [];
for i = 1:numel(q_values)
    mask = candidates(:,1) == q_values(i);
    group = candidates(mask, :);
    scores = local_peak_score(group, opts);
    [~, best_idx] = max(scores);
    kept = [kept; group(best_idx, :)]; %#ok<AGROW>
    reject_idx = true(size(group, 1), 1);
    reject_idx(best_idx) = false;
    if any(reject_idx)
        rejected = [rejected; local_reject_rows(group(reject_idx,:), label, 'duplicate_signed_q')]; %#ok<AGROW>
    end
end
[~, order] = sort(kept(:,1), 'ascend');
kept = kept(order, :);
end


function score = local_peak_score(rows, opts)
if isempty(rows)
    score = [];
    return
end
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
