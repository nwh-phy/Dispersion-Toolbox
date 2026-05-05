function result = qe_assign_peak_branches_by_windows(peaks, specs, opts)
%QE_ASSIGN_PEAK_BRANCHES_BY_WINDOWS Assign auto-fit peaks to editable branches.
%   RESULT = QE_ASSIGN_PEAK_BRANCHES_BY_WINDOWS(PEAKS, SPECS) assigns each
%   qe_auto_fit peak to at most one branch according to SPECS(i).energy_window_meV.
%   Overlapping windows are resolved by nearest window center.
%
%   PEAKS follows the qe_auto_fit convention:
%   [q, E, Gamma, R2, A, E_lo, E_hi, G_lo, G_hi, A_lo, A_hi, raw_height]

if nargin < 1 || isempty(peaks)
    peaks = zeros(0, 12);
end
if nargin < 2
    specs = [];
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

opts = local_apply_defaults(opts);
specs = local_normalize_specs(specs);
peaks = local_normalize_peaks(peaks);

n_specs = numel(specs);
branches = cell(n_specs, 1);
rejected = local_empty_rejected_table();
candidate_n = zeros(n_specs, 1);
assigned_n = zeros(n_specs, 1);

if isempty(peaks)
    result = local_result(branches, specs, rejected, candidate_n, assigned_n);
    return
end

[usable, rejected_quality] = local_filter_quality(peaks, opts);
rejected = [rejected; rejected_quality];

if isempty(usable)
    result = local_result(branches, specs, rejected, candidate_n, assigned_n);
    return
end

[assignment, candidate_n] = local_assign_to_nearest_window(usable, specs);
unassigned = assignment == 0;
if any(unassigned)
    rejected = [rejected; local_reject_rows(usable(unassigned,:), 'none', 'outside_windows')];
end

for i = 1:n_specs
    rows = usable(assignment == i, :);
    [branches{i}, duplicate_rejected] = local_keep_best_per_signed_q(rows, specs(i), opts);
    rejected = [rejected; duplicate_rejected]; %#ok<AGROW>
    assigned_n(i) = size(branches{i}, 1);
end

result = local_result(branches, specs, rejected, candidate_n, assigned_n);
end


function opts = local_apply_defaults(opts)
opts = local_set_default(opts, 'q_skip_Ainv', 0);
opts = local_set_default(opts, 'min_R2', -Inf);
opts = local_set_default(opts, 'max_gamma_ratio', Inf);
opts = local_set_default(opts, 'score_column', 12);
end


function s = local_set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end


function specs = local_normalize_specs(specs)
if isempty(specs)
    specs = [
        struct('name', 'Branch 1', 'energy_window_meV', [500 2100], 'enabled', true)
        struct('name', 'Branch 2', 'energy_window_meV', [1700 2600], 'enabled', true)
        struct('name', 'Branch 3', 'energy_window_meV', [2800 3800], 'enabled', true)
        ];
end

if ~isstruct(specs)
    error('qe_assign_peak_branches_by_windows:invalidSpecs', ...
        'Branch specs must be a struct array.');
end

for i = 1:numel(specs)
    if ~isfield(specs(i), 'name') || strlength(string(specs(i).name)) == 0
        specs(i).name = sprintf('Branch %d', i);
    end
    if ~isfield(specs(i), 'energy_window_meV')
        error('qe_assign_peak_branches_by_windows:missingWindow', ...
            'Branch spec %d is missing energy_window_meV.', i);
    end
    win = double(specs(i).energy_window_meV);
    win = sort(win(:).');
    if numel(win) ~= 2 || any(~isfinite(win)) || win(1) == win(2)
        error('qe_assign_peak_branches_by_windows:invalidWindow', ...
            'Branch spec %d must have a finite two-value energy window.', i);
    end
    specs(i).energy_window_meV = win;
    if ~isfield(specs(i), 'enabled') || isempty(specs(i).enabled)
        specs(i).enabled = true;
    end
end
end


function peaks = local_normalize_peaks(peaks)
if size(peaks, 2) == 0
    peaks = zeros(size(peaks, 1), 12);
elseif size(peaks, 2) < 12
    peaks(:, end+1:12) = NaN;
end
end


function [usable, rejected] = local_filter_quality(peaks, opts)
rejected = local_empty_rejected_table();
if isempty(peaks)
    usable = peaks;
    return
end

q = peaks(:,1);
energy = peaks(:,2);
gamma = peaks(:,3);
r2 = peaks(:,4);
gamma_ratio = gamma ./ max(energy, eps);

bad = ~isfinite(q) | ~isfinite(energy) | energy <= 0 | ...
    abs(q) < opts.q_skip_Ainv | r2 < opts.min_R2 | ...
    gamma_ratio > opts.max_gamma_ratio;

if any(bad)
    rejected = [rejected; local_reject_rows(peaks(bad,:), 'all', 'quality')];
end
usable = peaks(~bad, :);
end


function [assignment, candidate_n] = local_assign_to_nearest_window(peaks, specs)
n_peaks = size(peaks, 1);
n_specs = numel(specs);
assignment = zeros(n_peaks, 1);
best_score = inf(n_peaks, 1);
candidate_n = zeros(n_specs, 1);

energy = peaks(:,2);
for i = 1:n_specs
    if ~specs(i).enabled
        continue
    end
    win = specs(i).energy_window_meV;
    in_window = energy >= win(1) & energy <= win(2);
    candidate_n(i) = sum(in_window);
    center = mean(win);
    half_width = max(diff(win) / 2, eps);
    score = abs(energy - center) ./ half_width;
    better = in_window & score < best_score;
    assignment(better) = i;
    best_score(better) = score(better);
end
end


function [kept, rejected] = local_keep_best_per_signed_q(rows, spec, opts)
rejected = local_empty_rejected_table();
if isempty(rows)
    kept = rows;
    return
end

q_values = unique(rows(:,1));
kept = [];
for i = 1:numel(q_values)
    group = rows(rows(:,1) == q_values(i), :);
    scores = local_peak_score(group, opts);
    [~, best_idx] = max(scores);
    kept = [kept; group(best_idx, :)]; %#ok<AGROW>

    reject_idx = true(size(group, 1), 1);
    reject_idx(best_idx) = false;
    if any(reject_idx)
        rejected = [rejected; local_reject_rows(group(reject_idx,:), spec.name, 'duplicate_signed_q')]; %#ok<AGROW>
    end
end
kept = sortrows(kept, 1);
end


function score = local_peak_score(rows, opts)
r2 = rows(:,4);
damping = min(1, rows(:,2) ./ max(rows(:,3), 1));
amp = rows(:,5);
if opts.score_column <= size(rows, 2)
    raw_height = rows(:, opts.score_column);
    use_raw = isfinite(raw_height) & raw_height > 0;
    amp(use_raw) = raw_height(use_raw);
end
amp_scale = max(abs(amp(isfinite(amp))), eps);
score = max(r2, 0) .* max(damping, 0) .* max(amp ./ amp_scale, eps);
score(~isfinite(score)) = 0;
end


function result = local_result(branches, specs, rejected, candidate_n, assigned_n)
names = string({specs.name}).';
windows = vertcat(specs.energy_window_meV);
result = struct();
result.branches = branches;
result.specs = specs;
result.rejected = rejected;
result.summary = table((1:numel(specs)).', names, windows(:,1), windows(:,2), ...
    candidate_n(:), assigned_n(:), ...
    'VariableNames', {'branch_index', 'name', 'energy_min_meV', 'energy_max_meV', ...
    'candidate_n', 'assigned_n'});
end


function rejected = local_reject_rows(rows, branch, reason)
rejected = local_empty_rejected_table();
if isempty(rows)
    return
end
n = size(rows, 1);
rejected = table(rows(:,1), rows(:,2), repmat({char(string(branch))}, n, 1), ...
    repmat({reason}, n, 1), ...
    'VariableNames', {'q_Ainv', 'energy_meV', 'branch', 'reason'});
end


function rejected = local_empty_rejected_table()
rejected = table(zeros(0,1), zeros(0,1), cell(0,1), cell(0,1), ...
    'VariableNames', {'q_Ainv', 'energy_meV', 'branch', 'reason'});
end
