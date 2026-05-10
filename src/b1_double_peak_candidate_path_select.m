function selection = b1_double_peak_candidate_path_select(candidate_points, options)
%B1_DOUBLE_PEAK_CANDIDATE_PATH_SELECT Select a continuous B1 double-peak path.
%
% Candidate rows represent fitted lower/upper Lorentz pairs at one q/bin.
% Dynamic programming chooses one candidate per q while penalizing poor fits,
% edge picks, overbroad upper peaks, and discontinuous q-to-q jumps.

arguments
    candidate_points table
    options.sessionKey {mustBeTextScalar} = ""
    options.upperMediumJumpThresholdMeV (1,1) double = NaN
    options.largeJumpThresholdMeV (1,1) double = 250
    options.edgeEnergyMeV (1,2) double = [650 1900]
    options.upperMaxGammaOverE (1,1) double = 1.4
    options.upperMaxGammaMeV (1,1) double = 1600
    options.minPeakSeparationMeV (1,1) double = 10
    options.upperTrendMode {mustBeTextScalar} = "none"
    options.upperTrendAnchorQAbsAinv (1,1) double = 0.005
    options.upperTrendSmallQAbsAinv (1,1) double = 0.02
    options.upperTrendPlateauQAbsAinv (1,1) double = 0.06
    options.upperTrendAnchorMaxMeV (1,1) double = 700
    options.upperTrendLowUpperExemptMeV (1,2) double = [400 650]
    options.upperTrendDirectionToleranceMeV (1,1) double = 60
    options.upperTrendPlateauJumpThresholdMeV (1,1) double = 120
    options.upperTrendIsolatedJumpThresholdMeV (1,1) double = 120
    options.upperTrendAnchorPenalty (1,1) double = 250
    options.upperTrendWrongDirectionPenalty (1,1) double = 150
    options.upperTrendIsolatedJumpPenalty (1,1) double = 350
    options.upperTrendPlateauJumpPenalty (1,1) double = 150
end

if isempty(candidate_points) || height(candidate_points) == 0
    selection = struct('path_selection', local_empty_selection_table(), ...
        'selected_candidates', candidate_points);
    return
end
local_require_columns(candidate_points, {'candidate_id', 'q_Ainv', ...
    'lower_energy_meV', 'upper_energy_meV', 'R2'});

threshold = options.upperMediumJumpThresholdMeV;
if ~isfinite(threshold)
    threshold = local_default_medium_jump_threshold(options.sessionKey);
end

candidate_points = sortrows(candidate_points, {'q_Ainv', 'candidate_id'});
q_unique = unique(candidate_points.q_Ainv, 'stable');
nq = numel(q_unique);
groups = cell(nq, 1);
for qi = 1:nq
    groups{qi} = find(candidate_points.q_Ainv == q_unique(qi));
end

cost = cell(nq, 1);
prev = cell(nq, 1);
node_costs = cell(nq, 1);
for qi = 1:nq
    node_costs{qi} = local_node_cost(candidate_points(groups{qi}, :), ...
        options);
end

if local_use_rapidrise_plateau(options)
    [selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
        local_select_second_order(candidate_points, groups, node_costs, ...
        threshold, options);
else
    [selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
        local_select_first_order(candidate_points, groups, node_costs, ...
        threshold, options);
end

selected = candidate_points(selected_rows, :);
selection_tbl = local_selection_table(selected, node_cost, transition_cost, ...
    path_cost, threshold, options.largeJumpThresholdMeV, options, ...
    trend_penalty);
selection = struct('path_selection', selection_tbl, ...
    'selected_candidates', selected);
end


function local_require_columns(tbl, names)
missing = setdiff(names, tbl.Properties.VariableNames);
if ~isempty(missing)
    error('b1_double_peak_candidate_path_select:MissingColumn', ...
        'Candidate table missing required columns: %s', ...
        strjoin(missing, ', '));
end
end


function tf = local_use_rapidrise_plateau(options)
tf = strcmpi(char(string(options.upperTrendMode)), 'rapidrise_plateau');
end


function [selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
    local_select_first_order(candidate_points, groups, node_costs, ...
    threshold, options)
nq = numel(groups);
cost = cell(nq, 1);
prev = cell(nq, 1);
for qi = 1:nq
    idx = groups{qi};
    node = node_costs{qi};
    cost{qi} = inf(numel(idx), 1);
    prev{qi} = zeros(numel(idx), 1);
    if qi == 1
        cost{qi} = node;
        continue
    end
    prev_idx = groups{qi - 1};
    for ci = 1:numel(idx)
        trans = local_transition_cost(candidate_points(prev_idx, :), ...
            candidate_points(idx(ci), :), threshold, options);
        [best, best_pos] = min(cost{qi - 1} + trans);
        cost{qi}(ci) = node(ci) + best;
        prev{qi}(ci) = best_pos;
    end
end

[~, pos] = min(cost{nq});
selected_positions = zeros(nq, 1);
for qi = nq:-1:1
    selected_positions(qi) = pos;
    if qi > 1
        pos = prev{qi}(pos);
        if pos < 1
            pos = 1;
        end
    end
end
[selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
    local_cost_selected_path(candidate_points, groups, selected_positions, ...
    node_costs, threshold, options);
end


function [selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
    local_select_second_order(candidate_points, groups, node_costs, ...
    threshold, options)
nq = numel(groups);
if nq < 3
    [selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
        local_select_first_order(candidate_points, groups, node_costs, ...
        threshold, options);
    return
end

pair_cost = cell(nq, 1);
pair_prev = cell(nq, 1);
idx1 = groups{1};
idx2 = groups{2};
pair_cost{2} = inf(numel(idx1), numel(idx2));
pair_prev{2} = zeros(numel(idx1), numel(idx2));
for p1 = 1:numel(idx1)
    for p2 = 1:numel(idx2)
        trans = local_transition_cost(candidate_points(idx1(p1), :), ...
            candidate_points(idx2(p2), :), threshold, options);
        pair_cost{2}(p1, p2) = node_costs{1}(p1) + ...
            node_costs{2}(p2) + trans;
    end
end

for qi = 3:nq
    prevprev_idx = groups{qi - 2};
    prev_idx = groups{qi - 1};
    curr_idx = groups{qi};
    pair_cost{qi} = inf(numel(prev_idx), numel(curr_idx));
    pair_prev{qi} = zeros(numel(prev_idx), numel(curr_idx));
    for pp = 1:numel(prev_idx)
        for cc = 1:numel(curr_idx)
            previous = candidate_points(prev_idx(pp), :);
            current = candidate_points(curr_idx(cc), :);
            trans = local_transition_cost(previous, current, ...
                threshold, options);
            trend = zeros(numel(prevprev_idx), 1);
            for kk = 1:numel(prevprev_idx)
                trend(kk) = local_upper_trend_cost( ...
                    candidate_points(prevprev_idx(kk), :), previous, ...
                    current, options);
            end
            [best, best_pos] = min(pair_cost{qi - 1}(:, pp) + ...
                trans + trend);
            pair_cost{qi}(pp, cc) = node_costs{qi}(cc) + best;
            pair_prev{qi}(pp, cc) = best_pos;
        end
    end
end

[~, linear_pos] = min(pair_cost{nq}(:));
[prev_pos, curr_pos] = ind2sub(size(pair_cost{nq}), linear_pos);
selected_positions = zeros(nq, 1);
selected_positions(nq - 1) = prev_pos;
selected_positions(nq) = curr_pos;
for qi = nq:-1:3
    selected_positions(qi - 2) = pair_prev{qi}( ...
        selected_positions(qi - 1), selected_positions(qi));
    if selected_positions(qi - 2) < 1
        selected_positions(qi - 2) = 1;
    end
end

[selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
    local_cost_selected_path(candidate_points, groups, selected_positions, ...
    node_costs, threshold, options);
end


function [selected_rows, node_cost, transition_cost, path_cost, trend_penalty] = ...
    local_cost_selected_path(candidate_points, groups, selected_positions, ...
    node_costs, threshold, options)
nq = numel(groups);
selected_rows = zeros(nq, 1);
node_cost = NaN(nq, 1);
transition_cost = zeros(nq, 1);
trend_penalty = zeros(nq, 1);
path_cost = NaN(nq, 1);
for qi = 1:nq
    selected_rows(qi) = groups{qi}(selected_positions(qi));
    node_cost(qi) = node_costs{qi}(selected_positions(qi));
    if qi > 1
        previous = candidate_points(selected_rows(qi - 1), :);
        current = candidate_points(selected_rows(qi), :);
        transition_cost(qi) = local_transition_cost(previous, current, ...
            threshold, options);
    end
    if qi > 2
        trend_penalty(qi) = local_upper_trend_cost( ...
            candidate_points(selected_rows(qi - 2), :), ...
            candidate_points(selected_rows(qi - 1), :), ...
            candidate_points(selected_rows(qi), :), options);
    end
    if qi == 1
        path_cost(qi) = node_cost(qi);
    else
        path_cost(qi) = path_cost(qi - 1) + node_cost(qi) + ...
            transition_cost(qi) + trend_penalty(qi);
    end
end
end


function node = local_node_cost(tbl, options)
n = height(tbl);
node = inf(n, 1);
for i = 1:n
    lower = tbl.lower_energy_meV(i);
    upper = tbl.upper_energy_meV(i);
    if ~isfinite(lower) || ~isfinite(upper) || ...
            upper - lower < options.minPeakSeparationMeV
        continue
    end
    r2 = tbl.R2(i);
    if ~isfinite(r2)
        r2 = 0;
    end
    value = 100 * max(0, 1 - r2);
    edge = sort(options.edgeEnergyMeV);
    if (lower < edge(1) || lower > edge(2) || ...
            upper < edge(1) || upper > edge(2)) && ...
            ~local_small_q_low_upper_exempt(tbl, i, upper, options)
        value = value + 50;
    end
    if local_use_rapidrise_plateau(options)
        q_abs = local_candidate_q_abs(tbl, i);
        if q_abs <= options.upperTrendAnchorQAbsAinv && ...
                upper > options.upperTrendAnchorMaxMeV
            value = value + options.upperTrendAnchorPenalty;
        end
    end
    upper_gamma = local_numeric_column(tbl, 'upper_gamma_meV', i, NaN);
    gamma_over_E = upper_gamma ./ max(abs(upper), eps);
    if isfinite(gamma_over_E) && gamma_over_E > options.upperMaxGammaOverE
        value = value + 200;
    end
    if isfinite(upper_gamma) && upper_gamma > options.upperMaxGammaMeV
        value = value + 200;
    end
    node(i) = value;
end
end


function tf = local_small_q_low_upper_exempt(tbl, idx, upper, options)
tf = false;
if ~local_use_rapidrise_plateau(options)
    return
end
range = sort(options.upperTrendLowUpperExemptMeV);
q_abs = local_candidate_q_abs(tbl, idx);
tf = q_abs <= options.upperTrendSmallQAbsAinv && ...
    upper >= range(1) && upper <= range(2);
end


function trans = local_transition_cost(previous_tbl, current_tbl, ...
    medium_threshold, options)
trans = inf(height(previous_tbl), 1);
cur_lower = current_tbl.lower_energy_meV(1);
cur_upper = current_tbl.upper_energy_meV(1);
for i = 1:height(previous_tbl)
    lower_jump = abs(cur_lower - previous_tbl.lower_energy_meV(i));
    upper_jump = abs(cur_upper - previous_tbl.upper_energy_meV(i));
    value = 0;
    rapid_expected = local_expected_rapidrise_transition( ...
        previous_tbl(i, :), current_tbl, options);
    if upper_jump > medium_threshold && ~rapid_expected
        value = value + 50;
    end
    if local_plateau_transition(previous_tbl(i, :), current_tbl, options) && ...
            upper_jump > options.upperTrendPlateauJumpThresholdMeV
        value = value + options.upperTrendPlateauJumpPenalty;
    end
    if lower_jump > options.largeJumpThresholdMeV || ...
            (upper_jump > options.largeJumpThresholdMeV && ~rapid_expected)
        value = value + 300;
    end
    trans(i) = value;
end
end


function tf = local_expected_rapidrise_transition(previous_tbl, current_tbl, ...
    options)
tf = false;
if ~local_use_rapidrise_plateau(options)
    return
end
prev_q_abs = local_candidate_q_abs(previous_tbl, 1);
cur_q_abs = local_candidate_q_abs(current_tbl, 1);
if max(prev_q_abs, cur_q_abs) > options.upperTrendPlateauQAbsAinv
    return
end
delta_absq = cur_q_abs - prev_q_abs;
delta_upper = current_tbl.upper_energy_meV(1) - ...
    previous_tbl.upper_energy_meV(1);
tol = options.upperTrendDirectionToleranceMeV;
if abs(delta_absq) < eps
    tf = abs(delta_upper) <= options.largeJumpThresholdMeV;
elseif delta_absq > 0
    tf = delta_upper >= -tol;
else
    tf = delta_upper <= tol;
end
end


function tf = local_plateau_transition(previous_tbl, current_tbl, options)
tf = local_use_rapidrise_plateau(options) && ...
    local_candidate_q_abs(previous_tbl, 1) >= options.upperTrendPlateauQAbsAinv && ...
    local_candidate_q_abs(current_tbl, 1) >= options.upperTrendPlateauQAbsAinv;
end


function value = local_upper_trend_cost(prevprev_tbl, previous_tbl, ...
    current_tbl, options)
value = 0;
if ~local_use_rapidrise_plateau(options)
    return
end
q0 = local_candidate_q_abs(prevprev_tbl, 1);
q1 = local_candidate_q_abs(previous_tbl, 1);
q2 = local_candidate_q_abs(current_tbl, 1);
e0 = prevprev_tbl.upper_energy_meV(1);
e1 = previous_tbl.upper_energy_meV(1);
e2 = current_tbl.upper_energy_meV(1);
d1 = e1 - e0;
d2 = e2 - e1;
thr = options.upperTrendIsolatedJumpThresholdMeV;
if abs(d1) > thr && abs(d2) > thr && sign(d1) ~= sign(d2)
    is_center_valley = d1 < 0 && d2 > 0 && q1 <= options.upperTrendSmallQAbsAinv;
    if ~is_center_valley
        value = value + options.upperTrendIsolatedJumpPenalty;
    end
end

expected = sign(q2 - q1);
tol = options.upperTrendDirectionToleranceMeV;
if expected > 0 && d2 < -tol
    value = value + options.upperTrendWrongDirectionPenalty;
elseif expected < 0 && d2 > tol
    value = value + options.upperTrendWrongDirectionPenalty;
end

if min([q0 q1 q2]) >= options.upperTrendPlateauQAbsAinv && ...
        abs(d2) > options.upperTrendPlateauJumpThresholdMeV
    value = value + options.upperTrendPlateauJumpPenalty;
    if abs(d2) > options.largeJumpThresholdMeV
        value = value + 300;
    end
end
end


function value = local_numeric_column(tbl, name, idx, default_value)
value = default_value;
if ismember(name, tbl.Properties.VariableNames)
    value = tbl.(name)(idx);
end
end


function q_abs = local_candidate_q_abs(tbl, idx)
if ismember('q_abs_Ainv', tbl.Properties.VariableNames)
    q_abs = tbl.q_abs_Ainv(idx);
else
    q_abs = abs(tbl.q_Ainv(idx));
end
end


function tbl = local_selection_table(selected, node_cost, transition_cost, ...
    path_cost, medium_threshold, large_threshold, options, trend_penalty)
n = height(selected);
source = local_text_column(selected, 'candidate_source');
delta_upper = [NaN; diff(selected.upper_energy_meV)];
region = cell(n, 1);
for i = 1:n
    q_abs = local_candidate_q_abs(selected, i);
    if q_abs <= options.upperTrendSmallQAbsAinv
        region{i} = 'small_q_rapidrise';
    elseif q_abs <= options.upperTrendPlateauQAbsAinv
        region{i} = 'transition_to_plateau';
    else
        region{i} = 'plateau';
    end
end
mode = repmat({char(string(options.upperTrendMode))}, n, 1);
tbl = table(selected.candidate_id, selected.q_Ainv, ...
    selected.q_abs_Ainv, source, selected.lower_energy_meV, ...
    selected.upper_energy_meV, node_cost(:), transition_cost(:), ...
    trend_penalty(:), path_cost(:), delta_upper(:), region, mode, ...
    repmat(medium_threshold, n, 1), repmat(large_threshold, n, 1), ...
    'VariableNames', {'selected_candidate_id', 'q_Ainv', 'q_abs_Ainv', ...
    'candidate_source', 'lower_energy_meV', 'upper_energy_meV', ...
    'node_cost', 'transition_cost', 'upper_trend_penalty_meV', ...
    'path_cost', 'upper_delta_from_previous_meV', ...
    'upper_trend_region', 'upper_trend_mode', ...
    'upper_medium_jump_threshold_meV', 'large_jump_threshold_meV'});
end


function values = local_text_column(tbl, name)
if ismember(name, tbl.Properties.VariableNames)
    raw = tbl.(name);
else
    raw = repmat({''}, height(tbl), 1);
end
if isstring(raw)
    values = cellstr(raw);
elseif iscell(raw)
    values = raw;
else
    values = cellstr(string(raw));
end
end


function threshold = local_default_medium_jump_threshold(session_key)
text = lower(char(string(session_key)));
if contains(text, '20w')
    threshold = 120;
else
    threshold = 150;
end
end


function tbl = local_empty_selection_table()
names = {'selected_candidate_id', 'q_Ainv', 'q_abs_Ainv', ...
    'candidate_source', 'lower_energy_meV', 'upper_energy_meV', ...
    'node_cost', 'transition_cost', 'upper_trend_penalty_meV', ...
    'path_cost', 'upper_delta_from_previous_meV', ...
    'upper_trend_region', 'upper_trend_mode', ...
    'upper_medium_jump_threshold_meV', 'large_jump_threshold_meV'};
types = {'double', 'double', 'double', 'cell', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'cell', 'cell', ...
    'double', 'double'};
tbl = table('Size', [0 numel(names)], 'VariableTypes', types, ...
    'VariableNames', names);
end
