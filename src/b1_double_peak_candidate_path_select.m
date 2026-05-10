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
    idx = groups{qi};
    node = local_node_cost(candidate_points(idx, :), options);
    node_costs{qi} = node;
    cost{qi} = inf(numel(idx), 1);
    prev{qi} = zeros(numel(idx), 1);
    if qi == 1
        cost{qi} = node;
        continue
    end
    prev_idx = groups{qi - 1};
    for ci = 1:numel(idx)
        trans = local_transition_cost(candidate_points(prev_idx, :), ...
            candidate_points(idx(ci), :), threshold, ...
            options.largeJumpThresholdMeV);
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

selected_rows = zeros(nq, 1);
node_cost = NaN(nq, 1);
transition_cost = zeros(nq, 1);
path_cost = NaN(nq, 1);
for qi = 1:nq
    selected_rows(qi) = groups{qi}(selected_positions(qi));
    node_cost(qi) = node_costs{qi}(selected_positions(qi));
    path_cost(qi) = cost{qi}(selected_positions(qi));
    if qi > 1
        previous = candidate_points(selected_rows(qi - 1), :);
        current = candidate_points(selected_rows(qi), :);
        transition_cost(qi) = local_transition_cost(previous, current, ...
            threshold, options.largeJumpThresholdMeV);
    end
end

selected = candidate_points(selected_rows, :);
selection_tbl = local_selection_table(selected, node_cost, transition_cost, ...
    path_cost, threshold, options.largeJumpThresholdMeV);
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
    if lower < edge(1) || lower > edge(2) || ...
            upper < edge(1) || upper > edge(2)
        value = value + 50;
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


function trans = local_transition_cost(previous_tbl, current_tbl, ...
    medium_threshold, large_threshold)
trans = inf(height(previous_tbl), 1);
cur_lower = current_tbl.lower_energy_meV(1);
cur_upper = current_tbl.upper_energy_meV(1);
for i = 1:height(previous_tbl)
    lower_jump = abs(cur_lower - previous_tbl.lower_energy_meV(i));
    upper_jump = abs(cur_upper - previous_tbl.upper_energy_meV(i));
    value = 0;
    if upper_jump > medium_threshold
        value = value + 50;
    end
    if lower_jump > large_threshold || upper_jump > large_threshold
        value = value + 300;
    end
    trans(i) = value;
end
end


function value = local_numeric_column(tbl, name, idx, default_value)
value = default_value;
if ismember(name, tbl.Properties.VariableNames)
    value = tbl.(name)(idx);
end
end


function tbl = local_selection_table(selected, node_cost, transition_cost, ...
    path_cost, medium_threshold, large_threshold)
n = height(selected);
source = local_text_column(selected, 'candidate_source');
tbl = table(selected.candidate_id, selected.q_Ainv, ...
    selected.q_abs_Ainv, source, selected.lower_energy_meV, ...
    selected.upper_energy_meV, node_cost(:), transition_cost(:), ...
    path_cost(:), repmat(medium_threshold, n, 1), ...
    repmat(large_threshold, n, 1), ...
    'VariableNames', {'selected_candidate_id', 'q_Ainv', 'q_abs_Ainv', ...
    'candidate_source', 'lower_energy_meV', 'upper_energy_meV', ...
    'node_cost', 'transition_cost', 'path_cost', ...
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
    'node_cost', 'transition_cost', 'path_cost', ...
    'upper_medium_jump_threshold_meV', 'large_jump_threshold_meV'};
types = {'double', 'double', 'double', 'cell', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double'};
tbl = table('Size', [0 numel(names)], 'VariableTypes', types, ...
    'VariableNames', names);
end
