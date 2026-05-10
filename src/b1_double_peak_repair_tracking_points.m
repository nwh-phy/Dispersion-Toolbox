function out = b1_double_peak_repair_tracking_points(lower_points, upper_points, ...
    failures, options)
%B1_DOUBLE_PEAK_REPAIR_TRACKING_POINTS Repair logged B1 tracking outliers.
%
% This is a deterministic post-fit repair pass. It never removes points and
% never substitutes a single-peak result. Large energy jumps are replaced by
% same-branch interpolation from neighboring q points and every replacement is
% recorded in a repair log.

arguments
    lower_points table
    upper_points table
    failures table = table()
    options.largeJumpThresholdMeV (1,1) double = 250
end

threshold = max(1, double(options.largeJumpThresholdMeV));
lower_points = local_ensure_repair_columns(lower_points);
upper_points = local_ensure_repair_columns(upper_points);

[lower_points, lower_log] = local_repair_branch(lower_points, threshold);
[upper_points, upper_log] = local_repair_branch(upper_points, threshold);

if height(lower_points) == 0
    combined_points = upper_points;
elseif height(upper_points) == 0
    combined_points = lower_points;
else
    combined_points = [lower_points; upper_points];
end
if ~isempty(combined_points) && all(ismember({'q_Ainv', 'branch', ...
        'energy_meV'}, combined_points.Properties.VariableNames))
    combined_points = sortrows(combined_points, {'q_Ainv', 'branch', ...
        'energy_meV'});
end

out = struct();
out.lower_points = lower_points;
out.upper_points = upper_points;
out.combined_points = combined_points;
out.repair_log = [lower_log; upper_log];
out.exclusion_points = local_failure_exclusions(failures);
end


function points = local_ensure_repair_columns(points)
if isempty(points) || ~istable(points)
    points = local_empty_points_like();
    return
end
n = height(points);
if ~ismember('repair_source', points.Properties.VariableNames)
    points.repair_source = repmat({'raw_fit'}, n, 1);
end
if ~ismember('original_energy_meV', points.Properties.VariableNames)
    if ismember('energy_meV', points.Properties.VariableNames)
        points.original_energy_meV = points.energy_meV;
    else
        points.original_energy_meV = NaN(n, 1);
    end
end
if ~ismember('repair_detail', points.Properties.VariableNames)
    points.repair_detail = repmat({''}, n, 1);
end
end


function [points, repair_log] = local_repair_branch(points, threshold)
repair_log = local_empty_repair_log_table();
if height(points) < 3 || ~all(ismember({'q_Ainv', 'energy_meV'}, ...
        points.Properties.VariableNames))
    return
end

[~, order] = sort(points.q_Ainv);
points = points(order, :);
for pass = 1:2
    changed = false;
    energy = points.energy_meV;
    for i = 1:height(points)
        [should_repair, repaired_energy, detail] = local_repair_candidate( ...
            points.q_Ainv, energy, i, threshold);
        if ~should_repair
            continue
        end
        original_energy = points.energy_meV(i);
        points.energy_meV(i) = repaired_energy;
        points.original_energy_meV(i) = original_energy;
        points.repair_source{i} = 'jump_repair_interpolated_energy';
        points.repair_detail{i} = detail;
        points = local_shift_energy_ci(points, i, repaired_energy - original_energy);
        repair_log = [repair_log; local_repair_log_row(points(i, :), ...
            original_energy, repaired_energy, detail)]; %#ok<AGROW>
        changed = true;
    end
    if ~changed
        break
    end
end
end


function [tf, repaired_energy, detail] = local_repair_candidate(q, energy, idx, ...
    threshold)
tf = false;
repaired_energy = NaN;
detail = '';
n = numel(energy);
if ~isfinite(energy(idx))
    return
end
if idx > 1 && idx < n && isfinite(energy(idx - 1)) && ...
        isfinite(energy(idx + 1))
    left_jump = abs(energy(idx) - energy(idx - 1));
    right_jump = abs(energy(idx) - energy(idx + 1));
    neighbor_jump = abs(energy(idx + 1) - energy(idx - 1));
    if left_jump > threshold && right_jump > threshold && ...
            neighbor_jump <= threshold
        repaired_energy = interp1(q([idx - 1, idx + 1]), ...
            energy([idx - 1, idx + 1]), q(idx), 'linear');
        tf = isfinite(repaired_energy);
        detail = sprintf(['internal point replaced by same-branch ', ...
            'linear interpolation; jumps %.3g/%.3g meV'], ...
            left_jump, right_jump);
        return
    end
end
if idx == 1 && n >= 3 && isfinite(energy(2)) && isfinite(energy(3))
    first_jump = abs(energy(1) - energy(2));
    if first_jump > threshold && abs(energy(2) - energy(3)) <= threshold
        repaired_energy = interp1(q(2:3), energy(2:3), q(1), ...
            'linear', 'extrap');
        tf = isfinite(repaired_energy);
        detail = sprintf('endpoint extrapolated from next two points; jump %.3g meV', ...
            first_jump);
    end
elseif idx == n && n >= 3 && isfinite(energy(n - 1)) && ...
        isfinite(energy(n - 2))
    last_jump = abs(energy(n) - energy(n - 1));
    if last_jump > threshold && abs(energy(n - 1) - energy(n - 2)) <= threshold
        repaired_energy = interp1(q((n - 2):(n - 1)), ...
            energy((n - 2):(n - 1)), q(n), 'linear', 'extrap');
        tf = isfinite(repaired_energy);
        detail = sprintf('endpoint extrapolated from previous two points; jump %.3g meV', ...
            last_jump);
    end
end
end


function points = local_shift_energy_ci(points, idx, delta)
for name = {'E_ci_lo', 'E_ci_hi'}
    field = name{1};
    if ismember(field, points.Properties.VariableNames) && isfinite(delta)
        points.(field)(idx) = points.(field)(idx) + delta;
    end
end
end


function exclusions = local_failure_exclusions(failures)
exclusions = local_empty_exclusions_table();
if isempty(failures) || ~istable(failures) || height(failures) == 0
    return
end
for i = 1:height(failures)
    branch = NaN;
    branch_label = '';
    exclusions = [exclusions; table(failures.q_Ainv(i), ...
        failures.q_abs_Ainv(i), branch, {branch_label}, ...
        {'unrepaired_failure'}, NaN, NaN, failures.source_mode(i), ...
        failures.source_q_count(i), failures.source_q_Ainv(i), ...
        failures.source_q_index(i), failures.status(i), ...
        'VariableNames', exclusions.Properties.VariableNames)]; %#ok<AGROW>
end
end


function row = local_repair_log_row(point, original_energy, repaired_energy, ...
    detail)
row = table(point.q_Ainv(1), point.q_abs_Ainv(1), point.branch(1), ...
    point.branch_label(1), point.repair_source(1), original_energy, ...
    repaired_energy, point.source_mode(1), point.source_q_count(1), ...
    point.source_q_Ainv(1), point.source_q_index(1), {detail}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'branch', ...
    'branch_label', 'repair_source', 'original_energy_meV', ...
    'repaired_energy_meV', 'source_mode', 'source_q_count', ...
    'source_q_Ainv', 'source_q_index', 'detail'});
end


function tbl = local_empty_repair_log_table()
tbl = table('Size', [0 12], ...
    'VariableTypes', {'double', 'double', 'double', 'cell', 'cell', ...
    'double', 'double', 'cell', 'double', 'cell', 'cell', 'cell'}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'branch', ...
    'branch_label', 'repair_source', 'original_energy_meV', ...
    'repaired_energy_meV', 'source_mode', 'source_q_count', ...
    'source_q_Ainv', 'source_q_index', 'detail'});
end


function tbl = local_empty_exclusions_table()
tbl = local_empty_repair_log_table();
end


function tbl = local_empty_points_like()
tbl = table('Size', [0 3], ...
    'VariableTypes', {'cell', 'double', 'cell'}, ...
    'VariableNames', {'repair_source', 'original_energy_meV', ...
    'repair_detail'});
end
