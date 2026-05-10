function metrics = b1_double_peak_upper_stability_metrics(points, session_key, options)
%B1_DOUBLE_PEAK_UPPER_STABILITY_METRICS Score upper-branch continuity.
%
% This diagnostic is intentionally separate from the physical fit score. It
% catches medium-scale upper-branch wobble and overbroad Lorentz components
% that can pass the older >250 meV jump metric.

arguments
    points table
    session_key {mustBeTextScalar} = ""
    options.mediumJumpThresholdMeV (1,1) double = NaN
    options.largeJumpThresholdMeV (1,1) double = 180
    options.maxGammaOverE (1,1) double = 1.4
    options.maxGammaMeV (1,1) double = 1600
    options.lowR2Threshold (1,1) double = 0.50
end

threshold = options.mediumJumpThresholdMeV;
if ~isfinite(threshold)
    threshold = local_default_medium_jump_threshold(session_key);
end

n_points = height(points);
medium_jump_count = 0;
large_jump_count = 0;
max_abs_jump = NaN;
overbroad_count = 0;
low_r2_count = 0;
median_gamma_over_E = NaN;
median_R2 = NaN;

if n_points >= 2 && all(ismember({'q_Ainv', 'energy_meV'}, ...
        points.Properties.VariableNames))
    [~, order] = sort(double(points.q_Ainv));
    energy = double(points.energy_meV(order));
    delta = abs(diff(energy));
    medium_jump_count = sum(delta > threshold, 'omitnan');
    large_jump_count = sum(delta > options.largeJumpThresholdMeV, 'omitnan');
    max_abs_jump = max(delta, [], 'omitnan');
end

if n_points >= 1 && all(ismember({'energy_meV', 'gamma_meV'}, ...
        points.Properties.VariableNames))
    energy = double(points.energy_meV(:));
    gamma = double(points.gamma_meV(:));
    gamma_over_E = gamma ./ max(abs(energy), eps);
    bad_ratio = isfinite(gamma_over_E) & ...
        gamma_over_E > options.maxGammaOverE;
    bad_abs = isfinite(gamma) & gamma > options.maxGammaMeV;
    overbroad_count = sum(bad_ratio | bad_abs);
    median_gamma_over_E = median(gamma_over_E(isfinite(gamma_over_E)), ...
        'omitnan');
end

if n_points >= 1 && ismember('R2', points.Properties.VariableNames)
    r2 = double(points.R2(:));
    low_r2_count = sum(isfinite(r2) & r2 < options.lowR2Threshold);
    median_R2 = median(r2(isfinite(r2)), 'omitnan');
end

stability_score = 6 * medium_jump_count + 20 * large_jump_count + ...
    5 * overbroad_count + 2 * low_r2_count;

metrics = table({char(string(session_key))}, n_points, threshold, ...
    options.largeJumpThresholdMeV, medium_jump_count, large_jump_count, ...
    max_abs_jump, overbroad_count, options.maxGammaOverE, ...
    options.maxGammaMeV, low_r2_count, options.lowR2Threshold, ...
    median_gamma_over_E, median_R2, stability_score, ...
    'VariableNames', {'session_key', 'n_upper_points', ...
    'upper_medium_jump_threshold_meV', ...
    'upper_large_jump_threshold_meV', 'upper_medium_jump_count', ...
    'upper_large_jump_count', 'upper_max_abs_jump_meV', ...
    'upper_overbroad_count', 'upper_max_gamma_over_E', ...
    'upper_max_gamma_meV', 'upper_low_r2_count', ...
    'upper_low_r2_threshold', 'upper_median_gamma_over_E', ...
    'upper_median_R2', 'upper_stability_score'});
end


function threshold = local_default_medium_jump_threshold(session_key)
text = lower(char(string(session_key)));
if contains(text, '20w')
    threshold = 120;
else
    threshold = 150;
end
end
