function result = fit_loss_function(energy_meV, spectrum, options)
%FIT_LOSS_FUNCTION  Multi-peak spectral fitting with pluggable peak models.
%
%   Fits the spectrum to:
%     S(E) = B0 * E^(-alpha) + sum_i { peak_model(E0_i, width_i, A_i, E) }
%   or, if pre_subtracted=true (background already removed):
%     S(E) = C + sum_i { peak_model(E0_i, width_i, A_i, E) }
%
%   Name-value options:
%       E_min, E_max        - fitting window (meV)
%       max_peaks           - max peaks to detect (default: 3)
%       min_prominence      - peak prominence threshold (default: 0.15)
%       smooth_width        - smoothing window for peak detection (default: 25)
%       peak_model          - 'lorentz' (default), 'gaussian', 'voigt', 'damped_ho'
%       pre_subtracted      - true if background already subtracted (default: false)
%
%   See also: peak_models, propagate_seed_peaks, fit_quasi2d_plasmon

arguments
    energy_meV (:,1) double
    spectrum   (:,1) double
    options.E_min          (1,1) double = 50
    options.E_max          (1,1) double = Inf
    options.max_peaks      (1,1) double {mustBePositive, mustBeInteger} = 3
    options.min_prominence (1,1) double = 0.15
    options.smooth_width   (1,1) double {mustBePositive} = 25
    options.initial_guesses (:,1) double = []
    options.peak_model     (1,:) char = 'lorentz'
    options.pre_subtracted (1,1) logical = false
    options.bootstrap_ci_samples (1,1) double = NaN
end

%% Load peak model definition
pk_model = peak_models(options.peak_model);
n_peak_params = pk_model.n_params;
if isnan(options.bootstrap_ci_samples)
    if strcmpi(options.peak_model, 'fano')
        options.bootstrap_ci_samples = 25;
    else
        options.bootstrap_ci_samples = 0;
    end
elseif ~isfinite(options.bootstrap_ci_samples) || options.bootstrap_ci_samples < 0 || ...
        options.bootstrap_ci_samples ~= floor(options.bootstrap_ci_samples)
    error('fit_loss_function:InvalidBootstrapSamples', ...
        'bootstrap_ci_samples must be a nonnegative integer.');
end

%% Prepare data
E = double(energy_meV(:));
S = double(spectrum(:));

if isinf(options.E_max)
    options.E_max = max(E);
end
mask = E >= options.E_min & E <= options.E_max & isfinite(S);
E_w = E(mask);
S_w = S(mask);

if numel(E_w) < 10
    error('fit_loss_function:InsufficientData', ...
        'Need at least 10 data points in the fitting window.');
end

% Normalize
S_scale = max(abs(S_w));
if S_scale == 0; S_scale = 1; end
S_n = S_w / S_scale;

%% Estimate detection background/baseline from the data edges
n_pts = numel(E_w);
n_edge = max(5, round(n_pts * 0.12));
edge_idx = [1:n_edge, (n_pts-n_edge+1):n_pts]';

if options.pre_subtracted
    % Background was already removed upstream. Use only a constant baseline
    % for peak detection so initialization does not subtract a second model.
    C_init = local_robust_constant_baseline(S_n, edge_idx);
    S_detrended = max(S_n - C_init, 0);
else
    % Fit log(S) ~ -alpha*log(E) + log(B0) using the first and last portions.
    E_edge = E_w(edge_idx);
    S_edge = S_n(edge_idx);
    valid = S_edge > 0 & E_edge > 0;

    if sum(valid) >= 4
        P = polyfit(log(E_edge(valid)), log(S_edge(valid)), 1);
        alpha_init = -P(1);
        B0_init = exp(P(2));
    else
        alpha_init = 1.0;
        B0_init = S_n(1) * E_w(1);
    end
    alpha_init = max(0.1, min(alpha_init, 5));

    bg_est = B0_init * E_w .^ (-alpha_init);
    bg_est = min(bg_est, S_n * 0.95);  % don't over-subtract
    S_detrended = max(S_n - bg_est, 0);
end

smooth_spec = smoothdata(S_detrended, 'gaussian', options.smooth_width);
max_val = max(smooth_spec);

if ~isempty(options.initial_guesses)
    %% Manual initial guesses — skip findpeaks
    locs = options.initial_guesses(:);
    n_peaks = min(numel(locs), options.max_peaks);
    locs = locs(1:n_peaks);
    widths = ones(n_peaks, 1) * 300;
    pks = zeros(n_peaks, 1);
    for ig = 1:n_peaks
        [~, nearest] = min(abs(E_w - locs(ig)));
        pks(ig) = max(S_detrended(nearest), 0.01);
    end

elseif max_val > 0
    %% Auto peak detection
    min_prom = max_val * options.min_prominence;
    [pks, locs, widths] = findpeaks(smooth_spec, E_w, ...
        'MinPeakProminence', min_prom, ...
        'MinPeakDistance', 150, ...
        'SortStr', 'descend');
    n_peaks = min(numel(pks), options.max_peaks);

    if n_peaks == 0
        [~, max_idx] = max(smooth_spec);
        locs = E_w(max_idx);
        pks = smooth_spec(max_idx);
        widths = 300;
        n_peaks = 1;
    else
        pks = pks(1:n_peaks);
        locs = locs(1:n_peaks);
        widths = widths(1:n_peaks);
    end

else
    %% No signal — use global max
    [~, max_idx] = max(S_n);
    locs = E_w(max_idx);
    pks = S_n(max_idx);
    widths = 300;
    n_peaks = 1;
end

% Sort ascending by energy
[locs, si] = sort(locs);
pks = pks(si);
widths = widths(si);

%% Build model parameters
% pre_subtracted mode: [C, peak1 params..., peakN params...]
% simultaneous mode:   [B0, alpha, peak1 params..., peakN params...]
if options.pre_subtracted
    n_bg_params = 1;  % single constant offset
else
    n_bg_params = 2;  % B0, alpha
end
n_params = n_bg_params + n_peak_params * n_peaks;
p0 = zeros(n_params, 1);
lb = zeros(n_params, 1);
ub = zeros(n_params, 1);

if options.pre_subtracted
    p0(1) = C_init;    lb(1) = -Inf;   ub(1) = Inf;  % allow negative offset
else
    % Background: B0, alpha
    B0_floor = max(S_n(1), mean(S_n(1:min(5,end)))) * 0.01;
    p0(1) = B0_init;      lb(1) = B0_floor;  ub(1) = Inf;
    p0(2) = alpha_init;   lb(2) = 0.01;      ub(2) = 4;
end

% Peaks — use model-specific guesses and bounds
for i = 1:n_peaks
    base = n_bg_params + n_peak_params*(i-1);
    idx = base + (1:n_peak_params);

    % Get initial guess from model
    [~, nearest_idx] = min(abs(E_w - locs(i)));
    pk_detrended = max(S_detrended(nearest_idx), pks(i));
    p_guess = pk_model.guess_fn(locs(i), widths(i), pk_detrended);
    p0(idx) = p_guess(:);

    % Get bounds from model
    bnd = pk_model.bounds_fn(locs(i), widths(i), options.E_min, options.E_max);
    lb(idx) = bnd.lb(:);
    ub(idx) = bnd.ub(:);
end

model = @(p, E_in) local_composite_model(p, E_in, n_peaks, pk_model, n_bg_params, n_peak_params);

%% Fit
has_jacobian = false;
J_fit = [];
fit_opts = [];
try
    fit_opts = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 15000, ...
        'MaxIterations', 3000, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12);
    [p_fit, ~, residual_vec, ~, ~, ~, J_fit] = lsqcurvefit(model, p0, E_w, S_n, lb, ub, fit_opts);
    has_jacobian = true;
catch
    cost = @(p) sum((model(p, E_w) - S_n).^2);
    fmin_opts = optimset('Display', 'off', 'MaxFunEvals', 15000, ...
        'MaxIter', 3000, 'TolFun', 1e-12, 'TolX', 1e-12);
    p_fit = fminsearch(cost, p0, fmin_opts);
end

%% Compute confidence intervals from Jacobian
N_data = numel(E_w);
param_ci_all = NaN(n_params, 1);  % 95% half-width for each param
if has_jacobian && ~isempty(J_fit)
    try
        J_full = full(J_fit);
        MSE = sum(residual_vec.^2) / max(N_data - n_params, 1);
        % Covariance matrix: C = inv(J'*J) * MSE
        JtJ = J_full' * J_full;
        % Use pseudo-inverse for numerical stability
        C = pinv(JtJ) * MSE;
        param_se = sqrt(abs(diag(C)));
        % 95% CI half-width: t-value ≈ 1.96 for large N
        param_ci_all = 1.96 * param_se;
    catch
        % If CI computation fails, leave as NaN
    end
end

%% Extract results
if options.pre_subtracted
    B0_fit = 0;
    alpha_fit = 0;
else
    B0_fit = abs(p_fit(1));
    alpha_fit = abs(p_fit(2));
end

omega_p_vals = zeros(n_peaks, 1);
gamma_vals = zeros(n_peaks, 1);
amp_vals = zeros(n_peaks, 1);
peak_param_vals = zeros(n_peaks, n_peak_params);
peak_param_ci_lo = NaN(n_peaks, n_peak_params);
peak_param_ci_hi = NaN(n_peaks, n_peak_params);
omega_p_ci = NaN(n_peaks, 2);  % [lower, upper]
gamma_ci = NaN(n_peaks, 2);
amplitude_ci = NaN(n_peaks, 2);

for i = 1:n_peaks
    base = n_bg_params + n_peak_params*(i-1);
    idx = base + (1:n_peak_params);
    params = p_fit(idx);
    params(1:min(3, n_peak_params)) = abs(params(1:min(3, n_peak_params)));
    peak_param_vals(i, :) = params(:).';
    omega_p_vals(i) = peak_param_vals(i, 1);
    gamma_vals(i) = peak_param_vals(i, 2);
    amp_vals(i) = peak_param_vals(i, 3);

    % Extract CI for each peak parameter
    for j = 1:n_peak_params
        ci = param_ci_all(idx(j));
        peak_param_ci_lo(i, j) = peak_param_vals(i, j) - ci;
        peak_param_ci_hi(i, j) = peak_param_vals(i, j) + ci;
    end
    omega_p_ci(i, :) = [peak_param_ci_lo(i, 1), peak_param_ci_hi(i, 1)];
    gamma_ci(i, :)   = [peak_param_ci_lo(i, 2), peak_param_ci_hi(i, 2)];
    amplitude_ci(i,:) = [peak_param_ci_lo(i, 3), peak_param_ci_hi(i, 3)];
end

[omega_p_vals, si] = sort(omega_p_vals);
gamma_vals = gamma_vals(si);
amp_vals = amp_vals(si);
peak_param_vals = peak_param_vals(si, :);
peak_param_ci_lo = peak_param_ci_lo(si, :);
peak_param_ci_hi = peak_param_ci_hi(si, :);
omega_p_ci = omega_p_ci(si, :);
gamma_ci = gamma_ci(si, :);
amplitude_ci = amplitude_ci(si, :);

% Post-fit filter: remove peaks with amplitude < 10% of max
max_amp = max(amp_vals);
keep = amp_vals >= 0.1 * max_amp;
if any(keep)
    omega_p_vals = omega_p_vals(keep);
    gamma_vals = gamma_vals(keep);
    amp_vals = amp_vals(keep);
    peak_param_vals = peak_param_vals(keep, :);
    peak_param_ci_lo = peak_param_ci_lo(keep, :);
    peak_param_ci_hi = peak_param_ci_hi(keep, :);
    omega_p_ci = omega_p_ci(keep, :);
    gamma_ci = gamma_ci(keep, :);
    amplitude_ci = amplitude_ci(keep, :);
    n_peaks = sum(keep);
end

%% Fit quality
S_pred = model(p_fit, E_w);
SS_res = sum((S_n - S_pred).^2);
SS_tot = sum((S_n - mean(S_n)).^2);
R_sq = 1 - SS_res / max(SS_tot, eps);

%% Generate smooth curves
E_fine = linspace(min(E_w), max(E_w), 500)';

% Background curve
if options.pre_subtracted
    C_fit = p_fit(1);
    bg_curve = ones(size(E_fine)) * C_fit * S_scale;
else
    bg_curve = B0_fit * E_fine .^ (-alpha_fit) * S_scale;
end

% Total fit and individual peaks
curve_total = bg_curve;
peak_curves = cell(n_peaks, 1);
peak_heights = zeros(n_peaks, 1);
apex_energy_vals = NaN(n_peaks, 1);
apex_offset_vals = NaN(n_peaks, 1);

for i = 1:n_peaks
    single = local_eval_peak(pk_model, peak_param_vals(i, :), E_fine);
    single = single * S_scale;
    peak_curves{i} = single;
    [peak_heights(i), apex_idx] = max(single);
    apex_energy_vals(i) = E_fine(apex_idx);
    apex_offset_vals(i) = apex_energy_vals(i) - omega_p_vals(i);
    curve_total = curve_total + single;
end

gamma_ratio_vals = gamma_vals ./ max(omega_p_vals, eps);
[peak_valid_vals, peak_quality_notes] = local_peak_quality( ...
    omega_p_vals, gamma_ratio_vals, apex_offset_vals);
[apex_energy_ci, apex_ci_method] = local_apex_energy_ci( ...
    options.bootstrap_ci_samples, has_jacobian, fit_opts, model, p_fit, lb, ub, ...
    E_w, S_n, n_bg_params, n_peak_params, pk_model, peak_param_vals, ...
    omega_p_vals, omega_p_ci, apex_energy_vals, E_fine);

%% Build result
result = struct();
result.n_peaks = n_peaks;
result.omega_p = omega_p_vals;
result.gamma = gamma_vals;
result.amplitude = amp_vals * S_scale;
result.omega_p_ci = omega_p_ci;           % [n_peaks × 2] 95% CI
result.gamma_ci = gamma_ci;               % [n_peaks × 2] 95% CI
result.amplitude_ci = amplitude_ci * S_scale;  % [n_peaks × 2] 95% CI
if options.pre_subtracted
    result.B0 = 0;
    result.alpha = 0;
    result.offset = p_fit(1) * S_scale;
else
    result.B0 = B0_fit * S_scale;
    result.alpha = alpha_fit;
    result.offset = 0;
end
result.pre_subtracted = options.pre_subtracted;
result.R_squared = R_sq;
result.energy_fit = E_fine;
result.curve_fit = curve_total;
result.bg_curve = bg_curve;
result.peak_curves = peak_curves;
result.peak_heights = peak_heights;
result.apex_energy_meV = apex_energy_vals;
result.apex_energy_ci = apex_energy_ci;
result.apex_ci_method = apex_ci_method;
result.apex_offset_meV = apex_offset_vals;
result.gamma_ratio = gamma_ratio_vals;
result.peak_valid = peak_valid_vals;
result.peak_quality_notes = peak_quality_notes;
peak_param_result_vals = peak_param_vals;
peak_param_result_ci_lo = peak_param_ci_lo;
peak_param_result_ci_hi = peak_param_ci_hi;
if n_peak_params >= 3
    peak_param_result_vals(:, 3) = peak_param_result_vals(:, 3) * S_scale;
    peak_param_result_ci_lo(:, 3) = peak_param_result_ci_lo(:, 3) * S_scale;
    peak_param_result_ci_hi(:, 3) = peak_param_result_ci_hi(:, 3) * S_scale;
end
result.peak_param_names = pk_model.param_names;
result.peak_param_values = peak_param_result_vals;
result.peak_param_ci_lower = peak_param_result_ci_lo;
result.peak_param_ci_upper = peak_param_result_ci_hi;
if strcmpi(options.peak_model, 'fano') && n_peak_params >= 4
    result.fano_q = peak_param_vals(:, 4);
    result.fano_q_ci = [peak_param_ci_lo(:, 4), peak_param_ci_hi(:, 4)];
end
result.energy_data = E_w;
result.spectrum_data = S_w;
result.residuals = (S_n - S_pred) * S_scale;
result.peak_model_name = pk_model.name;

end


function S = local_composite_model(p, E, n_peaks, pk_model, n_bg_params, n_peak_params)
    % Composite model: background + sum of peaks
    if n_bg_params == 1
        % Pre-subtracted mode: S(E) = C + sum { peaks }
        S = ones(size(E)) * p(1);
    else
        % Simultaneous mode: S(E) = B0 * E^(-alpha) + sum { peaks }
        B0 = abs(p(1));
        alpha = abs(p(2));
        S = B0 * E .^ (-alpha);
    end

    for i = 1:n_peaks
        base = n_bg_params + n_peak_params*(i-1);
        params = p(base + (1:n_peak_params));
        params(1:min(3, n_peak_params)) = abs(params(1:min(3, n_peak_params)));
        S = S + local_eval_peak(pk_model, params, E);
    end
end


function S = local_eval_peak(pk_model, params, E)
    args = num2cell(params(:).');
    S = pk_model.model_fn(args{:}, E);
end


function C = local_robust_constant_baseline(S, edge_idx)
    S = S(:);
    finite = isfinite(S);
    S_finite = S(finite);
    if isempty(S_finite)
        C = 0;
        return
    end

    edge_idx = edge_idx(edge_idx >= 1 & edge_idx <= numel(S));
    edge_vals = S(edge_idx);
    edge_vals = edge_vals(isfinite(edge_vals));

    sorted_S = sort(S_finite);
    n_tail = max(1, round(numel(sorted_S) * 0.20));
    lower_tail = sorted_S(1:n_tail);

    baseline_vals = [edge_vals(:); lower_tail(:)];
    baseline_vals = baseline_vals(isfinite(baseline_vals));
    if isempty(baseline_vals)
        C = 0;
    else
        C = median(baseline_vals);
    end
end


function [apex_ci, method] = local_apex_energy_ci( ...
    n_samples, has_jacobian, fit_opts, model, p_fit, lb, ub, E_w, S_n, ...
    n_bg_params, n_peak_params, pk_model, peak_param_vals, omega_p_vals, ...
    omega_p_ci, apex_energy_vals, E_fine)

    n_targets = numel(apex_energy_vals);
    apex_ci = NaN(n_targets, 2);
    boot_valid = false(n_targets, 1);
    method = 'fallback';

    if n_samples > 0 && has_jacobian && ~isempty(fit_opts) && n_targets > 0
        [boot_ci, boot_valid] = local_bootstrap_apex_ci( ...
            n_samples, fit_opts, model, p_fit, lb, ub, E_w, S_n, ...
            n_bg_params, n_peak_params, pk_model, peak_param_vals, E_fine);
        apex_ci(boot_valid, :) = boot_ci(boot_valid, :);
        if any(boot_valid)
            method = 'bootstrap';
        end
    end

    apex_ci = local_fill_apex_ci(apex_ci, apex_energy_vals, omega_p_vals, omega_p_ci, E_fine);
    if any(boot_valid) && any(~boot_valid)
        method = 'bootstrap+fallback';
    end
end


function [boot_ci, valid] = local_bootstrap_apex_ci( ...
    n_samples, fit_opts, model, p_fit, lb, ub, E_w, S_n, ...
    n_bg_params, n_peak_params, pk_model, target_peak_params, E_fine)

    n_targets = size(target_peak_params, 1);
    boot_ci = NaN(n_targets, 2);
    valid = false(n_targets, 1);
    if n_targets == 0
        return
    end

    base_pred = model(p_fit, E_w);
    residual_pool = S_n - base_pred;
    residual_pool = residual_pool(isfinite(residual_pool));
    if isempty(residual_pool)
        return
    end

    try
        boot_opts = optimoptions(fit_opts, ...
            'Display', 'off', ...
            'MaxIterations', 200, ...
            'MaxFunctionEvaluations', 2000);
    catch
        boot_opts = fit_opts;
    end

    stream = RandStream('mt19937ar', 'Seed', 5489);
    boot_apex = NaN(n_samples, n_targets);
    n_boot_peaks = max(0, floor((numel(p_fit) - n_bg_params) / n_peak_params));

    for s = 1:n_samples
        sample_idx = randi(stream, numel(residual_pool), numel(E_w), 1);
        boot_y = base_pred + residual_pool(sample_idx);
        try
            p_boot = lsqcurvefit(model, p_fit, E_w, boot_y, lb, ub, boot_opts);
        catch
            continue
        end

        boot_params = local_extract_peak_params(p_boot, n_bg_params, n_peak_params, n_boot_peaks);
        boot_apex(s, :) = local_match_bootstrap_apex( ...
            boot_params, target_peak_params, pk_model, E_fine);
    end

    min_valid = max(5, ceil(0.25 * n_samples));
    for i = 1:n_targets
        vals = boot_apex(:, i);
        vals = vals(isfinite(vals));
        if numel(vals) < min_valid
            continue
        end

        boot_ci(i, :) = [local_percentile(vals, 2.5), local_percentile(vals, 97.5)];
        valid(i) = all(isfinite(boot_ci(i, :)));
    end
end


function peak_params = local_extract_peak_params(p, n_bg_params, n_peak_params, n_peaks)
    peak_params = NaN(n_peaks, n_peak_params);
    for i = 1:n_peaks
        base = n_bg_params + n_peak_params * (i - 1);
        params = p(base + (1:n_peak_params));
        params(1:min(3, n_peak_params)) = abs(params(1:min(3, n_peak_params)));
        peak_params(i, :) = params(:).';
    end
end


function matched_apex = local_match_bootstrap_apex(boot_params, target_peak_params, pk_model, E_fine)
    n_targets = size(target_peak_params, 1);
    matched_apex = NaN(1, n_targets);
    if isempty(boot_params)
        return
    end

    available = true(size(boot_params, 1), 1);
    for i = 1:n_targets
        distances = abs(boot_params(:, 1) - target_peak_params(i, 1));
        distances(~available) = Inf;
        [best_dist, best_idx] = min(distances);
        if ~isfinite(best_dist)
            continue
        end

        matched_apex(i) = local_param_apex_energy(pk_model, boot_params(best_idx, :), E_fine);
        available(best_idx) = false;
    end
end


function apex = local_param_apex_energy(pk_model, params, E_fine)
    apex = NaN;
    try
        single = local_eval_peak(pk_model, params, E_fine);
        [~, apex_idx] = max(single);
        apex = E_fine(apex_idx);
    catch
    end

    if ~isfinite(apex) && ~isempty(params)
        apex = params(1);
    end
end


function apex_ci = local_fill_apex_ci(apex_ci, apex_energy_vals, omega_p_vals, omega_p_ci, E_fine)
    min_half_width = local_min_energy_half_width(E_fine);
    for i = 1:numel(apex_energy_vals)
        center = apex_energy_vals(i);
        if ~isfinite(center) && i <= numel(omega_p_vals)
            center = omega_p_vals(i);
        end
        if ~isfinite(center)
            continue
        end

        lo = apex_ci(i, 1);
        hi = apex_ci(i, 2);
        if ~isfinite(lo) || ~isfinite(hi) || lo >= center || hi <= center
            half_width = local_fallback_half_width(i, omega_p_vals, omega_p_ci, min_half_width);
            apex_ci(i, :) = [center - half_width, center + half_width];
            continue
        end

        lo = min(lo, center - min_half_width);
        hi = max(hi, center + min_half_width);
        apex_ci(i, :) = [lo, hi];
    end
end


function half_width = local_fallback_half_width(i, omega_p_vals, omega_p_ci, min_half_width)
    half_width = NaN;
    if i <= size(omega_p_ci, 1) && i <= numel(omega_p_vals)
        bounds = omega_p_ci(i, :);
        if all(isfinite(bounds)) && isfinite(omega_p_vals(i))
            half_width = max(abs(bounds - omega_p_vals(i)));
        end
    end
    if ~isfinite(half_width) || half_width <= 0
        half_width = min_half_width;
    end
    half_width = max(half_width, min_half_width);
end


function half_width = local_min_energy_half_width(E_fine)
    E_fine = E_fine(isfinite(E_fine));
    if numel(E_fine) >= 2
        step = median(diff(sort(E_fine)));
    else
        step = NaN;
    end
    if ~isfinite(step) || step <= 0
        step = 2;
    end
    half_width = max(step / 2, 1);
end


function value = local_percentile(vals, pct)
    vals = sort(vals(:));
    n = numel(vals);
    if n == 0
        value = NaN;
        return
    end
    if n == 1
        value = vals(1);
        return
    end

    pos = 1 + (n - 1) * pct / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        value = vals(lo);
    else
        value = vals(lo) + (pos - lo) * (vals(hi) - vals(lo));
    end
end


function [is_valid, notes] = local_peak_quality(omega_p, gamma_ratio, apex_offset)
    omega_p = omega_p(:);
    gamma_ratio = gamma_ratio(:);
    apex_offset = apex_offset(:);

    n = numel(omega_p);
    is_valid = true(n, 1);
    notes = repmat({''}, n, 1);
    max_offset = max(150 * ones(n, 1), 0.25 * max(omega_p, eps));

    for i = 1:n
        reasons = {};
        if ~isfinite(gamma_ratio(i)) || gamma_ratio(i) > 2.0
            is_valid(i) = false;
            reasons{end+1} = sprintf('Gamma/E0=%.2g', gamma_ratio(i)); %#ok<AGROW>
        end
        if ~isfinite(apex_offset(i)) || abs(apex_offset(i)) > max_offset(i)
            is_valid(i) = false;
            reasons{end+1} = sprintf('apex offset=%.0f meV', apex_offset(i)); %#ok<AGROW>
        end

        if isempty(reasons)
            notes{i} = 'ok';
        else
            notes{i} = strjoin(reasons, '; ');
        end
    end
end
