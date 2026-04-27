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
end

%% Load peak model definition
pk_model = peak_models(options.peak_model);

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
% pre_subtracted mode: [C, Ep1, G1, A1, ..., EpN, GN, AN]  (1 + 3*N params)
% simultaneous mode:   [B0, alpha, Ep1, G1, A1, ..., EpN, GN, AN]  (2 + 3*N params)
if options.pre_subtracted
    n_bg_params = 1;  % single constant offset
else
    n_bg_params = 2;  % B0, alpha
end
n_params = n_bg_params + 3 * n_peaks;
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
    base = n_bg_params + 3*(i-1);

    % Get initial guess from model
    [~, nearest_idx] = min(abs(E_w - locs(i)));
    pk_detrended = max(S_detrended(nearest_idx), pks(i));
    p_guess = pk_model.guess_fn(locs(i), widths(i), pk_detrended);
    p0(base+1) = p_guess(1);
    p0(base+2) = p_guess(2);
    p0(base+3) = p_guess(3);

    % Get bounds from model
    bnd = pk_model.bounds_fn(locs(i), widths(i), options.E_min, options.E_max);
    lb(base+1) = bnd.lb(1);
    lb(base+2) = bnd.lb(2);
    lb(base+3) = bnd.lb(3);
    ub(base+1) = bnd.ub(1);
    ub(base+2) = bnd.ub(2);
    ub(base+3) = bnd.ub(3);
end

model = @(p, E_in) local_composite_model(p, E_in, n_peaks, pk_model.model_fn, n_bg_params);

%% Fit
has_jacobian = false;
J_fit = [];
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
omega_p_ci = NaN(n_peaks, 2);  % [lower, upper]
gamma_ci = NaN(n_peaks, 2);
amplitude_ci = NaN(n_peaks, 2);

for i = 1:n_peaks
    base = n_bg_params + 3*(i-1);
    omega_p_vals(i) = abs(p_fit(base+1));
    gamma_vals(i) = abs(p_fit(base+2));
    amp_vals(i) = abs(p_fit(base+3));

    % Extract CI for each peak parameter
    ci_E0 = param_ci_all(base+1);
    ci_Gm = param_ci_all(base+2);
    ci_Am = param_ci_all(base+3);
    omega_p_ci(i, :) = [omega_p_vals(i) - ci_E0, omega_p_vals(i) + ci_E0];
    gamma_ci(i, :)   = [gamma_vals(i) - ci_Gm, gamma_vals(i) + ci_Gm];
    amplitude_ci(i,:) = [amp_vals(i) - ci_Am, amp_vals(i) + ci_Am];
end

[omega_p_vals, si] = sort(omega_p_vals);
gamma_vals = gamma_vals(si);
amp_vals = amp_vals(si);
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

for i = 1:n_peaks
    Ep = omega_p_vals(i);
    Gm = gamma_vals(i);
    Am = amp_vals(i);
    single = pk_model.model_fn(Ep, Gm, Am, E_fine);
    single = single * S_scale;
    peak_curves{i} = single;
    peak_heights(i) = max(single);
    curve_total = curve_total + single;
end

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
result.energy_data = E_w;
result.spectrum_data = S_w;
result.residuals = (S_n - S_pred) * S_scale;
result.peak_model_name = pk_model.name;

end


function S = local_composite_model(p, E, n_peaks, peak_fn, n_bg_params)
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
        base = n_bg_params + 3*(i-1);
        E0 = abs(p(base+1));
        width = abs(p(base+2));
        A = abs(p(base+3));
        S = S + peak_fn(E0, width, A, E);
    end
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
