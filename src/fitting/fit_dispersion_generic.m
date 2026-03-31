function result = fit_dispersion_generic(q_Ainv, energy_meV, options)
%FIT_DISPERSION_GENERIC  Universal dispersion E(q) fitter with pluggable models.
%
%   result = fit_dispersion_generic(q_Ainv, energy_meV)
%   result = fit_dispersion_generic(q_Ainv, energy_meV, 'model', 'acoustic_linear')
%
%   Fits any registered dispersion model from dispersion_models.m to
%   experimental (q, E) data. Computes group velocity numerically.
%
%   Name-value options:
%       model       - Model name (default: 'quasi2d_plasmon')
%       confidence  - Weights for fitting (same size as q_Ainv)
%       epsilon_s   - Substrate dielectric constant (default: 1)
%       n_fit_pts   - Number of points for smooth fit curve (default: 200)
%
%   Output struct fields:
%       params        - fitted parameter vector
%       param_names   - cell array of parameter names
%       model_name    - display name of the model
%       model_label   - annotation string from label_fn
%       q_fit, E_fit  - smooth fitted dispersion curve
%       q_data, E_data - input data
%       residuals_meV - fit residuals
%       R_squared     - coefficient of determination
%       RMSE_meV      - root mean square error
%       group_velocity_q, group_velocity - v_g(q) computed numerically
%       params_ci     - [n_params × 2] 95% CI from Jacobian (NaN if unavailable)
%
%   See also: dispersion_models, fit_quasi2d_plasmon

arguments
    q_Ainv     (:,1) double
    energy_meV (:,1) double
    options.model      (1,:) char = 'quasi2d_plasmon'
    options.confidence (:,1) double = ones(size(q_Ainv))
    options.epsilon_s  (1,1) double {mustBePositive} = 1
    options.n_fit_pts  (1,1) double {mustBePositive, mustBeInteger} = 200
end

%% Load dispersion model
dm = dispersion_models(options.model);

%% Prepare data
q_raw = double(q_Ainv(:));
q = abs(q_raw);
E = double(energy_meV(:));
w = double(options.confidence(:));
if numel(w) ~= numel(q)
    w = ones(size(q));
end

% Remove invalid points
valid = isfinite(q) & isfinite(E) & q > 0 & E > 0 & isfinite(w);
q_raw = q_raw(valid);
q = q(valid);
E = E(valid);
w = w(valid);

if numel(q) < 2
    error('fit_dispersion_generic:InsufficientData', ...
        'Need at least 2 valid data points for fitting.');
end

%% Initial guess and bounds
p0 = dm.guess_fn(q, E);
bnd = dm.bounds_fn(q, E);
lb = bnd.lb;
ub = bnd.ub;

% Ensure column vectors for consistency
p0 = p0(:)';

%% Fit
model_fn = dm.model_fn;
weighted_model = @(p, q_in) sqrt(w) .* model_fn(p, q_in);
weighted_data = sqrt(w) .* E;

has_jacobian = false;
J_fit = [];

try
    fit_opts = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 5000, ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12);
    [p_fit, ~, residual_vec, ~, ~, ~, J_fit] = ...
        lsqcurvefit(weighted_model, p0, q, weighted_data, lb, ub, fit_opts);
    has_jacobian = true;
catch
    cost = @(p) sum(w .* (model_fn(abs(p), q) - E).^2);
    fmin_opts = optimset('Display', 'off', 'MaxFunEvals', 5000, ...
        'MaxIter', 1000, 'TolFun', 1e-12, 'TolX', 1e-12);
    p_fit = abs(fminsearch(cost, p0, fmin_opts));
end

p_fit = abs(p_fit);  % ensure positive parameters

%% Confidence intervals from Jacobian
n_params = numel(p_fit);
params_ci = NaN(n_params, 2);
if has_jacobian && ~isempty(J_fit)
    try
        J_full = full(J_fit);
        N_data = numel(q);
        MSE = sum(residual_vec.^2) / max(N_data - n_params, 1);
        C = pinv(J_full' * J_full) * MSE;
        param_se = sqrt(abs(diag(C)));
        ci_half = 1.96 * param_se;
        params_ci(:, 1) = p_fit(:) - ci_half;
        params_ci(:, 2) = p_fit(:) + ci_half;
    catch
    end
end

%% Compute fit curve (handle signed q)
q_max_abs = max(abs(q_raw)) * 1.2;
q_has_neg = any(q_raw < 0);
q_has_pos = any(q_raw > 0);

if q_has_neg && q_has_pos
    q_fit_abs = linspace(0, q_max_abs, options.n_fit_pts)';
    q_fit = [-flipud(q_fit_abs(2:end)); q_fit_abs];
    E_fit_abs = model_fn(p_fit, q_fit_abs);
    E_fit = [flipud(E_fit_abs(2:end)); E_fit_abs];
elseif q_has_neg
    q_fit_abs = linspace(0, q_max_abs, options.n_fit_pts)';
    q_fit = -flipud(q_fit_abs);
    E_fit = flipud(model_fn(p_fit, q_fit_abs));
else
    q_fit = linspace(0, q_max_abs, options.n_fit_pts)';
    E_fit = model_fn(p_fit, q_fit);
end

%% Residuals and R²
E_pred = model_fn(p_fit, q);
residuals = E - E_pred;
SS_res = sum(w .* residuals.^2);
SS_tot = sum(w .* (E - mean(E)).^2);
R_squared = 1 - SS_res / max(SS_tot, eps);
RMSE = sqrt(mean(residuals.^2));

%% Numerical group velocity: v_g = dE/dq
dq = 1e-4;  % Å⁻¹
q_vg = linspace(max(0.001, min(abs(q_raw))), max(abs(q_raw)) * 1.1, options.n_fit_pts)';
E_plus = model_fn(p_fit, q_vg + dq);
E_minus = model_fn(p_fit, q_vg - dq);
vg = (E_plus - E_minus) / (2 * dq);  % meV·Å

%% Derived quantities (model-specific)
% For quasi2d_plasmon, also compute E_flat, rho0, q_c
extra = struct();
if strcmp(options.model, 'quasi2d_plasmon') && numel(p_fit) >= 2
    A_fit = p_fit(1);
    rho0_fit = p_fit(2);
    eps_bg = (1 + options.epsilon_s) / 2;
    extra.A = A_fit;
    extra.rho0 = rho0_fit;
    extra.E_flat_meV = sqrt(A_fit / rho0_fit);
    extra.q_c_Ainv = eps_bg / rho0_fit;
    extra.epsilon_s = options.epsilon_s;
end

%% Build result struct
result = struct();
result.params = p_fit;
result.param_names = dm.param_names;
result.params_ci = params_ci;
result.model_name = dm.name;
result.model_label = dm.label_fn(p_fit);
result.q_fit = q_fit;
result.E_fit = E_fit;
result.q_data = q_raw;
result.E_data = E;
result.weights = w;
result.residuals_meV = residuals;
result.R_squared = R_squared;
result.RMSE_meV = RMSE;
result.group_velocity_q = q_vg;
result.group_velocity = vg;

% Merge extra fields (model-specific derived quantities)
fn = fieldnames(extra);
for i = 1:numel(fn)
    result.(fn{i}) = extra.(fn{i});
end

fprintf('  Dispersion fit [%s]:\n', dm.name);
for i = 1:numel(p_fit)
    if ~isnan(params_ci(i,1))
        fprintf('    %s = %.4g  [%.4g, %.4g]\n', dm.param_names{i}, ...
            p_fit(i), params_ci(i,1), params_ci(i,2));
    else
        fprintf('    %s = %.4g\n', dm.param_names{i}, p_fit(i));
    end
end
fprintf('    R²   = %.4f\n', R_squared);
fprintf('    RMSE = %.1f meV\n', RMSE);

end
