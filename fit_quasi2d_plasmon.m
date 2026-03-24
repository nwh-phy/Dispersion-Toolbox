function result = fit_quasi2d_plasmon(q_Ainv, energy_meV, options)
%FIT_QUASI2D_PLASMON  Fit the quasi-2D plasmon dispersion model.
%
%   result = fit_quasi2d_plasmon(q_Ainv, energy_meV)
%   result = fit_quasi2d_plasmon(q_Ainv, energy_meV, confidence=w, ...)
%
%   Fits the da Jornada et al. (Nat. Commun. 2020) model for
%   dispersionless plasmons in atomically thin quasi-2D metals:
%
%       E(q) = sqrt( A * |q| / ( (1+eps_s)/2 + rho0*|q| ) )
%
%   where:
%       A    = Drude weight parameter  [meV^2 * Angstrom]
%       rho0 = interband screening length  [Angstrom]
%       eps_s = substrate dielectric constant (default 1, suspended)
%
%   Reference values for monolayer 2H-TaS2:
%       rho0 ~ 25-28 Angstrom (Da Jornada 2020, Do et al. 2025)
%       E_flat ~ 1 eV in the dispersionless region
%
%   Inputs:
%       q_Ainv      - wave vector (1/Angstrom), will use |q|
%       energy_meV  - plasmon energy (meV)
%
%   Name-value options:
%       confidence  - weights for fitting (same size as q_Ainv)
%       epsilon_s   - substrate dielectric constant (default: 1)
%       rho0_init   - initial guess for rho0 (default: 25 Angstrom)
%       n_fit_pts   - number of points for smooth fit curve (default: 200)
%
%   Output struct fields:
%       A, rho0       - fitted parameters
%       epsilon_s     - substrate dielectric constant used
%       E_flat_meV    - asymptotic plasmon energy sqrt(A/rho0)
%       q_c_Ainv      - critical wavevector (1+eps_s)/(2*rho0)
%       q_fit, E_fit  - smooth fitted dispersion curve
%       q_data, E_data - input data used for fitting
%       residuals_meV - fit residuals
%       R_squared     - coefficient of determination
%       RMSE_meV      - root mean square error
%       group_velocity_q, group_velocity - v_g(q) in Angstrom*meV/hbar
%       model_label   - descriptive string for plot legend
%
%   See also: track_plasmon_dispersion, track_plasmon_ridge

arguments
    q_Ainv (:,1) double
    energy_meV (:,1) double
    options.confidence (:,1) double = ones(size(q_Ainv))
    options.epsilon_s (1,1) double {mustBePositive} = 1
    options.rho0_init (1,1) double {mustBePositive} = 25
    options.n_fit_pts (1,1) double {mustBePositive, mustBeInteger} = 200
end

%% Validate and prepare data
assert(numel(q_Ainv) == numel(energy_meV), ...
    'fit_quasi2d_plasmon:SizeMismatch', ...
    'q_Ainv and energy_meV must have the same number of elements.');

q_raw = double(q_Ainv(:));           % keep signed q
q = abs(q_raw);                      % fit uses |q|
E = double(energy_meV(:));
w = double(options.confidence(:));
if numel(w) ~= numel(q)
    w = ones(size(q));
end

% Remove NaN / Inf / zero-q points
valid = isfinite(q) & isfinite(E) & q > 0 & E > 0 & isfinite(w);
q_raw = q_raw(valid);
q = q(valid);
E = E(valid);
w = w(valid);

if numel(q) < 2
    error('fit_quasi2d_plasmon:InsufficientData', ...
        'Need at least 2 valid data points for fitting.');
end

eps_s = double(options.epsilon_s);
eps_bg = (1 + eps_s) / 2;  % background dielectric from substrate

%% Define model: E(q) = sqrt( A * q / (eps_bg + rho0*q) )
%  Parameters: p = [A, rho0]
model_func = @(p, q_in) sqrt( abs(p(1)) .* q_in ./ (eps_bg + abs(p(2)) .* q_in) );

%% Initial guess
rho0_guess = options.rho0_init;
% Estimate A from the data: at large q, E ~ sqrt(A/rho0)
% At small q, E ~ sqrt(A*q/eps_bg)
E_max = max(E);
A_guess = E_max^2 * rho0_guess;  % from E_flat = sqrt(A/rho0)

p0 = [A_guess, rho0_guess];

%% Fit using lsqcurvefit if available, else fminsearch
lb = [0, 0.1];              % lower bounds
ub = [Inf, 500];             % upper bounds (rho0 up to 500 Å)

try
    % Weighted least squares objective
    weighted_model = @(p, q_in) sqrt(w) .* model_func(p, q_in);
    weighted_data = sqrt(w) .* E;

    fit_options = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 5000, ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12);

    p_fit = lsqcurvefit(weighted_model, p0, q, weighted_data, lb, ub, fit_options);
catch
    % Fallback to fminsearch (no bounds)
    cost = @(p) sum(w .* (model_func(abs(p), q) - E).^2);
    fmin_opts = optimset('Display', 'off', 'MaxFunEvals', 5000, 'MaxIter', 1000, ...
        'TolFun', 1e-12, 'TolX', 1e-12);
    p_fit = abs(fminsearch(cost, p0, fmin_opts));
end

A_fit = abs(p_fit(1));
rho0_fit = abs(p_fit(2));

%% Compute fit curve — span from min(q_raw) to max(q_raw) including negatives
q_max_abs = max(abs(q_raw)) * 1.2;
q_has_neg = any(q_raw < 0);
q_has_pos = any(q_raw > 0);

if q_has_neg && q_has_pos
    % Both sides: generate symmetric curve
    q_fit_abs = linspace(0, q_max_abs, options.n_fit_pts)';
    q_fit = [-flipud(q_fit_abs(2:end)); q_fit_abs];
    E_fit_abs = model_func([A_fit, rho0_fit], q_fit_abs);
    E_fit = [flipud(E_fit_abs(2:end)); E_fit_abs];
elseif q_has_neg
    q_fit_abs = linspace(0, q_max_abs, options.n_fit_pts)';
    q_fit = -flipud(q_fit_abs);
    E_fit = flipud(model_func([A_fit, rho0_fit], q_fit_abs));
else
    q_fit = linspace(0, q_max_abs, options.n_fit_pts)';
    E_fit = model_func([A_fit, rho0_fit], q_fit);
end

%% Compute residuals and R²
E_pred = model_func([A_fit, rho0_fit], q);
residuals = E - E_pred;
SS_res = sum(w .* residuals.^2);
SS_tot = sum(w .* (E - mean(E)).^2);
R_squared = 1 - SS_res / max(SS_tot, eps);
RMSE = sqrt(mean(residuals.^2));

%% Derived quantities
E_flat = sqrt(A_fit / rho0_fit);       % asymptotic energy (meV)
q_c = eps_bg / rho0_fit;               % critical wavevector (1/Å)

%% Group velocity: v_g = dE/dq
%  E = sqrt(A*q/(eps_bg + rho0*q))
%  dE/dq = A*eps_bg / (2 * (eps_bg + rho0*q)^2 * E)
q_vg = linspace(max(q_fit(2), 1e-4), max(q) * 1.2, options.n_fit_pts)';
E_vg = model_func([A_fit, rho0_fit], q_vg);
vg = A_fit * eps_bg ./ (2 * (eps_bg + rho0_fit * q_vg).^2 .* max(E_vg, eps));

%% Build result struct
result = struct();
result.A = A_fit;
result.rho0 = rho0_fit;
result.epsilon_s = eps_s;
result.E_flat_meV = E_flat;
result.q_c_Ainv = q_c;
result.q_fit = q_fit;
result.E_fit = E_fit;
result.q_data = q_raw;                % signed q data
result.E_data = E;
result.weights = w;
result.residuals_meV = residuals;
result.R_squared = R_squared;
result.RMSE_meV = RMSE;
result.group_velocity_q = q_vg;
result.group_velocity = vg;
result.model_label = sprintf( ...
    '\\omega_p = \\surd(A q / (%.1f + \\rho_0 q)),  \\rho_0=%.1f \\AA,  E_{flat}=%.0f meV', ...
    eps_bg, rho0_fit, E_flat);
result.model_equation = 'E(q) = sqrt(A*|q| / ((1+eps_s)/2 + rho0*|q|))';
result.reference = 'da Jornada et al., Nat. Commun. 11, 1013 (2020)';

fprintf('  Quasi-2D plasmon fit:\n');
fprintf('    A          = %.2f  meV^2·Å\n', A_fit);
fprintf('    rho0       = %.2f  Å\n', rho0_fit);
fprintf('    E_flat     = %.1f  meV\n', E_flat);
fprintf('    q_c        = %.4f  Å^{-1}\n', q_c);
fprintf('    R²         = %.4f\n', R_squared);
fprintf('    RMSE       = %.1f  meV\n', RMSE);
end
