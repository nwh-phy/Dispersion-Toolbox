function model = dispersion_models(name)
%DISPERSION_MODELS  Registry of dispersion model definitions for E(q) fitting.
%
%   model = dispersion_models('quasi2d_plasmon')
%   model = dispersion_models('acoustic_linear')
%   model = dispersion_models('optical_constant')
%   model = dispersion_models('optical_quadratic')
%
%   Each model struct contains:
%       name        - display name
%       n_params    - number of fit parameters
%       param_names - cell array of parameter names + units
%       model_fn    - @(p, q) → E(q)  dispersion function
%       guess_fn    - @(q, E) → p0    initial guess from data
%       bounds_fn   - @(q, E) → struct with lb, ub
%       label_fn    - @(p_fit) → string   annotation label
%
%   Usage:
%       m = dispersion_models('acoustic_linear');
%       E = m.model_fn([v_s], q);
%
%   See also: fit_dispersion_generic, fit_quasi2d_plasmon, peak_models

arguments
    name (1,:) char {mustBeMember(name, ...
        {'quasi2d_plasmon', 'acoustic_linear', ...
         'optical_constant', 'optical_quadratic'})} = 'quasi2d_plasmon'
end

switch name
    case 'quasi2d_plasmon'
        model = quasi2d_plasmon_model();
    case 'acoustic_linear'
        model = acoustic_linear_model();
    case 'optical_constant'
        model = optical_constant_model();
    case 'optical_quadratic'
        model = optical_quadratic_model();
end

end


%% ═══════════════════════════════════════════════════════════════
%  Quasi-2D Plasmon: E = sqrt(A*|q| / (eps_bg + rho0*|q|))
%  da Jornada et al., Nat. Commun. 11, 1013 (2020)
%  p = [A, rho0], with eps_bg = 1 (adjustable via options)
%  ═══════════════════════════════════════════════════════════════
function m = quasi2d_plasmon_model()
    m.name = 'Quasi-2D Plasmon';
    m.n_params = 2;
    m.param_names = {'A (meV²·Å)', 'ρ₀ (Å)'};

    eps_bg = 1;  % default: suspended, (1+eps_s)/2 with eps_s=1

    m.model_fn = @(p, q) sqrt(abs(p(1)) .* abs(q) ./ (eps_bg + abs(p(2)) .* abs(q)));

    m.guess_fn = @(q, E) quasi2d_guess(q, E);

    m.bounds_fn = @(q, E) struct('lb', [0, 0.1], 'ub', [Inf, 500]);

    m.label_fn = @(p) sprintf('ρ₀=%.1f Å, E_{flat}=%.0f meV', ...
        abs(p(2)), sqrt(abs(p(1)) / max(abs(p(2)), 0.1)));
end

function p0 = quasi2d_guess(q, E)
    rho0 = 25;
    E_max = max(E);
    A = E_max^2 * rho0;
    p0 = [A, rho0];
end


%% ═══════════════════════════════════════════════════════════════
%  Acoustic Linear: E = v_s * |q|
%  p = [v_s]  (meV·Å, i.e. speed of sound in these units)
%  ═══════════════════════════════════════════════════════════════
function m = acoustic_linear_model()
    m.name = 'Acoustic Linear';
    m.n_params = 1;
    m.param_names = {'v_s (meV·Å)'};

    m.model_fn = @(p, q) abs(p(1)) .* abs(q);

    m.guess_fn = @(q, E) acoustic_guess(q, E);

    m.bounds_fn = @(q, E) struct('lb', 0, 'ub', Inf);

    m.label_fn = @(p) sprintf('v_s = %.1f meV·Å', abs(p(1)));
end

function p0 = acoustic_guess(q, E)
    % v_s ≈ E/q from median
    valid = abs(q) > 0 & E > 0;
    if any(valid)
        ratios = E(valid) ./ abs(q(valid));
        p0 = median(ratios);
    else
        p0 = 1000;
    end
end


%% ═══════════════════════════════════════════════════════════════
%  Optical Constant: E = ω₀
%  p = [omega_0]  (meV)
%  ═══════════════════════════════════════════════════════════════
function m = optical_constant_model()
    m.name = 'Optical Constant';
    m.n_params = 1;
    m.param_names = {'ω₀ (meV)'};

    m.model_fn = @(p, q) abs(p(1)) .* ones(size(q));

    m.guess_fn = @(q, E) median(E);

    m.bounds_fn = @(q, E) struct('lb', 0, 'ub', max(E) * 2);

    m.label_fn = @(p) sprintf('ω₀ = %.0f meV', abs(p(1)));
end


%% ═══════════════════════════════════════════════════════════════
%  Optical Quadratic: E = ω₀ − β·q²
%  p = [omega_0, beta]  (meV, meV·Å²)
%  ═══════════════════════════════════════════════════════════════
function m = optical_quadratic_model()
    m.name = 'Optical Quadratic';
    m.n_params = 2;
    m.param_names = {'ω₀ (meV)', 'β (meV·Å²)'};

    m.model_fn = @(p, q) abs(p(1)) - abs(p(2)) .* q.^2;

    m.guess_fn = @(q, E) optical_quad_guess(q, E);

    m.bounds_fn = @(q, E) struct('lb', [0, 0], 'ub', [max(E)*2, Inf]);

    m.label_fn = @(p) sprintf('ω₀=%.0f meV, β=%.1f meV·Å²', abs(p(1)), abs(p(2)));
end

function p0 = optical_quad_guess(q, E)
    omega0 = max(E);
    % Estimate beta from the curvature
    q_max = max(abs(q));
    E_min = min(E);
    if q_max > 0
        beta = (omega0 - E_min) / q_max^2;
    else
        beta = 1000;
    end
    p0 = [omega0, beta];
end
