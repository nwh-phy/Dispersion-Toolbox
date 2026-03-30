function model = peak_models(name)
%PEAK_MODELS  Registry of peak model definitions for spectral fitting.
%
%   model = peak_models('lorentz')
%   model = peak_models('gaussian')
%   model = peak_models('voigt')
%   model = peak_models('damped_ho')
%
%   Each model struct contains:
%       name        - display name
%       n_params    - parameters per peak (always 3: E0, width, amplitude)
%       param_names - {'E0', 'Gamma/Sigma', 'A'}
%       model_fn    - @(E0, width, A, E) → S(E)  single peak contribution
%       guess_fn    - @(E0_est, width_est, amp_est) → [p1, p2, p3]
%       bounds_fn   - @(E0_est, E_min, E_max) → struct with lb, ub
%
%   Usage in fit_loss_function:
%       result = fit_loss_function(E, S, 'peak_model', 'voigt');
%
%   See also: fit_loss_function, propagate_seed_peaks

arguments
    name (1,:) char {mustBeMember(name, ...
        {'lorentz', 'gaussian', 'voigt', 'damped_ho'})} = 'lorentz'
end

switch name
    case 'lorentz'
        model = lorentz_model();
    case 'gaussian'
        model = gaussian_model();
    case 'voigt'
        model = voigt_model();
    case 'damped_ho'
        model = damped_ho_model();
end

end


%% ═══════════════════════════════════════════════════════════════
%  Drude-Lorentz: A * E * Gamma / ((E^2 - E0^2)^2 + E^2 * Gamma^2)
%  Standard model for plasmons and metallic excitations.
%  ═══════════════════════════════════════════════════════════════
function m = lorentz_model()
    m.name = 'Drude-Lorentz';
    m.n_params = 3;
    m.param_names = {'E0 (meV)', 'Gamma (meV)', 'A'};

    m.model_fn = @(E0, Gamma, A, E) ...
        A .* E .* Gamma ./ ((E.^2 - E0^2).^2 + E.^2 .* Gamma^2);

    m.guess_fn = @(E0_est, width_est, amp_est) ...
        [E0_est, max(width_est * 0.5, 50), amp_est * E0_est * max(width_est * 0.5, 50)];

    m.bounds_fn = @lorentz_bounds;
end

function b = lorentz_bounds(E0_est, width_est, E_min, E_max)
    b.lb = [max(100, E0_est - max(width_est*2, 400)), 10, 0];
    b.ub = [min(E_max, E0_est + max(width_est*2, 400)), 5000, Inf];
end


%% ═══════════════════════════════════════════════════════════════
%  Gaussian: A * exp(-(E - E0)^2 / (2*sigma^2))
%  Useful for narrow, resolution-limited peaks (phonons).
%  ═══════════════════════════════════════════════════════════════
function m = gaussian_model()
    m.name = 'Gaussian';
    m.n_params = 3;
    m.param_names = {'E0 (meV)', 'Sigma (meV)', 'A'};

    m.model_fn = @(E0, sigma, A, E) ...
        A .* exp(-(E - E0).^2 ./ (2 * sigma^2));

    m.guess_fn = @(E0_est, width_est, amp_est) ...
        [E0_est, max(width_est * 0.4, 5), amp_est];

    m.bounds_fn = @gaussian_bounds;
end

function b = gaussian_bounds(E0_est, width_est, E_min, E_max)
    b.lb = [max(E_min, E0_est - max(width_est*2, 200)), 1, 0];
    b.ub = [min(E_max, E0_est + max(width_est*2, 200)), 2000, Inf];
end


%% ═══════════════════════════════════════════════════════════════
%  Pseudo-Voigt: η * Lorentz + (1-η) * Gaussian
%  Approximates the Voigt profile without numerical convolution.
%  η is computed from Gamma and sigma via Thompson formula.
%  3 params per peak: E0, FWHM_total, A (η derived internally).
%  ═══════════════════════════════════════════════════════════════
function m = voigt_model()
    m.name = 'Pseudo-Voigt';
    m.n_params = 3;
    m.param_names = {'E0 (meV)', 'FWHM (meV)', 'A'};

    m.model_fn = @pseudo_voigt_peak;

    m.guess_fn = @(E0_est, width_est, amp_est) ...
        [E0_est, max(width_est * 0.8, 10), amp_est];

    m.bounds_fn = @voigt_bounds;
end

function S = pseudo_voigt_peak(E0, fwhm, A, E)
    % Pseudo-Voigt: fixed η=0.5 mix of Lorentz and Gaussian
    % (Could be extended to fit η, but 3-param keeps API uniform)
    eta = 0.5;
    sigma = fwhm / (2 * sqrt(2 * log(2)));  % Gaussian sigma from FWHM
    gamma_half = fwhm / 2;                   % Lorentz half-width

    L = (1/pi) * gamma_half ./ ((E - E0).^2 + gamma_half^2);
    G = (1 / (sigma * sqrt(2*pi))) .* exp(-(E - E0).^2 ./ (2 * sigma^2));

    S = A .* (eta * L + (1 - eta) * G);
end

function b = voigt_bounds(E0_est, width_est, E_min, E_max)
    b.lb = [max(E_min, E0_est - max(width_est*2, 200)), 2, 0];
    b.ub = [min(E_max, E0_est + max(width_est*2, 200)), 3000, Inf];
end


%% ═══════════════════════════════════════════════════════════════
%  Damped Harmonic Oscillator (Bose-weighted):
%    S(E) = A * n(E,T) * Gamma * E0 / ((E^2 - E0^2)^2 + E^2*Gamma^2)
%  where n(E,T) = 1/(exp(E/kT) - 1) is the Bose-Einstein factor.
%  Appropriate for phonon loss spectra at finite temperature.
%  Uses T=300K by default (encoded in the model).
%  ═══════════════════════════════════════════════════════════════
function m = damped_ho_model()
    m.name = 'Damped HO (Bose)';
    m.n_params = 3;
    m.param_names = {'E0 (meV)', 'Gamma (meV)', 'A'};

    m.model_fn = @dho_peak;

    m.guess_fn = @(E0_est, width_est, amp_est) ...
        [E0_est, max(width_est * 0.5, 5), amp_est * E0_est * max(width_est * 0.5, 5)];

    m.bounds_fn = @dho_bounds;
end

function S = dho_peak(E0, Gamma, A, E)
    kT = 25.85;  % meV at 300 K (Boltzmann constant * T)

    % Bose-Einstein factor (loss side: E > 0)
    bose = 1 ./ (exp(E ./ kT) - 1);
    bose(E <= 0) = 0;
    bose(~isfinite(bose)) = 0;

    % Damped harmonic oscillator spectral function
    chi_imag = A .* Gamma .* E ./ ((E.^2 - E0^2).^2 + E.^2 .* Gamma^2);

    S = (bose + 1) .* chi_imag;  % loss: proportional to (n+1)*χ''
end

function b = dho_bounds(E0_est, width_est, E_min, E_max)
    b.lb = [max(5, E0_est - max(width_est*2, 100)), 1, 0];
    b.ub = [min(E_max, E0_est + max(width_est*2, 100)), 1000, Inf];
end
