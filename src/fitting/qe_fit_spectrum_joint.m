function result = qe_fit_spectrum_joint(E, Y, q_Ainv, options)
%QE_FIT_SPECTRUM_JOINT  Joint ZLP fit of one spectrum, in the result format of fit_loss_function.
%
%   result = qe_fit_spectrum_joint(E, Y, q_Ainv, n_peaks=2, signal_window=[300 1800])
%
%   E, Y are the unprocessed spectrum of one q channel: counts, ZLP included,
%   no background subtracted, no normalisation or denoising. The fit is
%   qe_zlp_joint_fit: ZLP + loss peaks convolved with the fitted ZLP, data
%   from the gain side up to signal_window(2). Unless kinematic=false the loss
%   peaks carry the 2D kinematic prefactor of the channel (q_Ainv, width dq,
%   beam_kV; see qe_kinematic_prefactor), so E0 and the widths describe the
%   loss function itself. Unless weights='none' the residuals are weighted by
%   Poisson noise, sigma = sqrt(max(Y, 1)); chi2_red is then relative when the
%   counts are not single-electron counts.
%
%   The fields shared with fit_loss_function (n_peaks, omega_p, gamma,
%   apex_energy_meV, gamma_ratio, peak_valid, energy_fit, curve_fit, bg_curve,
%   peak_curves, R_squared, peak_model_name) let interactive_qe_browser draw
%   and accept the result as it does for the background + peak fit. The curves
%   are in the units of Y on E >= 50 meV up to signal_window(2); bg_curve is
%   the ZLP plus auxiliary peaks. apex_energy_meV is the maximum of each peak's
%   loss function (Drude-Lorentz 'lorentz': analytic apex; 'lorentz_symmetric':
%   E0). The full qe_zlp_joint_fit output is in result.joint.
%
%   See also: qe_zlp_joint_fit, qe_kinematic_prefactor, fit_loss_function

arguments
    E (:,1) double
    Y (:,1) double
    q_Ainv (1,1) double
    options.n_peaks (1,1) double {mustBeInteger, mustBePositive} = 1
    options.peak_model char {mustBeMember(options.peak_model, {'lorentz', 'lorentz_symmetric'})} = 'lorentz'
    options.signal_window (1,2) double = [300 1800]
    options.kinematic (1,1) logical = true
    options.beam_kV (1,1) double {mustBePositive} = 30
    options.dq (1,1) double {mustBeNonnegative} = 0
    options.weights char {mustBeMember(options.weights, {'poisson', 'none'})} = 'poisson'
    options.initial_E0 (1,:) double = []
    options.n_starts (1,1) double {mustBeInteger, mustBePositive} = 8
end

assert(numel(E) == numel(Y), 'qe_fit_spectrum_joint:Size', 'E and Y must have the same length');
win = sort(options.signal_window);
args = {'fit_window', [-Inf win(2)], 'signal_window', win, 'peak_model', options.peak_model, ...
    'n_peaks', options.n_peaks, 'n_starts', options.n_starts};
if strcmp(options.weights, 'poisson')
    args = [args, {'noise_sigma', sqrt(max(Y, 1))}];
end
guess = options.initial_E0(isfinite(options.initial_E0));
if numel(guess) == options.n_peaks
    args = [args, {'initial_E0', sort(guess)}];
end
kin = struct('used', false, 'beam_kV', options.beam_kV, 'q_Ainv', q_Ainv, 'dq', options.dq);
if options.kinematic
    assert(abs(q_Ainv) > 0 || options.dq > 0, 'qe_fit_spectrum_joint:ZeroQ', ...
        'the kinematic prefactor at q = 0 needs the channel width dq');
    Kfun = qe_kinematic_prefactor(abs(q_Ainv), dq=options.dq, form='2d', beam_kV=options.beam_kV);
    args = [args, {'prefactor', Kfun}];
    kin.used = true;
end

f = qe_zlp_joint_fit(E, Y, args{:});

np = options.n_peaks;
result = struct();
result.method = 'zlp_joint';
kin_label = '';
if kin.used, kin_label = ' + kinematic'; end
result.peak_model_name = sprintf('ZLP joint, %s%s', local_model_label(options.peak_model), kin_label);
result.kinematic = kin;
result.weights = options.weights;
result.signal_window = win;
result.success = f.success;
result.joint = f;
result.n_peaks = np;
if ~f.success
    [result.omega_p, result.gamma, result.amplitude, result.apex_energy_meV, result.gamma_ratio] = deal(nan(np, 1));
    result.peak_valid = false(np, 1);
    result.energy_fit = zeros(0, 1); result.curve_fit = zeros(0, 1); result.bg_curve = zeros(0, 1);
    result.peak_curves = repmat({zeros(0, 1)}, 1, np);
    result.R_squared = NaN; result.chi2_red = NaN; result.status = f.numerical_status;
    return
end

result.omega_p = f.parameters(:, 1);
result.gamma = f.parameters(:, 2);
result.amplitude = f.parameters(:, 3);
result.apex_energy_meV = local_apex(options.peak_model, result.omega_p, result.gamma);
result.gamma_ratio = result.gamma ./ result.omega_p;
result.peak_valid = ~any(f.peak_boundary, 2) & ~f.zero_amplitude(:);
result.status = f.numerical_status;

show = E >= 50 & E <= win(2);
result.energy_fit = E(show);
result.curve_fit = f.prediction(show);
result.bg_curve = f.background(show);
result.peak_curves = arrayfun(@(p) f.peaks(show, p), 1:np, 'UniformOutput', false);
sig = E >= win(1) & E <= win(2) & isfinite(Y);
res = Y(sig) - f.prediction(sig);
result.R_squared = 1 - sum(res .^ 2) / sum((Y(sig) - mean(Y(sig))) .^ 2);
result.chi2_red = f.chi2_red_signal;
end


function apex = local_apex(model, E0, G)
% Maximum of each peak's loss function. Drude-Lorentz E G / ((E^2 - E0^2)^2 + E^2 G^2):
% E_apex^2 = (c + sqrt(c^2 + 12 E0^4)) / 6 with c = 2 E0^2 - G^2. Symmetric Lorentzian: E0.
if strcmp(model, 'lorentz')
    c = 2 * E0 .^ 2 - G .^ 2;
    apex = sqrt((c + sqrt(c .^ 2 + 12 * E0 .^ 4)) / 6);
else
    apex = E0;
end
end


function label = local_model_label(model)
if strcmp(model, 'lorentz')
    label = 'DL';
else
    label = 'sym. Lorentz';
end
end
