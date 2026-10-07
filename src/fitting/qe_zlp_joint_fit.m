function fit = qe_zlp_joint_fit(E, Y, options)
%QE_ZLP_JOINT_FIT  Joint ZLP + loss-peak fit with no subtracted background.
%
%   fit = qe_zlp_joint_fit(E, Y, signal_window=[300 1800], n_peaks=2)
%
%   Model (method follows EELS-ZEST, R. Mao 2026, re-implemented for q-EELS):
%
%     Y(E) = sum_k a_k Z_k(E - c) + sum_j A_j [f_j (*) K](E) (+ b)
%
%   Z_k   asymmetric Pearson VII ZLP components, [1 + (x/s)^2/m]^(-m) with
%         separate (s, m) on the gain (x<0) and loss (x>=0) sides.
%   K     instrument kernel = fitted ZLP, truncated to +/-kernel_halfwidth and
%         normalised to unit area ('fitted_zlp'); 'core' uses the narrowest
%         component, 'none' skips the convolution.
%   f_j   peak_models shape evaluated on the loss side; the gain side is its
%         detailed-balance mirror exp(-|E|/kT) f_j(|E|). E0 and widths are
%         therefore intrinsic (deconvolved) and measured from the ZLP centre.
%
%   The gain side (E<c) carries no plasmon weight at E0 >> kT, so it fixes the
%   ZLP tail; a soft prior on log(sR/sL) and log(mR/mL) transfers that shape to
%   the loss side. Nonlinear parameters (c, s, m, E0, width) are optimised by
%   lsqnonlin; amplitudes are solved inside every evaluation as a nonnegative
%   least-squares problem (variable projection), optionally with the
%   inequality ZLP(E) <= Y(E) + margin on the loss side.
%
%   Loss features between the ZLP core and signal_window (phonons, low-energy
%   excitations) either get their own detailed-balance peaks (aux_windows,
%   one aux_model peak per row, reported in aux_parameters and counted as
%   background for the signal window) or are dropped from the data term
%   (exclude_windows). The overshoot inequality still applies there.
%
%   Residuals use a Poisson-like noise model sigma = g*sqrt(local level), with
%   g estimated from the high-pass scatter, so chi2_red ~ 1 for a good fit.
%
%   Requires a uniform energy grid inside fit_window.
%
%   See also: peak_models, qe_compare_component_models, fit_loss_function

arguments
 E (:,1) double
 Y (:,1) double
 options.fit_window (1,2) double = [-Inf Inf]
 options.signal_window (1,2) double = [300 1800]
 options.peak_model char {mustBeMember(options.peak_model,{'lorentz','lorentz_symmetric'})} = 'lorentz_symmetric'
 options.n_peaks (1,1) double {mustBeNonnegative,mustBeInteger} = 1
 options.initial_E0 (1,:) double = []
 options.initial_width (1,:) double = []
 options.n_zlp (1,1) double {mustBeMember(options.n_zlp,[1 2 3])} = 2
 options.core_halfwidth (1,1) double = NaN
 options.core_weight (1,1) double {mustBeNonnegative} = 0.05
 options.gain_weight (1,1) double {mustBeNonnegative} = 3
 options.prefit_loss_max (1,1) double = NaN
 options.asym_sigma (1,1) double {mustBePositive} = 0.15
 options.guard_weight (1,1) double {mustBeNonnegative} = 0
 options.overshoot (1,1) logical = true
 options.overshoot_sigma (1,1) double {mustBeNonnegative} = 3
 options.kernel_mode char {mustBeMember(options.kernel_mode,{'fitted_zlp','core','none'})} = 'fitted_zlp'
 options.kernel_halfwidth (1,1) double {mustBePositive} = 300
 options.temperature_K (1,1) double {mustBePositive} = 300
 options.include_constant (1,1) logical = false
 options.exclude_windows (:,2) double = zeros(0, 2)
 options.aux_windows (:,2) double = zeros(0, 2)
 options.aux_model char {mustBeMember(options.aux_model,{'lorentz','lorentz_symmetric'})} = 'lorentz'
 options.n_starts (1,1) double {mustBePositive,mustBeInteger} = 8
 options.seed (1,1) double = 20261007
 options.max_iterations (1,1) double {mustBePositive,mustBeInteger} = 400
end

%% Data on a uniform grid
keep = E >= options.fit_window(1) & E <= options.fit_window(2);
E = E(keep); Y = Y(keep);
assert(numel(E) >= 20 && all(diff(E) > 0), 'qe_zlp_joint_fit:Axis', 'Need >=20 increasing energies.');
dE = median(diff(E));
assert(max(abs(diff(E) - dE)) < 1e-6 * max(1, dE), 'qe_zlp_joint_fit:NonUniform', ...
    'Energy grid must be uniform inside fit_window.');
valid = isfinite(Y);
assert(any(E < 0 & valid) && any(E > 0 & valid), 'qe_zlp_joint_fit:Range', ...
    'fit_window must contain both gain (E<0) and loss (E>0) data.');
Yf = Y; Yf(~valid) = 0;

%% ZLP landmarks, noise model, weights
[zmax, imax] = max(Yf .* (abs(E) < 100));
c0 = E(imax);
fwhm0 = local_fwhm(E, Yf, imax, dE);
if isnan(options.core_halfwidth), options.core_halfwidth = 2.5 * fwhm0; end
if isnan(options.prefit_loss_max), options.prefit_loss_max = fwhm0; end
kT = 0.08617333262 * options.temperature_K;

local = movmedian(Yf, 7);
level_floor = max(1e-6 * zmax, prctile(local(valid & local > 0), 1));
level = max(local, level_floor);
hp = (Yf - movmean(Yf, 5)) ./ sqrt(level) / sqrt(0.8);
noise_region = valid & abs(E - c0) > options.core_halfwidth;
g = 1.4826 * median(abs(hp(noise_region) - median(hp(noise_region))));
if ~(isfinite(g) && g > 0), g = 1; end
sigma = g * sqrt(level);

core = 1 ./ (1 + exp((abs(E - c0) - options.core_halfwidth) / max(dE, 1)));
mult = options.core_weight * core + (1 - core);
mult(E < c0 - options.core_halfwidth) = mult(E < c0 - options.core_halfwidth) * options.gain_weight;
sqrtW = sqrt(mult) ./ sigma;
sqrtW(~valid) = 0;
for k = 1:size(options.exclude_windows, 1)
    sqrtW(E >= options.exclude_windows(k, 1) & E <= options.exclude_windows(k, 2)) = 0;
end

sig = E >= options.signal_window(1) & E <= options.signal_window(2);
guard = sig .* sqrt(options.guard_weight) ./ sigma;
guard(~valid) = 0;

margin = options.overshoot_sigma * sigma;
over_mask = valid & E > c0 + options.core_halfwidth;

%% Fixed problem description shared by both stages
P = struct('E', E, 'Y', Yf, 'valid', valid, 'dE', dE, 'sqrtW', sqrtW, 'guard', guard, ...
    'margin', margin, 'over_mask', over_mask, 'kT', kT, 'nz', options.n_zlp, 'opts', options);
P.models = {};
lags = (-round(options.kernel_halfwidth / dE):round(options.kernel_halfwidth / dE)).' * dE;
P.lags = lags;
P.E_ext = (E(1) - lags(end)) + (0:(numel(E) + numel(lags) - 2)).' * dE;

%% ZLP bounds and start (component 1 narrow core, then progressively wider)
nz = options.n_zlp;
s_core = max(fwhm0 / 2, dE / 2);
zl = zeros(1, 4 * nz); zu = zl; z0 = zl;
for k = 1:nz
    if k == 1
        s_lo = dE / 4; s_hi = 3 * fwhm0; s_init = s_core; m_init = 3;
    else
        s_lo = fwhm0 / 2; s_hi = 20 * fwhm0; s_init = s_core * 2.5^(k - 1); m_init = 1;
    end
    ii = 4 * (k - 1) + (1:4);
    zl(ii) = [log(s_lo) log(s_lo) log(0.55) log(0.55)];
    zu(ii) = [log(s_hi) log(s_hi) log(50) log(50)];
    z0(ii) = [log(s_init) log(s_init) log(m_init) log(m_init)];
end
u_lo = [c0 - 3 * dE, zl];
u_hi = [c0 + 3 * dE, zu];
u0 = [c0, z0];

lsq = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', options.max_iterations, ...
    'MaxFunctionEvaluations', 200 * options.max_iterations, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

%% Stage A: ZLP-only prefit on the gain side and the loss-side core
PA = P; PA.np = 0;
PA.sqrtW(E > c0 + options.prefit_loss_max) = 0;
PA.guard(:) = 0;
[uA, ~, ~, flagA] = lsqnonlin(@(u) local_residual(u, PA), u0, u_lo, u_hi, lsq);

%% Stage B: joint fit, multistart over peak seeds (aux peaks mid-window in starts 1-2, random after)
np = options.n_peaks;
aw = options.aux_windows;
na = size(aw, 1);
P.np = np + na;
P.models = [repmat({peak_models(options.peak_model)}, 1, np), repmat({peak_models(options.aux_model)}, 1, na)];
sw = options.signal_window;
p_lo = [repmat([sw(1) log(2 * dE)], 1, np), reshape([aw(:, 1) log(2 * dE) * ones(na, 1)].', 1, [])];
p_hi = [repmat([sw(2) log(5000)], 1, np), reshape([aw(:, 2) log(diff(aw, 1, 2))].', 1, [])];
aux0 = reshape([mean(aw, 2) log(diff(aw, 1, 2) / 2)].', 1, []);
lb = [u_lo p_lo]; ub = [u_hi p_hi];
stream = RandStream('mt19937ar', 'Seed', options.seed);
widths = [80 250 650 1500];
n_starts = options.n_starts;
if np == 0, n_starts = 1; end
cand = struct('start', cell(1, n_starts), 'u0', [], 'u', [], 'cost', Inf, 'exitflag', NaN, 'message', '');
for st = 1:n_starts
    pk0 = zeros(1, 2 * np);
    if np > 0
        if st == 1 && ~isempty(options.initial_E0)
            e0 = options.initial_E0(1:np);
            w0 = 300 * ones(1, np);
            if ~isempty(options.initial_width), w0 = options.initial_width(1:np); end
        elseif st <= 2
            e0 = linspace(sw(1), sw(2), np + 2); e0 = e0(2:end - 1);
            w0 = widths(1 + mod(st, 4)) * ones(1, np);
        else
            e0 = sort(sw(1) + diff(sw) * (0.05 + 0.9 * rand(stream, 1, np)));
            w0 = widths(1 + mod(st - 1, 4)) * (0.6 + 0.8 * rand(stream, 1, np));
        end
        pk0(1:2:end) = e0; pk0(2:2:end) = log(w0);
    end
    ax0 = aux0;
    if st > 2 && na > 0
        ax0(1:2:end) = aw(:, 1).' + diff(aw, 1, 2).' .* (0.1 + 0.8 * rand(stream, 1, na));
        ax0(2:2:end) = log(2 * dE) + rand(stream, 1, na) .* (log(diff(aw, 1, 2)).' - log(2 * dE));
    end
    x0 = min(max([uA pk0 ax0], lb + 1e-9), ub - 1e-9);
    cand(st).start = st; cand(st).u0 = x0;
    try
        [u, cost, ~, flag, out] = lsqnonlin(@(u) local_residual(u, P), x0, lb, ub, lsq);
        cand(st).u = u; cand(st).cost = cost; cand(st).exitflag = flag; cand(st).message = out.message;
    catch ME
        cand(st).message = [ME.identifier ': ' ME.message];
    end
end
ok = find([cand.exitflag] > 0 & isfinite([cand.cost]));
fit = struct('success', false, 'numerical_status', 'failed', 'candidates', cand, ...
    'prefit_exitflag', flagA, 'effective_options', options);
if isempty(ok), return; end
[~, kbest] = min([cand(ok).cost]);
best = ok(kbest);
u = cand(best).u;

%% Package
[~, M] = local_residual(u, P);
th = local_decode(u, P);
fit.success = true;
fit.selected_start = best;
fit.cost = cand(best).cost;
fit.energy_meV = E;
fit.observed = Y;
fit.sigma = sigma;
fit.weights = sqrtW .^ 2;
fit.prediction = M.prediction;
fit.residual = Y - M.prediction;
fit.zlp = M.zlp;
fit.zlp_components = M.zlp_components;
fit.constant = M.constant;
fit.aux_peaks = M.peaks(:, np + (1:na));
fit.background = M.zlp + M.constant + sum(fit.aux_peaks, 2);
[~, MA] = local_residual(uA, PA);
fit.prefit_zlp = MA.zlp;
fit.kernel = struct('lag_meV', lags, 'K', M.kernel, 'mode', options.kernel_mode);
fit.zlp_parameters = struct('center_meV', th.c, 'height', M.coef(1:nz).', ...
    'sL', th.sL, 'sR', th.sR, 'mL', th.mL, 'mR', th.mR, ...
    'fwhm_meV', local_fwhm(E, M.zlp, find(E >= th.c, 1), dE));
pars = nan(np + na, 3);
for j = 1:np + na
    pars(j, :) = [th.E0(j) th.W(j) M.coef(nz + j) / M.col_scale(j)];
end
[~, order] = sort(pars(1:np, 1));
fit.parameters = pars(order, :);
fit.parameter_names = {'E0_meV', 'width_meV', 'amplitude'};
fit.aux_parameters = pars(np + (1:na), :);
fit.peaks = M.peaks(:, order);
fit.peaks_intrinsic = M.peaks_intrinsic(:, order);

r = fit.residual ./ sigma;
dof = max(nnz(valid & sig) - (2 + 3 * (np + na)), 1);
fit.chi2_red_signal = sum(r(valid & sig) .^ 2) / dof;
fit.chi2_red_gain = mean(r(valid & E < c0 - options.core_halfwidth) .^ 2);
ys = Y(valid & sig);
fit.r2_signal = 1 - sum(fit.residual(valid & sig) .^ 2) / sum((ys - mean(ys)) .^ 2);
fit.zlp_fraction_signal = sum(M.zlp(valid & sig)) / sum(ys);
fit.background_fraction_signal = sum(fit.background(valid & sig)) / sum(ys);
fit.overshoot_active = any(M.zlp(over_mask) > Y(over_mask) + margin(over_mask) - 1e-9 * zmax);

span = ub - lb;
hit = abs(u - lb) < 1e-4 * span | abs(u - ub) < 1e-4 * span;
fit.boundary = hit;
pb = reshape(hit(numel(u_lo) + 1:end), 2, []).';
fit.peak_boundary = pb(order, :);
fit.aux_boundary = pb(np + (1:na), :);
fit.zero_amplitude = fit.parameters(:, 3).' <= 0;
fit.numerical_status = 'converged';
if any(fit.peak_boundary(:)), fit.numerical_status = 'boundary'; end
if any(fit.zero_amplitude), fit.numerical_status = 'zero_amplitude'; end
if np == 2 && abs(diff(fit.parameters(:, 1))) < dE, fit.numerical_status = 'collapse_sampling_scale'; end
end


function [r, M] = local_residual(u, P)
th = local_decode(u, P);
E = P.E; n = numel(E); nz = P.nz; np = P.np;
x = E - th.c;
Zc = zeros(n, nz);
for k = 1:nz
    Zc(:, k) = local_pearson(x, th.sL(k), th.sR(k), th.mL(k), th.mR(k));
end

% Peak columns, convolved with a first-pass kernel, then with the fitted ZLP.
F = zeros(numel(P.E_ext), np);
for j = 1:np
    F(:, j) = local_peak(P, th.E0(j), th.W(j), P.models{j});
end
xk = P.lags - th.c;
Kc = zeros(numel(P.lags), nz);
for k = 1:nz
    Kc(:, k) = local_pearson(xk, th.sL(k), th.sR(k), th.mL(k), th.mR(k));
end
switch P.opts.kernel_mode
    case 'none'
        K = double(abs(xk) == min(abs(xk)));
    otherwise
        K = Kc(:, 1);
end
[coef, A, col_scale] = local_solve(Zc, local_convolve(F, K / sum(K), P), P);
if strcmp(P.opts.kernel_mode, 'fitted_zlp') && np > 0 && any(coef(1:nz) > 0)
    K = Kc * coef(1:nz);
    [coef, A, col_scale] = local_solve(Zc, local_convolve(F, K / sum(K), P), P);
end

pred = A * coef;
r = P.sqrtW .* (pred - P.Y);
r = [r; P.guard .* (Zc * coef(1:nz))];
if isfinite(P.opts.asym_sigma)
    r = [r; (log(th.sR) - log(th.sL)).' / P.opts.asym_sigma; (log(th.mR) - log(th.mL)).' / P.opts.asym_sigma];
end

if nargout > 1
    M.coef = coef;
    M.col_scale = col_scale;
    M.prediction = pred;
    M.zlp_components = Zc .* coef(1:nz).';
    M.zlp = sum(M.zlp_components, 2);
    M.peaks = A(:, nz + (1:np)) .* coef(nz + (1:np)).';
    M.peaks_intrinsic = zeros(n, np);
    H = (numel(P.lags) - 1) / 2;
    inside = H + (1:n);
    for j = 1:np
        M.peaks_intrinsic(:, j) = F(inside, j) * coef(nz + j) / col_scale(j);
    end
    M.constant = 0;
    if P.opts.include_constant, M.constant = coef(end); end
    M.kernel = K / sum(K);
end
end


function [coef, A, col_scale] = local_solve(Zc, Pk, P)
% Nonnegative amplitudes; columns scaled to unit max for conditioning.
col_scale = max(abs(Pk), [], 1);
col_scale(col_scale == 0) = 1;
A = [Zc, Pk ./ col_scale];
if P.opts.include_constant, A = [A, ones(size(Zc, 1), 1)]; end
C = [P.sqrtW .* A; P.guard .* [Zc, zeros(size(Zc, 1), size(A, 2) - size(Zc, 2))]];
d = [P.sqrtW .* P.Y; zeros(size(Zc, 1), 1)];
cn = vecnorm(C); cn(cn == 0) = 1;
coef = lsqnonneg(C ./ cn, d) ./ cn.';
if P.opts.overshoot && any(P.over_mask)
    nz = size(Zc, 2);
    Ai = zeros(nnz(P.over_mask), size(A, 2));
    Ai(:, 1:nz) = Zc(P.over_mask, :);
    bi = P.Y(P.over_mask) + P.margin(P.over_mask);
    if any(Ai * coef > bi)
        o = optimoptions('lsqlin', 'Algorithm', 'active-set', 'Display', 'off');
        x = lsqlin(C ./ cn, d, Ai ./ cn, bi, [], [], zeros(size(coef)), [], zeros(size(coef)), o);
        coef = x ./ cn.';
    end
end
end


function Pk = local_convolve(F, K, P)
n = numel(P.E);
Pk = zeros(n, size(F, 2));
for j = 1:size(F, 2)
    Pk(:, j) = conv(F(:, j), K, 'valid');
end
end


function f = local_peak(P, E0, W, model)
% Loss side from peak_models; gain side is the detailed-balance mirror.
Ea = abs(P.E_ext);
f = model.model_fn(E0, W, 1, Ea);
neg = P.E_ext < 0;
f(neg) = f(neg) .* exp(-Ea(neg) / P.kT);
f(~isfinite(f)) = 0;
end


function y = local_pearson(x, sL, sR, mL, mR)
y = zeros(size(x));
l = x < 0;
y(l) = (1 + (x(l) / sL) .^ 2 / mL) .^ (-mL);
y(~l) = (1 + (x(~l) / sR) .^ 2 / mR) .^ (-mR);
end


function th = local_decode(u, P)
nz = P.nz;
z = reshape(u(2:1 + 4 * nz), 4, nz);
th.c = u(1);
th.sL = exp(z(1, :)); th.sR = exp(z(2, :));
th.mL = exp(z(3, :)); th.mR = exp(z(4, :));
pk = reshape(u(2 + 4 * nz:end), 2, []);
if isempty(pk), pk = zeros(2, 0); end
th.E0 = pk(1, :);
th.W = exp(pk(2, :));
end


function w = local_fwhm(E, Y, i0, dE)
% Interpolated full width at half maximum around index i0.
h = Y(i0) / 2;
l = find(Y(1:i0) < h, 1, 'last');
r = i0 - 1 + find(Y(i0:end) < h, 1, 'first');
if isempty(l) || isempty(r), w = 4 * dE; return; end
el = E(l) + (h - Y(l)) / (Y(l + 1) - Y(l)) * dE;
er = E(r - 1) + (Y(r - 1) - h) / (Y(r - 1) - Y(r)) * dE;
w = max(er - el, dE);
end
