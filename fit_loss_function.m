function result = fit_loss_function(energy_meV, spectrum, options)
%FIT_LOSS_FUNCTION  Multi-peak Drude-Lorentz + power-law background fit.
%
%   Fits the spectrum to:
%     S(E) = B0 * E^(-alpha) + sum_i { A_i * E * Gamma_i / ((E^2-Ep_i^2)^2 + E^2*Gamma_i^2) }
%
%   The power-law term models the ZLP tail explicitly, preventing the
%   Lorentz peaks from absorbing the background.
%
%   Name-value options:
%       E_min, E_max        - fitting window (meV)
%       max_peaks           - max peaks to detect (default: 3)
%       min_prominence      - peak prominence threshold (default: 0.15)
%       smooth_width        - smoothing window for peak detection (default: 25)

arguments
    energy_meV (:,1) double
    spectrum   (:,1) double
    options.E_min          (1,1) double = 50
    options.E_max          (1,1) double = Inf
    options.max_peaks      (1,1) double {mustBePositive, mustBeInteger} = 3
    options.min_prominence (1,1) double = 0.15
    options.smooth_width   (1,1) double {mustBePositive} = 25
    options.initial_guesses (:,1) double = []
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

%% Estimate power-law background from the data edges
% Fit log(S) ~ -alpha*log(E) + log(B0) using the first and last portions
n_pts = numel(E_w);
n_edge = max(5, round(n_pts * 0.12));
edge_idx = [1:n_edge, (n_pts-n_edge+1):n_pts]';
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

%% Peak detection: subtract estimated background, then findpeaks
bg_est = B0_init * E_w .^ (-alpha_init);
bg_est = min(bg_est, S_n * 0.95);  % don't over-subtract
S_detrended = max(S_n - bg_est, 0);

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

%% Build model: [B0, alpha, Ep1, Gamma1, A1, ..., EpN, GammaN, AN]
% Total params: 2 + 3*n_peaks
n_params = 2 + 3 * n_peaks;
p0 = zeros(n_params, 1);
lb = zeros(n_params, 1);
ub = zeros(n_params, 1);

% Background: B0, alpha
p0(1) = B0_init;      lb(1) = 0;     ub(1) = Inf;
p0(2) = alpha_init;   lb(2) = 0.01;  ub(2) = 5;

% Peaks
for i = 1:n_peaks
    base = 2 + 3*(i-1);
    
    p0(base+1) = locs(i);
    lb(base+1) = max(100, locs(i) - max(widths(i)*2, 400));
    ub(base+1) = min(options.E_max, locs(i) + max(widths(i)*2, 400));
    
    Gamma_guess = max(widths(i) * 0.5, 50);
    p0(base+2) = Gamma_guess;
    lb(base+2) = 10;
    ub(base+2) = 5000;
    
    % Amplitude guess: use detrended peak height
    [~, nearest_idx] = min(abs(E_w - locs(i)));
    pk_detrended = max(S_detrended(nearest_idx), pks(i));
    A_guess = pk_detrended * locs(i) * Gamma_guess;
    p0(base+3) = A_guess;
    lb(base+3) = 0;
    ub(base+3) = Inf;
end

model = @(p, E_in) local_model(p, E_in, n_peaks);

%% Fit
try
    fit_opts = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 15000, ...
        'MaxIterations', 3000, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12);
    p_fit = lsqcurvefit(model, p0, E_w, S_n, lb, ub, fit_opts);
catch
    cost = @(p) sum((model(p, E_w) - S_n).^2);
    fmin_opts = optimset('Display', 'off', 'MaxFunEvals', 15000, ...
        'MaxIter', 3000, 'TolFun', 1e-12, 'TolX', 1e-12);
    p_fit = fminsearch(cost, p0, fmin_opts);
end

%% Extract results
B0_fit = abs(p_fit(1));
alpha_fit = abs(p_fit(2));

omega_p_vals = zeros(n_peaks, 1);
gamma_vals = zeros(n_peaks, 1);
amp_vals = zeros(n_peaks, 1);

for i = 1:n_peaks
    base = 2 + 3*(i-1);
    omega_p_vals(i) = abs(p_fit(base+1));
    gamma_vals(i) = abs(p_fit(base+2));
    amp_vals(i) = abs(p_fit(base+3));
end

[omega_p_vals, si] = sort(omega_p_vals);
gamma_vals = gamma_vals(si);
amp_vals = amp_vals(si);

% Post-fit filter: remove peaks with amplitude < 10% of max
max_amp = max(amp_vals);
keep = amp_vals >= 0.1 * max_amp;
if any(keep)
    omega_p_vals = omega_p_vals(keep);
    gamma_vals = gamma_vals(keep);
    amp_vals = amp_vals(keep);
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
bg_curve = B0_fit * E_fine .^ (-alpha_fit) * S_scale;

% Total fit and individual peaks
curve_total = B0_fit * E_fine .^ (-alpha_fit) * S_scale;
peak_curves = cell(n_peaks, 1);
peak_heights = zeros(n_peaks, 1);

for i = 1:n_peaks
    Ep = omega_p_vals(i);
    Gm = gamma_vals(i);
    Am = amp_vals(i);
    single = Am .* E_fine .* Gm ./ ((E_fine.^2 - Ep^2).^2 + E_fine.^2 .* Gm^2);
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
result.B0 = B0_fit * S_scale;
result.alpha = alpha_fit;
result.R_squared = R_sq;
result.energy_fit = E_fine;
result.curve_fit = curve_total;
result.bg_curve = bg_curve;
result.peak_curves = peak_curves;
result.peak_heights = peak_heights;
result.energy_data = E_w;
result.spectrum_data = S_w;
result.residuals = (S_n - S_pred) * S_scale;

end


function S = local_model(p, E, n_peaks)
    % S(E) = B0 * E^(-alpha) + sum { Lorentz peaks }
    B0 = abs(p(1));
    alpha = abs(p(2));
    S = B0 * E .^ (-alpha);
    
    for i = 1:n_peaks
        base = 2 + 3*(i-1);
        Ep = abs(p(base+1));
        Gm = abs(p(base+2));
        Am = abs(p(base+3));
        S = S + Am .* E .* Gm ./ ((E.^2 - Ep^2).^2 + E.^2 .* Gm^2);
    end
end
