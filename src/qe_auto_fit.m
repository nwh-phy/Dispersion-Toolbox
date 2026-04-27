function results = qe_auto_fit(qe, qe_raw, opts)
%QE_AUTO_FIT  Automatic multi-peak Drude-Lorentz dispersion extraction.
%
%   results = qe_auto_fit(qe, qe_raw, opts)
%
%   Fits loss functions across all q-channels in the specified range,
%   extracts ω_p(q) and Γ(q), and separates into dispersion branches.
%
%   Inputs:
%     qe      — preprocessed qeData struct (possibly area-normalized)
%     qe_raw  — qeData without area normalization (for raw peak heights)
%     opts    — struct with fitting parameters:
%       .E_min         double — lower energy bound for fitting (meV)
%       .E_max         double — upper energy bound for fitting (meV)
%       .q_start       double — q range start
%       .q_end         double — q range end
%       .prominence    double — findpeaks min prominence
%       .smooth_width  double — pre-fit smoothing width
%       .max_peaks     double — max peaks per spectrum
%       .peak_model    char   — 'lorentz', 'gaussian', 'voigt', 'damped_ho'
%       .pre_subtracted logical — true when qe has already had background removed
%       .guesses       double — manual peak guesses (meV), [] for blind mode
%       .seed_idx      double — q-channel index for seed propagation
%       .max_shift     double — max shift between adjacent q-channels (meV)
%       .energy_mask   logical — mask for qe.energy_meV (from energy range)
%       .energy_axis   double  — vector of energies in the mask
%       .R2_threshold  double  — min R² to keep (default 0.3)
%       .verbose       logical — print progress (default false)
%       .progress_fn   function_handle — callback(fraction, message) or []
%
%   Output:
%     results struct with:
%       .all_peaks    — Nx12 matrix of fitted peaks
%       .branches     — cell array of per-branch subsets
%       .fit_details  — per-channel fit result structs
%       .n_success    — number of successfully fitted channels
%       .used_seed    — whether seed propagation was used
%
%   See also: fit_loss_function, propagate_seed_peaks, interactive_qe_browser

arguments
    qe      struct
    qe_raw  struct
    opts    struct
end

% --- Defaults ---
if ~isfield(opts, 'R2_threshold'), opts.R2_threshold = 0.3; end
if ~isfield(opts, 'verbose'), opts.verbose = false; end
if ~isfield(opts, 'progress_fn'), opts.progress_fn = []; end
if ~isfield(opts, 'pre_subtracted'), opts.pre_subtracted = false; end

mask = opts.energy_mask;
energy_axis = opts.energy_axis;
q_axis = qe.q_Ainv;
n_q = numel(q_axis);
q_start = opts.q_start;
q_end = opts.q_end;

% ═══════════ DUAL MODE: Seed vs Blind ═══════════
if ~isempty(opts.guesses)
    % ── SEED MODE: propagate from seed q-channel ──
    report_progress(opts, 0.1, sprintf("Seed propagation (%d guesses)...", numel(opts.guesses)));

    intensity = double(qe.intensity(mask, :));

    try
        prop = propagate_seed_peaks(intensity, energy_axis, q_axis(:), ...
            'seed_guesses', opts.guesses(:), ...
            'seed_idx', opts.seed_idx, ...
            'direction', 'both', ...
            'max_shift', opts.max_shift, ...
            'peak_model', opts.peak_model, ...
            'pre_subtracted', opts.pre_subtracted, ...
            'E_min', opts.E_min, 'E_max', opts.E_max, ...
            'smooth_width', opts.smooth_width, ...
            'verbose', false);
    catch ME
        error('qe_auto_fit:SeedFailed', 'Seed propagation failed: %s', ME.message);
    end

    if isempty(prop.peaks) || size(prop.peaks, 1) < 2
        error('qe_auto_fit:InsufficientPeaks', 'Seed propagation: insufficient peaks found');
    end

    all_peaks = prop.peaks;
    fit_details = prop.fit_details;
    n_success = prop.n_success;
    used_seed = true;

else
    % ── BLIND MODE: per-channel findpeaks ──
    all_peaks = [];
    fit_details = cell(n_q, 1);
    used_seed = false;

    n_in_range = sum(q_axis >= q_start & q_axis <= q_end);
    report_progress(opts, 0.05, sprintf("Fitting %d channels (blind)...", n_in_range));

    n_success = 0;
    for k = 1:n_q
        q_val = q_axis(k);
        if q_val < q_start || q_val > q_end
            continue
        end

        spectrum = double(qe.intensity(mask, k));
        if all(spectrum == 0) || all(isnan(spectrum))
            continue
        end

        try
            result = fit_loss_function(energy_axis, spectrum, ...
                'E_min', opts.E_min, 'E_max', opts.E_max, ...
                'min_prominence', opts.prominence, ...
                'smooth_width', opts.smooth_width, ...
                'max_peaks', opts.max_peaks, ...
                'initial_guesses', [], ...
                'peak_model', opts.peak_model, ...
                'pre_subtracted', opts.pre_subtracted);
            fit_details{k} = result;
            n_success = n_success + 1;

            % Raw (non-areanorm) peak heights at same positions
            raw_spectrum = double(qe_raw.intensity(mask, k));
            raw_peak_heights = NaN(result.n_peaks, 1);
            for p = 1:result.n_peaks
                raw_peak_heights(p) = measure_peak_height( ...
                    energy_axis, raw_spectrum, result.omega_p(p), result.gamma(p));
            end

            for p = 1:result.n_peaks
                if result.omega_p(p) < opts.E_min
                    continue
                end
                % [q, E, Γ, R², A, E_ci_lo, E_ci_hi, G_ci_lo, G_ci_hi, A_ci_lo, A_ci_hi, raw_h]
                all_peaks(end+1, :) = [ ...
                    q_val, ...
                    result.omega_p(p), ...
                    result.gamma(p), ...
                    result.R_squared, ...
                    result.amplitude(p), ...
                    result.omega_p_ci(p, 1), result.omega_p_ci(p, 2), ...
                    result.gamma_ci(p, 1), result.gamma_ci(p, 2), ...
                    result.amplitude_ci(p, 1), result.amplitude_ci(p, 2), ...
                    raw_peak_heights(p)]; %#ok<AGROW>
            end
        catch
        end

        if mod(k, 10) == 0
            report_progress(opts, k/n_q * 0.9, sprintf("Auto-fitting... %d%%", round(k/n_q*100)));
        end
    end
end

% ═══════════ Filter by R² ═══════════
if isempty(all_peaks) || size(all_peaks, 1) < 2
    error('qe_auto_fit:InsufficientPeaks', 'Auto-fit: insufficient peaks found');
end

good = all_peaks(:, 4) > opts.R2_threshold;
peaks = all_peaks(good, :);

if size(peaks, 1) < 2
    error('qe_auto_fit:InsufficientGoodPeaks', ...
        'Auto-fit: only %d good peaks (R²>%.2f)', size(peaks,1), opts.R2_threshold);
end

% ═══════════ Separate into branches ═══════════
if used_seed && size(peaks, 2) >= 6
    % Seed mode: use branch_id column
    unique_branches = unique(peaks(:, 6));
    n_branches = numel(unique_branches);
    branches = cell(n_branches, 1);
    for b = 1:n_branches
        b_id = unique_branches(b);
        branches{b} = peaks(peaks(:,6) == b_id, 1:5);
        [~, si] = sort(branches{b}(:,1));
        branches{b} = branches{b}(si, :);
    end
else
    % Blind mode: 1D gap-based clustering
    unique_q = unique(peaks(:,1));
    max_p = 1;
    for i = 1:numel(unique_q)
        max_p = max(max_p, sum(peaks(:,1) == unique_q(i)));
    end

    if max_p > 1
        energies = peaks(:,2);
        E_sorted = sort(energies);
        gaps = diff(E_sorted);
        [~, sort_gap_idx] = sort(gaps, 'descend');

        boundaries = sort(E_sorted(sort_gap_idx(1:max_p-1)));

        branches = cell(max_p, 1);
        for b = 1:max_p
            if b == 1
                bmask = peaks(:,2) <= boundaries(1);
            elseif b == max_p
                bmask = peaks(:,2) > boundaries(end);
            else
                bmask = peaks(:,2) > boundaries(b-1) & peaks(:,2) <= boundaries(b);
            end

            b_peaks = peaks(bmask, :);
            [~, si] = sort(b_peaks(:,1));
            branches{b} = b_peaks(si, :);
        end
    else
        branches = {peaks};
        [~, si] = sort(branches{1}(:,1));
        branches{1} = branches{1}(si, :);
    end

    branches = branches(~cellfun('isempty', branches));
end

report_progress(opts, 1.0, "Auto-fit complete");

% ═══════════ Return ═══════════
results = struct();
results.all_peaks = peaks;
results.branches = branches;
results.fit_details = fit_details;
results.n_success = n_success;
results.used_seed = used_seed;

end


%% ═══════════ Local helpers ═══════════

function report_progress(opts, fraction, message)
    if isa(opts.progress_fn, 'function_handle')
        opts.progress_fn(fraction, message);
    end
    if opts.verbose
        fprintf('[qe_auto_fit] %.0f%% — %s\n', fraction*100, message);
    end
end
