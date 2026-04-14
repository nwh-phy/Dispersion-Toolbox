function [qe_out, bg_diag] = qe_preprocess(qe_in, opts)
%QE_PREPROCESS  Apply EELS preprocessing pipeline to a qeData struct.
%
%   qe_out = qe_preprocess(qe_in, opts)
%   [qe_out, bg_diag] = qe_preprocess(qe_in, opts)
%
%   Applies a configurable preprocessing chain:
%     0. Cosmic-ray removal (optional, median-threshold spike detection)
%     1. Normalization  (ZLP Peak or Area — beam current correction)
%     2. Denoising      (Wiener2D, Savitzky-Golay, or BM3D)
%     3. Quasi-elastic background removal (ZLP tail subtraction)
%        Models: Power/ExpPoly3/Pearson/PearsonVII/Power2/Exp2/Auto
%     4. Deconvolution  (Lucy-Richardson)
%
%   Physical note (Do et al. 2025, Fung et al. 2020):
%     Step 3 removes the quasi-elastic ZLP tail via **subtraction**, not
%     division. This preserves the absolute intensity for A(q) analysis.
%     ZLP Peak normalization (step 1) corrects for beam current only.
%
%   Background subtraction opts fields (all have sensible defaults):
%     .bg_method      char     — model name (Power/ExpPoly3/Pearson/PearsonVII/Power2/Exp2/Auto)
%     .bg_win_lo      [lo hi]  — primary fit window in meV (default [50 300])
%     .bg_win_hi      [lo hi]  — secondary fit window in meV ([] = single window)
%     .bg_iterative   logical  — enable iterative refit (default false)
%     .bg_asym_penalty double  — negative-residual penalty factor (default 1)
%
%   Optional second output bg_diag is a struct array (one per q-channel):
%     .rsquare, .rmse, .bg_curve, .residuals, .h_param, .snr,
%     .energy_fit, .iterations, .selected_method, .candidate_methods,
%     .candidate_scores, .candidate_details, .linear_rmse,
%     .neg_fraction, .neg_area_fraction, .bg_fraction
%
%   See also: interactive_qe_browser, fit_loss_function

arguments
    qe_in   struct
    opts    struct
end

qe_out = qe_in;
bg_diag = [];

if isfield(opts, 'do_despike') && opts.do_despike
    qe_out = apply_cosmic_ray_removal(qe_out, opts);
end

if isfield(opts, 'do_normalize') && opts.do_normalize
    qe_out = apply_normalization(qe_out, opts);
end

if isfield(opts, 'do_denoise') && opts.do_denoise
    qe_out = apply_denoise(qe_out, opts);
end

if isfield(opts, 'do_bg_sub') && opts.do_bg_sub
    want_diag = nargout >= 2;
    [qe_out, bg_diag] = apply_bg_subtraction(qe_out, opts, want_diag);
end

if isfield(opts, 'do_deconv') && opts.do_deconv
    qe_out = apply_deconvolution(qe_out, opts);
end

end


%% ═══════════ Normalization ═══════════

function qe_out = apply_normalization(qe_in, opts)
    qe_out = qe_in;
    energy_axis = double(qe_in.energy_meV(:));
    norm_min = opts.norm_min;
    norm_max = opts.norm_max;
    intensity = double(qe_in.intensity);
    method = char(opts.norm_method);

    if strcmpi(method, 'ZLP Peak')
        % ZLP Peak normalization: divide by ZLP peak height.
        % Preserves q-dependent spectral weight for A(q) studies.
        %
        % CRITICAL: Always use fixed [-50, 50] meV window for ZLP Peak
        % normalization, regardless of user display range. The ZLP is
        % at E≈0 by definition; using the display range (e.g. [300, 3836])
        % would find the plasmon peak instead of the zero-loss peak.
        zlp_fixed_min = -50;
        zlp_fixed_max =  50;
        zlp_mask = energy_axis >= zlp_fixed_min & energy_axis <= zlp_fixed_max;
        if ~any(zlp_mask)
            [~, nearest] = min(abs(energy_axis));
            zlp_mask(nearest) = true;
        end
        for qi = 1:size(intensity, 2)
            zlp_peak = max(intensity(zlp_mask, qi));
            if isfinite(zlp_peak) && zlp_peak > eps
                intensity(:, qi) = intensity(:, qi) ./ zlp_peak;
            end
        end
    else
        % Area normalization: divide by integrated spectral weight.
        % WARNING: destroys q-dependent A(q) scaling information.
        norm_mask = energy_axis >= norm_min & energy_axis <= norm_max;
        if ~any(norm_mask)
            norm_mask = true(size(energy_axis));
        end
        norm_energy = energy_axis(norm_mask);
        for qi = 1:size(intensity, 2)
            window_spec = max(intensity(norm_mask, qi), 0);
            if numel(norm_energy) > 1
                area = trapz(norm_energy, window_spec);
            else
                area = sum(window_spec);
            end
            if isfinite(area) && area > eps
                intensity(:, qi) = intensity(:, qi) ./ area;
            end
        end
    end
    qe_out.intensity = intensity;
end


function [qe_out, bg_diag] = apply_bg_subtraction(qe_in, opts, want_diag)
    % Quasi-elastic background removal (Do et al. 2025).
    %
    % Fits a model to the ZLP tail (quasi-elastic scattering) in the
    % low-loss region and subtracts it, isolating the inelastic loss
    % features (plasmon peaks, interband transitions).
    %
    % Features:
    %   - Curve Fitting Toolbox fit() with full gof stats
    %   - Dual-window fitting: anchor before AND after loss features
    %   - Iterative refit: detect negative signal → increase weight → refit
    %   - Models: Power, Power2, ExpPoly3, Pearson, PearsonVII, Exp2
    %     (Exp2 recommended for low-loss; Fung et al. 2020)
    %   - Diagnostics: R², RMSE, h-parameter, SNR per q-channel
    %
    % Physical note:
    %   This is "quasi-elastic background removal" — the ZLP tail extends
    %   far beyond its FWHM and must be subtracted before peak fitting.
    %   Do NOT confuse with "normalization" (beam current correction).

    if nargin < 3, want_diag = false; end

    qe_out = qe_in;
    energy_axis = double(qe_in.energy_meV(:));
    intensity = double(qe_in.intensity);
    n_q = size(intensity, 2);
    method = char(opts.bg_method);
    candidate_methods = local_background_model_candidates(opts, method);

    % --- Defaults for new opts fields (backward compatible) ---
    if isfield(opts, 'bg_win_lo') && numel(opts.bg_win_lo) == 2
        win1 = opts.bg_win_lo;
    else
        win1 = [50, 300];
    end
    if isfield(opts, 'bg_win_hi') && ~isempty(opts.bg_win_hi) && numel(opts.bg_win_hi) == 2
        win2 = opts.bg_win_hi;
        use_dual = true;
    else
        win2 = [];
        use_dual = false;
    end
    if isfield(opts, 'bg_iterative') && opts.bg_iterative
        do_iterative = true;
    else
        do_iterative = false;
    end
    if isfield(opts, 'bg_asym_penalty') && opts.bg_asym_penalty > 1
        asym_penalty = opts.bg_asym_penalty;
    else
        asym_penalty = 1;
    end

    % --- Build fit mask ---
    fit_mask = energy_axis >= win1(1) & energy_axis <= win1(2);
    if use_dual
        fit_mask = fit_mask | (energy_axis >= win2(1) & energy_axis <= win2(2));
    end
    if nnz(fit_mask) < 5
        bg_diag = [];
        return
    end

    % --- Signal region mask (between windows, for iterative check) ---
    if use_dual
        signal_region = energy_axis > win1(2) & energy_axis < win2(1);
    else
        signal_region = energy_axis > win1(2);
    end

    e_fit = energy_axis(fit_mask);
    sub_mask = energy_axis > win1(1);
    e_sub = energy_axis(sub_mask);

    % --- Pre-allocate diagnostics ---
    if want_diag
        empty_diag = struct('rsquare', NaN, 'rmse', NaN, ...
            'bg_curve', zeros(size(energy_axis)), ...
            'residuals', zeros(nnz(fit_mask), 1), ...
            'h_param', NaN, 'snr', NaN, ...
            'energy_fit', e_fit, 'iterations', 1, ...
            'selected_method', char(method), ...
            'candidate_methods', {candidate_methods}, ...
            'candidate_scores', NaN(1, numel(candidate_methods)), ...
            'candidate_details', {struct([])}, ...
            'linear_rmse', NaN, ...
            'neg_fraction', NaN, 'neg_area_fraction', NaN, ...
            'bg_fraction', NaN);
        bg_diag = repmat(empty_diag, 1, n_q);
    else
        bg_diag = [];
    end

    % Suppress Curve Fitting Toolbox warnings during batch fitting
    warnState = warning;
    warning('off', 'curvefit:fit:noStartPoint');
    warning('off', 'curvefit:fit:iterationLimitReached');
    warning('off', 'curvefit:fit:nonDoubleYData');
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:illConditionedMatrix');
    cleanupObj = onCleanup(@() warning(warnState));

    % --- Selective q-channel processing for fast preview ---
    if isfield(opts, 'q_indices') && ~isempty(opts.q_indices)
        q_list = opts.q_indices(:)';
        q_list = q_list(q_list >= 1 & q_list <= n_q);
    else
        q_list = 1:n_q;
    end

    for qi = q_list
        spec = intensity(:, qi);
        s_fit = spec(fit_mask);
        s_fit_pos = max(s_fit, eps);

        try
            [bg, rsq, rmse_val, residuals, n_iter, selected_method, candidate_scores, candidate_details] = ...
                local_fit_background_dispatch( ...
                    e_fit, s_fit_pos, e_sub, method, candidate_methods, ...
                    asym_penalty, spec, sub_mask);

            % --- Iterative refit (Ruishi Qi style) ---
            if do_iterative && any(signal_region)
                bg_full = interp1(e_sub, bg, energy_axis, 'linear', 0);
                neg_in_signal = (spec - bg_full < 0) & signal_region;
                if any(neg_in_signal)
                    % Increase weight on negative-signal points and refit
                    wt = ones(size(e_fit));
                    neg_fit_overlap = fit_mask & neg_in_signal;
                    % Map back to fit indices
                    neg_idx_in_fit = neg_fit_overlap(fit_mask);
                    wt(neg_idx_in_fit) = 2;
                    [bg, rsq, rmse_val, residuals, ~] = fit_background( ...
                        e_fit, s_fit_pos, e_sub, selected_method, asym_penalty, wt);
                    bg = local_apply_background_safety_cap(spec(sub_mask), bg);
                    candidate_details = local_refresh_selected_candidate_details( ...
                        candidate_details, selected_method, spec(sub_mask), bg, rsq, rmse_val, residuals);
                    candidate_scores = [candidate_details.score];
                    n_iter = 2;
                end
            end

            bg = real(bg(:));
            bg(~isfinite(bg)) = 0;
            bg = max(bg, 0);

            subtracted_signal = spec(sub_mask) - bg;
            intensity(sub_mask, qi) = subtracted_signal;

            % --- Compute diagnostics ---
            if want_diag
                bg_diag(qi).rsquare = rsq;
                bg_diag(qi).rmse = rmse_val;
                bg_diag(qi).iterations = n_iter;
                bg_diag(qi).selected_method = char(selected_method);
                bg_diag(qi).candidate_methods = candidate_methods;
                bg_diag(qi).candidate_scores = candidate_scores;
                bg_diag(qi).candidate_details = candidate_details;

                selected_idx = find(strcmp(candidate_methods, selected_method), 1);
                if ~isempty(selected_idx) && selected_idx <= numel(candidate_details)
                    bg_diag(qi).linear_rmse = candidate_details(selected_idx).linear_rmse;
                end

                % Full background curve for visualization
                bg_full_vis = zeros(size(energy_axis));
                bg_full_vis(sub_mask) = bg;
                bg_diag(qi).bg_curve = bg_full_vis;
                bg_diag(qi).residuals = residuals;

                % h-parameter (Fung et al.): h = (I_B + var(I_B)) / I_B
                I_B = mean(bg, 'omitnan');
                var_I_B = var(bg, 'omitnan');
                if I_B > eps
                    bg_diag(qi).h_param = (I_B + var_I_B) / I_B;
                end

                % SNR (Fung et al.): SNR = I_k / sqrt(I_k + h * I_B)
                I_k = max(max(subtracted_signal), eps);
                h = bg_diag(qi).h_param;
                if isfinite(h) && I_B > eps
                    bg_diag(qi).snr = I_k / sqrt(abs(I_k) + h * I_B);
                end

                [neg_fraction, neg_area_fraction, bg_fraction] = ...
                    local_background_quality_metrics(spec(sub_mask), subtracted_signal, bg);
                bg_diag(qi).neg_fraction = neg_fraction;
                bg_diag(qi).neg_area_fraction = neg_area_fraction;
                bg_diag(qi).bg_fraction = bg_fraction;
            end
        catch
        end
    end
    qe_out.intensity = intensity;
end


function methods = local_background_model_candidates(opts, requested_method)
    supported = {'Power', 'ExpPoly3', 'Pearson', 'PearsonVII', 'Power2', 'Exp2'};

    if strcmpi(requested_method, 'Auto')
        if isfield(opts, 'bg_candidate_methods') && ~isempty(opts.bg_candidate_methods)
            raw_methods = opts.bg_candidate_methods;
        else
            raw_methods = supported;
        end
    else
        raw_methods = {requested_method};
    end

    methods = {};
    for k = 1:numel(raw_methods)
        method = char(string(raw_methods(k)));
        if any(strcmpi(method, supported)) && ~any(strcmpi(method, methods))
            methods{end+1} = method; %#ok<AGROW>
        end
    end

    if isempty(methods)
        methods = {'Power'};
    end
end


function [bg, rsq, rmse_val, residuals, n_iter, selected_method, candidate_scores, candidate_details] = ...
        local_fit_background_dispatch(e_fit, s_fit, e_sub, requested_method, ...
        candidate_methods, asym_penalty, full_spec, sub_mask)
    if ~strcmpi(requested_method, 'Auto') || numel(candidate_methods) == 1
        selected_method = candidate_methods{1};
        [bg, rsq, rmse_val, residuals, n_iter] = fit_background( ...
            e_fit, s_fit, e_sub, selected_method, asym_penalty);
        bg = local_apply_background_safety_cap(full_spec(sub_mask), bg);
        candidate_details = local_make_candidate_detail( ...
            selected_method, full_spec(sub_mask), bg, rsq, rmse_val, residuals);
        candidate_scores = candidate_details.score;
        return
    end

    candidate_details = repmat(local_empty_candidate_detail(), 1, numel(candidate_methods));
    candidate_scores = inf(1, numel(candidate_methods));
    best = struct('score', inf, 'method', candidate_methods{1}, ...
        'bg', zeros(size(e_sub)), 'rsq', NaN, 'rmse', NaN, ...
        'residuals', [], 'n_iter', 1);

    for ci = 1:numel(candidate_methods)
        method = candidate_methods{ci};
        try
            [bg_i, rsq_i, rmse_i, residuals_i, n_iter_i] = fit_background( ...
                e_fit, s_fit, e_sub, method, asym_penalty);
            bg_i = local_apply_background_safety_cap(full_spec(sub_mask), bg_i);
            detail_i = local_make_candidate_detail( ...
                method, full_spec(sub_mask), bg_i, rsq_i, rmse_i, residuals_i);
            score_i = detail_i.score;
        catch
            bg_i = zeros(size(e_sub));
            rsq_i = NaN;
            rmse_i = NaN;
            residuals_i = [];
            n_iter_i = 1;
            detail_i = local_empty_candidate_detail(method);
            score_i = inf;
        end

        candidate_details(ci) = detail_i;
        candidate_scores(ci) = score_i;
        if score_i < best.score
            best.score = score_i;
            best.method = method;
            best.bg = bg_i;
            best.rsq = rsq_i;
            best.rmse = rmse_i;
            best.residuals = residuals_i;
            best.n_iter = n_iter_i;
        end
    end

    bg = best.bg;
    rsq = best.rsq;
    rmse_val = best.rmse;
    residuals = best.residuals;
    n_iter = best.n_iter;
    selected_method = best.method;
end


function detail = local_make_candidate_detail(method, raw_signal, bg, rsq, rmse_val, residuals)
    raw_signal = raw_signal(:);
    bg = real(bg(:));
    bg(~isfinite(bg)) = 0;
    bg = max(bg, 0);
    subtracted_signal = raw_signal - bg;
    linear_residuals = raw_signal - bg;
    linear_rmse = sqrt(mean(linear_residuals.^2, 'omitnan'));
    [neg_fraction, neg_area_fraction, bg_fraction] = ...
        local_background_quality_metrics(raw_signal, subtracted_signal, bg);

    detail = struct();
    detail.method = char(method);
    detail.rsquare = rsq;
    detail.rmse = rmse_val;
    detail.linear_rmse = linear_rmse;
    detail.neg_fraction = neg_fraction;
    detail.neg_area_fraction = neg_area_fraction;
    detail.bg_fraction = bg_fraction;
    detail.score = local_background_score(linear_rmse, method, ...
        neg_fraction, neg_area_fraction, bg_fraction, raw_signal);
    detail.residual_count = numel(residuals);
end


function detail = local_empty_candidate_detail(method)
    if nargin < 1
        method = '';
    end
    detail = struct('method', char(method), 'rsquare', NaN, 'rmse', NaN, ...
        'linear_rmse', NaN, 'neg_fraction', NaN, 'neg_area_fraction', NaN, ...
        'bg_fraction', NaN, 'score', inf, 'residual_count', 0);
end


function candidate_details = local_refresh_selected_candidate_details( ...
        candidate_details, selected_method, raw_signal, bg, rsq, rmse_val, residuals)
    selected_idx = find(strcmp({candidate_details.method}, char(selected_method)), 1);
    if isempty(selected_idx)
        selected_idx = numel(candidate_details) + 1;
    end
    candidate_details(selected_idx) = local_make_candidate_detail( ...
        selected_method, raw_signal, bg, rsq, rmse_val, residuals);
end


function bg = local_apply_background_safety_cap(raw_signal, bg)
    raw_signal = raw_signal(:);
    bg = real(bg(:));
    bg(~isfinite(bg)) = 0;
    bg = max(bg, 0);
    local_signal = max(raw_signal, 0);
    local_smooth = smoothdata(local_signal, 'gaussian', ...
        max(5, round(numel(local_signal) * 0.05)));
    bg = min(bg, local_smooth * 0.9);
end


function score = local_background_score(linear_rmse, method, neg_fraction, neg_area_fraction, bg_fraction, raw_signal)
    if ~isfinite(linear_rmse)
        score = inf;
        return
    end

    raw_signal = raw_signal(:);
    scale = max(median(abs(raw_signal), 'omitnan'), eps);
    rmse_term = linear_rmse / scale;
    complexity_term = local_background_model_complexity(method);
    score = rmse_term + 2.0 * neg_area_fraction + 0.5 * neg_fraction + ...
        0.2 * max(bg_fraction - 0.95, 0) + complexity_term;
end


function [neg_fraction, neg_area_fraction, bg_fraction] = ...
        local_background_quality_metrics(raw_signal, subtracted_signal, bg)
    raw_signal = raw_signal(:);
    subtracted_signal = subtracted_signal(:);
    bg = bg(:);

    valid = isfinite(raw_signal) & isfinite(subtracted_signal) & isfinite(bg);
    if ~any(valid)
        neg_fraction = NaN;
        neg_area_fraction = NaN;
        bg_fraction = NaN;
        return
    end

    raw_signal = raw_signal(valid);
    subtracted_signal = subtracted_signal(valid);
    bg = bg(valid);

    neg_mask = subtracted_signal < 0;
    neg_fraction = nnz(neg_mask) / max(numel(subtracted_signal), 1);
    neg_area = sum(-subtracted_signal(neg_mask), 'omitnan');
    total_area = sum(abs(subtracted_signal), 'omitnan') + eps;
    neg_area_fraction = neg_area / total_area;

    raw_area = sum(max(raw_signal, 0), 'omitnan') + eps;
    bg_fraction = sum(max(bg, 0), 'omitnan') / raw_area;
end


function penalty = local_background_model_complexity(method)
    switch char(method)
        case 'Power'
            penalty = 0.00;
        case 'Power2'
            penalty = 0.02;
        case 'Exp2'
            penalty = 0.03;
        case 'ExpPoly3'
            penalty = 0.04;
        case 'Pearson'
            penalty = 0.05;
        case 'PearsonVII'
            penalty = 0.06;
        otherwise
            penalty = 0.08;
    end
end


function [bg, rsq, rmse_val, residuals, n_iter] = fit_background( ...
        e_fit, s_fit, e_sub, method, asym_penalty, weights)
    % Fit background model using Curve Fitting Toolbox.
    %
    % Returns fitted background on e_sub, plus goodness-of-fit metrics.

    if nargin < 6
        weights = ones(size(e_fit));
    end
    n_iter = 1;
    rsq = NaN;
    rmse_val = NaN;
    residuals = [];

    e_fit = e_fit(:);
    s_fit = s_fit(:);
    weights = weights(:);

    switch method
        case 'Power'
            valid = e_fit > 0 & s_fit > 0;
            if nnz(valid) < 3
                bg = zeros(size(e_sub));
                return
            end
            [cfun, gof] = fit(e_fit(valid), s_fit(valid), 'power1', ...
                'Weight', weights(valid));
            bg = cfun(e_sub);
            residuals = s_fit(valid) - cfun(e_fit(valid));

        case 'Power2'
            valid = e_fit > 0 & s_fit > 0;
            if nnz(valid) < 4
                bg = zeros(size(e_sub));
                return
            end
            ft = fittype('a*x^b + c', 'independent', 'x');
            p0 = [max(s_fit), -1, min(s_fit) * 0.1];
            try
                [cfun, gof] = fit(e_fit(valid), s_fit(valid), ft, ...
                    'StartPoint', p0, 'Weight', weights(valid), ...
                    'Lower', [0 -Inf -Inf], 'Upper', [Inf 0 Inf]);
            catch
                % Fallback to Power1
                [cfun, gof] = fit(e_fit(valid), s_fit(valid), 'power1', ...
                    'Weight', weights(valid));
            end
            bg = cfun(e_sub);
            residuals = s_fit(valid) - cfun(e_fit(valid));

        case 'ExpPoly3'
            valid = s_fit > 0;
            if nnz(valid) < 5
                bg = zeros(size(e_sub));
                return
            end
            % Center and scale energy to reduce conditioning
            e_center = mean(e_fit(valid));
            e_scale = max(std(e_fit(valid)), eps);
            e_norm = (e_fit(valid) - e_center) / e_scale;
            e_sub_norm = (e_sub - e_center) / e_scale;
            [cfun, gof] = fit(e_norm, log(s_fit(valid)), 'poly3', ...
                'Weight', weights(valid));
            bg = exp(cfun(e_sub_norm));
            residuals = log(s_fit(valid)) - cfun(e_norm);

        case 'Pearson'
            % Log-log quadratic (original behavior)
            valid = e_fit > 0 & s_fit > 0;
            if nnz(valid) < 4
                bg = zeros(size(e_sub));
                return
            end
            % Center and scale log-energy to reduce conditioning
            le = log(e_fit(valid));
            le_center = mean(le);
            le_scale = max(std(le), eps);
            le_norm = (le - le_center) / le_scale;
            le_sub_norm = (log(e_sub) - le_center) / le_scale;
            [cfun, gof] = fit(le_norm, log(s_fit(valid)), 'poly2', ...
                'Weight', weights(valid));
            bg = exp(cfun(le_sub_norm));
            residuals = log(s_fit(valid)) - cfun(le_norm);

        case 'PearsonVII'
            % Pearson VII function: physical ZLP tail model (Ruishi Qi)
            % I * w^(2m) / (w^2 + (2^(1/m)-1)*(2*E-2*t0)^2)^m
            if nnz(s_fit > 0) < 4
                bg = zeros(size(e_sub));
                return
            end
            init = [max(s_fit), 10, 1.5, 0.01];
            ub = [max(s_fit)*10, 50, 5, 1];
            lb = [max(s_fit)/2, 2, 0.1, -1];
            opts_lsq = optimset('Display', 'off');
            try
                pbest = lsqcurvefit(@pearson_vii_diff, init, ...
                    [e_fit'; s_fit'], zeros(size(s_fit')), ...
                    lb, ub, opts_lsq);
                bg = pearson_vii(pbest, e_sub(:)');
                bg = bg(:);
                % Manual R² calculation
                fitted_vals = pearson_vii(pbest, e_fit(:)');
                ss_res = sum((s_fit - fitted_vals(:)).^2);
                ss_tot = sum((s_fit - mean(s_fit)).^2);
                gof = struct('rsquare', 1 - ss_res/max(ss_tot, eps), ...
                    'rmse', sqrt(ss_res / max(numel(s_fit)-4, 1)));
                residuals = s_fit - fitted_vals(:);
            catch
                bg = zeros(size(e_sub));
                gof = struct('rsquare', NaN, 'rmse', NaN);
                residuals = [];
            end

        case 'Exp2'
            % Two-term exponential: a*exp(b*x) + c*exp(d*x)
            if nnz(s_fit > 0) < 5
                bg = zeros(size(e_sub));
                return
            end
            try
                p0_log = polyfit(real(log(e_fit(e_fit > 0))), ...
                    real(log(s_fit(e_fit > 0))), 1);
                stt = [exp(p0_log(2))*0.8, p0_log(1)*0.001, ...
                       exp(p0_log(2))*0.2, p0_log(1)*0.0005];
                [cfun, gof] = fit(e_fit, s_fit, 'exp2', ...
                    'StartPoint', stt, 'Weight', weights);
            catch
                % Fallback: let MATLAB auto-start
                try
                    [cfun, gof] = fit(e_fit, s_fit, 'exp2', ...
                        'Weight', weights);
                catch
                    bg = zeros(size(e_sub));
                    return
                end
            end
            bg = cfun(e_sub);
            residuals = s_fit - cfun(e_fit);

        otherwise
            bg = zeros(size(e_sub));
            return
    end

    % Extract gof metrics
    if exist('gof', 'var') && isstruct(gof)
        rsq = gof.rsquare;
        if isfield(gof, 'rmse')
            rmse_val = gof.rmse;
        end
    end

    % --- Apply asymmetric penalty post-correction ---
    % If penalty > 1 and background overshoots data, nudge down
    if asym_penalty > 1
        bg = bg(:);
        s_sub = interp1(e_fit, s_fit, e_sub, 'linear', 'extrap');
        overshoot = bg > s_sub(:);
        if any(overshoot)
            excess = bg(overshoot) - s_sub(overshoot);
            bg(overshoot) = s_sub(overshoot) + excess / asym_penalty;
        end
    end
end


function ydata = pearson_vii(para, ene)
    % Pearson VII profile (Ruishi Qi implementation)
    I = para(1); w = para(2); m = para(3); t0 = para(4);
    ydata = I * (w).^(2*m) ./ ((w).^2 + (2^(1/m)-1) * (2*ene - 2*t0).^2).^m;
end


function diff_val = pearson_vii_diff(para, ene_spec)
    % Custom loss function for Pearson VII: energy-weighted, asymmetric penalty
    ene = ene_spec(1,:);
    spec = ene_spec(2,:);
    ydata = pearson_vii(para, ene);
    diff_val = (spec - ydata) .* (abs(ene)).^0.4;
    diff_val(diff_val < 0) = diff_val(diff_val < 0) * 4;
end


%% ═══════════ Deconvolution ═══════════

function qe_out = apply_deconvolution(qe_in, opts)
    % Lucy-Richardson deconvolution using the ZLP (q≈0 channel) as PSF.
    qe_out = qe_in;
    energy_axis = double(qe_in.energy_meV(:));
    intensity = double(qe_in.intensity);
    n_q = size(intensity, 2);
    n_iter = round(opts.deconv_iter);

    % --- Extract PSF from q≈0 channel ---
    if isfield(qe_in, 'q_zero_index') && isfinite(qe_in.q_zero_index)
        q0_idx = round(qe_in.q_zero_index);
    else
        [~, q0_idx] = max(max(intensity, [], 1));
    end
    q0_idx = max(1, min(q0_idx, n_q));
    q0_range = max(1, q0_idx - 1):min(n_q, q0_idx + 1);
    ref_spectrum = mean(intensity(:, q0_range), 2, 'omitnan');
    ref_spectrum = max(ref_spectrum, 0);

    % --- Use ONLY the ZLP peak region as PSF ---
    [~, zlp_idx] = max(ref_spectrum);
    zlp_energy = energy_axis(zlp_idx);
    zlp_half_width = 100;  % meV
    zlp_mask = energy_axis >= (zlp_energy - zlp_half_width) & ...
               energy_axis <= (zlp_energy + zlp_half_width);
    psf = ref_spectrum(zlp_mask);
    psf = max(psf, 0);
    psf_sum = sum(psf);
    if psf_sum > 0
        psf = psf / psf_sum;
    else
        return
    end

    % --- Deconvolve each q-channel ---
    for qi = 1:n_q
        spec = intensity(:, qi);
        spec = max(spec, 0);
        if max(spec) <= 0
            continue
        end
        try
            intensity(:, qi) = deconvlucy(spec, psf, n_iter);
        catch
            intensity(:, qi) = lr_deconv(spec, psf, n_iter);
        end
    end
    qe_out.intensity = intensity;
end


%% ═══════════ Lucy-Richardson (fallback) ═══════════

function result = lr_deconv(signal, psf, n_iter)
    % Manual Lucy-Richardson deconvolution (no toolbox dependency)
    signal = double(signal(:));
    psf = double(psf(:));
    estimate = signal;
    psf_flip = flipud(psf);
    for iter = 1:n_iter %#ok<FXUP>
        blurred = conv(estimate, psf, 'same');
        blurred(blurred < eps) = eps;
        ratio = signal ./ blurred;
        correction = conv(ratio, psf_flip, 'same');
        estimate = estimate .* correction;
    end
    result = max(estimate, 0);
end


%% ═══════════ Cosmic Ray Removal ═══════════

function qe_out = apply_cosmic_ray_removal(qe_in, ~)
    % Remove cosmic ray spikes from the intensity map.
    %
    % Uses a median-filter based approach: any pixel whose value exceeds
    % the local median by more than N standard deviations is replaced
    % with the median value. This is the standard "sigma-clipping" method.

    qe_out = qe_in;
    intensity = double(qe_in.intensity);
    threshold_sigma = 5;  % Spike detection threshold (σ)

    % 2D median filter as reference
    med_img = medfilt2(intensity, [3 3], 'symmetric');

    % Estimate noise from median absolute deviation
    residuals = intensity - med_img;
    noise_est = median(abs(residuals(:))) / 0.6745;
    if noise_est < eps, noise_est = 1; end

    % Flag spikes: |residual| > threshold * noise
    spike_mask = abs(residuals) > threshold_sigma * noise_est;
    n_spikes = nnz(spike_mask);

    if n_spikes > 0
        intensity(spike_mask) = med_img(spike_mask);
        fprintf('  Cosmic ray removal: %d spikes removed (%.2f%% of pixels)\n', ...
            n_spikes, 100 * n_spikes / numel(intensity));
    end
    qe_out.intensity = intensity;
end


%% ═══════════ Denoising ═══════════

function qe_out = apply_denoise(qe_in, opts)
    % Denoise the q-E intensity map.
    %   Wiener2D — 2D adaptive Wiener filter
    %   SavGol   — 1D Savitzky-Golay per q-channel
    %   BM3D     — Block-matching and 3D filtering (requires Nion toolbox)
    qe_out = qe_in;
    intensity = double(qe_in.intensity);
    method = char(opts.denoise_method);

    switch method
        case 'Wiener2D'
            sigma_input = opts.denoise_sigma;
            if sigma_input <= 0
                noise_est = median(abs(intensity(:) - median(intensity(:)))) / 0.6745;
            else
                noise_est = sigma_input;
            end
            try
                denoised = wiener2(intensity, [3 5], noise_est^2);
            catch
                kernel = ones(3, 5) / 15;
                denoised = conv2(intensity, kernel, 'same');
            end
            qe_out.intensity = denoised;

        case 'SavGol'
            sg_order = round(opts.sg_order);
            sg_framelen = round(opts.sg_framelen);
            if mod(sg_framelen, 2) == 0
                sg_framelen = sg_framelen + 1;
            end
            sg_order = min(sg_order, sg_framelen - 1);
            n_q = size(intensity, 2);
            for qi = 1:n_q
                spec = intensity(:, qi);
                if numel(spec) >= sg_framelen
                    try
                        intensity(:, qi) = sgolayfilt(spec, sg_order, sg_framelen);
                    catch
                        intensity(:, qi) = movmean(spec, sg_framelen);
                    end
                end
            end
            qe_out.intensity = intensity;

        case 'BM3D'
            % BM3D denoising using BM3D_QRS from Nion toolbox.
            %   opts.denoise_sigma controls the noise overestimation factor:
            %     0  = auto (default, BMfactor=1)
            %     >0 = used as BMfactor (higher = stronger denoising)
            if isfield(opts, 'denoise_sigma') && opts.denoise_sigma > 0
                bmfactor = opts.denoise_sigma;
            else
                bmfactor = 1;
            end
            try
                denoised = BM3D_QRS(intensity, bmfactor);
                qe_out.intensity = denoised;
            catch e
                warning('qe_preprocess:bm3d', ...
                    'BM3D_QRS not available (%s). Falling back to Wiener2D.', e.message);
                noise_est = median(abs(intensity(:) - median(intensity(:)))) / 0.6745;
                try
                    qe_out.intensity = wiener2(intensity, [3 5], noise_est^2);
                catch
                    kernel = ones(3, 5) / 15;
                    qe_out.intensity = conv2(intensity, kernel, 'same');
                end
            end
    end
end
