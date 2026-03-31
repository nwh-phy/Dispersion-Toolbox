function qe_out = qe_preprocess(qe_in, opts)
%QE_PREPROCESS  Apply EELS preprocessing pipeline to a qeData struct.
%
%   qe_out = qe_preprocess(qe_in, opts)
%
%   Applies a configurable preprocessing chain:
%     1. Normalization  (ZLP Peak or Area)
%     2. Denoising      (Wiener2D or Savitzky-Golay)
%     3. Background subtraction (Power / ExpPoly3 / Pearson)
%     4. Deconvolution  (Lucy-Richardson)
%
%   opts fields:
%     .do_normalize   logical  — enable normalization step
%     .norm_method    char     — 'ZLP Peak' | 'Area'
%     .norm_min       double   — window low bound (meV)
%     .norm_max       double   — window high bound (meV)
%     .do_denoise     logical  — enable denoising step
%     .denoise_method char     — 'Wiener2D' | 'SavGol'
%     .denoise_sigma  double   — noise sigma (0 = auto-estimate)
%     .sg_order       double   — Savitzky-Golay polynomial order
%     .sg_framelen    double   — Savitzky-Golay frame length (must be odd)
%     .do_bg_sub      logical  — enable background subtraction
%     .bg_method      char     — 'Power' | 'ExpPoly3' | 'Pearson'
%     .do_deconv      logical  — enable deconvolution
%     .deconv_iter    double   — Lucy-Richardson iterations
%
%   See also: interactive_qe_browser, fit_loss_function

arguments
    qe_in   struct
    opts    struct
end

qe_out = qe_in;

if isfield(opts, 'do_normalize') && opts.do_normalize
    qe_out = apply_normalization(qe_out, opts);
end

if isfield(opts, 'do_denoise') && opts.do_denoise
    qe_out = apply_denoise(qe_out, opts);
end

if isfield(opts, 'do_bg_sub') && opts.do_bg_sub
    qe_out = apply_bg_subtraction(qe_out, opts);
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
        zlp_mask = energy_axis >= norm_min & energy_axis <= norm_max;
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


%% ═══════════ Background Subtraction ═══════════

function qe_out = apply_bg_subtraction(qe_in, opts)
    % EELS background subtraction with single low-energy window fitting.
    %
    % Uses ONLY the low-energy ZLP tail region [50, win1_hi] as anchor.
    % The high-energy window was removed because it can contain
    % interband transitions or plasmon tails, contaminating the
    % background estimate and causing over-subtraction.
    %
    % Safety cap: background is limited to at most 90% of the
    % smoothed local signal — prevents complete erasure at high |q|.
    qe_out = qe_in;
    energy_axis = double(qe_in.energy_meV(:));
    intensity = double(qe_in.intensity);
    method = char(opts.bg_method);

    % --- Define fit window ---
    % Low-energy window only: [fit_lo, win1_hi]
    % Between ZLP tail and the onset of loss features.
    fit_lo = 50;
    win1_hi = 300;  % conservative upper bound below typical Bi loss onset

    fit_mask = energy_axis >= fit_lo & energy_axis <= win1_hi;

    if nnz(fit_mask) < 5
        return
    end

    e_fit = energy_axis(fit_mask);

    % Only subtract for energies above fit_lo
    sub_mask = energy_axis > fit_lo;
    e_sub = energy_axis(sub_mask);

    for qi = 1:size(intensity, 2)
        spec = intensity(:, qi);
        s_fit = spec(fit_mask);
        s_fit_pos = max(s_fit, eps);

        try
            switch method
                case 'Power'
                    valid = e_fit > 0 & s_fit_pos > 0;
                    if nnz(valid) < 3, continue; end
                    p = polyfit(log(e_fit(valid)), log(s_fit_pos(valid)), 1);
                    bg = exp(polyval(p, log(e_sub)));

                case 'ExpPoly3'
                    valid = s_fit_pos > 0;
                    if nnz(valid) < 5, continue; end
                    p = polyfit(e_fit(valid), log(s_fit_pos(valid)), 3);
                    bg = exp(polyval(p, e_sub));

                case 'Pearson'
                    valid = e_fit > 0 & s_fit_pos > 0;
                    if nnz(valid) < 4, continue; end
                    p = polyfit(log(e_fit(valid)), log(s_fit_pos(valid)), 2);
                    bg = exp(polyval(p, log(e_sub)));

                otherwise
                    continue
            end

            bg = real(bg(:));
            bg(~isfinite(bg)) = 0;
            bg = max(bg, 0);

            % --- Safety cap: prevent over-subtraction ---
            local_signal = max(spec(sub_mask), 0);
            local_smooth = smoothdata(local_signal, 'gaussian', ...
                max(5, round(numel(local_signal) * 0.05)));
            bg = min(bg, local_smooth * 0.9);

            intensity(sub_mask, qi) = spec(sub_mask) - bg;
        catch
        end
    end
    qe_out.intensity = intensity;
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


%% ═══════════ Denoising ═══════════

function qe_out = apply_denoise(qe_in, opts)
    % Denoise the q-E intensity map.
    %   Wiener2D — 2D adaptive Wiener filter
    %   SavGol   — 1D Savitzky-Golay per q-channel
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
    end
end
