%% analyze_dispersion.m �?Automated plasmon dispersion analysis
%  Compares 10w and 20w defocus datasets with full preprocessing pipeline.
%
%  Pipeline:  load �?crop q/E �?denoise �?area-normalize �?BG subtract
%             �?deconvolve �?auto-fit Drude-Lorentz �?branch separation
%             �?da Jornada quasi-2D plasmon model fitting
%
%  Output:    analysis_results.mat, bi_dispersion_comparison.png
%
%  Usage:     Run from the project root directory.

clearvars; close all; clc;

src_dir = fileparts(mfilename('fullpath'));
if isempty(src_dir); src_dir = pwd; end
project_dir = fileparts(src_dir);  % go up from src/ to project root
cd(project_dir);

%% ══════════════════════════════════════════════════════════════�?%  Configuration
%  ══════════════════════════════════════════════════════════════�?cfg = struct();
cfg.q_range   = [-0.15, 0.15];    % Å⁻�?cfg.E_range   = [300, 3836];      % meV
cfg.zlp_window = [-50, 50];       % ZLP peak extraction window for normalization
cfg.E_min_fit = 300;              % fit window lower bound
cfg.E_max_fit = 3836;             % fit window upper bound
cfg.max_peaks = 2;                % Branch 1 = plasmon, Branch 2 = interband
cfg.min_prominence = 0.15;
cfg.smooth_width   = 25;
cfg.q_skip    = 0.005;            % skip |q| < this (ZLP dominated)
cfg.R2_threshold = 0.3;           % quality filter for auto-fit
cfg.max_gamma_ratio = 2.0;        % reject overdamped peaks: Γ/ω�?> this
cfg.denoise_method = 'Wiener2D';
cfg.deconv_iter    = 5;
cfg.smooth_gauss_w = 3;           % gaussian smoothing width for display

sessions = struct( ...
    'name',  {'10w defocus', '20w defocus'}, ...
    'path',  {fullfile(project_dir, '20260120 Bi', '590 PL2 10w 0.004 10sx300'), ...
              fullfile(project_dir, '20260120 Bi', 'no pl2 20w 0.004 10sx300 2film')}, ...
    'dq',    {0.005, 0.0025} ...
);

results = cell(numel(sessions), 1);

%% ══════════════════════════════════════════════════════════════�?%  Main analysis loop
%  ══════════════════════════════════════════════════════════════�?for s = 1:numel(sessions)
    fprintf('\n══════════════════════════════════════════════\n');
    fprintf('  Session %d: %s\n', s, sessions(s).name);
    fprintf('  Path: %s\n', sessions(s).path);
    fprintf('══════════════════════════════════════════════\n');

    %% 1. Load data
    dataset = load_qe_dataset(sessions(s).path, sessions(s).dq);
    qe = dataset.qe;
    fprintf('  Loaded: %s\n', dataset.label);
    fprintf('  Size: %d energy × %d q channels\n', size(qe.intensity));
    fprintf('  dq = %.4f Å⁻�?pixel\n', dataset.dq_Ainv);

    %% 2. Crop to q and E range
    q_axis = qe.q_Ainv;
    energy_meV = qe.energy_meV;
    q_mask = q_axis >= cfg.q_range(1) & q_axis <= cfg.q_range(2);
    E_mask = energy_meV >= cfg.E_range(1) & energy_meV <= cfg.E_range(2);

    intensity = double(qe.intensity(E_mask, q_mask));
    energy_axis = energy_meV(E_mask);
    q_cropped = q_axis(q_mask);
    n_q = numel(q_cropped);
    fprintf('  After crop: %d energy × %d q channels\n', size(intensity));

    % Also keep the full-energy intensity for deconvolution PSF extraction
    intensity_full_E = double(qe.intensity(:, q_mask));
    energy_full = energy_meV;

    %% 3. Denoise (Wiener2D on full-E data, then crop)
    fprintf('  Denoising (Wiener2D)...\n');
    noise_est = median(abs(intensity_full_E(:) - median(intensity_full_E(:)))) / 0.6745;
    try
        denoised_full = wiener2(intensity_full_E, [3 5], noise_est^2);
    catch
        kernel = ones(3, 5) / 15;
        denoised_full = conv2(intensity_full_E, kernel, 'same');
    end
    intensity = denoised_full(E_mask, :);
    intensity_full_E = denoised_full;

    %% 4. Deconvolution (Lucy-Richardson) �?BEFORE normalization
    %  PSF must be extracted from raw denoised data to represent the true
    %  instrument response, not a sample-normalized artifact.
    fprintf('  Deconvolution (Lucy-Richardson, %d iter)...\n', cfg.deconv_iter);

    % Extract PSF from q�? channel (full energy range, RAW denoised)
    [~, q0_idx] = min(abs(q_cropped));
    q0_range = max(1, q0_idx-1):min(n_q, q0_idx+1);
    ref_spectrum = mean(intensity_full_E(:, q0_range), 2, 'omitnan');
    ref_spectrum = max(ref_spectrum, 0);

    % Use only the ZLP peak region as PSF
    [~, zlp_idx] = max(ref_spectrum);
    zlp_energy = energy_full(zlp_idx);
    zlp_half_width = 100;  % meV
    zlp_mask_psf = energy_full >= (zlp_energy - zlp_half_width) & ...
                   energy_full <= (zlp_energy + zlp_half_width);
    psf = ref_spectrum(zlp_mask_psf);
    psf = max(psf, 0);
    psf_sum = sum(psf);

    if psf_sum > 0
        psf = psf / psf_sum;
        for qi = 1:n_q
            spec = intensity(:, qi);
            spec = max(spec, 0);
            if max(spec) <= 0; continue; end
            try
                intensity(:, qi) = deconvlucy(spec, psf, cfg.deconv_iter);
            catch
                % Manual Lucy-Richardson fallback
                estimate = spec;
                psf_flip = flipud(psf);
                for iter = 1:cfg.deconv_iter
                    blurred = conv(estimate, psf, 'same');
                    blurred(blurred < eps) = eps;
                    ratio = spec ./ blurred;
                    correction = conv(ratio, psf_flip, 'same');
                    estimate = estimate .* correction;
                end
                intensity(:, qi) = max(estimate, 0);
            end
        end
    end

    %% 5. ZLP peak normalization
    %  Divide by ZLP peak height to normalize for beam current variations
    %  while preserving q-dependent spectral weight.
    %  ZLP peak �?elastic scattering cross-section × beam current,
    %  which is a q-independent systematic factor.
    %  Unlike area normalization, this does NOT destroy A(q) scaling.
    fprintf('  ZLP peak normalization [%.0f, %.0f] meV window...\n', ...
        cfg.zlp_window(1), cfg.zlp_window(2));
    zlp_norm_mask = energy_full >= cfg.zlp_window(1) & ...
                    energy_full <= cfg.zlp_window(2);
    for qi = 1:n_q
        zlp_peak = max(intensity_full_E(zlp_norm_mask, qi));
        if isfinite(zlp_peak) && zlp_peak > eps
            intensity(:, qi) = intensity(:, qi) ./ zlp_peak;
        end
    end

    %% 6. Background subtraction (single low-energy window)
    %  Power-law fit in [50, 300] meV only �?this region is physically
    %  clean (below Bi plasmon/interband onset). The high-energy window
    %  was removed to avoid contamination from loss features.
    fprintf('  Background subtraction (single-window [50, 300] meV)...\n');
    bg_fit_lo = 50;
    bg_fit_hi = 300;
    bg_fit_mask = energy_axis >= bg_fit_lo & energy_axis <= bg_fit_hi;

    if nnz(bg_fit_mask) >= 5
        e_fit = energy_axis(bg_fit_mask);
        sub_mask = energy_axis > bg_fit_lo;
        e_sub = energy_axis(sub_mask);

        for qi = 1:n_q
            spec = intensity(:, qi);
            s_fit = spec(bg_fit_mask);
            s_fit_pos = max(s_fit, eps);

            try
                valid = e_fit > 0 & s_fit_pos > 0;
                if nnz(valid) < 3; continue; end
                p = polyfit(log(e_fit(valid)), log(s_fit_pos(valid)), 1);
                bg = exp(polyval(p, log(e_sub)));
                bg = real(bg(:));
                bg(~isfinite(bg)) = 0;
                bg = max(bg, 0);

                % Safety cap: prevent over-subtraction
                local_signal = max(spec(sub_mask), 0);
                local_smooth = smoothdata(local_signal, 'gaussian', ...
                    max(5, round(numel(local_signal) * 0.05)));
                bg = min(bg, local_smooth * 0.9);

                intensity(sub_mask, qi) = spec(sub_mask) - bg;
            catch
            end
        end
    end

    %% 7. Automated multi-peak Drude-Lorentz fitting per q-channel
    fprintf('  Auto-fitting dispersion (%d channels)...\n', n_q);
    all_peaks = [];   % [q, omega_p, gamma, R², amplitude]
    fit_details = cell(n_q, 1);
    n_success = 0;

    for k = 1:n_q
        q_val = q_cropped(k);

        % Skip channels near q=0
        if abs(q_val) < cfg.q_skip
            continue
        end

        spectrum = intensity(:, k);
        if all(spectrum == 0) || all(isnan(spectrum))
            continue
        end

        try
            result = fit_loss_function(energy_axis, spectrum, ...
                'E_min', cfg.E_min_fit, 'E_max', cfg.E_max_fit, ...
                'min_prominence', cfg.min_prominence, ...
                'smooth_width', cfg.smooth_width, ...
                'max_peaks', cfg.max_peaks, ...
                'pre_subtracted', true);
            fit_details{k} = result;
            n_success = n_success + 1;

            for p = 1:result.n_peaks
                if result.omega_p(p) < cfg.E_min_fit
                    continue
                end
                % Reject overdamped peaks (Γ/ω�?> threshold)
                if result.gamma(p) / result.omega_p(p) > cfg.max_gamma_ratio
                    continue
                end
                all_peaks(end+1, :) = [ ...
                    q_val, ...
                    result.omega_p(p), ...
                    result.gamma(p), ...
                    result.R_squared, ...
                    result.amplitude(p)]; %#ok<AGROW>
            end
        catch
        end

        if mod(k, 20) == 0
            fprintf('    Progress: %d/%d (%.0f%%)\n', k, n_q, k/n_q*100);
        end
    end

    fprintf('  Fitted %d/%d channels, found %d peaks\n', n_success, n_q, size(all_peaks, 1));

    %% 8. Filter by R² and separate into branches
    if isempty(all_peaks) || size(all_peaks, 1) < 2
        warning('Session %d: insufficient peaks found, skipping.', s);
        results{s} = struct('name', sessions(s).name, 'success', false);
        continue
    end

    good = all_peaks(:, 4) > cfg.R2_threshold;
    peaks = all_peaks(good, :);
    fprintf('  After R² filter (>%.2f): %d peaks\n', cfg.R2_threshold, size(peaks, 1));
    fprintf('  After Γ/ω�?filter (�?.1f): %d peaks (applied upstream)\n', ...
        cfg.max_gamma_ratio, size(peaks, 1));

    % Separate into branches using k-means clustering on peak energy.
    % Do this BEFORE ±q averaging to preserve the full energy distribution.
    n_clusters = cfg.max_peaks;
    energies = peaks(:, 2);

    % Manual k-means (no Statistics Toolbox dependency)
    e_min_k = min(energies);
    e_max_k = max(energies);
    centroids = linspace(e_min_k, e_max_k, n_clusters)';

    for km_iter = 1:50
        dists = abs(energies - centroids');
        [~, idx] = min(dists, [], 2);
        new_centroids = centroids;
        for kk = 1:n_clusters
            members = energies(idx == kk);
            if ~isempty(members)
                new_centroids(kk) = mean(members);
            end
        end
        if max(abs(new_centroids - centroids)) < 1
            centroids = new_centroids;
            break
        end
        centroids = new_centroids;
    end

    % Final assignment
    dists = abs(energies - centroids');
    [~, idx] = min(dists, [], 2);

    % Sort centroids by energy (branch 1 = lowest energy = plasmon)
    [~, sort_order] = sort(centroids);
    reorder = zeros(n_clusters, 1);
    for kk = 1:n_clusters
        reorder(sort_order(kk)) = kk;
    end
    idx = reorder(idx);

    % Build branch arrays, then ±q symmetrize within each branch
    branches = cell(n_clusters, 1);
    for kk = 1:n_clusters
        br_raw = peaks(idx == kk, :);
        if ~isempty(br_raw)
            branches{kk} = local_symmetrize_peaks(br_raw);
        end
    end
    branches = branches(~cellfun('isempty', branches));
    n_branches = numel(branches);
    fprintf('  Separated into %d branches (±q averaged per branch)\n', n_branches);

    %% 9. Fit da Jornada quasi-2D plasmon model per branch
    %  Use R²-weighted fitting: w = R² × min(1, ω�?Γ)
    fit_results_branch = cell(n_branches, 1);
    for b = 1:n_branches
        br = branches{b};
        fprintf('  Branch %d: %d pts, E range [%.0f, %.0f] meV\n', ...
            b, size(br, 1), min(br(:, 2)), max(br(:, 2)));

        if size(br, 1) >= 5
            try
                % Compute per-point weights: high R² + underdamped �?high weight
                w_R2 = br(:, 4);  % R² column
                w_damping = min(1, br(:, 2) ./ max(br(:, 3), 1));  % ω�?Γ capped at 1
                weights = w_R2 .* w_damping;
                weights = weights / max(weights);  % normalize to [0,1]

                disp_result = fit_quasi2d_plasmon(br(:, 1), br(:, 2), ...
                    'confidence', weights);
                fit_results_branch{b} = disp_result;
                fprintf('    �?ρ₀ = %.1f Å, E_flat = %.0f meV, R² = %.4f\n', ...
                    disp_result.rho0, disp_result.E_flat_meV, disp_result.R_squared);
            catch ME
                fprintf('    �?Fit failed: %s\n', ME.message);
            end
        else
            fprintf('    �?Too few points for dispersion fit\n');
        end
    end

    %% Store results
    res = struct();
    res.name = sessions(s).name;
    res.path = sessions(s).path;
    res.dq = sessions(s).dq;
    res.success = true;
    res.cfg = cfg;
    res.q_axis = q_cropped;
    res.energy_axis = energy_axis;
    res.intensity = intensity;
    res.all_peaks = all_peaks;
    res.peaks_filtered = peaks;
    res.branches = branches;
    res.n_branches = n_branches;
    res.fit_results = fit_results_branch;
    res.fit_details = fit_details;
    res.n_success = n_success;
    results{s} = res;
end

%% ══════════════════════════════════════════════════════════════�?%  Save results
%  ══════════════════════════════════════════════════════════════�?save_path = fullfile(project_dir, 'figures', 'analysis_results.mat');
save(save_path, 'results', 'cfg', 'sessions');
fprintf('\n�?Results saved to %s\n', save_path);

%% ══════════════════════════════════════════════════════════════�?%  Comparative figure
%  ══════════════════════════════════════════════════════════════�?fprintf('\n  Generating comparison figure...\n');
fig = figure('Name', 'Bi Plasmon Dispersion: 10w vs 20w', ...
    'Color', 'w', 'Position', [50 50 1600 1000], ...
    'NumberTitle', 'off');

branch_colors = [0.85 0.15 0.15; 0.15 0.35 0.85; 0.15 0.7 0.3; 0.6 0.2 0.8];

for s = 1:numel(sessions)
    res = results{s};
    if ~res.success; continue; end

    %% Row 1: q-E map with overlaid dispersion points
    ax1 = subplot(3, 2, s);

    % Apply the same asinh-stretch display as the GUI for contrast
    raw_map = res.intensity;
    finite_vals = raw_map(isfinite(raw_map));
    lower_lim = prctile(finite_vals, 0.5);
    upper_lim = prctile(finite_vals, 99.8);
    if ~isfinite(lower_lim) || ~isfinite(upper_lim) || upper_lim <= lower_lim
        lower_lim = min(finite_vals);
        upper_lim = max(finite_vals);
    end
    scale = max(upper_lim - lower_lim, eps);
    clipped = max(raw_map - lower_lim, 0);
    stretch = 8;
    display_map = asinh(stretch * clipped / scale) / asinh(stretch);

    imagesc(ax1, res.q_axis, res.energy_axis, display_map);
    set(ax1, 'YDir', 'normal');
    colormap(ax1, turbo);
    clim(ax1, [0 1]);
    hold(ax1, 'on');

    % Overlay branch points
    for b = 1:res.n_branches
        br = res.branches{b};
        col = branch_colors(min(b, size(branch_colors, 1)), :);
        marker = 'o';
        if b == 1
            % Filled for plasmon branch
            scatter(ax1, br(:, 1), br(:, 2), 20, col, 'filled', ...
                'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'w', 'LineWidth', 0.3);
        else
            % Open for interband
            scatter(ax1, br(:, 1), br(:, 2), 20, col, ...
                'MarkerFaceAlpha', 0.5, 'LineWidth', 1);
        end

        % Overlay fit curve
        if ~isempty(res.fit_results{b})
            fr = res.fit_results{b};
            plot(ax1, fr.q_fit, fr.E_fit, '-', 'Color', col, 'LineWidth', 2);
        end
    end
    hold(ax1, 'off');
    xlabel(ax1, 'q (Å^{-1})');
    ylabel(ax1, 'Energy (meV)');
    title(ax1, sprintf('q-E Map | %s', res.name), 'FontSize', 11);
    xlim(ax1, cfg.q_range);
    ylim(ax1, cfg.E_range);

    %% Row 2: Dispersion E(q) with fits
    ax2 = subplot(3, 2, 2 + s);
    hold(ax2, 'on');
    legend_entries = {};

    for b = 1:res.n_branches
        br = res.branches{b};
        col = branch_colors(min(b, size(branch_colors, 1)), :);

        scatter(ax2, br(:, 1), br(:, 2), 25, col, 'filled', ...
            'MarkerFaceAlpha', 0.6);

        if ~isempty(res.fit_results{b})
            fr = res.fit_results{b};
            plot(ax2, fr.q_fit, fr.E_fit, '-', 'Color', col, 'LineWidth', 2);
            legend_entries{end+1} = sprintf('B%d: \\rho_0=%.1fÅ  E_f=%.0fmeV', ...
                b, fr.rho0, fr.E_flat_meV); %#ok<AGROW>
        else
            legend_entries{end+1} = sprintf('B%d (no fit)', b); %#ok<AGROW>
        end
    end

    hold(ax2, 'off');
    xlabel(ax2, 'q (Å^{-1})');
    ylabel(ax2, 'Peak energy (meV)');
    title(ax2, sprintf('Dispersion | %s', res.name), 'FontSize', 11);
    legend(ax2, legend_entries, 'Location', 'best', 'FontSize', 8);
    grid(ax2, 'on');
    xlim(ax2, cfg.q_range);

    %% Row 3: Linewidth Γ(q)
    ax3 = subplot(3, 2, 4 + s);
    hold(ax3, 'on');

    for b = 1:res.n_branches
        br = res.branches{b};
        col = branch_colors(min(b, size(branch_colors, 1)), :);
        scatter(ax3, br(:, 1), br(:, 3), 25, col, 'filled', ...
            'MarkerFaceAlpha', 0.6, ...
            'DisplayName', sprintf('B%d Γ(q)', b));
    end

    hold(ax3, 'off');
    xlabel(ax3, 'q (Å^{-1})');
    ylabel(ax3, 'Linewidth Γ (meV)');
    title(ax3, sprintf('Γ(q) | %s', res.name), 'FontSize', 11);
    legend(ax3, 'Location', 'best', 'FontSize', 8);
    grid(ax3, 'on');
    xlim(ax3, cfg.q_range);
end

sgtitle(fig, 'Bi Thin Film Plasmon Dispersion �?10w vs 20w Defocus', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% Save figure
fig_path = fullfile(project_dir, 'figures', 'bi_dispersion_comparison.png');
exportgraphics(fig, fig_path, 'Resolution', 300);
fprintf('�?Figure saved to %s\n', fig_path);

%% ══════════════════════════════════════════════════════════════�?%  Summary table
%  ══════════════════════════════════════════════════════════════�?fprintf('\n════════════════════════════════════════════════════════\n');
fprintf('  SUMMARY: Bi Plasmon Dispersion Comparison\n');
fprintf('════════════════════════════════════════════════════════\n');
fprintf('%-16s | %-8s | %-8s | %-10s | %-10s | %-8s\n', ...
    'Session', 'Branch', 'N pts', 'ρ₀ (Å)', 'E_flat(meV)', 'R²');
fprintf('%-16s-+-%-8s-+-%-8s-+-%-10s-+-%-10s-+-%-8s\n', ...
    repmat('-',1,16), repmat('-',1,8), repmat('-',1,8), ...
    repmat('-',1,10), repmat('-',1,10), repmat('-',1,8));

for s = 1:numel(sessions)
    res = results{s};
    if ~res.success; continue; end
    for b = 1:res.n_branches
        br = res.branches{b};
        if ~isempty(res.fit_results{b})
            fr = res.fit_results{b};
            fprintf('%-16s | B%-7d | %-8d | %-10.1f | %-10.0f | %-8.4f\n', ...
                res.name, b, size(br, 1), fr.rho0, fr.E_flat_meV, fr.R_squared);
        else
            fprintf('%-16s | B%-7d | %-8d | %-10s | %-10s | %-8s\n', ...
                res.name, b, size(br, 1), 'N/A', 'N/A', 'N/A');
        end
    end
end
fprintf('════════════════════════════════════════════════════════\n');

fprintf('\n�?Analysis complete.\n');


%% ══════════════════════════════════════════════════════════════�?%  Helper: ±q symmetric averaging
%  ══════════════════════════════════════════════════════════════�?function peaks_sym = local_symmetrize_peaks(peaks)
%LOCAL_SYMMETRIZE_PEAKS  Average +q and -q peaks at matched |q|.
%   Exploits E(q)=E(-q) symmetry to reduce noise.
%   Input:  peaks = [q, omega_p, gamma, R², amplitude]
%   Output: peaks_sym = same format but using |q| and averaged values

    if isempty(peaks)
        peaks_sym = peaks;
        return
    end

    abs_q = abs(peaks(:, 1));
    unique_q = unique(round(abs_q, 5));  % group by |q| rounded to 5 decimals

    peaks_sym = zeros(0, size(peaks, 2));
    for i = 1:numel(unique_q)
        mask = abs(abs_q - unique_q(i)) < 1e-6;
        grp = peaks(mask, :);

        % Average all columns except q (use |q|)
        avg_row = mean(grp, 1);
        avg_row(1) = unique_q(i);  % use positive |q|

        peaks_sym(end+1, :) = avg_row; %#ok<AGROW>
    end
end
