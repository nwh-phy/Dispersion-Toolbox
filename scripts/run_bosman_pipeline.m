%% run_bosman_pipeline.m — Full Do et al. (2025) pipeline on Bi data
%  Phase 1: Preprocess (denoise + ZLP norm + quasi-elastic BG removal)
%  Phase 2: Drude-Lorentz fitting → extract ω_p(q), Γ(q), A(q)
%  Phase 3: Self-consistent I_kin / ε_inter decoupling
%  Phase 4: Publication figures (dispersion, I_kin log-log, loss map)
%
%  Outputs:  figures saved to output/bosman_pipeline/

clearvars; close all; clc;

%% Setup paths
project_dir = 'c:\Users\HP\Desktop\vibecoding\2D-metal-plasmon';
cd(project_dir);
run('startup.m');

output_dir = fullfile(project_dir, 'output', 'bosman_pipeline');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% Select dataset — use the 10w defocus (better signal)
data_path = fullfile(project_dir, '20260120 Bi', '590 PL2 10w 0.004 10sx300');
dq_Ainv = 0.005;

fprintf('\n══════════════════════════════════════════════\n');
fprintf('  Bosman/Do et al. Pipeline — Bi thin film\n');
fprintf('  Data: %s\n', data_path);
fprintf('══════════════════════════════════════════════\n\n');

%% ═══════════════════════════════════════════════════════════════
%  PHASE 1: Data Preprocessing
%  ═══════════════════════════════════════════════════════════════
fprintf('── PHASE 1: Preprocessing ──\n');

% 1.1 Load raw data
dataset = load_qe_dataset(data_path, dq_Ainv);
qe = dataset.qe;
fprintf('  Loaded: %d energy × %d q channels\n', size(qe.intensity));
fprintf('  Energy range: [%.0f, %.0f] meV\n', qe.energy_meV(1), qe.energy_meV(end));
fprintf('  q range: [%.4f, %.4f] Å⁻¹\n', qe.q_Ainv(1), qe.q_Ainv(end));

% 1.2 Configure preprocessing (Bosman standard)
pp = struct();
pp.do_despike    = false;
pp.do_normalize  = true;
pp.norm_method   = 'ZLP Peak';   % beam current correction (preserves A(q))
pp.norm_min      = -50;          % fixed ZLP window
pp.norm_max      = 50;
pp.do_denoise    = true;
pp.denoise_method = 'Wiener2D';  % fast and safe
pp.denoise_sigma = 0;            % auto-estimate
pp.do_bg_sub     = true;
pp.bg_method     = 'Power';      % power-law ZLP tail model
pp.bg_win_lo     = [50, 300];    % primary fit window (pre-plasmon)
pp.bg_win_hi     = [];           % single window (no dual)
pp.bg_iterative  = false;
pp.do_deconv     = false;        % skip deconv for now

% 1.3 Run preprocessing
fprintf('  Running: ZLP norm → Wiener2D → BG removal...\n');
[qe_pp, bg_diag] = qe_preprocess(qe, pp);
fprintf('  Done. BG diagnostics: %d channels\n', numel(bg_diag));

% Show median R² of BG fits
if ~isempty(bg_diag)
    rsq_vals = [bg_diag.rsquare];
    rsq_vals = rsq_vals(isfinite(rsq_vals));
    fprintf('  BG fit quality: median R² = %.4f\n', median(rsq_vals));
end

%% ═══════════════════════════════════════════════════════════════
%  PHASE 2: Single-spectrum Drude-Lorentz fitting
%  ═══════════════════════════════════════════════════════════════
fprintf('\n── PHASE 2: Dispersion fitting ──\n');

% Crop to analysis region
E_range = [300, 3836];
q_range = [-0.15, 0.15];
E_mask = qe_pp.energy_meV >= E_range(1) & qe_pp.energy_meV <= E_range(2);
q_mask = qe_pp.q_Ainv >= q_range(1) & qe_pp.q_Ainv <= q_range(2);

energy_axis = qe_pp.energy_meV(E_mask);
q_axis = qe_pp.q_Ainv(q_mask);
intensity = double(qe_pp.intensity(E_mask, q_mask));
n_q = numel(q_axis);

fprintf('  Analysis region: E=[%.0f, %.0f] meV, q=[%.3f, %.3f] Å⁻¹\n', ...
    E_range(1), E_range(2), q_range(1), q_range(2));
fprintf('  %d energy × %d q channels\n', numel(energy_axis), n_q);

% Auto-fit each q-channel
fit_cfg = struct();
fit_cfg.E_min = E_range(1);
fit_cfg.E_max = E_range(2);
fit_cfg.max_peaks = 2;
fit_cfg.min_prominence = 0.0005;  % ZLP-norm'd intensities are O(10⁻³)
fit_cfg.smooth_width = 25;
fit_cfg.q_skip = 0.005;
fit_cfg.R2_threshold = 0.20;
fit_cfg.max_gamma_ratio = 3.0;

all_peaks = [];
n_success = 0;

fprintf('  Fitting %d channels...\n', n_q);
for k = 1:n_q
    q_val = q_axis(k);
    if abs(q_val) < fit_cfg.q_skip, continue; end
    
    spec = intensity(:, k);
    if max(spec) <= 0, continue; end
    
    try
        result = fit_loss_function(energy_axis, spec, ...
            'E_min', fit_cfg.E_min, 'E_max', fit_cfg.E_max, ...
            'max_peaks', fit_cfg.max_peaks, ...
            'min_prominence', fit_cfg.min_prominence, ...
            'smooth_width', fit_cfg.smooth_width, ...
            'pre_subtracted', true);
        
        if result.R_squared >= fit_cfg.R2_threshold
            for p = 1:result.n_peaks
                omega_p = result.omega_p(p);
                gamma_p = result.gamma(p);
                amp_p   = result.amplitude(p);
                
                if gamma_p / omega_p > fit_cfg.max_gamma_ratio, continue; end
                
                raw_h = measure_peak_height( ...
                    energy_axis, spec, omega_p, gamma_p);
                
                row = [q_val, omega_p, gamma_p, result.R_squared, amp_p, ...
                       NaN, NaN, NaN, NaN, NaN, NaN, raw_h];
                all_peaks = [all_peaks; row]; %#ok<AGROW>
                n_success = n_success + 1;
            end
        end
    catch
    end
end

fprintf('  Fitted %d peaks from %d channels\n', size(all_peaks, 1), n_q);

% Branch separation (energy threshold — no toolbox needed)
if size(all_peaks, 1) >= 6
    % Simple separation: plasmon (< 1500 meV) vs interband (> 1500 meV)
    E_threshold = 1500;  % meV
    mask_lo = all_peaks(:,2) < E_threshold;
    mask_hi = all_peaks(:,2) >= E_threshold;
    
    branches = {};
    if any(mask_lo), branches{end+1} = all_peaks(mask_lo, :); end
    if any(mask_hi), branches{end+1} = all_peaks(mask_hi, :); end
    
    fprintf('  Branch separation: %d branches (threshold = %d meV)\n', numel(branches), E_threshold);
    for b = 1:numel(branches)
        fprintf('    B%d: %d peaks, <E>=%.0f meV, <Γ>=%.0f meV\n', ...
            b, size(branches{b}, 1), mean(branches{b}(:,2)), mean(branches{b}(:,3)));
    end
else
    branches = {all_peaks};
    fprintf('  Too few peaks for branch separation\n');
end

%% ═══════════════════════════════════════════════════════════════
%  PHASE 3: Self-consistent physics decoupling
%  ═══════════════════════════════════════════════════════════════
fprintf('\n── PHASE 3: Physics extraction (Do et al. 2025) ──\n');

phys = qe_physics_extract(branches, struct('branch_index', 1, ...
    'q_min_fit', 0.01, 'q_max_linear', 0.08, ...
    'epsilon_s', 1, 'symmetry_average', true));

fprintf('  Quasi-2D A = %.4g meV²·Å\n', phys.Drude_weight);
fprintf('  Keldysh ρ₀ = %.1f Å\n', phys.rho0);
fprintf('  I_kin slopes: low-q = %.2f (expect -3), high-q = %.2f (expect -2)\n', ...
    phys.I_kin_slope_lo, phys.I_kin_slope_hi);
if isfinite(phys.q_crossover)
    fprintf('  Dimensional crossover q* = %.4f Å⁻¹\n', phys.q_crossover);
end

%% ═══════════════════════════════════════════════════════════════
%  PHASE 4: Publication figures
%  ═══════════════════════════════════════════════════════════════
fprintf('\n── PHASE 4: Publication figures ──\n');

branch_colors = [0.1 0.5 0.9; 0.9 0.3 0.1; 0.2 0.7 0.3];

% ─── Figure 1: Dispersion relation ω_p(q) ───
fig1 = figure('Name', 'Dispersion', 'Color', 'w', 'Position', [50 200 700 500]);
hold on;
for b = 1:numel(branches)
    br = branches{b};
    col = branch_colors(mod(b-1, 3)+1, :);
    scatter(br(:,1), br(:,2), 40, col, 'filled', 'MarkerFaceAlpha', 0.7, ...
        'DisplayName', sprintf('Branch %d (<E>=%.0f meV)', b, mean(br(:,2))));
end
xlabel('q (Å^{-1})', 'FontSize', 13);
ylabel('\omega_p (meV)', 'FontSize', 13);
title('Plasmon Dispersion — Bi thin film', 'FontSize', 14);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

exportgraphics(fig1, fullfile(output_dir, 'fig1_dispersion.png'), 'Resolution', 300);
fprintf('  Saved: fig1_dispersion.png\n');

% ─── Figure 2: I_kin dimensional crossover (log-log) ───
fig2 = figure('Name', 'I_kin Crossover', 'Color', 'w', 'Position', [100 200 700 500]);
loglog(phys.q, phys.I_kin, 'ko', 'MarkerFaceColor', [0.1 0.5 0.9], ...
    'MarkerSize', 6, 'DisplayName', 'I_{kin}(q) = A(q) \times \epsilon_{inter}(q)');
hold on;

% Piecewise fits
q_med = median(phys.q);
q_lo = phys.q(phys.q <= q_med);
q_hi = phys.q(phys.q > q_med);
Ik_lo = phys.I_kin(phys.q <= q_med);
Ik_hi = phys.I_kin(phys.q > q_med);

if numel(q_lo) >= 2
    C_lo = exp(median(log(Ik_lo) - phys.I_kin_slope_lo * log(q_lo)));
    q_fit_lo = linspace(min(q_lo), max(q_lo), 100);
    loglog(q_fit_lo, C_lo * q_fit_lo.^phys.I_kin_slope_lo, '--', ...
        'Color', [0.8 0.2 0.2], 'LineWidth', 2, ...
        'DisplayName', sprintf('Low-q: q^{%.1f}', phys.I_kin_slope_lo));
end
if numel(q_hi) >= 2
    C_hi = exp(median(log(Ik_hi) - phys.I_kin_slope_hi * log(q_hi)));
    q_fit_hi = linspace(min(q_hi), max(q_hi), 100);
    loglog(q_fit_hi, C_hi * q_fit_hi.^phys.I_kin_slope_hi, '-.', ...
        'Color', [0.2 0.7 0.2], 'LineWidth', 2, ...
        'DisplayName', sprintf('High-q: q^{%.1f}', phys.I_kin_slope_hi));
end
if isfinite(phys.q_crossover)
    xline(phys.q_crossover, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
        'Label', sprintf('q* = %.3f', phys.q_crossover), ...
        'LabelHorizontalAlignment', 'left', ...
        'HandleVisibility', 'off');
end

xlabel('q (Å^{-1})', 'FontSize', 13);
ylabel('I_{kin}(q) (a.u.)', 'FontSize', 13);
title('EELS Prefactor — Dimensional Crossover', 'FontSize', 14);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

exportgraphics(fig2, fullfile(output_dir, 'fig2_ikin_crossover.png'), 'Resolution', 300);
fprintf('  Saved: fig2_ikin_crossover.png\n');

% ─── Figure 3: ε_inter + Keldysh fit ───
fig3 = figure('Name', 'Screening', 'Color', 'w', 'Position', [150 200 700 500]);
plot(phys.q, phys.epsilon_inter, 'ko', 'MarkerFaceColor', [0.1 0.5 0.9], ...
    'MarkerSize', 6, 'DisplayName', '\epsilon_{inter}(q) from \omega_p');
hold on;
q_fit = linspace(min(phys.q), max(phys.q), 200);
plot(q_fit, 1 + phys.rho0 * q_fit, 'r--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Keldysh: 1 + %.0f q', phys.rho0));
xlabel('q (Å^{-1})', 'FontSize', 13);
ylabel('\epsilon_{inter}(q)', 'FontSize', 13);
title(sprintf('Interband Screening (\\rho_0 = %.0f Å)', phys.rho0), 'FontSize', 14);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

exportgraphics(fig3, fullfile(output_dir, 'fig3_screening.png'), 'Resolution', 300);
fprintf('  Saved: fig3_screening.png\n');

% ─── Figure 4: Loss function map (side by side) ───
map_opts = struct();
map_opts.mode = 'experimental';
map_opts.E_range = E_range;
map_opts.q_range = q_range;
fig4 = qe_loss_map(qe_pp, phys, map_opts);

exportgraphics(fig4, fullfile(output_dir, 'fig4_loss_map.png'), 'Resolution', 300);
fprintf('  Saved: fig4_loss_map.png\n');

% ─── Figure 5: Γ & Amplitude dashboard ───
fig5 = qe_gamma_dashboard(all_peaks, branches);
exportgraphics(fig5, fullfile(output_dir, 'fig5_gamma_dashboard.png'), 'Resolution', 300);
fprintf('  Saved: fig5_gamma_dashboard.png\n');

% ─── Figure 6: Quality factor & damping ───
fig6 = figure('Name', 'Quality Factor', 'Color', 'w', 'Position', [250 200 700 500]);
subplot(1,2,1);
scatter(phys.q, phys.quality_factor, 40, [0.1 0.5 0.9], 'filled', 'MarkerFaceAlpha', 0.7);
xlabel('q (Å^{-1})'); ylabel('Q = \omega_p / \Gamma');
title('Plasmon Quality Factor'); grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1);

subplot(1,2,2);
scatter(phys.q, phys.gamma, 40, [0.9 0.3 0.1], 'filled', 'MarkerFaceAlpha', 0.7);
xlabel('q (Å^{-1})'); ylabel('\Gamma (meV)');
title('Damping Rate'); grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1);

sgtitle('Plasmon Lifetime Analysis — Bi', 'FontSize', 14);
exportgraphics(fig6, fullfile(output_dir, 'fig6_quality_factor.png'), 'Resolution', 300);
fprintf('  Saved: fig6_quality_factor.png\n');

%% Summary
fprintf('\n══════════════════════════════════════════════\n');
fprintf('  Pipeline complete! 6 figures saved to:\n');
fprintf('  %s\n', output_dir);
fprintf('══════════════════════════════════════════════\n');
