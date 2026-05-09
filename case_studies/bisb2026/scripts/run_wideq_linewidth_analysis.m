function output = run_wideq_linewidth_analysis(options)
%RUN_WIDEQ_LINEWIDTH_ANALYSIS Extend q coverage for B1/B3 linewidth checks.
%   This is a linewidth-only layer. It keeps the current Fano-apex peak
%   positions, then adds numerical Fano FWHM, single-peak Drude-Lorentz
%   comparison fits, and a pointwise uncertainty budget.

arguments
    options.reuseExisting (1,1) logical = true
end

script_path = mfilename('fullpath');
project_root = bisb_find_project_root(fileparts(script_path));
run(fullfile(project_root, 'startup.m'));

out_dir = fullfile(project_root, 'paper_results', 'wideq_linewidth_260506');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

session_outputs = local_load_or_run_wideq_sessions(project_root, options.reuseExisting);

gamma_summary = local_collect_gamma_summary(session_outputs);
writetable(gamma_summary, fullfile(out_dir, 'gamma_summary.csv'));

do_style_summary = local_build_do_style_summary(session_outputs, gamma_summary);
writetable(do_style_summary, fullfile(out_dir, 'do_style_linewidth_summary.csv'));

uncertainty_budget = local_build_uncertainty_budget(do_style_summary);
writetable(uncertainty_budget, fullfile(out_dir, 'linewidth_uncertainty_budget.csv'));

figure_paths = struct();
figure_paths.b1 = fullfile(out_dir, 'b1_linewidth_q.png');
figure_paths.b3 = fullfile(out_dir, 'b3_linewidth_q.png');
figure_paths.quality = fullfile(out_dir, 'b1_b3_quality_factor.png');
figure_paths.map = fullfile(out_dir, 'wideq_boundary_map.png');
figure_paths.fwhm_comparison = fullfile(out_dir, 'b1_b3_fwhm_lorentz_comparison.png');

local_plot_linewidth(gamma_summary, 1, figure_paths.b1);
local_plot_linewidth(gamma_summary, 3, figure_paths.b3);
local_plot_quality_factor(gamma_summary, figure_paths.quality);
local_plot_boundary_map(session_outputs, figure_paths.map);
local_plot_fwhm_lorentz_comparison(do_style_summary, figure_paths.fwhm_comparison);

addendum_path = fullfile(out_dir, 'wideq_linewidth_addendum.md');
local_write_addendum(addendum_path, gamma_summary, do_style_summary, ...
    uncertainty_budget);

output = struct();
output.output_dir = out_dir;
output.session_outputs = session_outputs;
output.gamma_summary = gamma_summary;
output.do_style_summary = do_style_summary;
output.uncertainty_budget = uncertainty_budget;
output.figure_paths = figure_paths;
output.addendum_path = addendum_path;

fprintf('Wide-q linewidth output directory: %s\n', out_dir);
end


function session_outputs = local_load_or_run_wideq_sessions(project_root, reuseExisting)
tags = { ...
    '590_gui_history_area_260506_wideq030', ...
    'n0_PL2_10w_gui_history_area_260506_wideq030', ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined_wideq030'};
mat_paths = cellfun(@(tag) fullfile(project_root, 'paper_results', tag, ...
    'analysis_results.mat'), tags, 'UniformOutput', false);

if reuseExisting && all(cellfun(@isfile, mat_paths))
    session_outputs = cell(1, numel(mat_paths));
    for i = 1:numel(mat_paths)
        loaded = load(mat_paths{i}, 'output');
        session_outputs{i} = loaded.output;
    end
    return
end

analysis_output = run_590_gui_history_area_analysis("all", ...
    qRangeOverride_Ainv=[-0.30 0.30], outputTagSuffix="_wideq030");
session_outputs = analysis_output.sessions;
end


function gamma_summary = local_collect_gamma_summary(session_outputs)
gamma_summary = table();
for i = 1:numel(session_outputs)
    out = session_outputs{i};
    for branch_id = [1 3]
        br = out.branches{branch_id};
        if isempty(br)
            continue
        end
        q = br(:, 1);
        energy = br(:, 2);
        gamma = br(:, 3);
        r2 = br(:, 4);
        if size(br, 2) >= 7
            ci_half = 0.5 * (br(:, 7) - br(:, 6));
        else
            ci_half = NaN(size(q));
        end
        gamma_over_E = gamma ./ max(energy, eps);
        quality_factor = energy ./ gamma;
        quality_factor(~isfinite(quality_factor) | gamma <= 0) = NaN;
        use_for_linewidth = isfinite(energy) & isfinite(gamma) & ...
            isfinite(gamma_over_E) & isfinite(ci_half) & ...
            gamma >= 50 & gamma_over_E > 0 & gamma_over_E <= 2.0 & ...
            r2 >= 0.30 & ci_half <= 500;
        if branch_id == 1
            use_for_linewidth = use_for_linewidth & gamma_over_E >= 0.5;
        end
        linewidth_status = repmat({'used'}, numel(q), 1);
        linewidth_status(~use_for_linewidth) = {'low_confidence'};

        tbl = table( ...
            repmat({out.session.key}, numel(q), 1), ...
            repmat({out.session.display_name}, numel(q), 1), ...
            repmat(branch_id, numel(q), 1), ...
            q, abs(q), energy, gamma, gamma_over_E, quality_factor, r2, ci_half, ...
            use_for_linewidth, linewidth_status, ...
            repmat(out.session.dq_Ainv, numel(q), 1), ...
            repmat({out.output_dir}, numel(q), 1), ...
            'VariableNames', {'session', 'session_label', 'branch', ...
            'q_Ainv', 'q_abs_Ainv', 'energy_meV', 'gamma_meV', ...
            'gamma_over_E', 'quality_factor', 'R2', ...
            'E_ci_half_meV', 'use_for_linewidth', 'linewidth_status', ...
            'q_resolution_Ainv', 'source_dir'});
        gamma_summary = [gamma_summary; tbl]; %#ok<AGROW>
    end
end
end


function do_tbl = local_build_do_style_summary(session_outputs, gamma_summary)
do_tbl = gamma_summary;
extra_names = {'fano_fwhm_meV', 'fano_fwhm_left_meV', ...
    'fano_fwhm_right_meV', 'fano_fwhm_apex_meV', ...
    'lorentz_E0_meV', 'lorentz_gamma_meV', 'lorentz_R2', ...
    'lorentz_gamma_over_E', 'fano_to_lorentz_width_ratio'};
for i = 1:numel(extra_names)
    do_tbl.(extra_names{i}) = NaN(height(do_tbl), 1);
end
do_tbl.fwhm_status = repmat({'not_evaluated'}, height(do_tbl), 1);
do_tbl.lorentz_status = repmat({'not_evaluated'}, height(do_tbl), 1);

for s = 1:numel(session_outputs)
    out = session_outputs{s};
    session_mask = strcmp(do_tbl.session, out.session.key);
    idx_rows = find(session_mask);
    for r = idx_rows(:)'
        if ~do_tbl.use_for_linewidth(r)
            do_tbl.fwhm_status{r} = 'skipped_low_confidence';
            do_tbl.lorentz_status{r} = 'skipped_low_confidence';
            continue
        end

        detail = local_fit_detail_for_q(out, do_tbl.q_Ainv(r));
        peak_idx = local_match_detail_peak(detail, do_tbl.energy_meV(r));
        if peak_idx > 0
            fwhm = measure_peak_fwhm(detail.energy_fit(:), ...
                detail.peak_curves{peak_idx}(:));
            do_tbl.fano_fwhm_meV(r) = fwhm.fwhm_meV;
            do_tbl.fano_fwhm_left_meV(r) = fwhm.left_meV;
            do_tbl.fano_fwhm_right_meV(r) = fwhm.right_meV;
            do_tbl.fano_fwhm_apex_meV(r) = fwhm.apex_meV;
            do_tbl.fwhm_status{r} = char(fwhm.status);
        else
            do_tbl.fwhm_status{r} = 'no_matching_fano_peak';
        end

        lorentz = local_lorentz_refit(detail, do_tbl.energy_meV(r), ...
            do_tbl.branch(r), do_tbl.fano_fwhm_meV(r));
        do_tbl.lorentz_E0_meV(r) = lorentz.E0_meV;
        do_tbl.lorentz_gamma_meV(r) = lorentz.gamma_meV;
        do_tbl.lorentz_R2(r) = lorentz.R2;
        do_tbl.lorentz_gamma_over_E(r) = lorentz.gamma_over_E;
        do_tbl.lorentz_status{r} = lorentz.status;
        if isfinite(do_tbl.fano_fwhm_meV(r)) && isfinite(lorentz.gamma_meV) ...
                && lorentz.gamma_meV > 0
            do_tbl.fano_to_lorentz_width_ratio(r) = ...
                do_tbl.fano_fwhm_meV(r) / lorentz.gamma_meV;
        end
    end
end
end


function detail = local_fit_detail_for_q(out, q_Ainv)
[~, qi] = min(abs(out.qe_pp.q_Ainv(:) - q_Ainv));
detail = [];
if isfield(out.fit_res, 'fit_details') && numel(out.fit_res.fit_details) >= qi
    detail = out.fit_res.fit_details{qi};
end
end


function peak_idx = local_match_detail_peak(detail, energy_meV)
peak_idx = 0;
if isempty(detail) || ~isfield(detail, 'apex_energy_meV') ...
        || isempty(detail.apex_energy_meV)
    return
end
[dist, idx] = min(abs(detail.apex_energy_meV(:) - energy_meV));
if isfinite(dist) && dist <= 250
    peak_idx = idx;
end
end


function lorentz = local_lorentz_refit(detail, energy_meV, branch_id, fano_fwhm_meV)
lorentz = struct('E0_meV', NaN, 'gamma_meV', NaN, 'R2', NaN, ...
    'gamma_over_E', NaN, 'status', 'not_evaluated');
if isempty(detail) || ~isfield(detail, 'energy_data') || isempty(detail.energy_data)
    lorentz.status = 'missing_fit_detail';
    return
end

half_window = local_lorentz_half_window(branch_id, energy_meV, fano_fwhm_meV);
E_min = max(min(detail.energy_data), energy_meV - half_window);
E_max = min(max(detail.energy_data), energy_meV + half_window);
if E_max - E_min < 250
    lorentz.status = 'window_too_narrow';
    return
end

try
    fit = fit_loss_function(detail.energy_data(:), detail.spectrum_data(:), ...
        E_min=E_min, E_max=E_max, max_peaks=1, ...
        min_prominence=0.02, smooth_width=9, ...
        initial_guesses=energy_meV, peak_model='lorentz', ...
        pre_subtracted=false, bootstrap_ci_samples=0);
catch
    lorentz.status = 'fit_failed';
    return
end

if fit.n_peaks < 1
    lorentz.status = 'no_peak';
    return
end

lorentz.E0_meV = fit.omega_p(1);
lorentz.gamma_meV = fit.gamma(1);
lorentz.R2 = fit.R_squared;
lorentz.gamma_over_E = fit.gamma(1) / max(fit.omega_p(1), eps);
if isfinite(lorentz.R2) && lorentz.R2 >= 0.30 && isfinite(lorentz.gamma_meV)
    lorentz.status = 'ok';
else
    lorentz.status = 'low_confidence';
end
end


function half_window = local_lorentz_half_window(branch_id, energy_meV, fano_fwhm_meV)
if branch_id == 1
    half_window = max(450, min(950, 0.75 * energy_meV));
else
    half_window = 520;
end
if isfinite(fano_fwhm_meV) && fano_fwhm_meV > 0
    half_window = max(half_window, min(1000, 0.75 * fano_fwhm_meV));
end
end


function budget = local_build_uncertainty_budget(do_tbl)
budget = do_tbl(:, {'session', 'session_label', 'branch', 'q_Ainv', ...
    'q_abs_Ainv', 'energy_meV', 'E_ci_half_meV', 'q_resolution_Ainv'});
budget.dispersion_slope_meV_per_Ainv = NaN(height(budget), 1);
budget.q_resolution_broadening_meV = NaN(height(budget), 1);
budget.instrument_resolution_meV = NaN(height(budget), 1);
budget.total_known_energy_uncertainty_meV = NaN(height(budget), 1);
budget.uncertainty_note = repmat({ ...
    'Fit CI and finite-q broadening are included; ZLP/instrument width is not yet propagated.'}, ...
    height(budget), 1);

sessions = unique(do_tbl.session, 'stable');
for s = 1:numel(sessions)
    for branch_id = [1 3]
        mask = strcmp(do_tbl.session, sessions{s}) & do_tbl.branch == branch_id;
        idx = find(mask);
        if numel(idx) < 3
            continue
        end
        [q_unique, e_unique] = local_mean_energy_by_abs_q( ...
            do_tbl.q_abs_Ainv(idx), do_tbl.energy_meV(idx));
        if numel(q_unique) < 3
            continue
        end
        denom = gradient(q_unique);
        denom(abs(denom) < eps) = NaN;
        slope = gradient(e_unique) ./ denom;
        interp_slope = interp1(q_unique, slope, do_tbl.q_abs_Ainv(idx), ...
            'linear', 'extrap');
        budget.dispersion_slope_meV_per_Ainv(idx) = interp_slope;
        budget.q_resolution_broadening_meV(idx) = ...
            abs(interp_slope) .* do_tbl.q_resolution_Ainv(idx) ./ 2;
    end
end
budget.total_known_energy_uncertainty_meV = hypot( ...
    budget.E_ci_half_meV, budget.q_resolution_broadening_meV);
end


function [q_unique, e_mean] = local_mean_energy_by_abs_q(q_abs, energy)
q_round = round(q_abs(:), 6);
[q_unique, ~, group_idx] = unique(q_round, 'stable');
e_mean = accumarray(group_idx, energy(:), [], @(x) mean(x, 'omitnan'));
[q_unique, order] = sort(q_unique);
e_mean = e_mean(order);
end


function local_plot_linewidth(gamma_summary, branch_id, out_path)
branch_tbl = gamma_summary(gamma_summary.branch == branch_id, :);
branch_tbl = branch_tbl(branch_tbl.use_for_linewidth, :);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 940 720]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
local_plot_session_points(ax1, branch_tbl, 'gamma_meV');
hold(ax1, 'off');
grid(ax1, 'on');
box(ax1, 'on');
xlim(ax1, [0 0.30]);
ylabel(ax1, '\Gamma_F_a_n_o parameter (meV)');
title(ax1, sprintf('B%d fitted Fano width parameter', branch_id));
legend(ax1, 'Location', 'best', 'FontSize', 8);

ax2 = nexttile;
hold(ax2, 'on');
local_plot_session_points(ax2, branch_tbl, 'gamma_over_E');
hold(ax2, 'off');
grid(ax2, 'on');
box(ax2, 'on');
xlim(ax2, [0 0.30]);
xlabel(ax2, '|q| (1/A)');
ylabel(ax2, '\Gamma_F_a_n_o/E');

exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function local_plot_quality_factor(gamma_summary, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 450]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for branch_id = [1 3]
    branch_tbl = gamma_summary(gamma_summary.branch == branch_id, :);
    branch_tbl = branch_tbl(branch_tbl.use_for_linewidth, :);
    ax = nexttile;
    hold(ax, 'on');
    local_plot_session_points(ax, branch_tbl, 'quality_factor');
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0 0.30]);
    xlabel(ax, '|q| (1/A)');
    ylabel(ax, 'Q = E/\Gamma_F_a_n_o');
    title(ax, sprintf('B%d quality factor from Fano parameter', branch_id));
    legend(ax, 'Location', 'best', 'FontSize', 8);
end

exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function local_plot_fwhm_lorentz_comparison(do_tbl, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 520]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for branch_id = [1 3]
    branch_tbl = do_tbl(do_tbl.branch == branch_id & do_tbl.use_for_linewidth, :);
    ax = nexttile;
    hold(ax, 'on');
    sessions = unique(branch_tbl.session, 'stable');
    for i = 1:numel(sessions)
        mask = strcmp(branch_tbl.session, sessions{i});
        sub = branch_tbl(mask, :);
        [q_sorted, order] = sort(sub.q_abs_Ainv);
        col = local_session_color(i);
        plot(ax, q_sorted, sub.fano_fwhm_meV(order), 'o', ...
            'LineStyle', 'none', 'MarkerSize', 4.5, 'LineWidth', 1.0, ...
            'Color', col, 'DisplayName', ...
            sprintf('%s Fano FWHM', char(string(sub.session_label{1}))));
        plot(ax, q_sorted, sub.lorentz_gamma_meV(order), 's', ...
            'LineStyle', 'none', 'MarkerSize', 4.5, 'LineWidth', 1.0, ...
            'Color', col, 'MarkerFaceColor', 'none', 'DisplayName', ...
            sprintf('%s Lorentz Gamma', char(string(sub.session_label{1}))));
    end
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0 0.30]);
    xlabel(ax, '|q| (1/A)');
    ylabel(ax, 'Width (meV)');
    title(ax, sprintf('B%d numerical FWHM vs Drude-Lorentz fit', branch_id));
    legend(ax, 'Location', 'best', 'FontSize', 6);
end

exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function local_plot_session_points(ax, tbl, value_name)
sessions = unique(tbl.session, 'stable');
markers = {'o', 's', '^'};
for i = 1:numel(sessions)
    mask = strcmp(tbl.session, sessions{i});
    sub = tbl(mask, :);
    [q_sorted, order] = sort(sub.q_abs_Ainv);
    y = sub.(value_name);
    y = y(order);
    plot(ax, q_sorted, y, markers{min(i, numel(markers))}, ...
        'LineStyle', 'none', 'MarkerSize', 5, 'LineWidth', 1.0, ...
        'Color', local_session_color(i), ...
        'DisplayName', char(string(sub.session_label{1})));
end
end


function local_plot_boundary_map(session_outputs, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1450 520]);
tiledlayout(fig, 1, numel(session_outputs), 'TileSpacing', 'compact', ...
    'Padding', 'compact');
for i = 1:numel(session_outputs)
    out = session_outputs{i};
    qe = out.qe_pp;
    snap = out.snap;
    q_lim = [-0.30 0.30];
    q_mask = qe.q_Ainv >= q_lim(1) & qe.q_Ainv <= q_lim(2);
    e_mask = qe.energy_meV >= min(snap.energyMin, snap.energyMax) & ...
        qe.energy_meV <= max(snap.energyMin, snap.energyMax);
    map = double(qe.intensity(e_mask, q_mask));

    ax = nexttile;
    imagesc(ax, qe.q_Ainv(q_mask), qe.energy_meV(e_mask), map);
    axis(ax, 'xy');
    colormap(ax, turbo);
    clim_vals = local_intensity_limits(map);
    if all(isfinite(clim_vals)) && clim_vals(1) < clim_vals(2)
        clim(ax, clim_vals);
    end
    hold(ax, 'on');
    local_overlay_branch(ax, out.branches{1}, local_branch_color(1), 'B1');
    local_overlay_branch(ax, out.branches{3}, local_branch_color(3), 'B3');
    hold(ax, 'off');
    xlim(ax, q_lim);
    ylim(ax, [min(snap.energyMin, snap.energyMax), max(snap.energyMin, snap.energyMax)]);
    xlabel(ax, 'q (1/A)');
    if i == 1
        ylabel(ax, 'Energy relative to ZLP (meV)');
    end
    title(ax, char(string(out.session.display_name)));
    grid(ax, 'on');
    legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'FontSize', 7);
end
exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function local_overlay_branch(ax, br, color, label_text)
if isempty(br)
    return
end
scatter(ax, br(:, 1), br(:, 2), 16, color, 'filled', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.6, ...
    'DisplayName', label_text);
end


function clim_vals = local_intensity_limits(map)
vals = double(map(:));
vals = vals(isfinite(vals));
if isempty(vals)
    clim_vals = [NaN NaN];
    return
end
vals = sort(vals);
lo_idx = max(1, round(0.02 * numel(vals)));
hi_idx = min(numel(vals), round(0.98 * numel(vals)));
clim_vals = [vals(lo_idx), vals(hi_idx)];
end


function color = local_session_color(index)
palette = [ ...
    0.090 0.350 0.750; ...
    0.850 0.325 0.098; ...
    0.000 0.500 0.250];
color = palette(min(index, size(palette, 1)), :);
end


function color = local_branch_color(branch_id)
palette = [ ...
    0.090 0.350 0.750; ...
    0.600 0.600 0.600; ...
    0.000 0.500 0.250];
color = palette(branch_id, :);
end


function local_write_addendum(addendum_path, gamma_summary, do_tbl, ...
    uncertainty_budget)
b1 = gamma_summary(gamma_summary.branch == 1, :);
b3 = gamma_summary(gamma_summary.branch == 3, :);
b1_do = do_tbl(do_tbl.branch == 1 & do_tbl.use_for_linewidth, :);
b3_do = do_tbl(do_tbl.branch == 3 & do_tbl.use_for_linewidth, :);

lines = {
    '# 宽 q 峰宽补充分析'
    ''
    '本补充页只讨论 B1/B3 的峰宽、相对峰宽和品质因子。动量范围扩展到 `[-0.30, 0.30] A^{-1}`，目的是检查更大动量下的可见性边界和阻尼增强线索；这里不把某一个 q 点定义为 Landau damping 临界点。'
    ''
    '## 1. 宽 q 图像中的可见性边界'
    ''
    '![wide q boundary map](wideq_boundary_map.png)'
    ''
    'B1 在较大 |q| 处更容易受到宽峰和低信噪比影响；20w 2film 的高 q 区域仍应作为低置信度边界处理。B3 整体更连续，对称性也更好。'
    ''
    '## 2. Fano 参数 Γ(q) 与 Q'
    ''
    '![B1 linewidth](b1_linewidth_q.png)'
    ''
    local_summary_sentence(b1, 'B1')
    ''
    '![B3 linewidth](b3_linewidth_q.png)'
    ''
    local_summary_sentence(b3, 'B3')
    ''
    '![B1 B3 quality factor](b1_b3_quality_factor.png)'
    ''
    'Q = E/Gamma 用于描述峰位能量和拟合宽度的相对尺度。若 Q 接近 1 或更低，该段更适合标为强阻尼或低置信度区域，而不是直接纳入主色散拟合。'
    ''
    '## 3. Do-style 对照：Fano 数值 FWHM 与 Drude-Lorentz Γ'
    ''
    '![FWHM Lorentz comparison](b1_b3_fwhm_lorentz_comparison.png)'
    ''
    local_do_style_sentence(b1_do, 'B1')
    ''
    local_do_style_sentence(b3_do, 'B3')
    ''
    '这张图的作用是把当前 Fano 峰位提取和 Do 等实验文献中常用的 Lorentz-Drude 阻尼参数放在同一张图里比较。若两者接近，说明峰形更接近单一阻尼振子；若两者差异很大，说明 Fano 非对称、连续谱耦合或局部背景对峰宽有显著贡献。'
    ''
    '## 4. 误差预算'
    ''
    local_uncertainty_sentence(uncertainty_budget)
    ''
    '当前误差预算已经包含峰位拟合 CI 和有限 q 分辨率导致的色散展宽估计；ZLP/仪器能量分辨率尚未传播进去，因此这些误差棒仍是当前处理口径下的已知部分，而不是完整实验误差。'
    ''
    '相关数据表：`gamma_summary.csv`、`do_style_linewidth_summary.csv`、`linewidth_uncertainty_budget.csv`。'
    };

fid = fopen(addendum_path, 'w', 'n', 'UTF-8');
if fid < 0
    error('run_wideq_linewidth_analysis:CannotWriteAddendum', ...
        'Cannot write %s', addendum_path);
end
cleanup = onCleanup(@() fclose(fid));
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
end


function sentence = local_summary_sentence(tbl, label)
if ismember('use_for_linewidth', tbl.Properties.VariableNames)
    tbl = tbl(tbl.use_for_linewidth, :);
end
if isempty(tbl)
    sentence = sprintf('%s 在当前宽 q 口径下没有可汇总的峰宽点。', label);
    return
end
q_max = max(tbl.q_abs_Ainv, [], 'omitnan');
gamma_med = median(tbl.gamma_meV, 'omitnan');
ratio_med = median(tbl.gamma_over_E, 'omitnan');
q_med = median(tbl.quality_factor, 'omitnan');
sentence = sprintf('%s 当前汇总到的最大 |q| 约为 %.3f A^{-1}，中位 Fano Γ 约 %.0f meV，中位 Γ/E 约 %.2f，中位 Q 约 %.2f。', ...
    label, q_max, gamma_med, ratio_med, q_med);
end


function sentence = local_do_style_sentence(tbl, label)
if isempty(tbl)
    sentence = sprintf('%s 暂无可用的 Do-style 峰宽对照点。', label);
    return
end
fwhm_med = median(tbl.fano_fwhm_meV, 'omitnan');
lorentz_med = median(tbl.lorentz_gamma_meV(strcmp(tbl.lorentz_status, 'ok')), ...
    'omitnan');
ratio_med = median(tbl.fano_to_lorentz_width_ratio, 'omitnan');
sentence = sprintf('%s 的中位 Fano 数值 FWHM 约 %.0f meV；Lorentz-Drude 对照拟合中位 Γ 约 %.0f meV；FWHM/Gamma 的中位比值约 %.2f。', ...
    label, fwhm_med, lorentz_med, ratio_med);
end


function sentence = local_uncertainty_sentence(tbl)
known = tbl.total_known_energy_uncertainty_meV;
known = known(isfinite(known));
if isempty(known)
    sentence = '当前没有可汇总的逐点误差预算。';
    return
end
sentence = sprintf('已知部分的总能量不确定度中位数约 %.1f meV，90%% 分位约 %.1f meV。', ...
    median(known, 'omitnan'), local_percentile(known, 90));
end


function value = local_percentile(vals, pct)
vals = sort(vals(:));
vals = vals(isfinite(vals));
if isempty(vals)
    value = NaN;
    return
end
pos = 1 + (numel(vals) - 1) * pct / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    value = vals(lo);
else
    value = vals(lo) + (pos - lo) * (vals(hi) - vals(lo));
end
end
