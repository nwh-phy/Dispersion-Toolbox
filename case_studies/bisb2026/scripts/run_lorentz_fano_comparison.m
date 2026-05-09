function output = run_lorentz_fano_comparison(options)
%RUN_LORENTZ_FANO_COMPARISON Compare Lorentz and Fano peak extraction.
%   Generates PPT-ready figures in q range [-0.15, 0.15] 1/A.

arguments
    options.reuseExisting (1,1) logical = true
end

script_path = mfilename('fullpath');
project_root = bisb_find_project_root(fileparts(script_path));
run(fullfile(project_root, 'startup.m'));

out_dir = fullfile(project_root, 'paper_results', ...
    'lorentz_fano_compare_260506');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fano_outputs = local_load_or_run_fano(project_root, options.reuseExisting);
lorentz_outputs = local_load_or_run_lorentz(project_root, options.reuseExisting);

comparison = table();
figure_rows = table();
for i = 1:numel(fano_outputs)
    fano = fano_outputs{i};
    lorentz = lorentz_outputs{i};
    key = fano.session.key;

    peak_path = fullfile(out_dir, ...
        sprintf('%s_peak_extraction_fano_vs_lorentz.png', key));
    gamma_path = fullfile(out_dir, ...
        sprintf('%s_gamma_fano_vs_lorentz.png', key));

    local_plot_peak_extraction(fano, lorentz, peak_path);
    local_plot_gamma_comparison(fano, lorentz, gamma_path);

    comparison = [comparison; local_compare_session(fano, lorentz)]; %#ok<AGROW>
    figure_rows = [figure_rows; table({key}, {peak_path}, {gamma_path}, ...
        'VariableNames', {'session', 'peak_extraction_figure', ...
        'gamma_comparison_figure'})]; %#ok<AGROW>
end

comparison_path = fullfile(out_dir, 'lorentz_fano_branch_comparison.csv');
writetable(comparison, comparison_path);

figure_index_path = fullfile(out_dir, 'figure_index.csv');
writetable(figure_rows, figure_index_path);

readme_path = fullfile(out_dir, 'lorentz_fano_compare_readme.md');
local_write_readme(readme_path, figure_rows, comparison_path);

output = struct();
output.output_dir = out_dir;
output.fano_outputs = fano_outputs;
output.lorentz_outputs = lorentz_outputs;
output.comparison = comparison;
output.figure_index = figure_rows;
output.comparison_path = comparison_path;
output.figure_index_path = figure_index_path;
output.readme_path = readme_path;

fprintf('Lorentz/Fano comparison output directory: %s\n', out_dir);
end


function outputs = local_load_or_run_fano(project_root, reuseExisting)
tags = strcat(local_base_tags(), '_fano_compare');
mat_paths = local_analysis_mat_paths(project_root, tags);
if reuseExisting && all(cellfun(@isfile, mat_paths))
    outputs = local_load_outputs(mat_paths);
    if local_outputs_have_session(outputs)
        return
    end
end
run_out = run_590_gui_history_area_analysis("all", ...
    qRangeOverride_Ainv=[-0.15 0.15], outputTagSuffix="_fano_compare", ...
    peakModelOverride="fano");
outputs = run_out.sessions;
end


function outputs = local_load_or_run_lorentz(project_root, reuseExisting)
tags = strcat(local_base_tags(), '_lorentz_compare');
mat_paths = local_analysis_mat_paths(project_root, tags);
if reuseExisting && all(cellfun(@isfile, mat_paths))
    outputs = local_load_outputs(mat_paths);
    if local_outputs_have_session(outputs)
        return
    end
end
run_out = run_590_gui_history_area_analysis("all", ...
    qRangeOverride_Ainv=[-0.15 0.15], outputTagSuffix="_lorentz_compare", ...
    peakModelOverride="lorentz");
outputs = run_out.sessions;
end


function tags = local_base_tags()
tags = { ...
    '590_gui_history_area_260506', ...
    'n0_PL2_10w_gui_history_area_260506', ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined'};
end


function paths = local_analysis_mat_paths(project_root, tags)
paths = cellfun(@(tag) fullfile(project_root, 'paper_results', tag, ...
    'analysis_results.mat'), tags, 'UniformOutput', false);
end


function outputs = local_load_outputs(mat_paths)
outputs = cell(1, numel(mat_paths));
for i = 1:numel(mat_paths)
    loaded = load(mat_paths{i}, 'output');
    outputs{i} = loaded.output;
end
end


function tf = local_outputs_have_session(outputs)
tf = true;
for i = 1:numel(outputs)
    tf = tf && isfield(outputs{i}, 'session') && isfield(outputs{i}.session, 'key');
end
end


function local_plot_peak_extraction(fano, lorentz, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1480 560]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
clim_vals = local_common_intensity_limits(fano.qe_pp, fano.snap);

local_plot_one_extraction_panel(nexttile, fano, 'Fano apex extraction', clim_vals);
local_plot_one_extraction_panel(nexttile, lorentz, 'Lorentz extraction', clim_vals);

sgtitle(fig, sprintf('%s | q range [-0.15, 0.15] 1/A', ...
    char(string(fano.session.display_name))));
exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function local_plot_one_extraction_panel(ax, out, title_text, clim_vals)
qe = out.qe_pp;
snap = out.snap;
q_lim = [-0.15 0.15];
q_mask = qe.q_Ainv >= q_lim(1) & qe.q_Ainv <= q_lim(2);
e_mask = qe.energy_meV >= min(snap.energyMin, snap.energyMax) & ...
    qe.energy_meV <= max(snap.energyMin, snap.energyMax);
map = double(qe.intensity(e_mask, q_mask));

imagesc(ax, qe.q_Ainv(q_mask), qe.energy_meV(e_mask), map);
axis(ax, 'xy');
colormap(ax, turbo);
if all(isfinite(clim_vals)) && clim_vals(1) < clim_vals(2)
    clim(ax, clim_vals);
end
hold(ax, 'on');
for b = 1:min(3, numel(out.branches))
    local_overlay_branch(ax, out.branches{b}, local_branch_color(b), ...
        sprintf('B%d', b), local_branch_marker(b));
end
hold(ax, 'off');
xlim(ax, q_lim);
ylim(ax, [min(snap.energyMin, snap.energyMax), max(snap.energyMin, snap.energyMax)]);
xlabel(ax, 'q (1/A)');
ylabel(ax, 'Energy relative to ZLP (meV)');
title(ax, title_text);
grid(ax, 'on');
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'FontSize', 8);
end


function local_overlay_branch(ax, br, color, label_text, marker)
if isempty(br)
    return
end
mask = br(:, 1) >= -0.15 & br(:, 1) <= 0.15;
br = br(mask, :);
if isempty(br)
    return
end
scatter(ax, br(:, 1), br(:, 2), 18, color, marker, 'filled', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'DisplayName', label_text);
end


function local_plot_gamma_comparison(fano, lorentz, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1350 760]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for b = 1:3
    ax = nexttile;
    hold(ax, 'on');
    local_plot_gamma_branch(ax, fano.branches, b, 'Fano gamma', ...
        local_branch_color(b), 'o');
    local_plot_gamma_branch(ax, lorentz.branches, b, 'Lorentz gamma', ...
        local_branch_color(b) * 0.75, 's');
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [-0.15 0.15]);
    ylabel(ax, sprintf('B%d Gamma (meV)', b));
    if b == 3
        xlabel(ax, 'q (1/A)');
    end
    title(ax, sprintf('B%d linewidth parameter comparison', b));
    legend(ax, 'Location', 'best', 'FontSize', 8);
end

sgtitle(fig, sprintf('%s | Lorentz vs Fano gamma', ...
    char(string(fano.session.display_name))));
exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function local_plot_gamma_branch(ax, branches, branch_id, label_text, color, marker)
if numel(branches) < branch_id || isempty(branches{branch_id})
    return
end
br = branches{branch_id};
mask = br(:, 1) >= -0.15 & br(:, 1) <= 0.15 & isfinite(br(:, 3));
br = br(mask, :);
if isempty(br)
    return
end
plot(ax, br(:, 1), br(:, 3), marker, 'LineStyle', 'none', ...
    'MarkerSize', 4.5, 'LineWidth', 1.0, 'Color', color, ...
    'DisplayName', label_text);
end


function comparison = local_compare_session(fano, lorentz)
comparison = table();
for b = 1:3
    fbr = local_branch_or_empty(fano.branches, b);
    lbr = local_branch_or_empty(lorentz.branches, b);
    [fmatch, lmatch] = local_pair_by_q(fbr, lbr);
    if isempty(fmatch)
        e_delta = NaN;
        g_delta = NaN;
    else
        e_delta = median(abs(fmatch(:, 2) - lmatch(:, 2)), 'omitnan');
        g_delta = median(abs(fmatch(:, 3) - lmatch(:, 3)), 'omitnan');
    end
    comparison = [comparison; table( ...
        {fano.session.key}, {fano.session.display_name}, b, ...
        size(fbr, 1), size(lbr, 1), size(fmatch, 1), ...
        median(fbr(:, 3), 'omitnan'), median(lbr(:, 3), 'omitnan'), ...
        e_delta, g_delta, ...
        'VariableNames', {'session', 'session_label', 'branch', ...
        'n_fano', 'n_lorentz', 'n_paired_q', ...
        'fano_gamma_median_meV', 'lorentz_gamma_median_meV', ...
        'median_abs_energy_delta_meV', ...
        'median_abs_gamma_delta_meV'})]; %#ok<AGROW>
end
end


function br = local_branch_or_empty(branches, branch_id)
if numel(branches) < branch_id || isempty(branches{branch_id})
    br = zeros(0, 12);
else
    br = branches{branch_id};
    br = br(br(:, 1) >= -0.15 & br(:, 1) <= 0.15, :);
end
end


function [fmatch, lmatch] = local_pair_by_q(fbr, lbr)
fmatch = [];
lmatch = [];
if isempty(fbr) || isempty(lbr)
    return
end
fq = round(fbr(:, 1), 6);
lq = round(lbr(:, 1), 6);
[common_q, fidx, lidx] = intersect(fq, lq, 'stable');
if isempty(common_q)
    return
end
fmatch = fbr(fidx, :);
lmatch = lbr(lidx, :);
end


function clim_vals = local_common_intensity_limits(qe, snap)
q_lim = [-0.15 0.15];
q_mask = qe.q_Ainv >= q_lim(1) & qe.q_Ainv <= q_lim(2);
e_mask = qe.energy_meV >= min(snap.energyMin, snap.energyMax) & ...
    qe.energy_meV <= max(snap.energyMin, snap.energyMax);
vals = double(qe.intensity(e_mask, q_mask));
vals = vals(isfinite(vals(:)));
if isempty(vals)
    clim_vals = [NaN NaN];
    return
end
vals = sort(vals);
lo_idx = max(1, round(0.02 * numel(vals)));
hi_idx = min(numel(vals), round(0.98 * numel(vals)));
clim_vals = [vals(lo_idx), vals(hi_idx)];
end


function color = local_branch_color(branch_id)
palette = [ ...
    0.090 0.350 0.750; ...
    0.500 0.500 0.500; ...
    0.000 0.500 0.250];
color = palette(branch_id, :);
end


function marker = local_branch_marker(branch_id)
markers = {'o', '^', 's'};
marker = markers{branch_id};
end


function local_write_readme(readme_path, figure_rows, comparison_path)
fid = fopen(readme_path, 'w', 'n', 'UTF-8');
if fid < 0
    error('run_lorentz_fano_comparison:CannotWriteReadme', ...
        'Cannot write %s', readme_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Lorentz vs Fano peak extraction comparison\n\n');
fprintf(fid, '- q range: `[-0.15, 0.15] 1/A`.\n');
fprintf(fid, '- Same area-normalized GUI-history settings; only the peak model is changed.\n');
fprintf(fid, '- Peak extraction figures compare branch positions on the physical q-E map.\n');
fprintf(fid, '- Gamma figures compare the fitted linewidth parameter for each branch.\n\n');
fprintf(fid, '## Figures\n\n');
for i = 1:height(figure_rows)
    fprintf(fid, '- `%s`\n', figure_rows.session{i});
    fprintf(fid, '  - `%s`\n', figure_rows.peak_extraction_figure{i});
    fprintf(fid, '  - `%s`\n', figure_rows.gamma_comparison_figure{i});
end
fprintf(fid, '\nSummary table: `%s`\n', comparison_path);
end
