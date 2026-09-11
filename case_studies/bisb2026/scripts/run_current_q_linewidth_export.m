function output = run_current_q_linewidth_export()
%RUN_CURRENT_Q_LINEWIDTH_EXPORT Export current-calibrated Gamma(q) figures.
%
% This reader-facing export uses the current branch point CSV files from the
% active no-BG area-normalized Fano-apex pipeline. It does not refit spectra.

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
out_dir = fullfile(project_root, 'paper_results', 'current_q_linewidth_260522');
if ~isfolder(out_dir)
    mkdir(out_dir);
end

sessions = local_sessions(project_root);
all_points = table();
for i = 1:numel(sessions)
    for branch = 1:3
        csv_path = fullfile(sessions(i).source_dir, ...
            sprintf('branch%d_points.csv', branch));
        if ~isfile(csv_path)
            error('run_current_q_linewidth_export:MissingBranchCsv', ...
                'Missing branch CSV: %s', csv_path);
        end
        tbl = readtable(csv_path);
        local_require_columns(tbl, {'q_Ainv', 'energy_meV', 'gamma_meV', ...
            'R2', 'E_ci_half_meV'});

        n = height(tbl);
        part = table();
        part.session = repmat(string(sessions(i).key), n, 1);
        part.session_label = repmat(string(sessions(i).label), n, 1);
        part.source_dir = repmat(string(sessions(i).source_dir), n, 1);
        part.branch = repmat(branch, n, 1);
        part.branch_label = repmat("B" + string(branch), n, 1);
        part.q_Ainv = tbl.q_Ainv;
        part.q_abs_Ainv = abs(tbl.q_Ainv);
        part.energy_meV = tbl.energy_meV;
        part.gamma_meV = tbl.gamma_meV;
        part.gamma_over_E = tbl.gamma_meV ./ max(tbl.energy_meV, eps);
        part.R2 = tbl.R2;
        part.E_ci_half_meV = tbl.E_ci_half_meV;
        all_points = [all_points; part]; %#ok<AGROW>
    end
end

summary = local_summary(all_points);
writetable(all_points, fullfile(out_dir, 'gamma_points_current_q.csv'));
writetable(summary, fullfile(out_dir, 'gamma_summary_current_q.csv'));

fig_branches_png = fullfile(out_dir, ...
    'fig_current_gamma_q_branches_by_dataset.png');
fig_branches_pdf = fullfile(out_dir, ...
    'fig_current_gamma_q_branches_by_dataset.pdf');
fig_b1_png = fullfile(out_dir, 'fig_current_b1_gamma_q_three_dataset.png');
fig_b1_pdf = fullfile(out_dir, 'fig_current_b1_gamma_q_three_dataset.pdf');
fig_b1_median_png = fullfile(out_dir, ...
    'fig_current_b1_gamma_q_paired_median.png');
fig_b1_median_pdf = fullfile(out_dir, ...
    'fig_current_b1_gamma_q_paired_median.pdf');

local_plot_branches_by_dataset(all_points, sessions, fig_branches_png, ...
    fig_branches_pdf);
local_plot_b1_by_dataset(all_points, sessions, fig_b1_png, fig_b1_pdf);
local_plot_b1_paired_median(all_points, sessions, fig_b1_median_png, ...
    fig_b1_median_pdf);
local_write_report(out_dir, sessions, summary, ...
    {fig_branches_png, fig_b1_png, fig_b1_median_png});

output = struct();
output.output_dir = out_dir;
output.points = all_points;
output.summary = summary;
output.figures = struct( ...
    'branches_png', fig_branches_png, ...
    'branches_pdf', fig_branches_pdf, ...
    'b1_png', fig_b1_png, ...
    'b1_pdf', fig_b1_pdf, ...
    'b1_median_png', fig_b1_median_png, ...
    'b1_median_pdf', fig_b1_median_pdf);

fprintf('Current-q linewidth output directory: %s\n', out_dir);
end


function sessions = local_sessions(project_root)
sessions = struct( ...
    'key', {}, 'label', {}, 'source_dir', {});

sessions(1).key = '10w_1film';
sessions(1).label = '10w 1film';
sessions(1).source_dir = fullfile(project_root, 'paper_results', ...
    '590_gui_history_area_260506');

sessions(2).key = '10w_1film_repeat';
sessions(2).label = '10w 1film repeat';
sessions(2).source_dir = fullfile(project_root, 'paper_results', ...
    'n0_PL2_10w_gui_history_area_260506');

sessions(3).key = '20w_2film';
sessions(3).label = '20w 2film';
sessions(3).source_dir = fullfile(project_root, 'paper_results', ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined');
end


function local_require_columns(tbl, names)
missing = setdiff(names, tbl.Properties.VariableNames);
if ~isempty(missing)
    error('run_current_q_linewidth_export:MissingColumn', ...
        'Missing columns: %s', strjoin(missing, ', '));
end
end


function summary = local_summary(all_points)
keys = unique(all_points(:, {'session', 'session_label', 'branch', ...
    'branch_label'}), 'rows', 'stable');
session = strings(height(keys), 1);
session_label = strings(height(keys), 1);
branch = zeros(height(keys), 1);
branch_label = strings(height(keys), 1);
n_points = zeros(height(keys), 1);
q_abs_min_Ainv = zeros(height(keys), 1);
q_abs_max_Ainv = zeros(height(keys), 1);
gamma_median_meV = zeros(height(keys), 1);
gamma_over_E_median = zeros(height(keys), 1);
E_ci_half_median_meV = zeros(height(keys), 1);
for i = 1:height(keys)
    mask = all_points.session == keys.session(i) & ...
        all_points.branch == keys.branch(i);
    sub = all_points(mask, :);
    session(i) = keys.session(i);
    session_label(i) = keys.session_label(i);
    branch(i) = keys.branch(i);
    branch_label(i) = keys.branch_label(i);
    n_points(i) = height(sub);
    q_abs_min_Ainv(i) = min(sub.q_abs_Ainv);
    q_abs_max_Ainv(i) = max(sub.q_abs_Ainv);
    gamma_median_meV(i) = median(sub.gamma_meV, 'omitnan');
    gamma_over_E_median(i) = median(sub.gamma_over_E, 'omitnan');
    E_ci_half_median_meV(i) = median(sub.E_ci_half_meV, 'omitnan');
end
summary = table(session, session_label, branch, branch_label, n_points, ...
    q_abs_min_Ainv, q_abs_max_Ainv, gamma_median_meV, ...
    gamma_over_E_median, E_ci_half_median_meV);
end


function local_plot_branches_by_dataset(all_points, sessions, out_png, out_pdf)
colors = local_branch_colors();
markers = {'o', 's', '^'};
fig = figure('Color', 'w', 'Position', [80 80 1500 820], ...
    'Visible', 'off');
t = tiledlayout(fig, 2, numel(sessions), 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(t, 'Current calibrated q: effective Fano linewidth by branch');

for i = 1:numel(sessions)
    session_mask = all_points.session == string(sessions(i).key);
    ax = nexttile(t, i);
    hold(ax, 'on');
    for branch = 1:3
        mask = session_mask & all_points.branch == branch;
        local_scatter(ax, all_points.q_abs_Ainv(mask), ...
            all_points.gamma_meV(mask), colors(branch, :), ...
            markers{branch}, sprintf('B%d', branch));
    end
    title(ax, sessions(i).label);
    xlabel(ax, '|q| (A^{-1})');
    if i == 1
        ylabel(ax, '\Gamma_{Fano} (meV)');
    end
    xlim(ax, [0 0.016]);
    ylim(ax, [0 2300]);
    grid(ax, 'on');
    box(ax, 'on');
    if i == 1
        legend(ax, 'Location', 'northwest');
    end

    ax = nexttile(t, numel(sessions) + i);
    hold(ax, 'on');
    for branch = 1:3
        mask = session_mask & all_points.branch == branch;
        local_scatter(ax, all_points.q_abs_Ainv(mask), ...
            all_points.gamma_over_E(mask), colors(branch, :), ...
            markers{branch}, sprintf('B%d', branch));
    end
    yline(ax, 1.0, '--', 'Q=1', 'Color', [0.45 0.45 0.45], ...
        'LabelHorizontalAlignment', 'left');
    xlabel(ax, '|q| (A^{-1})');
    if i == 1
        ylabel(ax, '\Gamma_{Fano}/E_p');
    end
    xlim(ax, [0 0.016]);
    ylim(ax, [0 2.1]);
    grid(ax, 'on');
    box(ax, 'on');
end

exportgraphics(fig, out_png, 'Resolution', 300);
exportgraphics(fig, out_pdf, 'ContentType', 'vector');
close(fig);
end


function local_plot_b1_by_dataset(all_points, sessions, out_png, out_pdf)
colors = lines(numel(sessions));
markers = {'o', 's', '^'};
fig = figure('Color', 'w', 'Position', [120 120 980 780], ...
    'Visible', 'off');
t = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(t, 'Current calibrated q: B1 effective linewidth');

for panel = 1:2
    ax = nexttile(t, panel);
    hold(ax, 'on');
    for i = 1:numel(sessions)
        mask = all_points.session == string(sessions(i).key) & ...
            all_points.branch == 1;
        if panel == 1
            y = all_points.gamma_meV(mask);
        else
            y = all_points.gamma_over_E(mask);
        end
        local_scatter(ax, all_points.q_abs_Ainv(mask), y, ...
            colors(i, :), markers{i}, sessions(i).label);
    end
    xlim(ax, [0 0.016]);
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '|q| (A^{-1})');
    if panel == 1
        ylabel(ax, '\Gamma_{Fano} (meV)');
        ylim(ax, [0 2300]);
        legend(ax, 'Location', 'southeast');
    else
        ylabel(ax, '\Gamma_{Fano}/E_p');
        ylim(ax, [0 2.1]);
        yline(ax, 1.0, '--', 'Q=1', 'Color', [0.45 0.45 0.45], ...
            'LabelHorizontalAlignment', 'left');
    end
end

exportgraphics(fig, out_png, 'Resolution', 300);
exportgraphics(fig, out_pdf, 'ContentType', 'vector');
close(fig);
end


function local_plot_b1_paired_median(all_points, sessions, out_png, out_pdf)
colors = lines(numel(sessions));
markers = {'o', 's', '^'};
fig = figure('Color', 'w', 'Position', [120 120 980 780], ...
    'Visible', 'off');
t = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(t, 'Current calibrated q: B1 linewidth, |q|-paired median');

for panel = 1:2
    ax = nexttile(t, panel);
    hold(ax, 'on');
    for i = 1:numel(sessions)
        mask = all_points.session == string(sessions(i).key) & ...
            all_points.branch == 1;
        q = all_points.q_abs_Ainv(mask);
        if panel == 1
            y = all_points.gamma_meV(mask);
        else
            y = all_points.gamma_over_E(mask);
        end
        pale = 0.70 + 0.30 .* colors(i, :);
        scatter(ax, q, y, 16, ...
            'Marker', markers{i}, ...
            'MarkerEdgeColor', pale, ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
        [q_med, y_med] = local_median_by_abs_q(q, y);
        plot(ax, q_med, y_med, '-', ...
            'Color', colors(i, :), ...
            'Marker', markers{i}, ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', 'w', ...
            'LineWidth', 1.6, ...
            'DisplayName', sessions(i).label);
    end
    xlim(ax, [0 0.016]);
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '|q| (A^{-1})');
    if panel == 1
        ylabel(ax, '\Gamma_{Fano} (meV)');
        ylim(ax, [0 2300]);
        legend(ax, 'Location', 'southeast');
    else
        ylabel(ax, '\Gamma_{Fano}/E_p');
        ylim(ax, [0 2.1]);
        yline(ax, 1.0, '--', 'Q=1', 'Color', [0.45 0.45 0.45], ...
            'LabelHorizontalAlignment', 'left');
    end
end

exportgraphics(fig, out_png, 'Resolution', 300);
exportgraphics(fig, out_pdf, 'ContentType', 'vector');
close(fig);
end


function [q_med, y_med] = local_median_by_abs_q(q, y)
q_key = round(q(:) .* 1e6) ./ 1e6;
[uq, ~, idx] = unique(q_key);
y_med = zeros(size(uq));
for i = 1:numel(uq)
    y_med(i) = median(y(idx == i), 'omitnan');
end
[q_med, order] = sort(uq);
y_med = y_med(order);
end


function local_scatter(ax, x, y, color, marker, name)
scatter(ax, x, y, 24, ...
    'Marker', marker, ...
    'MarkerEdgeColor', color, ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.0, ...
    'DisplayName', name);
end


function colors = local_branch_colors()
colors = [ ...
    0.00 0.35 0.75; ...
    0.85 0.33 0.10; ...
    0.12 0.55 0.25];
end


function local_write_report(out_dir, sessions, summary, figure_paths)
report_path = fullfile(out_dir, 'current_q_linewidth_report.md');
fid = fopen(report_path, 'w');
if fid < 0
    error('run_current_q_linewidth_export:CannotWriteReport', ...
        'Cannot write report: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Current calibrated q linewidth export\n\n');
fprintf(fid, 'Generated from current branch point CSV files. This export does not refit spectra.\n\n');
fprintf(fid, '## Sources\n\n');
for i = 1:numel(sessions)
    fprintf(fid, '- `%s`: `%s`\n', sessions(i).label, sessions(i).source_dir);
end
fprintf(fid, '\n## Figures\n\n');
for i = 1:numel(figure_paths)
    [~, name, ext] = fileparts(figure_paths{i});
    fprintf(fid, '- `%s%s`\n', name, ext);
end
fprintf(fid, '\n## Summary\n\n');
fprintf(fid, '| Dataset | Branch | N | |q|max (A^-1) | median Gamma (meV) | median Gamma/E | median CI half (meV) |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|\n');
for i = 1:height(summary)
    fprintf(fid, '| %s | %s | %d | %.5f | %.1f | %.3f | %.1f |\n', ...
        summary.session_label(i), summary.branch_label(i), ...
        summary.n_points(i), summary.q_abs_max_Ainv(i), ...
        summary.gamma_median_meV(i), summary.gamma_over_E_median(i), ...
        summary.E_ci_half_median_meV(i));
end
fprintf(fid, '\n## Manuscript boundary\n\n');
fprintf(fid, ['Gamma is an effective Fano linewidth from the same branch ' ...
    'tracking tables used for the current dispersion plots. It supports ' ...
    'a linewidth/damping boundary, not an independent microscopic lifetime ' ...
    'or absolute oscillator-strength conclusion.\n']);
end
