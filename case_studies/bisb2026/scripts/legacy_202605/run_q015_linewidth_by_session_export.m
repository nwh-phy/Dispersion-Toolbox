function output = run_q015_linewidth_by_session_export()
%RUN_Q015_LINEWIDTH_BY_SESSION_EXPORT Split linewidth plots by session.

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

source_dir = fullfile(project_root, 'paper_results', 'wideq_linewidth_260506');
out_dir = fullfile(project_root, 'paper_results', ...
    'linewidth_q015_by_session_260507');
if ~isfolder(out_dir)
    mkdir(out_dir);
end

gamma_path = fullfile(source_dir, 'gamma_summary.csv');
do_path = fullfile(source_dir, 'do_style_linewidth_summary.csv');
if ~isfile(gamma_path) || ~isfile(do_path)
    error('run_q015_linewidth_by_session_export:MissingInput', ...
        'Missing linewidth CSV inputs in %s.', source_dir);
end

gamma_tbl = readtable(gamma_path);
do_tbl = readtable(do_path);
local_require_columns(gamma_tbl, {'session', 'session_label', 'branch', ...
    'q_abs_Ainv', 'gamma_meV', 'gamma_over_E', 'use_for_linewidth'});
local_require_columns(do_tbl, {'session', 'session_label', 'branch', ...
    'q_abs_Ainv', 'fano_fwhm_meV', 'lorentz_gamma_meV', ...
    'use_for_linewidth', 'fwhm_status'});

q_max_Ainv = 0.15;
gamma_q015 = gamma_tbl(gamma_tbl.q_abs_Ainv <= q_max_Ainv & ...
    local_truthy(gamma_tbl.use_for_linewidth), :);
do_q015 = do_tbl(do_tbl.q_abs_Ainv <= q_max_Ainv & ...
    local_truthy(do_tbl.use_for_linewidth), :);

gamma_filtered_csv = fullfile(out_dir, 'gamma_summary_q015.csv');
do_filtered_csv = fullfile(out_dir, 'do_style_linewidth_summary_q015.csv');
writetable(gamma_q015, gamma_filtered_csv);
writetable(do_q015, do_filtered_csv);

sessions = unique(string(gamma_q015.session), 'stable');
figure_rows = table();
branch_figure_rows = table();

for i = 1:numel(sessions)
    session_key = sessions(i);
    session_gamma = gamma_q015(string(gamma_q015.session) == session_key, :);
    session_do = do_q015(string(do_q015.session) == session_key, :);
    if isempty(session_gamma)
        continue;
    end

    label = char(string(session_gamma.session_label(1)));
    file_key = char(session_key);
    png_path = fullfile(out_dir, sprintf('%s_linewidth_q015.png', file_key));
    pdf_path = fullfile(out_dir, sprintf('%s_linewidth_q015.pdf', file_key));

    local_plot_session_linewidth(session_gamma, session_do, label, ...
        q_max_Ainv, png_path, pdf_path);

    figure_rows = [figure_rows; table({file_key}, {label}, ...
        {png_path}, {pdf_path}, ...
        'VariableNames', {'session', 'session_label', ...
        'png_path', 'pdf_path'})]; %#ok<AGROW>

    for branch_id = [1, 3]
        branch_png_path = fullfile(out_dir, ...
            sprintf('%s_B%d_linewidth_q015.png', file_key, branch_id));
        branch_pdf_path = fullfile(out_dir, ...
            sprintf('%s_B%d_linewidth_q015.pdf', file_key, branch_id));
        local_plot_one_branch_linewidth(session_gamma, session_do, label, ...
            branch_id, q_max_Ainv, branch_png_path, branch_pdf_path);

        branch_figure_rows = [branch_figure_rows; table({file_key}, ...
            {label}, branch_id, {branch_png_path}, {branch_pdf_path}, ...
            'VariableNames', {'session', 'session_label', ...
            'branch', 'png_path', 'pdf_path'})]; %#ok<AGROW>
    end
end

figure_index_csv = fullfile(out_dir, 'figure_index.csv');
writetable(figure_rows, figure_index_csv);
branch_figure_index_csv = fullfile(out_dir, 'figure_index_by_branch.csv');
writetable(branch_figure_rows, branch_figure_index_csv);

output = struct();
output.output_dir = out_dir;
output.figure_index = figure_rows;
output.branch_figure_index = branch_figure_rows;
output.figure_index_csv = figure_index_csv;
output.branch_figure_index_csv = branch_figure_index_csv;
output.gamma_filtered_csv = gamma_filtered_csv;
output.do_filtered_csv = do_filtered_csv;
output.q_max_Ainv = q_max_Ainv;

fprintf('q<=%.2f linewidth-by-session export complete.\n', q_max_Ainv);
fprintf('  Output directory: %s\n', out_dir);
fprintf('  Figures: %d\n', height(figure_rows));
fprintf('  Branch-split figures: %d\n', height(branch_figure_rows));
end


function local_require_columns(tbl, required)
names = tbl.Properties.VariableNames;
for i = 1:numel(required)
    if ~ismember(required{i}, names)
        error('run_q015_linewidth_by_session_export:MissingColumn', ...
            'Input table is missing required column "%s".', required{i});
    end
end
end


function mask = local_truthy(values)
if islogical(values)
    mask = values;
elseif isnumeric(values)
    mask = values ~= 0;
elseif iscell(values)
    mask = strcmpi(string(values), "true") | strcmp(string(values), "1");
else
    mask = strcmpi(string(values), "true") | strcmp(string(values), "1");
end
mask = mask(:);
end


function local_plot_session_linewidth(gamma_tbl, do_tbl, session_label, ...
        q_max_Ainv, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1180 780]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

branches = [1, 3];
for i = 1:numel(branches)
    branch_id = branches(i);
    ax = nexttile(i);
    local_plot_gamma_panel(ax, gamma_tbl, branch_id, q_max_Ainv);

    ax = nexttile(i + 2);
    local_plot_fwhm_panel(ax, do_tbl, branch_id, q_max_Ainv);
end

sgtitle(fig, sprintf('%s | linewidth analysis, |q| <= %.2f A^{-1}', ...
    session_label, q_max_Ainv), 'FontWeight', 'bold');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_one_branch_linewidth(gamma_tbl, do_tbl, session_label, ...
        branch_id, q_max_Ainv, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 780 760]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
local_plot_gamma_panel(ax, gamma_tbl, branch_id, q_max_Ainv);

ax = nexttile;
local_plot_fwhm_panel(ax, do_tbl, branch_id, q_max_Ainv);

sgtitle(fig, sprintf('%s | B%d linewidth, |q| <= %.2f A^{-1}', ...
    session_label, branch_id, q_max_Ainv), 'FontWeight', 'bold');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_gamma_panel(ax, tbl, branch_id, q_max_Ainv)
sub = tbl(tbl.branch == branch_id, :);
sub = sortrows(sub, 'q_abs_Ainv');
col = local_branch_color(branch_id);
hold(ax, 'on');
if isempty(sub)
    text(ax, 0.5, 0.5, 'No retained points', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
else
    scatter(ax, sub.q_abs_Ainv, sub.gamma_meV, 22, col, 'filled', ...
        'MarkerFaceAlpha', 0.85, 'DisplayName', '\Gamma_F_a_n_o');
end
hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [0, q_max_Ainv]);
xlabel(ax, '|q| (A^{-1})');
ylabel(ax, '\Gamma_F_a_n_o (meV)');
title(ax, sprintf('B%d fitted width parameter', branch_id));
end


function local_plot_fwhm_panel(ax, tbl, branch_id, q_max_Ainv)
sub = tbl(tbl.branch == branch_id, :);
if ismember('fwhm_status', sub.Properties.VariableNames)
    sub = sub(strcmp(string(sub.fwhm_status), "ok"), :);
end
sub = sub(isfinite(sub.fano_fwhm_meV) | isfinite(sub.lorentz_gamma_meV), :);
sub = sortrows(sub, 'q_abs_Ainv');
col = local_branch_color(branch_id);
hold(ax, 'on');
if isempty(sub)
    text(ax, 0.5, 0.5, 'No retained FWHM points', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
else
    scatter(ax, sub.q_abs_Ainv, sub.fano_fwhm_meV, 22, col, 'filled', ...
        'MarkerFaceAlpha', 0.85, 'DisplayName', 'Fano FWHM');
    scatter(ax, sub.q_abs_Ainv, sub.lorentz_gamma_meV, 28, ...
        'Marker', 's', 'MarkerEdgeColor', col, ...
        'MarkerFaceColor', 'none', 'LineWidth', 0.9, ...
        'DisplayName', 'Lorentz \Gamma');
end
hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [0, q_max_Ainv]);
xlabel(ax, '|q| (A^{-1})');
ylabel(ax, 'Width (meV)');
title(ax, sprintf('B%d FWHM / Lorentz-width check', branch_id));
legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
end


function col = local_branch_color(branch_id)
switch branch_id
    case 1
        col = [0.120, 0.470, 0.900];
    case 3
        col = [0.930, 0.280, 0.300];
    otherwise
        col = [0.2, 0.2, 0.2];
end
end
