function output = run_thesis_pipeline(cfg, run_opts)
%RUN_THESIS_PIPELINE Execute the reproducible thesis q-EELS analysis chain.
%   OUTPUT = RUN_THESIS_PIPELINE() runs the baseline thesis pipeline for the
%   registered sessions and writes provenance-controlled outputs under
%   paper_results/thesis_qeels_pipeline/baseline/.

if nargin < 1 || isempty(cfg)
    cfg = thesis_config();
end
if nargin < 2 || isempty(run_opts)
    run_opts = struct();
end
run_opts = local_apply_run_defaults(run_opts, cfg);

sessions = thesis_sessions(cfg);
if isfield(run_opts, 'session_names') && ~isempty(run_opts.session_names)
    keep = ismember({sessions.name}, run_opts.session_names);
    sessions = sessions(keep);
else
    sessions = sessions([sessions.enabled]);
end

if run_opts.write_outputs && ~exist(run_opts.output_dir, 'dir')
    mkdir(run_opts.output_dir);
end

results = repmat(local_empty_session_result(), 1, numel(sessions));
for si = 1:numel(sessions)
    results(si) = local_run_one_session(sessions(si), cfg, run_opts);
end

manifest = build_thesis_manifest(cfg, sessions, results, struct());
output = struct();
output.cfg = cfg;
output.sessions = sessions;
output.results = results;
output.manifest = manifest;
output.output_dir = run_opts.output_dir;

if run_opts.write_outputs
    save(fullfile(run_opts.output_dir, 'analysis_results.mat'), 'output', 'cfg', 'sessions', 'results', 'manifest');
    local_write_json(fullfile(run_opts.output_dir, 'analysis_manifest.json'), manifest);
    local_write_summary_csv(fullfile(run_opts.output_dir, 'branch_summary.csv'), results);
    writetable(summarize_thesis_results(results, cfg.pipeline.version), ...
        fullfile(run_opts.output_dir, 'thesis_summary.csv'));
end
end


function run_opts = local_apply_run_defaults(run_opts, cfg)
if ~isfield(run_opts, 'write_outputs'), run_opts.write_outputs = true; end
if ~isfield(run_opts, 'output_dir') || isempty(run_opts.output_dir)
    run_opts.output_dir = cfg.output.baseline_dir;
end
if ~isfield(run_opts, 'session_names'), run_opts.session_names = {}; end
if ~isfield(run_opts, 'verbose'), run_opts.verbose = true; end
end


function result = local_run_one_session(session, cfg, run_opts)
result = local_empty_session_result();
result.session_name = session.name;
result.session_role = session.role;
result.data_path = session.path;
result.dq_Ainv = session.dq_Ainv;
result.success = false;

try
    if run_opts.verbose
        fprintf('\n[thesis] Session %s\n', session.name);
        fprintf('  path: %s\n', session.path);
    end
    dataset = load_qe_dataset(session.path, session.dq_Ainv);
    qe = dataset.qe;

    [qe_pp, bg_diag] = qe_preprocess(qe, cfg.preprocess);
    fit_opts = cfg.fit;
    fit_opts.energy_mask = qe_pp.energy_meV >= cfg.analysis.E_range_meV(1) & ...
                           qe_pp.energy_meV <= cfg.analysis.E_range_meV(2);
    fit_opts.energy_axis = qe_pp.energy_meV(fit_opts.energy_mask);
    fit_opts.q_start = cfg.analysis.q_range_Ainv(1);
    fit_opts.q_end = cfg.analysis.q_range_Ainv(2);
    fit_opts.verbose = false;

    fit_res = qe_auto_fit(qe_pp, qe_pp, fit_opts);
    in_q = abs(fit_res.all_peaks(:,1)) >= cfg.analysis.q_skip_Ainv & ...
           fit_res.all_peaks(:,1) >= cfg.analysis.q_range_Ainv(1) & ...
           fit_res.all_peaks(:,1) <= cfg.analysis.q_range_Ainv(2);
    fit_res.all_peaks = fit_res.all_peaks(in_q, :);

    branches = assign_qe_branches(fit_res.all_peaks, cfg.branch);
    branches = symmetrize_qe_branches(branches, cfg.symmetry);
    model_results = fit_thesis_branch_models(branches, cfg);

    result.success = true;
    result.dataset_label = dataset.label;
    result.q_zero_index = qe.q_zero_index;
    result.energy_range_meV = [min(qe.energy_meV) max(qe.energy_meV)];
    result.q_range_Ainv = [min(qe.q_Ainv) max(qe.q_Ainv)];
    result.n_fit_success = fit_res.n_success;
    result.n_raw_peaks = size(fit_res.all_peaks, 1);
    result.branch_summary = branches.summary;
    result.branches = branches;
    result.models = model_results;
    result.bg_diag = bg_diag;
    result.fit_details = fit_res.fit_details;

    if run_opts.write_outputs
        local_export_session_outputs(result, run_opts.output_dir);
    end
catch ME
    result.success = false;
    result.error = ME.message;
    if run_opts.verbose
        fprintf('  FAILED: %s\n', ME.message);
    end
end
end


function result = local_empty_session_result()
result = struct();
result.session_name = '';
result.session_role = '';
result.data_path = '';
result.dq_Ainv = NaN;
result.success = false;
result.error = '';
result.dataset_label = '';
result.q_zero_index = NaN;
result.energy_range_meV = [NaN NaN];
result.q_range_Ainv = [NaN NaN];
result.n_fit_success = 0;
result.n_raw_peaks = 0;
result.branch_summary = struct('low_n', 0, 'high_n', 0, 'rejected_n', 0);
result.branches = struct();
result.models = struct();
result.bg_diag = struct([]);
result.fit_details = {};
end


function local_export_session_outputs(result, output_dir)
session_dir = fullfile(output_dir, result.session_name);
if ~exist(session_dir, 'dir')
    mkdir(session_dir);
end
local_write_peak_csv(fullfile(session_dir, 'low_branch_signed.csv'), result.branches.low.accepted);
local_write_peak_csv(fullfile(session_dir, 'high_branch_signed.csv'), result.branches.high.accepted);
local_write_peak_csv(fullfile(session_dir, 'low_branch_symmetrized.csv'), result.branches.low.symmetrized);
local_write_peak_csv(fullfile(session_dir, 'high_branch_symmetrized.csv'), result.branches.high.symmetrized);
if ~isempty(result.branches.rejected)
    writetable(result.branches.rejected, fullfile(session_dir, 'rejected_peaks.csv'));
end
save(fullfile(session_dir, 'session_result.mat'), 'result');
end


function local_write_peak_csv(path, peaks)
cols = {'q_Ainv','energy_meV','gamma_meV','R2','amplitude_fit', ...
        'E_ci_lo','E_ci_hi','gamma_ci_lo','gamma_ci_hi', ...
        'A_ci_lo','A_ci_hi','raw_height'};
if isempty(peaks)
    tbl = cell2table(cell(0, numel(cols)), 'VariableNames', cols);
else
    if size(peaks,2) < numel(cols)
        peaks(:, end+1:numel(cols)) = NaN;
    end
    tbl = array2table(peaks(:,1:numel(cols)), 'VariableNames', cols);
end
writetable(tbl, path);
end


function local_write_json(path, data)
fid = fopen(path, 'w');
if fid < 0
    error('run_thesis_pipeline:CannotWriteJson', 'Cannot write %s', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
end


function local_write_summary_csv(path, results)
rows = cell(numel(results), 8);
for i = 1:numel(results)
    rows{i,1} = results(i).session_name;
    rows{i,2} = results(i).success;
    rows{i,3} = results(i).branch_summary.low_n;
    rows{i,4} = results(i).branch_summary.high_n;
    rows{i,5} = results(i).branch_summary.rejected_n;
    rows{i,6} = results(i).n_raw_peaks;
    rows{i,7} = results(i).n_fit_success;
    rows{i,8} = results(i).error;
end
tbl = cell2table(rows, 'VariableNames', {'session','success','low_n','high_n','rejected_n','raw_peak_n','fit_success_n','error'});
writetable(tbl, path);
end
