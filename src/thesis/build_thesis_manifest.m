function manifest = build_thesis_manifest(cfg, sessions, results, meta)
%BUILD_THESIS_MANIFEST Create a provenance record for thesis q-EELS outputs.

if nargin < 1 || isempty(cfg)
    cfg = thesis_config();
end
if nargin < 2
    sessions = thesis_sessions(cfg);
end
if nargin < 3
    results = struct([]);
end
if nargin < 4 || isempty(meta)
    meta = struct();
end

manifest = struct();
manifest.generated_at = char(datetime('now', 'TimeZone', 'local', 'Format', "yyyyMMdd'T'HHmmss"));
manifest.pipeline = cfg.pipeline;
manifest.project_root = cfg.project_root;
manifest.analysis = cfg.analysis;
manifest.preprocess = cfg.preprocess;
manifest.fit = cfg.fit;
manifest.branch = cfg.branch;
manifest.symmetry = cfg.symmetry;
manifest.models = cfg.models;
manifest.sessions = local_session_manifest(sessions);
manifest.results = local_result_manifest(results);
manifest.git = local_git_meta(cfg.project_root, meta);
manifest.caution = struct( ...
    'sample_naming', cfg.notes.naming_caution, ...
    'physics_interpretation', cfg.notes.physics_caution, ...
    'runtime_scope', ['Manifest records the pipeline/configuration provenance. ', ...
                      'It does not by itself validate physical mode assignment.']);
end


function out = local_session_manifest(sessions)
out = sessions;
for i = 1:numel(out)
    if isfield(out(i), 'path')
        out(i).path_exists = exist(out(i).path, 'dir') == 7 || exist(out(i).path, 'file') == 2;
    end
end
end


function out = local_result_manifest(results)
if isempty(results)
    out = struct([]);
    return
end

template = struct( ...
    'session_name', '', ...
    'session_role', '', ...
    'data_path', '', ...
    'dq_Ainv', NaN, ...
    'success', false, ...
    'error', '', ...
    'dataset_label', '', ...
    'q_zero_index', NaN, ...
    'energy_range_meV', [NaN NaN], ...
    'q_range_Ainv', [NaN NaN], ...
    'n_fit_success', 0, ...
    'n_raw_peaks', 0, ...
    'branch_summary', struct('low_n', 0, 'high_n', 0, 'rejected_n', 0), ...
    'models', struct());
out = repmat(template, 1, numel(results));

for i = 1:numel(results)
    r = results(i);
    out(i).session_name = local_get_field(r, 'session_name', '');
    out(i).session_role = local_get_field(r, 'session_role', '');
    out(i).data_path = local_get_field(r, 'data_path', '');
    out(i).dq_Ainv = local_get_field(r, 'dq_Ainv', NaN);
    out(i).success = logical(local_get_field(r, 'success', false));
    out(i).error = local_get_field(r, 'error', '');
    out(i).dataset_label = local_get_field(r, 'dataset_label', '');
    out(i).q_zero_index = local_get_field(r, 'q_zero_index', NaN);
    out(i).energy_range_meV = local_get_vector_field(r, 'energy_range_meV', [NaN NaN]);
    out(i).q_range_Ainv = local_get_vector_field(r, 'q_range_Ainv', [NaN NaN]);
    out(i).n_fit_success = local_get_field(r, 'n_fit_success', 0);
    out(i).n_raw_peaks = local_get_field(r, 'n_raw_peaks', 0);
    out(i).branch_summary = local_branch_summary(local_get_field(r, 'branch_summary', struct()));
    if isfield(r, 'models')
        out(i).models = local_model_manifest(r.models);
    end
end
end


function summary = local_branch_summary(in)
summary = struct( ...
    'low_n', local_get_field(in, 'low_n', 0), ...
    'high_n', local_get_field(in, 'high_n', 0), ...
    'rejected_n', local_get_field(in, 'rejected_n', 0));
end


function models = local_model_manifest(in)
models = struct();
if isfield(in, 'low')
    models.low = local_one_model_manifest(in.low);
end
if isfield(in, 'high')
    models.high = local_one_model_manifest(in.high);
end
end


function out = local_one_model_manifest(in)
out = struct( ...
    'branch', local_get_field(in, 'branch', ''), ...
    'model', local_get_field(in, 'model', ''), ...
    'success', logical(local_get_field(in, 'success', false)), ...
    'error', local_get_field(in, 'error', ''), ...
    'fit', struct());
if isfield(in, 'fit')
    out.fit = local_fit_manifest(in.fit);
end
end


function out = local_fit_manifest(fit)
out = struct();
scalar_fields = {'R_squared', 'RMSE_meV', 'rho0', 'E_flat_meV', ...
    'q_c_Ainv', 'epsilon_s', 'A'};
for i = 1:numel(scalar_fields)
    name = scalar_fields{i};
    if isfield(fit, name)
        out.(name) = local_get_field(fit, name, NaN);
    end
end
if isfield(fit, 'params')
    out.params = local_get_vector_field(fit, 'params', []);
end
if isfield(fit, 'param_names')
    out.param_names = fit.param_names;
end
if isfield(fit, 'params_ci')
    out.params_ci = local_get_vector_field(fit, 'params_ci', []);
end
if isfield(fit, 'model_name')
    out.model_name = local_get_field(fit, 'model_name', '');
end
if isfield(fit, 'model_label')
    out.model_label = local_get_field(fit, 'model_label', '');
end
end


function value = local_get_field(s, field_name, default_value)
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    value = s.(field_name);
else
    value = default_value;
end
if isnumeric(value) || islogical(value)
    value = value(1);
elseif isstring(value)
    value = char(value(1));
elseif ischar(value)
    value = char(value);
elseif isstruct(value) && isstruct(default_value)
    % Keep small nested structs such as branch_summary; bulky result structs
    % are handled explicitly by local_result_manifest.
else
    value = default_value;
end
end


function value = local_get_vector_field(s, field_name, default_value)
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name)) && isnumeric(s.(field_name))
    value = s.(field_name);
else
    value = default_value;
end
end


function git = local_git_meta(root_dir, meta)
git = struct('commit', '', 'branch', '', 'status_short', '');
if isfield(meta, 'git_commit')
    git.commit = meta.git_commit;
end
if isfield(meta, 'git_branch')
    git.branch = meta.git_branch;
end
if isfield(meta, 'git_status_short')
    git.status_short = meta.git_status_short;
end

if isempty(git.commit)
    [ok, txt] = system(sprintf('git -C "%s" rev-parse HEAD', root_dir));
    if ok == 0
        git.commit = strtrim(txt);
    end
end
if isempty(git.branch)
    [ok, txt] = system(sprintf('git -C "%s" branch --show-current', root_dir));
    if ok == 0
        git.branch = strtrim(txt);
    end
end
if isempty(git.status_short)
    [ok, txt] = system(sprintf('git -C "%s" status --short', root_dir));
    if ok == 0
        git.status_short = strtrim(txt);
    end
end
end
