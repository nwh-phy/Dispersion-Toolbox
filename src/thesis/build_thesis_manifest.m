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
manifest.results = results;
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
