function paths = b1_forward_paths(cfg, run_root)
%B1_FORWARD_PATHS  Resolve explicit files or preserve the historical positional API.
if nargin < 2, run_root = []; end
if isfield(cfg, 'paths')
    paths = cfg.paths;
    assert(isempty(run_root) || strcmp(run_root, paths.run_root), ...
        'b1_forward_paths:Root', 'run_root differs from cfg.paths.run_root.');
elseif isfield(cfg, 'session')
    paths = b1_forward_config(cfg.session, run_root);
else
    paths = b1_forward_config([], run_root);
end
end
