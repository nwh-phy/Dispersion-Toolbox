function root = b1_simple_root()
% Project root of Dispersion-Toolbox, with src/ and the in-use BiSb scripts on the path
% (this folder, b1_binned/ and b1_lowq/; not the frozen legacy_* folders).
% The data and paper_results/ are gitignored and exist only in the main checkout,
% so the analysis scripts have to run there, not from a separate worktree.
here = fileparts(mfilename('fullpath'));
root = bisb_find_project_root(here);
addpath(genpath(fullfile(root, 'src')));
addpath(here, fullfile(here, 'b1_binned'), fullfile(here, 'b1_lowq'));
end
