function root = b1_simple_root()
% Project root of Dispersion-Toolbox, with src/ and the BiSb scripts on the path.
% The data and paper_results/ are gitignored and exist only in the main checkout,
% so the b1_simple_* scripts have to run there, not from a separate worktree.
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'case_studies', 'bisb2026', 'scripts'));
end
