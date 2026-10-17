function root = b1_simple_root()
% Main checkout of Dispersion-Toolbox. The b1_simple_* scripts live on the
% simple-pipeline-20261017 worktree, but the operators, the P2b configuration and
% the data are uncommitted files of the main checkout, so they are read from there.
root = '/Users/nwh/Documents/Dispersion-Toolbox';
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'case_studies', 'bisb2026', 'scripts'));
end
