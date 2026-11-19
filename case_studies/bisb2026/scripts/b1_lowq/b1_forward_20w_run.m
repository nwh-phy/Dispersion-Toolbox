function b1_forward_20w_run(out, stage)
%B1_FORWARD_20W_RUN  Run the P1 stages serially in one MATLAB process.
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
cfg = b1_forward_20w_config();
if nargin < 2, stage = 'all'; end
assert(isfolder(out) && startsWith(out,fullfile(cfg.paths.p2b_dir,'dispersion_20w_')), ...
    'Expected a new P1 output folder.');
tasks = {@() b1_forward_20w_inputs(cfg,out), @() b1_forward_20w_fit(cfg,out), ...
    @() b1_forward_20w_compare(cfg,out)};
names = {'S1','S2','S3'};
t = tic;
issues = {};
selected = find(strcmpi(stage,'all') | strcmpi(stage,names));
assert(~isempty(selected),'Unknown P1 stage.');
for k = selected
    try
        tk = tic;
        tasks{k}();
        fprintf('%s finished in %.2f s\n',names{k},toc(tk));
    catch err
        fprintf(2,'%s failed: %s\n',names{k},getReport(err,'extended'));
        issues{end+1} = names{k}; %#ok<AGROW>
    end
end
fprintf('Total %.2f s\n',toc(t));
assert(isempty(issues),'Failed P1 stages: %s.',strjoin(issues,', '));
end
