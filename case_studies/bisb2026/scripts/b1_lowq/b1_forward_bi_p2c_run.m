function b1_forward_bi_p2c_run(phase, run, ids)
%B1_FORWARD_BI_P2C_RUN  Sequenced phases with per-acquisition error isolation.
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
base = fullfile(root,'paper_results','b1_forward_bi_20261119','p2c');
run = fullfile(base,run);
assert(isfile(fullfile(run,'started_utc.txt')),'Missing P2c starting record.');
started = strtrim(fileread(fullfile(run,'started_utc.txt')));
deadline = datetime(started,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''','TimeZone','UTC')+hours(4)+minutes(50);
if nargin<3, ids = {'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'}; end
tasks = cellfun(@(id) b1_forward_bi_p2c_config(id,run),ids,'UniformOutput',false);
nFail = 0;
for k = 1:numel(tasks)
    if datetime('now','TimeZone','UTC')>=deadline
        fprintf(2,'Wall-clock cutoff. Unstarted: %s\n',strjoin(ids(k:end),', '));
        break
    end
    try
        P = tasks{k};
        out = fullfile(run,phase,P.id);
        assert(~isfolder(out),'Output already exists: %s.',out);
        mkdir(out);
        t = tic;
        fprintf('START %s %s\n',phase,P.id);
        switch phase
            case {'odd','odd_verified'}
                b1_forward_bi_p2c_odd(P,out);
            case {'asymmetry','asymmetry_verified'}
                b1_forward_bi_p2c_asymfit(P,out);
            otherwise
                assert(startsWith(phase,'low_'),'Unknown P2c phase %s.',phase);
                assert(isfile(fullfile(run,'step1_conclusion.log')), ...
                    'Complete and log step 1 before executing the low-window stage.');
                b1_forward_bi_p2c_low(P,out,phase,deadline);
        end
        fprintf('DONE %s %s %.2f s\n',phase,P.id,toc(t));
    catch err
        nFail = nFail+1;
        fprintf(2,'FAILED %s %s %s\n',phase,ids{k},getReport(err,'extended','hyperlinks','off'));
    end
end
assert(nFail==0,'%d task(s) failed.',nFail);
end
