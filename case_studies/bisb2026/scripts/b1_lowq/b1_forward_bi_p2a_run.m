function b1_forward_bi_p2a_run(phase,diagnostic_folder,ids)
%B1_FORWARD_BI_P2A_RUN  Run P2a for four datasets serially in one batch.
arguments
    phase (1,:) char {mustBeMember(phase,{'response','diagnostics'})}
    diagnostic_folder (1,:) char = 'b1_forward_bi_checks'
    ids (1,:) cell = {'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'}
end
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
assert(~isempty(regexp(diagnostic_folder,'^b1_forward_bi_[a-z0-9_]+$','once')),'Invalid diagnostic folder.');
failed = {};
for task = ids
    id = task{1};
    cfg = b1_forward_bi_p2a_config(id);
    t = tic;
    try
        fprintf('\n%s %s %s\n',id,phase,char(datetime('now')));
        if strcmp(phase,'response')
            b1_p2a_inputs([],cfg.input);
            b1_p2a_response([],cfg.response);
        else
            cfg.diagnostic_output_dir = fullfile(cfg.input.output_dir,diagnostic_folder);
            b1_forward_bi_p2a_check(cfg);
        end
        fprintf('DONE %s %s %.3f s\n',id,phase,toc(t));
    catch err
        failed{end+1} = id; %#ok<AGROW>
        fprintf(2,'FAIL %s %s\n%s\n',id,phase,getReport(err,'extended'));
    end
end
assert(isempty(failed),'P2a phase %s failed for %s.',phase,strjoin(failed,', '));
end
