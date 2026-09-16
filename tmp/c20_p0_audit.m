cd(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(pwd,'src'))); addpath(genpath(fullfile(pwd,'lib')));
out = fullfile(pwd,'paper_results','b1_components_v2','p0_20260912_18af79b6');
assert(~isfolder(out),'Audit directory exists'); mkdir(out);
diary(fullfile(out,'environment_and_inputs.txt'));
disp(version); disp(ver); which load_qe_dataset -all; which run_thesis_pipeline -all;
s = thesis_sessions();
tags = {'590_gui_history_area_260506','no_PL2_20w_2film_gui_history_area_260506_highq_refined','n0_PL2_10w_gui_history_area_260506'};
for i=1:numel(s)
 disp(s(i));
 for name={'eq3D.mat','eq3D_processed.mat','op_history_260506.mat'}
  p=fullfile(s(i).path,name{1}); if ~isfile(p), continue; end
  fprintf('\nFILE %s\n',p); disp(whos('-file',p));
  if startsWith(name{1},'eq3D')
   a=load(p); disp(fieldnames(a)); fprintf('a3 %s e [%g %g] de %g\n',mat2str(size(a.a3)),min(a.e),max(a.e),median(diff(a.e)));
   if isfield(a,'q'), fprintf('q [%g %g] dq %g\n',min(a.q),max(a.q),median(diff(a.q))); end
   if isfield(a,'import_provenance'), disp(a.import_provenance); end
  else
   a=load(p); disp(fieldnames(a)); disp(a.opHistory{end});
  end
 end
 p=fullfile(pwd,'paper_results',tags{i},'analysis_results.mat'); a=load(p); fprintf('\nHISTORICAL %s\n',p); disp(fieldnames(a)); disp(fieldnames(a.output));
 disp(a.output.qe_pp.dq_Ainv); disp(a.output.qe_pp.q_zero_index); fprintf('saved q [%g %g] step %.8g\n',min(a.output.qe_pp.q_Ainv),max(a.output.qe_pp.q_Ainv),median(diff(a.output.qe_pp.q_Ainv)));
 t=readtable(fullfile(pwd,'paper_results',tags{i},'branch1_points.csv')); disp(t.Properties.VariableNames); disp(t(1:min(3,height(t)),:));
end
diary off;
diary(fullfile(out,'test_baseline.txt'));
r=runtests('tests'); save(fullfile(out,'test_baseline.mat'),'r'); writetable(table(r),fullfile(out,'test_baseline.csv'));
r_case=runtests('case_studies/bisb2026/tests/test_b1_double_peak_binning_workflow.m'); save(fullfile(out,'test_case.mat'),'r_case'); writetable(table(r_case),fullfile(out,'test_case.csv'));
diary off;
