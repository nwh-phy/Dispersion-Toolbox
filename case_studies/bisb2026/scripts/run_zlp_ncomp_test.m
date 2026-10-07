function out = run_zlp_ncomp_test()
% Low-q ZLP check: do more Pearson components fit the gain side to the
% noise level at |q| = 0.0025 (ZLP ~0.6-1e6 counts), and does B1 then stop
% moving with the sub-window model? R2 (|q| = 0.0075) is the control.
% A1 spectra from the v7 mirror-q run, DL n=2 B1 peaks. Parent outputs stay read-only.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p7=fullfile(root,'paper_results','b1_mirror_q_v7','20261007T074250Z');
cfg=struct('labels',{{'R1','R1m','R2'}},'targets',[-.0025 .0025 .0075],'seed_offset',[1 1 2], ...
 'window',[300 1800],'model','lorentz','n_peaks',2,'J_starts',12,'seed',20260912, ...
 'n_zlp',[2 3 4],'core_halfwidths',[NaN 40]);
variants={'aux1',[30 300];'aux2f',[30 80;80 300]}; cfg.variants=variants;
id=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmss''Z'''));
out=fullfile(root,'paper_results','zlp_ncomp_test',[id '_A1']); assert(~isfolder(out));
mkdir(out); diary(fullfile(out,'execution.log')); dc=onCleanup(@()diary('off'));
disp(['RUN_DIR=' out]);
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root)); cfg.git_head=strtrim(head);
v7=load(fullfile(p7,'mirror_fits.mat')); keys=cellfun(@(a)a.key,v7.A,'UniformOutput',false);
rows={}; R=struct([]);
for r=1:numel(cfg.labels)
 d7=v7.A{strcmp(keys,[cfg.labels{r} '_' cfg.model])}; E=d7.E1(:); Y=d7.Y1(:);
 seed=cfg.seed+cfg.seed_offset(r);
 for nz=cfg.n_zlp
  for ch=cfg.core_halfwidths
   for v=1:size(variants,1)
    t0=tic; f=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
     'peak_model',cfg.model,'n_peaks',cfg.n_peaks,'n_starts',cfg.J_starts,'seed',seed, ...
     'aux_windows',variants{v,2},'core_halfwidth',ch,'n_zlp',nz); secs=toc(t0);
    p=nan(2,2); cl=NaN; cg=NaN; cs=NaN;
    if f.success
     p=f.parameters(:,1:2); cg=f.chi2_red_gain; cs=f.chi2_red_signal;
     x=f.energy_meV-f.zlp_parameters.center_meV; lo=x>=30&x<=300&isfinite(f.residual);
     cl=mean((f.residual(lo)./f.sigma(lo)).^2);
    end
    rows(end+1,:)={cfg.labels{r},cfg.targets(r),nz,f.effective_options.core_halfwidth,variants{v,1}, ...
     f.numerical_status,cg,cl,cs,p(1,1),p(1,2),p(2,1),p(2,2),secs}; %#ok<AGROW>
    fprintf('%s nz=%d core=%.0f %s: chi2 gain %.2f low %.2f signal %.2f | E0 %.1f / %.1f | %.0f s\n', ...
     cfg.labels{r},nz,f.effective_options.core_halfwidth,variants{v,1},cg,cl,cs,p(1,1),p(2,1),secs);
    R(end+1).label=cfg.labels{r}; R(end).n_zlp=nz; R(end).core=ch; R(end).variant=variants{v,1}; R(end).fit=f; %#ok<AGROW>
   end
  end
 end
end
T=cell2table(rows,'VariableNames',{'region','q_Ainv','n_zlp','core_halfwidth','variant','status', ...
 'chi2_gain','chi2_low','chi2_signal','E0_1','W_1','E0_2','W_2','seconds'});
writetable(T,fullfile(out,'ncomp_test.csv')); save(fullfile(out,'ncomp_test.mat'),'cfg','R','T','-v7.3');
disp(T);
end
