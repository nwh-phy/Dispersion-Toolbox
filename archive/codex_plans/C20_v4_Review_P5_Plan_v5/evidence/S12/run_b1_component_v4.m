function out = run_b1_component_v4()
% Scoped reconstruction; parent outputs stay read-only.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
p3=fullfile(root,'paper_results','b1_components_v3','20260912T035735540Z');
cfg=struct('targets',[-.0025 .0075 .0125],'Ns',[1 3 5],'window',[300 1800], ...
 'legacy_starts',24,'independent_starts',24,'extra_n2_starts',12,'seed',20260912, ...
 'backgrounds',{{'power_law','power_law_plus_nonnegative_constant'}}, ...
 'models',{{'lorentz_symmetric','lorentz'}},'scope','P4_repair_and_engineering_P5');
sources=[dir(fullfile(root,'case_studies','bisb2026','scripts','c20_v4*.m')); ...
 dir(fullfile(root,'case_studies','bisb2026','scripts','run_b1_component_v4.m')); dir(fullfile(root,'tests','C20V4ContractsTest.m'))];
paths=string(fullfile({sources.folder},{sources.name}));
paths=[paths,string(fullfile(root,{'src/fitting/qe_component_mapping.m','src/fitting/qe_component_prediction.m', ...
 'src/fitting/qe_compare_component_models.m','src/fitting/peak_models.m','src/fitting/measure_peak_fwhm.m', ...
 'src/qe_centered_bins.m','src/qe_prepare_count_bins.m','src/qe_zlp_integer_align.m','src/io/read_npy.m', ...
 'src/io/load_raw_session.m','src/io/make_qe_struct.m','src/io/resolve_q_crop_bounds.m', ...
 'src/io/infer_qe_dq_Ainv.m','src/thesis/thesis_sessions.m','tests/test_b1_component_pilot_v2.m', ...
 'case_studies/bisb2026/scripts/bisb_find_project_root.m'}))];
hashes=strings(size(paths)); for k=1:numel(paths), hashes(k)=c20_v4_io('hash',paths(k)); end
codehash=c20_v4_io('digest',char(join(paths+hashes,'|'))); cfghash=c20_v4_io('digest',jsonencode(cfg));
id=[char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z''')) '_' cfghash(1:8) '_' codehash(1:8)];
out=fullfile(root,'paper_results','b1_components_v4',id); assert(~isfolder(out));
for name={'provenance','audit','590_PL2_10w','tests','validation','590_PL2_10w/figures','checkpoints'}, mkdir(fullfile(out,name{1})); end
diary(fullfile(out,'provenance','execution.log')); dc=onCleanup(@()diary('off')); disp(['RUN_DIR=' out]);
c20_v4_io('json',fullfile(out,'provenance','config.json'),cfg); save(fullfile(out,'provenance','config.mat'),'cfg');
c20_v4_io('json',fullfile(out,'provenance','parents.json'),struct('v2',p2,'v3',p3,'code_hash',codehash));
parent_before={c20_v4_io('inventory',p2),c20_v4_io('inventory',p3)};
save(fullfile(out,'provenance','parent_hashes_before.mat'),'parent_before');
for k=1:numel(paths)
 rel=extractAfter(paths(k),[root filesep]); dst=fullfile(out,'provenance','source_snapshot',rel);
 if ~isfolder(fileparts(dst)), mkdir(fileparts(dst)); end
 copyfile(paths(k),dst);
end
writetable(table(paths',hashes','VariableNames',{'path','sha256'}),fullfile(out,'provenance','source_hashes.csv'));
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root)); [~,status]=system(sprintf('git -C "%s" status --short',root));
[~,diff]=system(sprintf('git -C "%s" diff',root)); [~,staged]=system(sprintf('git -C "%s" diff --cached',root));
c20_v4_io('json',fullfile(out,'provenance','repo.json'),struct('head',head,'status',status,'unstaged',diff,'staged',staged));
env=evalc('disp(version); disp(ver); which qe_compare_component_models -all; which read_npy -all; which qe_centered_bins -all; which peak_models -all;');
c20_v4_io('text',fullfile(out,'provenance','environment.txt'),env);
r=[runtests(fullfile(root,'tests','C20V4ContractsTest.m')),runtests(fullfile(root,'tests','test_b1_component_pilot_v2.m'))];
writetable(table(r),fullfile(out,'tests','results.csv')); save(fullfile(out,'tests','results.mat'),'r');
assert(all([r.Passed]),'Tests failed; execution stopped'); disp('STAGE tests passed');
s=load(fullfile(p2,'590_PL2_10w','fit_details.mat')); parentstats=c20_v4_export_fits(s.details,fullfile(out,'audit'));
assert(parentstats.models==108&&parentstats.candidates==2592&&parentstats.selected_boundary_models==98&&parentstats.boundary_flag_mismatches==0);
disp(parentstats);
l=load(fullfile(p2,'590_PL2_10w','L1_minimal.mat')); qe=l.d.qe;
bins=qe_centered_bins(qe,cfg.targets,cfg.Ns); assert(numel(bins)==9&&all([bins.valid]));
E=qe.energy_meV; members=unique([bins.source_channel]); member_spectra=qe.intensity(:,ismember(qe.source_channel,members));
save(fullfile(out,'590_PL2_10w','centered_binned_spectra.mat'),'bins','E','members','member_spectra','-v7');
bt=struct2table(rmfield(bins,{'sum','mean','variance_mean','valid_mask','indices','source_channel','q_members'}));
bt.native_members=string(arrayfun(@(b)mat2str(b.source_channel),bins,'UniformOutput',false)).';
writetable(bt,fullfile(out,'590_PL2_10w','centered_bins.csv'));
for j=1:3
 f=figure('Visible','off','Position',[100 100 1000 550]); hold on;
 for k=1:3, b=bins((j-1)*3+k); plot(E,b.mean,'DisplayName',sprintf('N=%d',b.N)); end
 xlim([300 1800]); legend('show'); xlabel('Energy (meV)'); ylabel('Sequence-summed q mean'); title(sprintf('Same-center q=%+.4f; overlapping data across N',cfg.targets(j)));
 exportgraphics(f,fullfile(out,'590_PL2_10w','figures',sprintf('same_center_R%d.png',j))); close(f);
end
legacy={}; independent={}; background={};
for mode=1:3
 dlist={};
 for j=1:3
  b=bins((j-1)*3+2); u=struct('bin_size_requested',3,'source_q_count',3,'q_Ainv',b.q_Ainv, ...
   'q_left',b.q_left,'q_right',b.q_right,'source_channel',b.source_channel);
  for mi=1:2
   key=sprintf('R%d_%s_mode%d',j,cfg.models{mi},mode); policy='legacy'; bg='power_law'; extra=0;
   if mode>1, policy='independent'; extra=cfg.extra_n2_starts; end
   if mode==3, bg=cfg.backgrounds{2}; end
   fits=qe_compare_component_models(E,b.mean,peak_model=cfg.models{mi},seed=cfg.seed+j, ...
    n_starts=24,start_policy=policy,extra_starts=extra,baseline_mode=bg);
   d=struct('key',key,'unit',u,'fits',fits,'role',mode,'bin',b,'data_hash',c20_v4_io('digest',jsonencode(b.mean)));
   save(fullfile(out,'checkpoints',[key '.mat']),'d','-v7'); dlist{end+1}=d; %#ok<AGROW>
   c20_v4_plot_fits(d,fullfile(out,'590_PL2_10w','figures',[key '.png']));
   disp(['FIT ' key]);
  end
 end
 if mode==1
  legacy=dlist; details=legacy; stats=c20_v4_export_fits(details,fullfile(out,'590_PL2_10w'));
  assert(stats.models==12&&stats.components==18&&stats.candidates==288);
  save(fullfile(out,'590_PL2_10w','fit_details.mat'),'details','-v7');
 elseif mode==2
  independent=dlist; c20_v4_export_fits(dlist,fullfile(out,'590_PL2_10w','independent_starts'));
 else
  background=dlist; c20_v4_export_fits(dlist,fullfile(out,'590_PL2_10w','background_B1'));
 end
 save(fullfile(out,'590_PL2_10w','solver_background_comparison.mat'),'legacy','independent','background','-v7');
end
manifest=jsondecode(fileread(fullfile(p2,'input_manifest.resolved.yaml')));
rec=manifest.sessions(strcmp({manifest.sessions.session_id},'590_PL2_10w')); assert(isscalar(rec));
frames=c20_v4_frames(rec,qe,bins,fullfile(out,'590_PL2_10w')); %#ok<NASGU>
parent_after={c20_v4_io('inventory',p2),c20_v4_io('inventory',p3)};
assert(isequal(parent_before,parent_after),'Parent contents changed'); save(fullfile(out,'provenance','parent_hashes_after.mat'),'parent_after');
c20_v4_io('json',fullfile(out,'stage_status.json'),struct('parent_mapping','passed', ...
 'legacy_reconstruction','12 models saved','centered_bins','9 valid','A1','tested_applied', ...
 'sequence_semantics','camera sequence; fixed position unknown', ...
 'numerical_checks','witness verified; inspect optimized violations','background_B1','12 models saved', ...
 'model_adequacy','not_assessed','identifiability','not_assessed','physical_origin','not_assessed','P5','not_run_yet'));
disp(['CORE_COMPLETE=' out]);
end
