function out = run_b1_data_constrained_pilot_v5()
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
parent=fullfile(root,'paper_results','b1_components_v4','20260912T045417369Z_2f64ff5c_4e9fdb89');
assert(strcmp(c20_v4_io('hash',fullfile(parent,'review_packet.zip')),'77dfe0c21292212469a0063406cbe279972f4e2817226e5466d5659a9fe9388c'));
cfg=struct('scope','data_constrained_pilot','window',[300 1800],'reference_window',[300 1800], ...
 'A1_starts',24,'A1_extra_n2',12,'member_starts',12,'seed',20260912, ...
 'models',{{'lorentz_symmetric','lorentz'}},'member_alignment','A0', ...
 'member_assumptions','linear center at five q; region-shared width/r; per-member nonnegative B and A');
files=[dir(fullfile(root,'case_studies','bisb2026','scripts','c20_v5*.m')); dir(fullfile(root,'case_studies','bisb2026','scripts','run_b1_data_constrained_pilot_v5.m')); ...
 dir(fullfile(root,'src','fitting','qe_*member*.m')); dir(fullfile(root,'src','fitting','qe_area_scopes.m')); dir(fullfile(root,'tests','C20V5Test.m'))];
paths=string(fullfile({files.folder},{files.name}));
paths=[paths,string(fullfile(root,{'src/fitting/qe_compare_component_models.m','src/fitting/qe_component_mapping.m', ...
 'src/fitting/qe_component_prediction.m','src/fitting/peak_models.m','src/fitting/measure_peak_fwhm.m', ...
 'src/qe_centered_bins.m','src/qe_prepare_count_bins.m','src/qe_zlp_integer_align.m','src/io/read_npy.m', ...
 'case_studies/bisb2026/scripts/c20_v4_export_fits.m','case_studies/bisb2026/scripts/c20_v4_io.m', ...
 'case_studies/bisb2026/scripts/c20_v4_plot_fits.m','tests/C20V4ContractsTest.m','tests/test_b1_component_pilot_v2.m'}))];
hashes=strings(size(paths)); for k=1:numel(paths), hashes(k)=c20_v4_io('hash',paths(k)); end
codehash=c20_v4_io('digest',char(join(paths+hashes,'|'))); cfghash=c20_v4_io('digest',jsonencode(cfg));
id=[char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z''')) '_' cfghash(1:8) '_' codehash(1:8)];
out=fullfile(root,'paper_results','b1_components_v5',id); assert(~isfolder(out));
for name={'provenance','audit','A0_A1','A0_A1/figures','member_models','sequence_models','q_center_diagnostics','validation','tests'}, mkdir(fullfile(out,name{1})); end
diary(fullfile(out,'provenance','execution.log')); dc=onCleanup(@()diary('off')); disp(['RUN_DIR=' out]);
c20_v4_io('json',fullfile(out,'config_resolved.json'),cfg);
c20_v4_io('json',fullfile(out,'parent_packet_reference.json'),struct('parent',parent,'zip_sha256',c20_v4_io('hash',fullfile(parent,'review_packet.zip'))));
c20_v4_io('json',fullfile(out,'provenance','acquisition_user_confirmation.json'),struct( ...
 'date','2026-09-12','source','direct user confirmation in current task', ...
 'same_region',true,'continuous_acquisition',true,'scan_during_sequence',false,'region_change',false,'beam_adjustment',false, ...
 'stationarity','not implied by same-region acquisition','independent_frames','not established'));
pp=jsondecode(fileread(fullfile(parent,'provenance','parents.json')));
parents={pp.v2,pp.v3,parent}; parent_hashes=cellfun(@(p)c20_v4_io('inventory',p),parents,'UniformOutput',false);
save(fullfile(out,'provenance','parent_hashes_before.mat'),'parents','parent_hashes');
for k=1:numel(paths)
 rel=extractAfter(paths(k),[root filesep]); dst=fullfile(out,'provenance','source_snapshot',rel);
 if ~isfolder(fileparts(dst)), mkdir(fileparts(dst)); end
 copyfile(paths(k),dst);
end
writetable(table(paths',hashes','VariableNames',{'path','sha256'}),fullfile(out,'provenance','source_hashes.csv'));
[~,status]=system(sprintf('git -C "%s" status --short',root)); [~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root));
[~,diff]=system(sprintf('git -C "%s" diff',root)); [~,staged]=system(sprintf('git -C "%s" diff --cached',root));
c20_v4_io('json',fullfile(out,'provenance','git.json'),struct('status',status,'head',head,'unstaged',diff,'staged',staged));
c20_v4_io('text',fullfile(out,'provenance','environment.txt'),evalc('disp(version); which qe_fit_member_models -all; which peak_models -all; which qe_area_scopes -all;'));
manifest=readtable(fullfile(parent,'review_packet','FILE_MANIFEST.csv'),TextType='string');
used={'590_PL2_10w/sequence_block_spectra.mat','590_PL2_10w/solver_background_comparison.mat','590_PL2_10w/centered_binned_spectra.mat'};
for k=1:numel(used)
 ix=replace(manifest.path,char(92),'/')==used{k}; assert(nnz(ix)==1);
 assert(strcmpi(c20_v4_io('hash',fullfile(parent,used{k})),manifest.sha256(ix)));
end
writetable(manifest(ismember(replace(manifest.path,char(92),'/'),string(used)),:),fullfile(out,'provenance','input_hashes.csv'));
r=[runtests(fullfile(root,'tests','C20V5Test.m')),runtests(fullfile(root,'tests','C20V4ContractsTest.m')),runtests(fullfile(root,'tests','test_b1_component_pilot_v2.m'))];
save(fullfile(out,'tests','results.mat'),'r'); writetable(table(r),fullfile(out,'tests','results.csv')); assert(all([r.Passed]));
old=load(fullfile(pp.v2,'590_PL2_10w','fit_details.mat')); areas={}; area_records={};
for d=old.details
 d=d{1};
 for f=d.fits
  model=peak_models(f.peak_model); curves=zeros(size(f.components));
  for j=1:f.n_components, p=f.parameters(j,:); curves(:,j)=model.model_fn(p(1),p(2),p(3),f.energy_meV); end
  a=qe_area_scopes(f.energy_meV,curves,cfg.reference_window);
  for j=1:f.n_components
   areas(end+1,:)={string(d.key),f.n_components,j,string(f.peak_model),f.parameters(j,3),a.area_fit_window(j),a.area_reference_window(j), ...
    string(mat2str(a.fit_window_meV)),string(mat2str(a.reference_window_meV)),string(mat2str(a.reference_actual_support_meV)),a.reference_window_fully_observed,string(a.integration_method)}; %#ok<AGROW>
  end
  area_records{end+1}=struct('key',d.key,'n',f.n_components,'parameters',f.parameters,'model',f.peak_model,'E',f.energy_meV,'components',curves,'areas',a); %#ok<AGROW>
 end
end
writetable(cell2table(areas,'VariableNames',{'key','n_components','component','peak_model','native_A','area_fit_window','area_reference_window', ...
 'fit_window_meV','reference_window_meV','actual_reference_support','reference_window_fully_observed','integration_method'}),fullfile(out,'audit','area_scope_corrected.csv'));
save(fullfile(out,'audit','area_arrays.mat'),'area_records','-v7'); disp('AREA_SCOPE_CORRECTED_NO_REFIT');
a=load(fullfile(parent,'590_PL2_10w','sequence_block_spectra.mat')); frames=a.frames;
a=load(fullfile(parent,'590_PL2_10w','solver_background_comparison.mat')); A0=a.independent;
save(fullfile(out,'A0_A1','A0_existing_mode2.mat'),'A0','-v7'); A1={}; paircheck={};
for j=1:3
 for mi=1:2
  base=A0{2*(j-1)+mi}; f0=base.fits(1); E=frames.E; ix=E>=cfg.window(1)&E<=cfg.window(2);
  Y0=sum(frames.sequence_bin_A0(:,:,j),2); Y1=sum(frames.sequence_bin_A1(:,:,j),2);
  assert(isequal(E(ix),f0.energy_meV)); err=max(abs(Y0(ix)-f0.observed)); assert(err<1e-8);
  fits=qe_compare_component_models(E,Y1,energy_window=cfg.window,peak_model=cfg.models{mi},n_starts=cfg.A1_starts, ...
   extra_starts=cfg.A1_extra_n2,start_policy='independent',seed=cfg.seed+j);
  d=struct('key',sprintf('A1_R%d_%s',j,cfg.models{mi}),'unit',base.unit,'fits',fits, ...
   'A0',Y0,'A1',Y1,'E',E,'input_hash',c20_v4_io('digest',jsonencode(Y1)),'code_hash',codehash,'config_hash',cfghash);
  save(fullfile(out,'A0_A1',[d.key '.mat']),'d','-v7'); A1{end+1}=d; %#ok<AGROW>
  paircheck(end+1,:)={j,string(cfg.models{mi}),err}; %#ok<AGROW>
  c20_v4_plot_fits(d,fullfile(out,'A0_A1','figures',[d.key '.png'])); disp(['FIT ' d.key]);
 end
end
save(fullfile(out,'A0_A1','fit_details.mat'),'A1','-v7'); c20_v4_export_fits(A1,fullfile(out,'A0_A1'));
writetable(cell2table(paircheck,'VariableNames',{'region','model','A0_input_error'}),fullfile(out,'A0_A1','input_pair_checks.csv'));
allmembers={};
for j=1:3
 idx=(j-1)*5+(1:5); X=squeeze(sum(frames.member_spectra_A0(frames.support,:,idx),2)); q=frames.q_members(idx); members=frames.native_members(idx);
 for mi=1:2
  fits=qe_fit_member_models(frames.E,q,X,energy_window=cfg.window,peak_model=cfg.models{mi},n_starts=cfg.member_starts, ...
   seed=cfg.seed+j,initial_fits=A0{2*(j-1)+mi}.fits);
  d=struct('key',sprintf('R%d_%s',j,cfg.models{mi}),'fits',fits,'members',members,'region',j, ...
   'input_hash',c20_v4_io('digest',jsonencode(X)),'code_hash',codehash,'config_hash',cfghash,'alignment','A0');
  c20_v5_export_members(d,fullfile(out,'member_models',d.key)); allmembers{end+1}=d; %#ok<AGROW>
  disp(['MEMBER_FIT ' d.key]);
 end
end
save(fullfile(out,'member_models','all_member_fits.mat'),'allmembers','-v7');
save(fullfile(out,'sequence_models','sequence_inputs.mat'),'frames','-v7');
for k=1:numel(parents), assert(isequal(parent_hashes{k},c20_v4_io('inventory',parents{k})),'Parent changed'); end
c20_v4_io('json',fullfile(out,'stage_status.json'),struct('area','corrected_without_refit','A1','12 models saved','members','12 M1/M2 saved','q_center','not_run_yet','controls','not_run_yet','inference','not_calibrated'));
disp(['CORE_COMPLETE=' out]);
end
