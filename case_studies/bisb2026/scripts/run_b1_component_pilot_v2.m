function out = run_b1_component_pilot_v2(options)
%RUN_B1_COMPONENT_PILOT_V2 Read-only source audit and P0-P3 component pilot.
% Does not execute historical wrappers, physical fits, or update old indexes.
arguments
 options.mode string {mustBeMember(options.mode,["pilot","validate_inputs_only"])} = "pilot"
 options.output_root string = ""
end
root=bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(genpath(fullfile(root,'lib')));
cfg=b1_component_config_v2();
[~,head]=system('git rev-parse HEAD'); head=strtrim(head);
cfg_hash=digest(uint8(jsonencode(cfg)));
if strlength(options.output_root)==0
 id=[char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z''')) '_' cfg_hash(1:8) '_' head(1:8)];
 out=fullfile(root,'paper_results','b1_components_v2',id);
else, out=char(options.output_root); end
assert(~isfolder(out)&&~isfile(out),'run_b1_component_pilot_v2:Exists','Output exists; refusing overwrite');
mkdir(out); mkdir(fullfile(out,'logs')); mkdir(fullfile(out,'tests'));
diary(fullfile(out,'logs','execution.txt')); diary_cleanup=onCleanup(@()diary('off'));
fprintf('RUN_DIR=%s\n',out);
writejson(fullfile(out,'config_resolved.json'),cfg); save(fullfile(out,'config_resolved.mat'),'cfg');
[~,status]=system('git status --short --branch'); [~,difftext]=system('git diff');
writejson(fullfile(out,'repo_state.json'),struct('head',head,'status',status,'diff',difftext));
env=evalc('disp(version); disp(struct2table(ver)); which load_raw_session -all; which load_qe_dataset -all; which run_thesis_pipeline -all');
writefile(fullfile(out,'environment.txt'),env);
source_files={'src/io/load_raw_session.m','src/io/load_qe_dataset.m','src/io/read_npy.m', ...
 'src/io/make_qe_struct.m','src/qe_prepare_count_bins.m','src/fitting/peak_models.m', ...
 'src/fitting/qe_compare_component_models.m','src/fitting/measure_peak_fwhm.m', ...
 'src/thesis/thesis_sessions.m','case_studies/bisb2026/scripts/b1_component_config_v2.m', ...
 'case_studies/bisb2026/scripts/run_b1_component_pilot_v2.m','tests/test_b1_component_pilot_v2.m'};
code_manifest=table();
for fi=1:numel(source_files)
 source=fullfile(root,source_files{fi}); target=fullfile(out,'source_snapshot',source_files{fi});
 parent=fileparts(target); if ~isfolder(parent), mkdir(parent); end
 copyfile(source,target);
 code_manifest=[code_manifest;table(string(source_files{fi}),string(filehash(source)), ...
  'VariableNames',{'path','sha256'})]; %#ok<AGROW>
end
writetable(code_manifest,fullfile(out,'source_code_hashes.csv'));
r=runtests(fullfile(root,'tests','test_b1_component_pilot_v2.m'));
save(fullfile(out,'tests','pilot_tests.mat'),'r');
writetable(table(r),fullfile(out,'tests','pilot_tests.csv'));
assert(all([r.Passed]),'Pilot tests failed; inputs and fits not executed');
migration(root,out);
tags={'590_gui_history_area_260506','no_PL2_20w_2film_gui_history_area_260506_highq_refined','n0_PL2_10w_gui_history_area_260506'};
manifest=struct('schema','JSON-compatible YAML 1.2','sessions',struct([]));
pilot=[];
for si=1:numel(cfg.sessions)
 s=cfg.sessions(si); dest=fullfile(out,s.name); mkdir(dest);
 raw=dir(fullfile(s.path,'*.npy')); assert(numel(raw)==1,'Expected exactly one registered raw NPY');
 rawpath=fullfile(raw.folder,raw.name); jsonpath=replace(rawpath,'.npy','.json');
 meta=jsondecode(fileread(jsonpath));
 assert(strcmp(meta.spatial_calibrations(end).units,'eV') && meta.spatial_calibrations(end).scale>0,'Unknown energy units');
 rec=struct('session_id',s.name,'registered_dq_Ainv',s.dq_Ainv,'files',struct([]));
 paths={rawpath,jsonpath,fullfile(s.path,'eq3D.mat'),fullfile(s.path,'eq3D_processed.mat'), ...
 fullfile(s.path,'op_history_260506.mat'),fullfile(root,'paper_results',tags{si},'analysis_results.mat'), ...
 fullfile(root,'paper_results',tags{si},'branch1_points.csv')};
 for pi=1:numel(paths)
  item=struct('path',paths{pi},'exists',isfile(paths{pi}),'sha256','','fields',{{}});
  if item.exists
   item.sha256=filehash(paths{pi});
   if endsWith(paths{pi},'.mat'), w=whos('-file',paths{pi}); item.fields={w.name}; end
  end
  if pi==1, rec.files=item; else, rec.files(pi)=item; end
 end
 historical=load(paths{6},'output'); h=historical.output;
 points=readtable(paths{7});
 rec.historical=struct('q_dq',median(diff(h.qe_pp.q_Ainv)), ...
 'q_zero_index',h.qe_pp.q_zero_index,'q_range',[min(h.qe_pp.q_Ainv),max(h.qe_pp.q_Ainv)], ...
 'fields',{fieldnames(h)},'csv_fields',{points.Properties.VariableNames}, ...
 'csv_q_range',[min(points.q_Ainv),max(points.q_Ainv)], ...
 'csv_on_registered_grid',all(abs(points.q_Ainv/s.dq_Ainv-round(points.q_Ainv/s.dq_Ainv))<1e-7));
 writejson(fullfile(dest,'historical_snapshot.json'),struct('preprocess',h.preprocess_opts,'snap',h.snap));
 % Full native q axis; metadata dimensions checked before importer use.
 nq=meta.metadata.hardware_source.camera_processing_parameters.sensor_dimensions(1);
 d=load_raw_session(rawpath,q_crop=[1 nq],dq_Ainv=s.dq_Ainv, ...
     write_cache=false,align_zlp=cfg.align_zlp,show_progress=false,invalid_count_threshold=cfg.invalid_count_threshold);
 qe=d.qe;
 assert(size(qe.intensity,2)==nq && abs(median(diff(qe.energy_meV))-1000*meta.spatial_calibrations(end).scale)<1e-9,'Raw axis identity mismatch');
 rec.actual=struct('read_file',rawpath,'level','L1','shape',size(qe.intensity), ...
  'energy_range_meV',[min(qe.energy_meV) max(qe.energy_meV)],'dE_meV',median(diff(qe.energy_meV)), ...
  'dq_Ainv',median(diff(qe.q_Ainv)),'q_zero_native_channel',qe.source_channel(qe.q_zero_index), ...
  'q_range',[min(qe.q_Ainv) max(qe.q_Ainv)],'source_channel',qe.source_channel, ...
  'transform','sum sequence frames; transpose to E x q; global ZLP energy zero; native q zero via make_qe_struct; no per-q alignment', ...
  'legacy_to_native_column_mapping','not established; no legacy arrays rescaled or used as pilot input', ...
  'noise','corrected detector counts; covariance unknown; unweighted descriptive LS; variance not assumed Poisson', ...
  'response','independent reference ZLP unavailable; no intrinsic width/lifetime', ...
  'mask','q_skip 0.0005 plus entire q channel exclusion for any saturation sentinel/nonfinite; wider empirical beam mask not established', ...
  'invalid_native_q',d.invalid_native_q);
 [~,zlp_by_q]=max(qe.intensity,[],1);
 qc=struct('frame',d.frame_diagnostics,'zlp_energy_meV_by_q',qe.energy_meV(zlp_by_q), ...
  'nonfinite_per_q',sum(~isfinite(qe.intensity),1),'negative_samples',nnz(qe.intensity<0));
 writejson(fullfile(dest,'raw_quality.json'),qc);
 save(fullfile(dest,'L1_minimal.mat'),'d','-v7.3');
 Ns=cfg.native_N; if s.dq_Ainv==.00025, Ns=[Ns cfg.matched_20w_N]; end
 allbins=cell(size(Ns)); binrows=table();
 for ni=1:numel(Ns)
  b=qe_prepare_count_bins(qe,Ns(ni),q_range_Ainv=cfg.q_range_Ainv, ...
      q_skip_Ainv=cfg.q_skip_Ainv,processing_level="L1");
  allbins{ni}=b;
  for k=1:numel(b.units)
   u=b.units(k); row=table(Ns(ni),k,u.q_Ainv,u.q_left,u.q_right,u.q_width,u.source_q_count,u.partial, ...
    string(u.source_q_index),string(u.source_q_Ainv),string(mat2str(u.source_channel)), ...
    'VariableNames',{'N_requested','bin_id','q_Ainv','q_left','q_right','q_width','source_q_count','partial','source_q_index','source_q_Ainv','native_channels'});
   binrows=[binrows;row]; %#ok<AGROW>
  end
 end
 writetable(binrows,fullfile(dest,'bins.csv')); save(fullfile(dest,'binned_spectra.mat'),'allbins','Ns','-v7.3');
 rec.bin_count=height(binrows);
 if si==1, manifest.sessions=rec; else, manifest.sessions(si)=rec; end
 if si==1, pilot=struct('qe',qe,'bins',{allbins},'dest',dest); end
 fprintf('P0-P2 input audited: %s; %d bins\n',s.name,height(binrows));
end
writejson(fullfile(out,'input_manifest.resolved.yaml'),manifest);
writejson(fullfile(out,'data_lineage.json'),manifest);
if options.mode=="pilot", fitpilot(pilot,cfg); end
% Re-hash every input after execution to verify source read-only behavior.
unchanged=true;
for si=1:numel(manifest.sessions)
 for item=manifest.sessions(si).files
  if item.exists, unchanged=unchanged && strcmp(filehash(item.path),item.sha256); end
 end
end
assert(unchanged,'Source hash changed during run');
writefile(fullfile(out,'source_integrity.txt'),'All listed existing source hashes unchanged after run.');
writefile(fullfile(out,'DECISIONS.md'),sprintf(['P0-P3 pilot only. Bi:Sb approximately 70:30; complete MoS2-metal-MoS2 stack.\n' ...
 'No physical fit, trend prior, Fano, weak-peak deletion, jump repair or historical reference reconstruction.\n' ...
 'New L1 uses native frame sums without per-q alignment. Frame drift and q-dependent ZLP offsets are saved, not corrected.\n' ...
 'Unknown detector covariance: variance_sum/mean are NaN; member scatter is separately named.\n' ...
 'All fit classifications remain unresolved_or_invalid pending P4-P5. No confidence intervals or intrinsic lifetimes.\n' ...
 'Three sessions were input-audited and binned; only 590 was fitted. This is not frozen cross-session validation.\n']));
fprintf('COMPLETED %s\n',out);
if options.mode=="pilot", write_b1_component_pilot_report_v2(out); end
end

function fitpilot(pilot,cfg)
dest=pilot.dest; mkdir(fullfile(dest,'figures')); comparisons=table(); parameters=table(); candidates=table(); selected=table(); details={};
for ni=1:numel(cfg.native_N)
 b=pilot.bins{ni}; qs=[b.units.q_Ainv];
 for ri=1:3
  allowed=find(sign(qs)==sign(cfg.representative_q_Ainv(ri)) & ~[b.units.partial]);
  assert(~isempty(allowed),'No valid full representative bin');
  [~,j]=min(abs(qs(allowed)-cfg.representative_q_Ainv(ri))); bi=allowed(j); u=b.units(bi);
  selected=[selected;table(cfg.native_N(ni),ri,bi,cfg.representative_q_Ainv(ri),u.q_Ainv,string(u.source_q_index), ...
   'VariableNames',{'N','representative','bin_id','target_q','actual_q','members'})]; %#ok<AGROW>
  f=figure('Visible','off','Position',[100 100 950 520]);
  plot(pilot.qe.energy_meV,pilot.qe.intensity(:,u.q_indices),'Color',[.7 .7 .7]); hold on;
  plot(b.energy_meV,b.mean(:,bi),'k','LineWidth',1.5); xlim([300 2100]);
  xlabel('Energy loss (meV)'); ylabel('Corrected detector counts');
  title(sprintf('L1 unaligned members + mean; q=%+.4f; N=%d; width=%.4g A^{-1}',u.q_Ainv,u.source_q_count,u.q_width));
  exportgraphics(f,fullfile(dest,'figures',sprintf('members_N%d_R%d.png',cfg.native_N(ni),ri))); close(f);
  for wi=1:size(cfg.energy_windows_meV,1)
   for mi=1:numel(cfg.models)
    fits=qe_compare_component_models(b.energy_meV,b.mean(:,bi), ...
     energy_window=cfg.energy_windows_meV(wi,:),peak_model=cfg.models{mi},n_starts=cfg.n_starts,seed=cfg.seed);
    key=sprintf('N%d_R%d_W%d_%s',cfg.native_N(ni),ri,cfg.energy_windows_meV(wi,2),cfg.models{mi});
    details{end+1}=struct('key',key,'unit',u,'fits',fits); %#ok<AGROW>
    for n=1:2
     z=fits(n); cnd=z.candidates;
     boundary=false; condition=NaN;
     if z.success, boundary=any(cnd(z.selected_start).boundary); condition=cnd(z.selected_start).jacobian_condition; end
     comparisons=[comparisons;table(string(key),n,z.success,string(z.numerical_status),z.sse,z.normalized_sse,z.selected_start,boundary,condition, ...
      string(z.scientific_status),'VariableNames',{'key','n','success','numerical_status','sse','normalized_sse','selected_start','boundary','jacobian_condition','scientific_status'})]; %#ok<AGROW>
     for c=cnd
      candidates=[candidates;table(string(key),n,c.start,c.exitflag,c.objective,any(c.boundary),c.jacobian_condition,string(mat2str(c.p0)),string(mat2str(c.p)),string(c.message), ...
       'VariableNames',{'key','n','start','exitflag','normalized_sse','boundary','jacobian_condition','p0_scaled','p_scaled','message'})]; %#ok<AGROW>
     end
     if z.success
      for p=1:n
       pars=z.parameters(p,:); diag=z.fwhm_diagnostic(p);
       parameters=[parameters;table(string(key),n,p,u.q_Ainv,pars(1),pars(2),pars(3),diag.apex_meV,diag.fwhm_meV,trapz(z.energy_meV,z.components(:,p)),string(diag.status), ...
        'VariableNames',{'key','n','component','q_Ainv','E0_meV','native_width_meV','native_amplitude','apex_meV','finite_window_floor_FWHM_meV','finite_window_area','fwhm_status'})]; %#ok<AGROW>
      end
     end
    end
    f=figure('Visible','off','Position',[100 100 1100 780]); tiledlayout(2,2);
    for n=1:2
     z=fits(n); nexttile(n); plot(z.energy_meV,z.observed,'Color',[.65 .65 .65]); hold on;
     plot(z.energy_meV,z.prediction,'k','LineWidth',1.4); plot(z.energy_meV,z.background,'--'); plot(z.energy_meV,z.components);
     title(sprintf('n=%d | %s',n,z.numerical_status),'Interpreter','none'); ylabel('Detector counts (mean)');
     nexttile(n+2); plot(z.energy_meV,z.residual); yline(0,'k:'); xlabel('Energy loss (meV)'); ylabel('Observed - model');
    end
    sgtitle(sprintf('L1 | %s | q=%+.4f N=%d width=%.4g | power-law BG | unresolved',key,u.q_Ainv,u.source_q_count,u.q_width),'Interpreter','none');
    exportgraphics(f,fullfile(dest,'figures',[key '.png'])); close(f);
    fprintf('FIT %s n1=%s n2=%s\n',key,fits(1).numerical_status,fits(2).numerical_status);
   end
  end
 end
end
writetable(selected,fullfile(dest,'representative_selection.csv'));
writetable(comparisons,fullfile(dest,'model_comparison.csv'));
writetable(comparisons(comparisons.n==1,:),fullfile(dest,'single_component_fits.csv'));
writetable(candidates(candidates.n==2,:),fullfile(dest,'double_component_candidates.csv'));
writetable(candidates,fullfile(dest,'all_multistart_candidates.csv'));
writetable(parameters,fullfile(dest,'component_parameters.csv'));
writetable(comparisons(~comparisons.success | comparisons.boundary,:),fullfile(dest,'fit_failures.csv'));
save(fullfile(dest,'fit_details.mat'),'details','-v7.3');
end

function migration(root,out)
files={'case_studies/bisb2026/scripts/run_b1_double_peak_binning_analysis.m','src/b1_double_peak_binning_extract.m', ...
 'src/b1_peak_evidence_audit_classify.m'}; rows=table();
for k=1:numel(files)
 lines=splitlines(string(fileread(fullfile(root,files{k}))));
 for i=1:numel(lines)
  if ~isempty(regexpi(lines(i),'(qRange|q_range|q_skip|lowQNoBin|low_q_no_bin|highQForceBin|high_q_force_bin|fitDenoiseQ|fit_denoise_q|tracking.*q.*Ainv|trend.*q.*Ainv|reference.*q|max_q_gap|small.*q)'))
   action="inactive in pilot; historical value preserved";
   if contains(lines(i),{'q_range','qRange'}), action="pilot explicit [-0.015,0.015]"; end
   if contains(lines(i),'q_skip'), action="pilot explicit 0.0005"; end
   if contains(lines(i),{'lowQNoBin','low_q_no_bin'}), action="fixed-N ignores adaptive protection; center mask retained"; end
   rows=[rows;table(string(files{k}),i,strtrim(lines(i)),action,string('native channel offsets via registry dq; no legacy point rescaling'), ...
    'VariableNames',{'file','line','legacy_rule','pilot_action','channel_mapping_policy'})]; %#ok<AGROW>
  end
 end
end
writetable(rows,fullfile(out,'q_rule_migration.csv'));
end
function writejson(path,value), writefile(path,jsonencode(value,PrettyPrint=true)); end
function writefile(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',value);
end
function h=filehash(path)
md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(path,'rb'); assert(fid>=0); c=onCleanup(@()fclose(fid));
while ~feof(fid), x=fread(fid,8*1024*1024,'*uint8'); md.update(typecast(x,'int8')); end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function h=digest(bytes)
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(bytes,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
