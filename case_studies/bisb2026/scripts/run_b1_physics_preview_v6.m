function out=run_b1_physics_preview_v6()
% Frozen sparse physics preview, not another estimator-development pipeline.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
p5=fullfile(root,'paper_results','b1_components_v5','20260912T055651568Z_c1dfab1b_fcd805b3');
cfg=struct('window',[300 1800],'N',3,'starts',24,'extra_n2',12,'seed',20260912, ...
 'targets590',[-.0125 -.0095 -.0075 -.0045 -.0025 .0025 .0045 .0075 .0095 .0125], ...
 'targetsRepeat',[-.0025 .0025 .0075 .0125],'models',{{'lorentz_symmetric','lorentz'}},'frozen_reference_abs_q',.0075);
assert(strcmp(c20_v4_io('hash',fullfile(p5,'delivery_final','review_packet.zip')),'cddcfa5ff41bbd6263e2d5988c0ce37c282f7f1e36dcabf25cc0993f9aaa515a'));
source={mfilename('fullpath'),which('qe_compare_component_models'),which('qe_centered_bins'),which('qe_zlp_integer_align'), ...
 which('qe_spectral_centroids'),which('qe_area_scopes'),which('peak_models'),which('read_npy'), ...
 which('c20_v4_io'),which('c20_v4_export_fits'),which('c20_v4_plot_fits'),which('thesis_sessions'),which('qe_prepare_count_bins'),which('qe_component_mapping'),which('measure_peak_fwhm')}; source{1}=[source{1} '.m'];
hashes=cellfun(@(p)c20_v4_io('hash',p),source,'UniformOutput',false);
codehash=c20_v4_io('digest',strjoin(hashes,'|')); cfghash=c20_v4_io('digest',jsonencode(cfg));
out=fullfile(root,'paper_results','b1_physics_preview',[char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z''')) '_' cfghash(1:8) '_' codehash(1:8)]);
assert(~isfolder(out)); mkdir(out); mkdir(fullfile(out,'appendix')); mkdir(fullfile(out,'figures'));
diary(fullfile(out,'appendix','execution.log')); cleanup=onCleanup(@()diary('off')); disp(['RUN_DIR=' out]);
c20_v4_io('json',fullfile(out,'config_resolved.json'),cfg);
before={c20_v4_io('inventory',p2),c20_v4_io('inventory',p5)};
save(fullfile(out,'appendix','parent_hashes_before.mat'),'before','p2','p5');
for k=1:numel(source), [~,name,ext]=fileparts(source{k}); copyfile(source{k},fullfile(out,'appendix',[name ext])); end
writetable(table(string(source)',string(hashes)','VariableNames',{'source','sha256'}),fullfile(out,'appendix','source_hashes.csv'));
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root)); [~,status]=system(sprintf('git -C "%s" status --short',root));
c20_v4_io('json',fullfile(out,'appendix','repo.json'),struct('head',head,'status',status));
r=runtests(fullfile(root,'tests','C20PhysicsPreviewTest.m')); writetable(table(r),fullfile(out,'appendix','tests.csv')); assert(all([r.Passed]));
manifest=jsondecode(fileread(fullfile(p2,'input_manifest.resolved.yaml'))); registry=thesis_sessions();
old=load(fullfile(p5,'A0_A1','fit_details.mat')); reuse=old.A1;
f=load(fullfile(p5,'sequence_models','sequence_inputs.mat')); saved590=f.frames;
session_ids={'590_PL2_10w','n0_PL2_10w_repeat'}; products=cell(1,2);
for si=1:2
 sid=session_ids{si}; reg=registry(strcmp({registry.name},sid)); rec=manifest.sessions(strcmp({manifest.sessions.session_id},sid));
 targets=cfg.targets590; if si==2, targets=cfg.targetsRepeat; end
 rawfile=rec.files(find(endsWith(string({rec.files.path}),'.npy'),1)); jsonfile=rec.files(find(endsWith(string({rec.files.path}),'.json'),1));
 if ~isfile(rawfile.path), c20_v4_io('text',fullfile(out,[sid '_missing.txt']),'Native input missing; no fallback fit'); continue; end
 assert(strcmpi(c20_v4_io('hash',rawfile.path),rawfile.sha256)); assert(strcmpi(c20_v4_io('hash',jsonfile.path),jsonfile.sha256));
 assert(strcmp(fileparts(rawfile.path),reg.path));
 l=load(fullfile(p2,sid,'L1_minimal.mat')); qe=l.d.qe;
 assert(abs(qe.dq_Ainv-reg.dq_Ainv)<1e-10);
 raw=read_npy(rawfile.path); E=qe.energy_meV; T=size(raw,1); nq=size(raw,2);
 reference_idx=find(abs(qe.q_Ainv)<=.001 & ~ismember(qe.source_channel,rec.actual.invalid_native_q));
 profiles=squeeze(sum(double(raw(:,reference_idx,:)),2)).';
 if si==1
  align=saved590.alignment; support=saved590.support;
 else
  align=qe_zlp_integer_align(E,profiles,profiles,[-100 100]); support=align.support; align=rmfield(align,'aligned');
 end
 assert(all(align.valid),'Invalid ZLP frame requires explicitly matched A0 subset');
 A1=nan(numel(support),nq); A0=nan(numel(E),nq); invalid=false(1,nq);
 for q=1:nq
  x=double(squeeze(raw(:,q,:))).'; invalid(q)=any(~isfinite(x)|x>=double(intmax('uint32')),'all');
  if invalid(q), continue; end
  A0(:,q)=sum(x,2); shifted=zeros(numel(support),1);
  for t=1:T, shifted=shifted+x(support+align.measured_offset_pixels(t),t); end
  A1(:,q)=shifted;
 end
 clear raw
 assert(isequal(find(invalid),rec.actual.invalid_native_q(:).'));
 assert(max(abs(A0(:,~invalid)-qe.intensity(:,~invalid)),[],'all')==0);
 qe1=qe; qe1.energy_meV=E(support); qe1.intensity=A1; bins=qe_centered_bins(qe1,targets,cfg.N);
 dest=fullfile(out,sid); mkdir(dest); mkdir(fullfile(dest,'spectra'));
 map=struct('session',sid,'E',E(support),'q',qe.q_Ainv,'A1',A1,'A0',A0(support,:), ...
  'source_channel',qe.source_channel,'alignment',align,'support',support,'raw_sha256',rawfile.sha256, ...
  'json_sha256',jsonfile.sha256,'invalid_q',invalid,'q_zero_native',rec.actual.q_zero_native_channel, ...
  'exposure_metadata',jsondecode(fileread(jsonfile.path)));
 save(fullfile(dest,'A1_full_q_map.mat'),'map','bins','-v7'); details={}; provenance={};
 for bi=1:numel(bins)
  b=bins(bi); if ~b.valid, continue; end
  u=struct('bin_size_requested',3,'source_q_count',3,'q_Ainv',b.q_Ainv,'q_left',b.q_left,'q_right',b.q_right,'source_channel',b.source_channel);
  for mi=1:2
   reused=false;
   if si==1
    for oi=1:numel(reuse)
     od=reuse{oi};
     if abs(od.unit.q_Ainv-b.q_Ainv)<1e-10&&strcmp(od.fits(1).peak_model,cfg.models{mi})
      fitmask=map.E>=cfg.window(1)&map.E<=cfg.window(2); opt=od.fits(1).effective_options;
      assert(max(abs(b.mean(fitmask)-od.fits(1).observed))<1e-8);
      assert(isequal(map.E(fitmask),od.fits(1).energy_meV)&&opt.n_starts==cfg.starts&&opt.extra_starts==cfg.extra_n2&&strcmp(opt.start_policy,'independent'));
      fits=od.fits; reused=true; break
     end
    end
   end
   if ~reused
    fits=qe_compare_component_models(map.E,b.mean,energy_window=cfg.window,peak_model=cfg.models{mi}, ...
     n_starts=cfg.starts,extra_starts=cfg.extra_n2,start_policy='independent',baseline_mode='power_law',seed=cfg.seed+round(mean(b.source_channel)));
   end
   key=sprintf('q%+.5f_%s',b.q_Ainv,cfg.models{mi}); d=struct('key',key,'unit',u,'fits',fits,'reused_v5',reused,'raw_sha256',rawfile.sha256,'config_hash',cfghash);
   save(fullfile(dest,'spectra',[key '.mat']),'d','-v7'); details{end+1}=d; %#ok<AGROW>
   provenance(end+1,:)={string(key),b.q_Ainv,string(cfg.models{mi}),reused}; %#ok<AGROW>
   c20_v4_plot_fits(d,fullfile(dest,'spectra',[key '.png'])); disp([sid ' ' key ' reused=' num2str(reused)]);
  end
 end
 save(fullfile(dest,'fit_details.mat'),'details','-v7'); c20_v4_export_fits(details,dest);
 writetable(cell2table(provenance,'VariableNames',{'key','q','model','reused_v5'}),fullfile(dest,'fit_origin.csv'));
 products{si}=struct('session',sid,'details',{details},'bins',bins,'map',map);
 assert(strcmpi(c20_v4_io('hash',rawfile.path),rawfile.sha256));
end
save(fullfile(out,'appendix','products.mat'),'products','cfg','-v7');
assert(isequal(before{1},c20_v4_io('inventory',p2))&&isequal(before{2},c20_v4_io('inventory',p5)));
disp(['PHYSICS_DATA_COMPLETE=' out]);
end
