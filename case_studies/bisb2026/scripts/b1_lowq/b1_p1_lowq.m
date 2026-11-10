function out = b1_p1_lowq(run_root, session, qlist)
% P1 low-q fits (settings: b1_p1_config) on single native channels, both
% signs, plus the old-style Stage A (ZLP only, gain side + core) for the
% prefit comparison and a frame bootstrap. The baseline fit of the same
% channel (cfg.baseline_stage) is loaded for the same-data comparison.
% Output: <run_root>/p1/<session>_<fingerprint>/q_<q>.mat and config.json;
% the fingerprint covers the configuration and the code of every file used.
if nargin<2, session='590_PL2_10w'; end
if nargin<3, qlist=[0.0005 0.001 0.002 0.005]; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_p1_config();
cf=fullfile(run_root,'stage1',sprintf('channels_%s.mat',session)); S=load(cf,'C'); C=S.C; mf=matfile(cf);
cfg.session=session; cfg.dq=C.dq; cfg.sigma_probe=C.sigma_probe;
files={fullfile(root,'src','fitting','qe_zlp_joint_fit.m'),fullfile(root,'src','fitting','qe_kinematic_prefactor.m'), ...
 fullfile(root,'src','fitting','peak_models.m'),[mfilename('fullpath') '.m'],which('b1_p1_config')};
cfg.code=struct('file',{},'sha256',{});
for k=1:numel(files), [~,nm,ex]=fileparts(files{k}); cfg.code(k)=struct('file',[nm ex],'sha256',b1_sha256({'file',files{k}})); end
fp=b1_sha256(cfg);
out=fullfile(run_root,'p1',sprintf('%s_%s',session,fp(1:8))); if ~isfolder(out), mkdir(out); end
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root));
meta=struct('fingerprint',fp,'git_head',strtrim(head),'matlab',version,'channels_file',cf,'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')),'cfg',cfg);
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
opts={'fit_window',cfg.fit_window,'signal_window',cfg.signal_window,'holdout_windows',cfg.holdout_windows, ...
 'peak_model',cfg.peak_model,'n_peaks',cfg.n_peaks,'peak_E0_bounds',cfg.peak_E0_bounds,'peak_width_max',cfg.peak_width_max, ...
 'aux_windows',cfg.aux_windows,'aux_width_max',cfg.aux_width_max,'aux_model',cfg.aux_model,'n_zlp',cfg.n_zlp, ...
 'asym_sigma',cfg.asym_sigma,'gain_weight',cfg.gain_weight,'core_weight',cfg.core_weight,'overshoot',cfg.overshoot, ...
 'overshoot_sigma',cfg.overshoot_sigma,'prefit_aux',cfg.prefit_aux,'prefit_loss_max',cfg.prefit_loss_max, ...
 'kernel_mode',cfg.kernel_mode,'kernel_halfwidth',cfg.kernel_halfwidth,'prefactor_floor_meV',cfg.prefactor_floor, ...
 'prefactor_on_aux',cfg.prefactor_on_aux};
for qa=qlist
 for sgn=[-1 1]
  qv=sgn*qa; file=fullfile(out,sprintf('q_%+.4f.mat',qv));
  if isfile(file), L=load(file,'R'); assert(strcmp(L.R.fingerprint,fp),'stale %s',file); fprintf('skip %s\n',file); continue; end
  t0=tic;
  k=find(abs(C.q-qv)<1e-7); assert(numel(k)==1,'channel %g missing',qv);
  Fb=double(mf.frames(:,:,k)); E=C.E; Y=sum(Fb,2); sg=noise_sigma(Fb);
  [Kfun,kinfo]=qe_kinematic_prefactor(C.q(k),dq=C.dq,form=cfg.kinematic_form,h_perp=cfg.h_perp,sigma_probe=C.sigma_probe,beam_kV=cfg.beam_kV);
  % baseline fit of the same channel; same data and noise model
  Lb=load(fullfile(run_root,cfg.baseline_stage,session,sprintf('bin_%+.4f_N1.mat',qv)),'R'); mb=Lb.R.fits(1);
  w=E>=cfg.fit_window(1)&E<=cfg.fit_window(2);
  assert(isequal(mb.curves.E(:),E(w)) && max(abs(mb.curves.observed(:)-Y(w)))==0, 'baseline data differ at q=%g',qv);
  dsig=max(abs(mb.curves.sigma(:)-sg(w))./sg(w));
  S0=cfg.start_E0; W0=cfg.start_width;
  if cfg.add_baseline_start, S0=[S0; mb.parameters(:,1).']; W0=[W0; mb.parameters(:,2).']; end
  seed=cfg.seed+round(1e5*abs(qv))+(qv>0);
  f=qe_zlp_joint_fit(E,Y,opts{:},'prefactor',Kfun,'noise_sigma',sg,'start_E0',S0,'start_width',W0, ...
   'n_starts',size(S0,1)+cfg.n_random_starts,'seed',seed);
  % old-style Stage A: ZLP only, gain side + core (original defaults, no holdout)
  fA0=qe_zlp_joint_fit(E,Y,'fit_window',cfg.fit_window,'signal_window',[300 1800],'n_peaks',0,'n_zlp',cfg.n_zlp, ...
   'noise_sigma',sg,'kernel_mode','fitted_zlp');
  % frame bootstrap, warm start from the optimum plus a few random starts
  T=size(Fb,2); rs=RandStream('mt19937ar','Seed',seed+7); B=struct('parameters',{},'aux_parameters',{},'status',{},'cost',{});
  for b=1:cfg.n_boot*f.success
   cnt=accumarray(randi(rs,T,T,1),1,[T 1]); Ys=Fb*cnt;
   fb=qe_zlp_joint_fit(E,Ys,opts{:},'prefactor',Kfun,'noise_sigma',sg,'start_E0',f.parameters(:,1).','start_width',f.parameters(:,2).', ...
    'n_starts',1+cfg.boot_random_starts,'seed',seed+100+b,'warm_u',f.u);
   if fb.success, B(end+1)=struct('parameters',fb.parameters,'aux_parameters',fb.aux_parameters,'status',fb.numerical_status,'cost',fb.cost); end %#ok<AGROW>
  end
  R=struct('fingerprint',fp,'session',session,'q',qv,'members',C.q(k),'q_mean',kinfo.q_mean,'sigma_probe',C.sigma_probe,'dq',C.dq, ...
   'E',E(w),'Y',Y(w),'sigma',sg(w),'baseline_sigma_maxreldiff',dsig,'p1',pack(f),'stageA_old',fA0.stageA, ...
   'baseline',struct('stage',cfg.baseline_stage,'parameters',mb.parameters,'aux_parameters',mb.aux_parameters,'status',mb.status, ...
   'cost',mb.cost,'chi2_signal',mb.chi2_signal,'curves',mb.curves),'boot',B,'seconds',toc(t0));
  save(file,'R','-v7');
  fprintf('%s q=%+.4f: P1 %s E0 %s W %s aux %s | baseline E0 %s | stageA chi2(gain) old %.2f new %.2f | %d boot | %.0f s\n', ...
   session,qv,f.numerical_status,mat2str(round(f.parameters(:,1)')),mat2str(round(f.parameters(:,2)')),mat2str(round(f.aux_parameters(:,1:2))), ...
   mat2str(round(mb.parameters(:,1)')),fA0.stageA.chi2_red_data,f.stageA.chi2_red_data,numel(B),R.seconds);
 end
end
end

function p=pack(f)
p=struct('success',f.success,'status',f.numerical_status);
if ~f.success, return; end
c=[f.candidates.cost];
p=struct('success',true,'status',f.numerical_status,'cost',f.cost,'u',f.u,'selected_start',f.selected_start,'start_costs',c, ...
 'parameters',f.parameters,'aux_parameters',f.aux_parameters,'peak_boundary',f.peak_boundary,'aux_boundary',f.aux_boundary, ...
 'zlp_parameters',f.zlp_parameters,'chi2_red_signal',f.chi2_red_signal,'chi2_red_gain',f.chi2_red_gain, ...
 'overshoot_active',f.overshoot_active,'holdout_mask',f.holdout_mask,'data_mask',f.data_mask,'overshoot_mask',f.overshoot_mask, ...
 'kernel',f.kernel,'stageA',f.stageA, ...
 'curves',struct('E',f.energy_meV,'observed',f.observed,'sigma',f.sigma,'prediction',f.prediction,'zlp',f.zlp, ...
  'peaks',f.peaks,'aux_peaks',f.aux_peaks,'background',f.background,'loss_function',f.loss_function,'prefactor',f.prefactor_grid));
end

function sg=noise_sigma(Fb)
% Same noise model as b1_adaptq_fit_worker / b1_fullq_build.
s=sum(Fb,1); m=mean(Fb,2); Fc=Fb.*(mean(s)./s);
ratio=var(Fc,0,2)./max(m,eps); g=movmedian(ratio,21,1,'omitnan'); g(~isfinite(g)|g<=0)=1;
sg=sqrt(g.*max(sum(Fb,2),1));
end
