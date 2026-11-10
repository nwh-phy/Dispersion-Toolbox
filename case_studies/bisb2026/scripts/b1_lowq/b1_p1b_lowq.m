function out = b1_p1b_lowq(run_root, session, qlist)
% P1b low-q fits (settings: b1_p1b_config) on the P1 channels: for each kernel
% version (KA, KB) a main fit with the 250 meV prefactor floor, fixed-parameter
% forward comparisons (floor 0 vs 250; refined vs default quadrature) and a
% short fit with floor 0. Noise from the common-beam reference (b1_p1b_diag).
% Channels listed in cfg.perturbation also get the holdout perturbation test.
% Output: <run_root>/p1b/<session>_<fingerprint>/q_<q>.mat and config.json.
if nargin<2, session='590_PL2_10w'; end
if nargin<3, qlist=[0.0005 0.001 0.002 0.005]; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_p1b_config();
cf=fullfile(run_root,'stage1',sprintf('channels_%s.mat',session)); S=load(cf,'C'); C=S.C; mf=matfile(cf);
Rf=load(fullfile(run_root,'p1b',['diag_' session],'common_ref.mat'),'ref'); ref=Rf.ref;
cfg.session=session; cfg.dq=C.dq; cfg.sigma_probe=C.sigma_probe;
files={fullfile(root,'src','fitting','qe_zlp_joint_fit.m'),fullfile(root,'src','fitting','qe_kinematic_prefactor.m'), ...
 fullfile(root,'src','fitting','peak_models.m'),[mfilename('fullpath') '.m'],which('b1_p1b_config'),which('b1_p1b_noise')};
cfg.code=struct('file',{},'sha256',{});
for k=1:numel(files), [~,nm,ex]=fileparts(files{k}); cfg.code(k)=struct('file',[nm ex],'sha256',b1_sha256({'file',files{k}})); end
cfg.common_ref_sha256=b1_sha256(typecast(ref(:),'uint8'));
fp=b1_sha256(cfg);
out=fullfile(run_root,'p1b',sprintf('%s_%s',session,fp(1:8))); if ~isfolder(out), mkdir(out); end
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root));
meta=struct('fingerprint',fp,'git_head',strtrim(head),'matlab',version,'channels_file',cf,'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')),'cfg',cfg);
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
base={'fit_window',cfg.fit_window,'signal_window',cfg.signal_window,'holdout_windows',cfg.holdout_windows, ...
 'peak_model',cfg.peak_model,'n_peaks',cfg.n_peaks,'peak_E0_bounds',cfg.peak_E0_bounds,'peak_width_max',cfg.peak_width_max, ...
 'aux_windows',cfg.aux_windows,'aux_width_max',cfg.aux_width_max,'aux_model',cfg.aux_model,'n_zlp',cfg.n_zlp, ...
 'asym_sigma',cfg.asym_sigma,'gain_weight',cfg.gain_weight,'core_weight',cfg.core_weight,'overshoot',cfg.overshoot, ...
 'overshoot_sigma',cfg.overshoot_sigma,'zlp_mode',cfg.zlp_mode,'kernel_iterations',cfg.kernel_iterations, ...
 'kernel_halfwidth',cfg.kernel_halfwidth,'prefit_aux',cfg.prefit_aux,'prefit_overshoot_local',cfg.prefit_overshoot_local, ...
 'prefactor_on_aux',cfg.prefactor_on_aux};
for qa=qlist
 for sgn=[-1 1]
  qv=sgn*qa; file=fullfile(out,sprintf('q_%+.4f.mat',qv));
  if isfile(file), L=load(file,'R'); assert(strcmp(L.R.fingerprint,fp),'stale %s',file); fprintf('skip %s\n',file); continue; end
  t0=tic;
  k=find(abs(C.q-qv)<1e-7); assert(numel(k)==1,'channel %g missing',qv);
  Fb=double(mf.frames(:,:,k)); E=C.E; Y=sum(Fb,2);
  sg_fd=b1_p1b_noise(Fb,E,ref,cfg.holdout_windows,'first_diff');
  sg_com=b1_p1b_noise(Fb,E,ref,cfg.holdout_windows,'common');
  sg_ch=b1_p1b_noise(Fb,E,ref,cfg.holdout_windows,'channel');
  w=E>=cfg.fit_window(1)&E<=cfg.fit_window(2);
  ql=cfg.quadrature_lo; qh=cfg.quadrature_hi;
  [Klo,kinfo]=qe_kinematic_prefactor(C.q(k),dq=C.dq,form=cfg.kinematic_form,h_perp=cfg.h_perp,sigma_probe=C.sigma_probe,beam_kV=cfg.beam_kV, ...
   n_sub=ql.n_sub,n_probe=ql.n_probe,probe_extent=ql.probe_extent,q_round=ql.q_round);
  Khi=qe_kinematic_prefactor(C.q(k),dq=C.dq,form=cfg.kinematic_form,h_perp=cfg.h_perp,sigma_probe=C.sigma_probe,beam_kV=cfg.beam_kV, ...
   n_sub=qh.n_sub,n_probe=qh.n_probe,probe_extent=qh.probe_extent,q_round=qh.q_round);
  L1=load(fullfile(run_root,cfg.p1_dir,sprintf('q_%+.4f.mat',qv)),'R'); R1=L1.R;
  assert(max(abs(R1.Y(:)-Y(w)))==0,'P1 data differ at q=%g',qv);
  S0=cfg.start_E0; W0=cfg.start_width;
  if cfg.add_p1_start, S0=[S0; R1.p1.parameters(:,1).']; W0=[W0; R1.p1.parameters(:,2).']; end
  seed=cfg.seed+round(1e5*abs(qv))+(qv>0);
  F=struct('kernel',{},'f250',{},'f0',{},'forward',{},'quadrature_ok',{},'f0_quadrature',{});
  for kv=1:numel(cfg.kernels)
   Kv=cfg.kernels(kv); kopt={'prefit_gain_min',Kv.prefit_gain_min,'prefit_loss_max',Kv.prefit_loss_max};
   f250=qe_zlp_joint_fit(E,Y,base{:},kopt{:},'prefactor',Klo,'prefactor_floor_meV',250,'noise_sigma',sg_fd, ...
    'start_E0',S0,'start_width',W0,'n_starts',size(S0,1)+cfg.n_random_starts,'seed',seed+10*kv);
   fw=forward_terms(f250,Klo,Khi);
   s=sg_fd(w); chk=(E(w)>=-180 & E(w)<=1800);
   qmax=max(abs(fw.d_quad0(chk))./s(chk)); quad_ok=qmax<cfg.quadrature_tol_sigma;
   K0=Klo; qname='lo'; if ~quad_ok, K0=Khi; qname='hi'; end
   f0=qe_zlp_joint_fit(E,Y,base{:},kopt{:},'prefactor',K0,'prefactor_floor_meV',0,'noise_sigma',sg_fd, ...
    'start_E0',[S0(1:size(cfg.start_E0,1),:); f250.parameters(:,1).'],'start_width',[W0(1:size(cfg.start_width,1),:); f250.parameters(:,2).'], ...
    'n_starts',size(cfg.start_E0,1)+1+cfg.n_random_starts_short,'seed',seed+10*kv+5,'warm_u',f250.u);
   fw.quad_max_sigma_floor0=qmax; fw.quad_max_sigma_floor250=max(abs(fw.d_quad250(chk))./s(chk));
   fw.floor_max_sigma=max(abs(fw.d_floor(chk))./s(chk));
   F(kv)=struct('kernel',Kv.name,'f250',pack(f250),'f0',pack(f0),'forward',fw,'quadrature_ok',quad_ok,'f0_quadrature',qname);
  end
  pert=struct([]);
  for pp=cfg.perturbation
   if abs(pp.q-qv)<1e-9, pert=perturbation_test(E,Fb,ref,cfg,base,Klo,S0,W0,seed,F(2).f250); end
  end
  R=struct('fingerprint',fp,'session',session,'q',qv,'members',C.q(k),'q_mean',kinfo.q_mean,'sigma_probe',C.sigma_probe,'dq',C.dq, ...
   'E',E(w),'Y',Y(w),'sigma',sg_fd(w),'sigma_common',sg_com(w),'sigma_channel',sg_ch(w),'fits',F,'perturbation',pert, ...
   'p1',struct('prediction',R1.p1.curves.prediction,'parameters',R1.p1.parameters,'zlp',R1.p1.curves.zlp,'aux_peaks',R1.p1.curves.aux_peaks,'peaks',R1.p1.curves.peaks), ...
   'baseline',struct('prediction',R1.baseline.curves.prediction,'parameters',R1.baseline.parameters,'zlp',R1.baseline.curves.zlp, ...
   'aux_peaks',R1.baseline.curves.aux_peaks,'peaks',R1.baseline.curves.peaks),'seconds',toc(t0));
  save(file,'R','-v7');
  for kv=1:numel(F)
   a=F(kv).f250; b=F(kv).f0;
   fprintf('%s q=%+.4f %s: f250 %s E0 %s W %s c %.1f | floor fwd max %.1f sigma, quad(0) %.3f sigma -> f0[%s] %s E0 %s W %s\n',session,qv,F(kv).kernel, ...
    a.status,mat2str(round(a.parameters(:,1)')),mat2str(round(a.parameters(:,2)')),a.zlp_parameters.center_meV, ...
    F(kv).forward.floor_max_sigma,F(kv).forward.quad_max_sigma_floor0,F(kv).f0_quadrature,b.status,mat2str(round(b.parameters(:,1)')),mat2str(round(b.parameters(:,2)')));
  end
  if ~isempty(pert), fprintf('  perturbation test: %s\n',pert.summary); end
  fprintf('  %.0f s\n',R.seconds);
 end
end
end

function fw=forward_terms(f,Klo,Khi)
% B1 contribution at the fitted parameters (amplitudes fixed) with different
% prefactors, rebuilt exactly as in the model: DL on the extended axis
% (detailed balance below zero), times the prefactor, convolved with the fit's
% kernel. d_floor = floor 0 - floor 250; d_quad0/250 = refined - default
% quadrature at floor 0 / 250. check = rebuilt(floor 250) vs the fit's B1 sum.
E=f.energy_meV; dE=median(diff(E)); lag=f.kernel.lag_meV; Kn=f.kernel.K/sum(f.kernel.K);
Ex=(E(1)-lag(end))+(0:(numel(E)+numel(lag)-2)).'*dE; Ea=abs(Ex);
m=peak_models('lorentz'); kT=0.08617333262*300; L=zeros(size(Ex));
for j=1:size(f.parameters,1)
 l=f.parameters(j,3)*m.model_fn(f.parameters(j,1),f.parameters(j,2),1,Ea); l(Ex<0)=l(Ex<0).*exp(-Ea(Ex<0)/kT); l(~isfinite(l))=0; L=L+l;
end
cv=@(K) conv(L.*K,Kn,'valid');
y250=cv(Klo(max(Ea,250))); y0=cv(Klo(Ea)); yh0=cv(Khi(Ea)); yh250=cv(Khi(max(Ea,250)));
fw=struct('E',E,'b1_floor250',y250,'b1_floor0',y0,'d_floor',y0-y250,'d_quad0',yh0-y0,'d_quad250',yh250-y250, ...
 'rebuild_check',max(abs(y250-sum(f.peaks,2)))/max(abs(sum(f.peaks,2))));
end

function p=perturbation_test(E,Fb,ref,cfg,base,Klo,S0,W0,seed,ref_fit)
% Change the frames inside the holdout band only; sigma outside, Stage A,
% kernel and all fitted parameters must not change (KB, floor 250).
h=E>=cfg.holdout_windows(1)&E<=cfg.holdout_windows(2); P=cfg.perturbation(1);
Fp=Fb; Fp(h,:)=P.scale*Fb(h,:)+P.offset_counts/size(Fb,2);
s0=b1_p1b_noise(Fb,E,ref,cfg.holdout_windows,'first_diff'); s1=b1_p1b_noise(Fp,E,ref,cfg.holdout_windows,'first_diff');
Kv=cfg.kernels(2); kopt={'prefit_gain_min',Kv.prefit_gain_min,'prefit_loss_max',Kv.prefit_loss_max};
fp=qe_zlp_joint_fit(E,sum(Fp,2),base{:},kopt{:},'prefactor',Klo,'prefactor_floor_meV',250,'noise_sigma',s1, ...
 'start_E0',S0,'start_width',W0,'n_starts',size(S0,1)+cfg.n_random_starts,'seed',seed+20);
same=@(a,b) isequal(a,b);
p=struct('sigma_outside_identical',same(s0(~h),s1(~h)),'parameters_identical',same(fp.parameters,ref_fit.parameters), ...
 'aux_identical',same(fp.aux_parameters,ref_fit.aux_parameters),'centre_identical',same(fp.zlp_parameters.center_meV,ref_fit.zlp_parameters.center_meV), ...
 'stageA_identical',same(fp.stageA.u,ref_fit.stageA.u),'kernel_identical',same(fp.kernel.K,ref_fit.kernel.K), ...
 'max_abs_dparam',max(abs(fp.parameters(:)-ref_fit.parameters(:))));
p.summary=sprintf('sigma %d, params %d, aux %d, centre %d, stageA %d, kernel %d (1 = identical); max|dparam| %.3g', ...
 p.sigma_outside_identical,p.parameters_identical,p.aux_identical,p.centre_identical,p.stageA_identical,p.kernel_identical,p.max_abs_dparam);
end

function p=pack(f)
p=struct('success',f.success,'status',f.numerical_status);
if ~f.success, p.stageA=f.stageA; return; end
p=struct('success',true,'status',f.numerical_status,'cost',f.cost,'u',f.u,'selected_start',f.selected_start,'start_costs',[f.candidates.cost], ...
 'parameters',f.parameters,'aux_parameters',f.aux_parameters,'peak_boundary',f.peak_boundary,'aux_boundary',f.aux_boundary, ...
 'zlp_parameters',f.zlp_parameters,'chi2_red_signal',f.chi2_red_signal,'chi2_red_gain',f.chi2_red_gain,'overshoot_active',f.overshoot_active, ...
 'holdout_mask',f.holdout_mask,'data_mask',f.data_mask,'overshoot_mask',f.overshoot_mask,'kernel',f.kernel,'stageA',f.stageA, ...
 'curves',struct('E',f.energy_meV,'observed',f.observed,'sigma',f.sigma,'prediction',f.prediction,'zlp',f.zlp, ...
  'peaks',f.peaks,'aux_peaks',f.aux_peaks,'background',f.background,'loss_function',f.loss_function,'prefactor',f.prefactor_grid));
end
