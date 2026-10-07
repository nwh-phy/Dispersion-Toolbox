function b1_fullq_fit_worker(run_root, session, bin_idx)
% B1 full-q task, stages 2/4: every fit configuration plus the frame
% bootstrap of the main configuration for the given bins of one session.
% One file per bin under <run_root>/stage2/<session>/ (existing files are
% skipped, so a run can resume). Bins are processed outward in |q| per sign,
% and each configuration also starts from the previous bin's optimum.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_fullq_config();
src=fullfile(run_root,'stage1',sprintf('frames_%s.mat',session));
S=load(src,'D'); D=S.D; mf=matfile(src);
out=fullfile(run_root,'stage2',session); if ~isfolder(out), mkdir(out); end
q=D.q_center(bin_idx); [~,o]=sortrows([sign(q(:)) abs(q(:))]); bin_idx=bin_idx(o);
C=configurations(cfg); warm=struct();
for b=bin_idx(:).'
 qb=D.q_center(b); file=fullfile(out,sprintf('bin_%+.4f.mat',qb)); side=sprintf('s%d',sign(qb)>0);
 if isfile(file)
  L=load(file,'R'); warm=remember(warm,side,L.R.fits); fprintf('skip %s\n',file); continue
 end
 E=D.E; Y=double(D.Ysum(:,b)); sg=double(D.sigma(:,b)); mq=D.bins(b).q_members;
 R=struct('session',session,'q',qb,'members',mq,'source_channel',D.bins(b).source_channel, ...
  'sigma_probe',D.sigma_probe,'cfg',cfg,'fits',struct([]),'boot',struct([]));
 t_bin=tic;
 for c=1:numel(C)
  if C(c).h_variant && abs(qb)>cfg.h_variant_qmax, continue; end
  [Kfun,kinfo]=qe_kinematic_prefactor(mq,form=C(c).form,h_perp=C(c).h,sigma_probe=D.sigma_probe,beam_kV=cfg.beam_kV);
  wu=[]; if isfield(warm,side) && isfield(warm.(side),C(c).key), wu=warm.(side).(C(c).key); end
  t0=tic;
  f=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
   'peak_model',C(c).model,'n_peaks',C(c).n,'n_starts',cfg.n_starts_other+(cfg.n_starts-cfg.n_starts_other)*C(c).main,'seed',cfg.seed+b, ...
   'aux_windows',C(c).aux,'n_zlp',cfg.n_zlp,'prefactor',Kfun,'prefactor_floor_meV',cfg.prefactor_floor,'prefactor_on_aux',cfg.prefactor_on_aux,'noise_sigma',sg,'warm_u',wu);
  s=compact(f,C(c),kinfo,toc(t0));
  if C(c).main || (C(c).n==1 && C(c).primary_aux && strcmp(C(c).form,'2d') && ~C(c).h_variant) || ...
    (strcmp(C(c).model,'lorentz_symmetric') && C(c).n==2 && C(c).primary_aux && strcmp(C(c).form,'2d') && ~C(c).h_variant)
   s.curves=curves(f);
  end
  if isempty(R.fits), R.fits=s; else, R.fits(end+1)=s; end
  if C(c).main, main_fit=f; Kmain=Kfun; end
 end
 % Frame bootstrap of the main configuration (warm start from its optimum).
 Fb=double(mf.frames(:,:,b)); T=size(Fb,2); rs=RandStream('mt19937ar','Seed',cfg.seed+1000*b);
 nboot=cfg.n_boot*main_fit.success;
 for k=1:nboot
  cnt=accumarray(randi(rs,T,T,1),1,[T 1]); Ys=Fb*cnt;
  fb=qe_zlp_joint_fit(E,Ys,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
   'peak_model',C(1).model,'n_peaks',C(1).n,'n_starts',cfg.boot_starts,'seed',cfg.seed+b+k, ...
   'aux_windows',C(1).aux,'n_zlp',cfg.n_zlp,'prefactor',Kmain,'prefactor_floor_meV',cfg.prefactor_floor,'prefactor_on_aux',cfg.prefactor_on_aux,'noise_sigma',sg,'warm_u',main_fit.u);
  sb=compact(fb,C(1),struct(),0); sb=rmfield(sb,{'u'});
  if isempty(R.boot), R.boot=sb; else, R.boot(end+1)=sb; end
 end
 R.seconds=toc(t_bin);
 save(file,'R','-v7');
 warm=remember(warm,side,R.fits);
 m=R.fits(1); fprintf('%s q=%+.4f done in %.0f s: main %s E0 %.0f/%.0f W %.0f/%.0f chi2 sig %.2f gain %.2f\n', ...
  session,qb,R.seconds,m.status,m.parameters(1,1),m.parameters(2,1),m.parameters(1,2),m.parameters(2,2),m.chi2_signal,m.chi2_gain);
end
end

function C=configurations(cfg)
% First entry is the main configuration: DL, n=2, main aux, 2D, h_main.
C=struct('key',{},'model',{},'n',{},'aux_name',{},'aux',{},'form',{},'h',{},'main',{},'h_variant',{},'primary_aux',{});
for mi=1:numel(cfg.models)
 for n=[2 1]
  for a=1:size(cfg.aux,1)
   for k=1:numel(cfg.kinematic)
    C(end+1)=entry(cfg.models{mi},n,cfg.aux(a,:),cfg.kinematic{k},cfg.h_main,false,a==1); %#ok<AGROW>
   end
  end
 end
end
C(1).main=true;
for h=cfg.h_variants
 C(end+1)=entry(cfg.models{1},2,cfg.aux(1,:),cfg.kinematic{1},h,true,true); %#ok<AGROW>
end
end

function e=entry(model,n,aux,form,h,hv,primary)
key=matlab.lang.makeValidName(sprintf('%s_n%d_%s_%s_h%g',model,n,aux{1},form,1e4*h));
e=struct('key',key,'model',model,'n',n,'aux_name',aux{1},'aux',aux{2},'form',form,'h',h, ...
 'main',false,'h_variant',hv,'primary_aux',primary);
end

function s=compact(f,c,kinfo,secs)
s=struct('key',c.key,'model',c.model,'n',c.n,'aux',c.aux_name,'form',c.form,'h',c.h,'main',c.main, ...
 'success',f.success,'status',f.numerical_status,'cost',NaN,'chi2_signal',NaN,'chi2_gain',NaN,'chi2_low',NaN, ...
 'parameters',nan(c.n,3),'aux_parameters',nan(size(c.aux,1),3),'peak_boundary',false(c.n,2), ...
 'loss_area',nan(1,c.n),'loss_fsum',nan(1,c.n),'background_fraction',NaN,'zlp',struct(),'u',[], ...
 'q_mean',NaN,'seconds',secs,'curves',struct());
if isfield(kinfo,'q_mean'), s.q_mean=kinfo.q_mean; end
if ~f.success, return; end
s.cost=f.cost; s.chi2_signal=f.chi2_red_signal; s.chi2_gain=f.chi2_red_gain;
x=f.energy_meV-f.zlp_parameters.center_meV; lo=x>=30&x<=300&isfinite(f.residual);
s.chi2_low=mean((f.residual(lo)./f.sigma(lo)).^2);
s.parameters=f.parameters; s.aux_parameters=f.aux_parameters; s.peak_boundary=f.peak_boundary;
Ep=f.energy_meV>0; Ev=f.energy_meV(Ep);
s.loss_area=trapz(Ev,f.loss_function(Ep,:)); s.loss_fsum=trapz(Ev,Ev.*f.loss_function(Ep,:));
s.background_fraction=f.background_fraction_signal; s.zlp=f.zlp_parameters; s.u=f.u;
end

function c=curves(f)
c=struct('E',f.energy_meV,'observed',f.observed,'sigma',f.sigma,'prediction',f.prediction,'zlp',f.zlp, ...
 'peaks',f.peaks,'aux_peaks',f.aux_peaks,'loss_function',f.loss_function,'prefactor',f.prefactor_grid);
end

function warm=remember(warm,side,fits)
for k=1:numel(fits)
 if fits(k).success && ~isempty(fits(k).u), warm.(side).(fits(k).key)=fits(k).u; end
end
end
