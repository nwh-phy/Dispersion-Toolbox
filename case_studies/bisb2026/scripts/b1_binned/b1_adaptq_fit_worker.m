function b1_adaptq_fit_worker(run_root, session, sides, part, nparts, variant)
% B1 adaptive-binning task, stage 2: fits of every bin of b1_adaptq_config
% for one session. Bins are fitted outward in |q| per sign, each
% configuration also warm-started from the previous bin's optimum (as
% b1_fullq_fit_worker). One file per bin under <run_root>/stage2/<session>/;
% existing files are skipped only when they carry the same configuration
% fingerprint (cfg_hash); a stale or unfingerprinted file stops the run.
% part/nparts splits the bins of a side.
% variant '2d_h5': Bosman-style main configuration (2D kinematic factor,
% h = 0.0005; user decision 2026-10-08), written to <run_root>/stage2_2d_h5/
% and warm-started from the 2D h = 0.002 fit of the same bin in stage2.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_adaptq_config();
if nargin<6, variant=''; end
stage='stage2';
if strcmp(variant,'2d_h5')
 cfg.configs={'lorentz',2,'2d'; 'lorentz',1,'2d'}; cfg.h_main=0.0005; cfg.h_variant_qmax=-Inf; stage='stage2_2d_h5';
end
src=fullfile(run_root,'stage1',sprintf('channels_%s.mat',session));
S=load(src,'C'); C=S.C; mf=matfile(src);
out=fullfile(run_root,stage,session); if ~isfolder(out), mkdir(out); end
cfg.dq=C.dq; cfg_hash=b1_sha256(cfg);
fid=fopen(fullfile(out,sprintf('config_%s.json',cfg_hash(1:12))),'w'); fprintf(fid,'%s',jsonencode(cfg,'PrettyPrint',true)); fclose(fid);
bins=adaptq_bins(cfg,C.q);
if nargin<4, part=1; nparts=1; end
for side=sides(:).'
 bs=bins([bins.sign]==side); [~,o]=sort([bs.abs_q]); bs=bs(o);
 n=numel(bs); edges=round(linspace(0,n,nparts+1)); bs=bs(edges(part)+1:edges(part+1));
 warm=struct();
 for b=bs
  file=fullfile(out,sprintf('bin_%+.4f_N%d.mat',b.q_center,b.N));
  if isfile(file)
   L=load(file,'R');
   assert(isfield(L.R,'cfg_hash') && strcmp(L.R.cfg_hash,cfg_hash), 'b1_adaptq_fit_worker:Stale', ...
    'existing %s was made with another configuration; move it or use a new stage', file);
   warm=remember(warm,L.R.fits); fprintf('skip %s\n',file); continue
  end
  Fb=zeros(numel(C.E),C.T);
  for k=b.channels, Fb=Fb+double(mf.frames(:,:,k)); end
  E=C.E; Y=sum(Fb,2); sg=noise_sigma(Fb,E); mq=C.q(b.channels);
  R=struct('session',session,'q',b.q_center,'N',b.N,'members',mq,'sigma_probe',C.sigma_probe, ...
   'cfg',cfg,'cfg_hash',cfg_hash,'fits',struct([]),'boot',struct([]));
  if ~isempty(variant)
   ref=fullfile(run_root,'stage2',session,sprintf('bin_%+.4f_N%d.mat',b.q_center,b.N));
   if isfile(ref), L=load(ref,'R'); warm=seed(warm,L.R.fits,cfg.h_main); end
  end
  t_bin=tic; main_fit=[];
  K=configurations(cfg,abs(b.q_center));
  for c=1:numel(K)
   [Kfun,kinfo]=qe_kinematic_prefactor(mq,dq=C.dq,form=K(c).form,h_perp=K(c).h,sigma_probe=C.sigma_probe,beam_kV=cfg.beam_kV);
   wu=[]; if isfield(warm,K(c).key), wu=warm.(K(c).key); end
   t0=tic;
   f=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
    'peak_model',K(c).model,'n_peaks',K(c).n,'n_starts',cfg.n_starts_other+(cfg.n_starts-cfg.n_starts_other)*K(c).main, ...
    'seed',cfg.seed+round(1e4*abs(b.q_center))+c,'aux_windows',cfg.aux{2},'n_zlp',cfg.n_zlp,'prefactor',Kfun, ...
    'prefactor_floor_meV',cfg.prefactor_floor,'prefactor_on_aux',cfg.prefactor_on_aux,'noise_sigma',sg,'warm_u',wu);
   s=compact(f,K(c),kinfo,toc(t0));
   if f.success && ~K(c).h_variant, s.curves=curves(f); end
   if isempty(R.fits), R.fits=s; else, R.fits(end+1)=s; end
   if K(c).main, main_fit=f; Kmain=Kfun; end
  end
  % Frame bootstrap of the main configuration.
  T=size(Fb,2); rs=RandStream('mt19937ar','Seed',cfg.seed+round(1e5*b.q_center)+7);
  for k=1:cfg.n_boot*(~isempty(main_fit)&&main_fit.success)
   cnt=accumarray(randi(rs,T,T,1),1,[T 1]); Ys=Fb*cnt;
   fb=qe_zlp_joint_fit(E,Ys,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
    'peak_model',K(1).model,'n_peaks',K(1).n,'n_starts',cfg.boot_starts,'seed',cfg.seed+k, ...
    'aux_windows',cfg.aux{2},'n_zlp',cfg.n_zlp,'prefactor',Kmain,'prefactor_floor_meV',cfg.prefactor_floor, ...
    'prefactor_on_aux',cfg.prefactor_on_aux,'noise_sigma',sg,'warm_u',main_fit.u);
   sb=compact(fb,K(1),struct(),0); sb=rmfield(sb,{'u'});
   if isempty(R.boot), R.boot=sb; else, R.boot(end+1)=sb; end
  end
  R.seconds=toc(t_bin); save(file,'R','-v7'); warm=remember(warm,R.fits);
  m=R.fits(1); fprintf('%s q=%+.4f N=%d done in %.0f s: main %s E0 %s W %s chi2 sig %.2f gain %.2f\n', ...
   session,b.q_center,b.N,R.seconds,m.status,mat2str(round(m.parameters(:,1)')),mat2str(round(m.parameters(:,2)')),m.chi2_signal,m.chi2_gain);
 end
end
end

function bins=adaptq_bins(cfg,q)
% Non-overlapping bins of N adjacent native channels of this data set (its own
% dq), same layout on both sides; a layout row covers lo < |q| <= hi (lo = 0
% keeps the lowest native channel).
bins=struct('sign',{},'abs_q',{},'q_center',{},'N',{},'channels',{});
tol=1e-9;
for side=[-1 1]
 for r=1:size(cfg.layout,1)
  lo=cfg.layout(r,1); hi=cfg.layout(r,2); N=cfg.layout(r,3);
  ch=find(sign(q)==side & abs(q)>lo+tol & abs(q)<=hi+tol); [~,o]=sort(abs(q(ch))); ch=ch(o);
  for j=1:N:numel(ch)-N+1
   cj=ch(j:j+N-1); members=q(cj);
   bins(end+1)=struct('sign',side,'abs_q',mean(abs(members)),'q_center',mean(members),'N',N,'channels',cj(:).'); %#ok<AGROW>
  end
 end
end
end

function sg=noise_sigma(Fb,E)
% As b1_fullq_build: per-frame variance/mean after removing the common-mode
% intensity, moving median along E, sigma^2 = g(E) * counts.
s=sum(Fb,1); m=mean(Fb,2); Fc=Fb.*(mean(s)./s);
ratio=var(Fc,0,2)./max(m,eps); g=movmedian(ratio,21,1,'omitnan'); g(~isfinite(g)|g<=0)=1;
sg=sqrt(g.*max(sum(Fb,2),1));
end

function K=configurations(cfg,aq)
K=struct('key',{},'model',{},'n',{},'form',{},'h',{},'main',{},'h_variant',{});
for i=1:size(cfg.configs,1)
 K(end+1)=entry(cfg.configs{i,1},cfg.configs{i,2},cfg.configs{i,3},cfg.h_main,false); %#ok<AGROW>
end
K(1).main=true;
if aq<=cfg.h_variant_qmax+1e-9
 for h=cfg.h_variants, K(end+1)=entry(cfg.configs{1,1},cfg.configs{1,2},cfg.configs{1,3},h,true); end %#ok<AGROW>
end
end

function e=entry(model,n,form,h,hv)
key=matlab.lang.makeValidName(sprintf('%s_n%d_%s_h%g',model,n,form,1e4*h));
e=struct('key',key,'model',model,'n',n,'form',form,'h',h,'main',false,'h_variant',hv);
end

function s=compact(f,c,kinfo,secs)
s=struct('key',c.key,'model',c.model,'n',c.n,'form',c.form,'h',c.h,'main',c.main,'h_variant',c.h_variant, ...
 'success',f.success,'status',f.numerical_status,'cost',NaN,'chi2_signal',NaN,'chi2_gain',NaN, ...
 'parameters',nan(c.n,3),'aux_parameters',[],'peak_boundary',false(c.n,2), ...
 'loss_area',nan(1,c.n),'background_fraction',NaN,'u',[],'q_mean',NaN,'seconds',secs,'curves',struct());
if isfield(kinfo,'q_mean'), s.q_mean=kinfo.q_mean; end
if ~f.success, return; end
s.cost=f.cost; s.chi2_signal=f.chi2_red_signal; s.chi2_gain=f.chi2_red_gain;
s.parameters=f.parameters; s.aux_parameters=f.aux_parameters; s.peak_boundary=f.peak_boundary;
Ep=f.energy_meV>0; s.loss_area=trapz(f.energy_meV(Ep),f.loss_function(Ep,:));
s.background_fraction=f.background_fraction_signal; s.u=f.u;
end

function c=curves(f)
c=struct('E',f.energy_meV,'observed',f.observed,'sigma',f.sigma,'prediction',f.prediction,'zlp',f.zlp, ...
 'peaks',f.peaks,'aux_peaks',f.aux_peaks,'background',f.background,'loss_function',f.loss_function,'prefactor',f.prefactor_grid);
end

function warm=remember(warm,fits)
for k=1:numel(fits)
 if fits(k).success && ~isempty(fits(k).u), warm.(fits(k).key)=fits(k).u; end
end
end

function warm=seed(warm,fits,h)
% Warm start for the h variant from the 2D fits of the same bin (other h).
for k=1:numel(fits)
 f=fits(k);
 if f.success && ~isempty(f.u) && strcmp(f.form,'2d') && ~f.h_variant
  warm.(matlab.lang.makeValidName(sprintf('%s_n%d_%s_h%g',f.model,f.n,f.form,1e4*h)))=f.u;
 end
end
end
