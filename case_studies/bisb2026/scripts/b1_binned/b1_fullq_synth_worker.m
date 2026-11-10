function b1_fullq_synth_worker(run_root, session, q_list)
% B1 full-q task, stage 3: two-branch synthetic spectra. The main-configuration
% fit of each selected bin is the truth; noise has the measured level and the
% measured lag-1 correlation along E (AR(1)); each replica is refit with the
% main configuration (warm start at the truth plus random starts).
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_fullq_config();
S=load(fullfile(run_root,'stage1',sprintf('frames_%s.mat',session)),'D'); D=S.D;
rho=median(D.noise_rho_loss(1,:),'omitnan'); if ~isfinite(rho), rho=0; end
out=fullfile(run_root,'stage3',session); if ~isfolder(out), mkdir(out); end
for qb=q_list(:).'
 file=fullfile(out,sprintf('synth_%+.4f.mat',qb));
 if isfile(file), fprintf('skip %s\n',file); continue; end
 L=load(fullfile(run_root,'stage2',session,sprintf('bin_%+.4f.mat',qb)),'R'); R=L.R;
 m=R.fits([R.fits.main]); c=m.curves; E=c.E; truth=c.prediction; sg=c.sigma;
 Kfun=qe_kinematic_prefactor(R.members,form=m.form,h_perp=m.h,sigma_probe=R.sigma_probe,beam_kV=cfg.beam_kV);
 rs=RandStream('mt19937ar','Seed',cfg.seed+round(1e5*abs(qb))+(qb>0));
 Ssyn=struct('k',{},'parameters',{},'aux_parameters',{},'status',{},'chi2_signal',{},'loss_area',{});
 for k=1:cfg.n_synth
  z=randn(rs,numel(E),1); e=z;
  for i=2:numel(e), e(i)=rho*e(i-1)+sqrt(1-rho^2)*z(i); end
  Ys=truth+sg.*e;
  f=qe_zlp_joint_fit(E,Ys,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
   'peak_model',m.model,'n_peaks',m.n,'n_starts',6,'seed',cfg.seed+k, ...
   'aux_windows',cfg.aux{1,2},'n_zlp',cfg.n_zlp,'prefactor',Kfun,'prefactor_floor_meV',cfg.prefactor_floor,'prefactor_on_aux',cfg.prefactor_on_aux,'noise_sigma',sg,'warm_u',m.u);
  p=nan(m.n,3); a=nan(size(cfg.aux{1,2},1),3); la=nan(1,m.n); cs=NaN;
  if f.success
   p=f.parameters; a=f.aux_parameters; cs=f.chi2_red_signal;
   Ep=f.energy_meV>0; la=trapz(f.energy_meV(Ep),f.loss_function(Ep,:));
  end
  Ssyn(end+1)=struct('k',k,'parameters',p,'aux_parameters',a,'status',f.numerical_status, ...
   'chi2_signal',cs,'loss_area',la); %#ok<AGROW>
 end
 Z=struct('session',session,'q',qb,'truth',m.parameters,'truth_loss_area',m.loss_area,'rho',rho,'replicas',Ssyn);
 save(file,'Z','-v7');
 P=cat(3,Ssyn.parameters); dp=P(:,1:2,:)-m.parameters(:,1:2);
 fprintf('%s q=%+.4f synth: bias E0 %s W %s | sd E0 %s W %s meV\n',session,qb, ...
  mat2str(round(mean(dp(:,1,:),3,'omitnan')')),mat2str(round(mean(dp(:,2,:),3,'omitnan')')), ...
  mat2str(round(std(P(:,1,:),0,3,'omitnan')')),mat2str(round(std(P(:,2,:),0,3,'omitnan')')));
end
end
