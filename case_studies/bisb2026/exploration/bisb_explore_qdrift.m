function bisb_explore_qdrift(out_dir)
% Controlled late-stage diagnostic: only q bin transport and its noise change.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
run(fullfile(root,'startup.m'));sessions={'590','n0','20w'};rows=struct([]);
for s=1:numel(sessions)
 sid=sessions{s};folder=fullfile(out_dir,sid);D=load(fullfile(folder,'prepared.mat'));
 Q=load(fullfile(folder,'qdrift.mat'));dest=fullfile(folder,'qdrift_fits');if ~isfolder(dest),mkdir(dest);end
 for j=1:numel(Q.q)
  b=Q.base_bins(j);
  for variant={'nominal','train_two'}
   v=variant{1};L=load(fullfile(folder,'fit',sprintf('bin%02d_%s.mat',b,v)),'result');base=L.result;
   opts=base.effective_options;opts.warm_u=base.u;opts.n_starts=4;
   opts.prefactor=qe_kinematic_prefactor(D.members{b}(:),dq=D.dq,h_perp=D.h,sigma_probe=D.sigma_probe,n_sub=53,n_probe=9);
   F=Q.frames(:,:,j);opts.noise_sigma=Q.sigma(:,j);
   if strcmp(v,'train_two'),F=F(:,1:end-50);opts.noise_sigma=Q.train_sigma(:,j);end
   args=namedargs2cell(opts);f=qe_zlp_joint_fit(Q.E,sum(F,2),args{:});
   row=struct('session',sid,'q_Ainv',Q.q(j),'variant',v,'success',f.success,'status',f.numerical_status, ...
    'E01',NaN,'E02',NaN,'apex1',NaN,'apex2',NaN,'chi2',NaN,'predictive_chi2',NaN);
   if f.success
    row.chi2=f.chi2_red_signal;
    for p=1:2
     e=f.parameters(p,1);g=f.parameters(p,2);c=2*e^2-g^2;
     row.(sprintf('E0%d',p))=e;row.(sprintf('apex%d',p))=sqrt((c+sqrt(c^2+12*e^4))/6);
    end
    if strcmp(v,'train_two')
     scale=sum(Q.beam(end-49:end))/sum(Q.beam(1:end-50));
     y=sum(Q.frames(:,end-49:end,j),2);sg=Q.test_sigma(:,j);[~,idx]=ismember(f.energy_meV,Q.E);
     mask=f.energy_meV>=300 & f.energy_meV<=1800;r=(y(idx)-scale*f.prediction)./sg(idx);
     row.predictive_chi2=mean(r(mask).^2);
    end
   end
   f.effective_options.prefactor=[];save(fullfile(dest,sprintf('bin%02d_%s.mat',b,v)),'f','row','-v7');
   if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
   fprintf('QDRIFT %s q=%+.4f %s %s E0=%.2f %.2f predictive=%.3f\n',sid,Q.q(j),v,row.status,row.E01,row.E02,row.predictive_chi2);
  end
 end
 writetable(struct2table(rows),fullfile(out_dir,'summary','qdrift_fits.csv'));
end
end
