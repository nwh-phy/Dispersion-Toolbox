function bisb_explore_synthetic(folder)
% Independent Python-generated single/two-peak response, high/low counts.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
run(fullfile(root,'startup.m'));D=load(fullfile(folder,'synthetic.mat'));
dest=fullfile(folder,'fits');if ~isfolder(dest),mkdir(dest);end
rows=struct([]);
for j=1:size(D.Y,2)
 K=qe_kinematic_prefactor(D.q(j),dq=.0005,h_perp=.0066,n_sub=105);
 for np=1:2
  file=fullfile(dest,sprintf('case%02d_n%d.mat',j,np));
  if isfile(file)
   L=load(file,'row');row=L.row;
  else
   tic_fit=tic;
   f=qe_zlp_joint_fit(D.E,D.Y(:,j),fit_window=[-300 1800],signal_window=[300 1800], ...
    n_peaks=np,peak_model='lorentz',prefactor=K,prefactor_floor_meV=250, ...
    noise_sigma=sqrt(.82*max(D.mu(:,j),1)),n_starts=6,seed=20261023+j,max_iterations=300);
   row=struct('case_id',j,'n_true',D.n_true(j),'q_Ainv',D.q(j),'height',D.height(j), ...
    'replicate',D.replicate(j),'n_fit',np,'success',f.success,'status',f.numerical_status, ...
    'E01',NaN,'E02',NaN,'apex1',NaN,'apex2',NaN,'chi2',NaN,'boundary',true,'seconds',toc(tic_fit));
   if f.success
    row.chi2=f.chi2_red_signal;row.boundary=any(f.peak_boundary(:))||any(f.zero_amplitude);
    for k=1:np
     e=f.parameters(k,1);g=f.parameters(k,2);c=2*e^2-g^2;
     row.(sprintf('E0%d',k))=e;row.(sprintf('apex%d',k))=sqrt((c+sqrt(c^2+12*e^4))/6);
    end
   end
   f.effective_options.prefactor=[];save(file,'f','row','-v7');
  end
  if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
  fprintf('SYNTH case%d truth%d height%d fit%d %s E0=[%.1f %.1f] chi2=%.4g\n', ...
   j,D.n_true(j),D.height(j),np,row.status,row.E01,row.E02,row.chi2);
 end
 writetable(struct2table(rows),fullfile(folder,'synthetic_fits.csv'));
end
end
