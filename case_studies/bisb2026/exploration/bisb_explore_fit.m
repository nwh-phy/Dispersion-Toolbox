function bisb_explore_fit(out_dir, stage, session_id)
% Existing joint fitter, explicit variants, resumable per-fit checkpoints.
% No historical derived MAT files, GUI changes, or implicit data paths.
arguments
 out_dir (1,:) char
 stage (1,:) char {mustBeMember(stage,{'check','fit','bootstrap'})} = 'fit'
 session_id (1,:) char = 'all'
end
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
run(fullfile(root,'startup.m'));
sessions={'590','n0','20w'};
if ~strcmp(session_id,'all'), sessions={session_id}; end
rows=struct([]);
for s=1:numel(sessions)
 sid=sessions{s}; folder=fullfile(out_dir,sid); D=load(fullfile(folder,'prepared.mat'));
 if strcmp(stage,'check')
  al=qe_zlp_integer_align(D.E_raw,D.reference,reshape(D.reference,size(D.reference,1),[],1),[-100 100]);
  assert(isequal(al.support,D.align_support),'explore:Support','Python/MATLAB supports differ');
  assert(isequal(al.valid,logical(D.align_valid)),'explore:Valid','Python/MATLAB frame masks differ');
  assert(isequal(al.measured_offset_pixels(al.valid),D.align_offsets(al.valid)),'explore:Shifts','Python/MATLAB shifts differ');
  fprintf('ALIGNMENT_PARITY %s PASS\n',sid);
  for q=[.0025 .005 .015]
   K=qe_kinematic_prefactor(q,dq=D.dq,h_perp=D.h,sigma_probe=D.sigma_probe,n_sub=53,n_probe=9);
   Kf=qe_kinematic_prefactor(q,dq=D.dq,h_perp=D.h,sigma_probe=D.sigma_probe,n_sub=105,n_probe=17);
   e=(300:25:1800)'; err=max(abs(K(e)./Kf(e)-1));
   fprintf('QUADRATURE %s q=%.4f max_relative=%.5g\n',sid,q,err);
   assert(err<.005,'explore:Quadrature','Acceptance quadrature differs by >=0.5 percent');
  end
  continue
 end
 dest=fullfile(folder,stage);if ~isfolder(dest),mkdir(dest);end
 signature=struct('prepared',file_sha(fullfile(folder,'prepared.mat')), ...
  'runner',file_sha([mfilename('fullpath') '.m']), ...
  'fitter',file_sha(which('qe_zlp_joint_fit')),'kinematic',file_sha(which('qe_kinematic_prefactor')));
 stamp=fullfile(dest,'signature.mat');
 if isfile(stamp)
  previous=load(stamp,'signature');
  assert(isequal(previous.signature,signature),'explore:StaleCache','Inputs or code changed; use a new output directory');
 else
  save(stamp,'signature');
 end
 if strcmp(stage,'bootstrap')
  do_bootstrap(D,dest,folder,sid);
  continue
 end
 variants={'nominal','one_peak','small_h','no_acceptance','no_kinematic', ...
           'poisson_weights','unaligned','symmetric','upper_2100','train_two','train_one'};
 for b=1:numel(D.q)
  base=[];
  for vi=1:numel(variants)
   v=variants{vi};
   if D.native(b) && vi>1,continue;end
   file=fullfile(dest,sprintf('bin%02d_%s.mat',b,v));
   if isfile(file)
    L=load(file,'result','row');f=L.result;row=L.row;
   else
    tic_fit=tic;
    [E,Y,sg,args,settings]=problem(D,b,v);
    if ~isempty(base) && base.success && settings.n_peaks==2 && ~strcmp(v,'symmetric')
     args=[args,{'warm_u',base.u}]; %#ok<AGROW>
    end
    f=qe_zlp_joint_fit(E,Y,args{:},noise_sigma=sg);
    row=compact(f,sid,b,D.q(b),D.width(b),D.native(b),v,toc(tic_fit));
    row.predictive_chi2=NaN;row.predictive_bias=NaN;
    if startsWith(v,'train_') && f.success
     scale=sum(D.beam(end-49:end))/sum(D.beam(1:end-50));
     observed=sum(double(D.frames(:,end-49:end,b)),2);
     sigma=D.test_sigma(:,b);[~,idx]=ismember(f.energy_meV,E);
     mask=f.energy_meV>=300 & f.energy_meV<=1800;
     r=(observed(idx)-scale*f.prediction)./sigma(idx);
     row.predictive_chi2=mean(r(mask).^2);row.predictive_bias=mean(r(mask));
    end
    result=f; result.effective_options.prefactor=[];
    save(file,'result','row','settings','-v7');
   end
   if strcmp(v,'nominal'),base=f;end
   if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
   fprintf('%s bin %02d q=%+.4f %s: %s E0=%s chi2=%.3g %.1fs\n', ...
    sid,b,D.q(b),v,row.status,mat2str([row.E01 row.E02],5),row.chi2,row.seconds);
  end
  writetable(struct2table(rows),fullfile(out_dir,['fits_' session_id '.csv']));
 end
end
end


function h=file_sha(path)
fid=fopen(path,'rb');assert(fid>=0,'explore:File','Cannot read %s',path);
cleanup=onCleanup(@() fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,8*1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end


function [E,Y,sg,args,settings]=problem(D,b,v)
E=D.E(:);F=double(D.frames(:,:,b));sg=D.sigma(:,b);
np=2;model='lorentz';h=D.h;probe=D.sigma_probe;kin=true;window=[300 1800];starts=4;
if strcmp(v,'nominal'),starts=8;end
if strcmp(v,'one_peak')||strcmp(v,'train_one'),np=1;starts=8;end
if strcmp(v,'small_h'),h=.0005;end
if strcmp(v,'no_acceptance'),h=0;probe=0;end
if strcmp(v,'no_kinematic'),kin=false;end
if strcmp(v,'symmetric'),model='lorentz_symmetric';end
if strcmp(v,'upper_2100'),window=[300 2100];end
if strcmp(v,'unaligned'),F=double(D.unaligned_frames(:,:,b));sg=D.unaligned_sigma(:,b);end
if startsWith(v,'train_'),F=F(:,1:end-50);sg=D.train_sigma(:,b);end
Y=sum(F,2);
if strcmp(v,'poisson_weights'),sg=sqrt(D.counts_per_electron*max(Y,1));end
args={'fit_window',[-300 window(2)],'signal_window',window,'peak_model',model,'n_peaks',np, ...
 'n_starts',starts,'seed',20261023+b,'aux_windows',[30 80;80 300],'prefactor_floor_meV',250, ...
 'prefactor_on_aux',false,'max_iterations',300};
if kin
 K=qe_kinematic_prefactor(D.members{b}(:),dq=D.dq,form='2d',beam_kV=D.beam_kV, ...
  h_perp=h,sigma_probe=probe,n_sub=53,n_probe=9);
 args=[args,{'prefactor',K}];
end
settings=struct('n_peaks',np,'peak_model',model,'h_Ainv',h,'probe_sigma_Ainv',probe, ...
 'kinematic',kin,'window',window,'prefactor_floor_meV',250,'native',D.native(b), ...
 'q_members',D.members{b},'n_frames',size(F,2),'noise',v);
end


function row=compact(f,sid,b,q,width,native,v,seconds)
row=struct('session',sid,'bin',b,'q_Ainv',q,'width_Ainv',width,'native',logical(native), ...
 'variant',v,'status',f.numerical_status,'success',f.success,'seconds',seconds, ...
 'E01',NaN,'E02',NaN,'Gamma1',NaN,'Gamma2',NaN,'apex1',NaN,'apex2',NaN, ...
 'amplitude1',NaN,'amplitude2',NaN,'chi2',NaN,'R2',NaN,'boundary1',true,'boundary2',true);
if ~f.success,return;end
row.chi2=f.chi2_red_signal;row.R2=f.r2_signal;
for p=1:size(f.parameters,1)
 e=f.parameters(p,1);g=f.parameters(p,2);c=2*e^2-g^2;
 ap=sqrt((c+sqrt(c^2+12*e^4))/6);
 if strcmp(f.effective_options.peak_model,'lorentz_symmetric'),ap=e;end
 row.(sprintf('E0%d',p))=e;row.(sprintf('Gamma%d',p))=g;row.(sprintf('apex%d',p))=ap;
 row.(sprintf('amplitude%d',p))=f.parameters(p,3);
 row.(sprintf('boundary%d',p))=any(f.peak_boundary(p,:))||f.zero_amplitude(p);
end
end


function do_bootstrap(D,dest,folder,sid)
% Shared draws across all chosen q spectra preserve same-frame covariance.
% IID and length-10 circular blocks are distinct uncertainty diagnostics.
select=find(~D.native & (abs(abs(D.q)-.005)<1e-9 | abs(D.q-.015)<1e-9));
nb=128;
rows=struct([]);
for block=[1 10]
 rs=RandStream('mt19937ar','Seed',20261023+block);
 T=size(D.frames,2);draws=zeros(T,nb);
 for r=1:nb
  starts=randi(rs,T,ceil(T/block),1)-1;
  ix=mod(starts+(0:block-1),T)+1;ix=reshape(ix.',[],1);draws(:,r)=ix(1:T);
 end
 for b=select
  L=load(fullfile(folder,'fit',sprintf('bin%02d_nominal.mat',b)),'result');base=L.result;
  if ~base.success,continue;end
  [E,~,sg,args,settings]=problem(D,b,'bootstrap');
  k=find(strcmp(args(1:2:end),'n_starts'));args{2*k}=1;
  args=[args,{'warm_u',base.u}]; %#ok<AGROW>
  F=double(D.frames(:,:,b));
  for r=1:nb
   file=fullfile(dest,sprintf('bin%02d_block%02d_rep%03d.mat',b,block,r));
   if isfile(file)
    L=load(file,'row');row=L.row;
   else
    tic_fit=tic;Y=sum(F(:,draws(:,r)),2);
    % Conditional noise gain, adjusted for each resampled mean count level.
    sg_boot=sg.*sqrt(max(Y,1)./max(sum(F,2),1));
    f=qe_zlp_joint_fit(E,Y,args{:},noise_sigma=sg_boot);
    row=compact(f,sid,b,D.q(b),D.width(b),false,'bootstrap',toc(tic_fit));
    row.block_length=block;row.replicate=r;
    draw=draws(:,r);
    save(file,'row','draw','settings','-v7');
   end
   if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
   if mod(r,16)==0,fprintf('BOOT %s bin%02d block%d %d/%d\n',sid,b,block,r,nb);end
  end
  writetable(struct2table(rows),fullfile(folder,'bootstrap.csv'));
 end
end
end
