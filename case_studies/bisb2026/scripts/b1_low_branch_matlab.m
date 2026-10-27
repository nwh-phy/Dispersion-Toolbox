function b1_low_branch_matlab(stage,cfg,out_root,session_id,part,nparts)
% Raw frames -> native/high-q 3/5 bins -> existing MATLAB joint fitter.
% All computation, figures and verification are native MATLAB. Stages are
% resumable. No derived arrays or parameters from the Python exploration.
arguments
 stage (1,:) char {mustBeMember(stage,{'prepare','pilot','fit','bootstrap','report','verify'})}
 cfg struct
 out_root (1,:) char
 session_id (1,:) char = 'all'
 part (1,1) double = 1
 nparts (1,1) double = 1
end
root=bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
if ~isfolder(out_root),mkdir(out_root);end
sessions=cfg.session;if ~strcmp(session_id,'all'),sessions=sessions(strcmp({sessions.id},session_id));end
assert(~isempty(sessions),'b1_low_branch_matlab:Session','Unknown session');
if strcmp(stage,'report'),report(out_root,cfg);return;end
for s=sessions
 folder=fullfile(out_root,s.id);if ~isfolder(folder),mkdir(folder);end
 file=fullfile(folder,'prepared.mat');
 if strcmp(stage,'prepare'),prepare(s,cfg,file);continue;end
 L=load(file,'D');D=L.D;
 if strcmp(stage,'verify'),verify_data(D,s,cfg,folder);continue;end
 dest=fullfile(folder,'fit');if ~isfolder(dest),mkdir(dest);end
 stamp=struct('prepared_sha256',b1_sha256({'file',file}),'options',cfg.fit, ...
  'fitter_sha256',b1_sha256({'file',which('qe_zlp_joint_fit')}), ...
  'kinematic_sha256',b1_sha256({'file',which('qe_kinematic_prefactor')}));
 ix=1:numel(D.units);
 if strcmp(stage,'pilot')
  [~,k]=min(abs([D.units.q_Ainv]-.001));ix=k;
 elseif strcmp(stage,'bootstrap')
  ix=find([D.units.N]==1 & any(abs(abs([D.units.q_Ainv].')-cfg.bootstrap.targets)<1e-10,2).');
 else
  ix=ix(mod(ix-1,nparts)==part-1);
 end
 for b=ix
  target=fullfile(dest,sprintf('bin_%03d.mat',b));
  if strcmp(stage,'bootstrap'),bootstrap(D,b,cfg,folder,target);continue;end
  if isfile(target)
   R=load(target,'stamp');assert(isequaln(R.stamp,stamp),'b1_low_branch_matlab:Stale','Use a new output folder');continue
  end
  tic_fit=tic;[E,Y,sg,args,K]=problem(D,b,cfg);fits=cell(1,2);
  for n=1:2
   fits{n}=qe_zlp_joint_fit(E,Y,args{:},n_peaks=n,noise_sigma=sg,prefactor=K);
   fits{n}.effective_options.prefactor=[];
  end
  rows=[compact(fits{1},D,b,1);compact(fits{2},D,b,2)];
  save(target,'fits','rows','stamp','-v7');
  fprintf('MATLAB %s %d/%d q=%+.5f N=%d %.1fs E0=%s chi2=%.3g test=%.3g\n', ...
   s.id,b,numel(D.units),D.units(b).q_Ainv,D.units(b).N,toc(tic_fit), ...
   mat2str([rows(2).low_E0_meV rows(2).high_E0_meV],5),rows(2).signal_msr,rows(2).test_msr);
 end
end
end

function prepare(s,cfg,file)
source=fullfile(cfg.data_root,'20260120 BiSb',s.folder,s.file);
signature=struct('source',source,'source_sha256',s.sha256,'session',s,'preparation',cfg.preparation, ...
 'alignment_code',b1_sha256({'file',which('qe_zlp_integer_align')}), ...
 'binning_code',b1_sha256({'file',which('qe_prepare_count_bins')}));
if isfile(file)
 L=load(file,'D');assert(isequaln(L.D.signature,signature),'b1_low_branch_matlab:Stale','Preparation changed');return
end
assert(strcmp(b1_sha256({'file',source}),s.sha256),'b1_low_branch_matlab:Source','Raw source hash changed');
raw=read_npy(source);assert(ndims(raw)==3&&size(raw,2)==512&&size(raw,3)==1028,'Unexpected raw axes');
meta=jsondecode(fileread(strrep(source,'.npy','.json')));ax=meta.spatial_calibrations;
eidx=find(contains(string({ax.units}),'eV'),1,'last');
Eraw=(ax(eidx).offset+(0:size(raw,3)-1)*ax(eidx).scale)'*1000;
q=((1:size(raw,2))-s.zero)*s.dq;refq=abs(q)<=cfg.preparation.reference_q+1e-12;
ref=squeeze(sum(double(raw(:,refq,:)),2)).';
search=find(Eraw>=cfg.preparation.raw_zlp_search(1)&Eraw<=cfg.preparation.raw_zlp_search(2));
[~,p]=max(sum(ref(search,:),2));origin=Eraw(search(p));
selected=find(abs(q)<=cfg.preparation.q_max+.002+1e-12);
X=permute(raw(:,selected,:),[3 1 2]);clear raw
assert(all(isfinite(X(:)))&&all(X(:)<double(intmax('uint32'))),'Invalid detector samples');
al=qe_zlp_integer_align(Eraw-origin,ref,X,cfg.preparation.alignment_window);clear X
assert(all(al.valid),'b1_low_branch_matlab:Frames','Invalid ZLP frames require review');
E=Eraw(al.support)-Eraw(al.reference_pixel);T=numel(al.valid);Q=q(selected);
F=al.aligned;al=rmfield(al,'aligned');
beam=zeros(1,T);
for t=1:T,beam(t)=sum(ref(al.support(abs(E)<=100)+al.measured_offset_pixels(t),t));end
assert(all(beam>0),'Nonpositive beam reference');
elastic=squeeze(sum(F(abs(E)<=cfg.preparation.q_centroid_energy,:,:),1));
central=abs(Q)<=cfg.preparation.q_centroid_range+1e-12;
mu=sum(elastic(:,central).*Q(central),2)./sum(elastic(:,central),2);
drift=movmedian(mu-median(mu),cfg.preparation.q_median_frames).';
if ~cfg.preparation.do_q_drift,drift(:)=0;end
profile=sum(elastic,1);[~,p]=max(profile);jj=p+(-2:2);
w=profile(jj)/sum(profile(jj));x=Q(jj);probe=sqrt(max(sum(w.*(x-sum(w.*x)).^2)-s.dq^2/12,0));
qe=struct('q_Ainv',Q,'dq_Ainv',s.dq,'intensity',ones(1,numel(Q)), ...
 'source_channel',selected,'energy_meV',0);
low=qe_prepare_count_bins(qe,1,q_range_Ainv=[-cfg.preparation.native_max cfg.preparation.native_max], ...
 q_skip_Ainv=cfg.preparation.q_min);
units=low.units;
for N=cfg.preparation.high_channels
 last=cfg.preparation.native_max+floor((cfg.preparation.q_max-cfg.preparation.native_max)/s.dq/N+1e-8)*N*s.dq;
 high=qe_prepare_count_bins(qe,N,q_range_Ainv=[-last last], ...
  q_skip_Ainv=cfg.preparation.native_max+s.dq/2);
 assert(~any([high.units.partial]),'Partial outer q bin');
 units=[units high.units]; %#ok<AGROW>
end
keep=E>=cfg.preparation.energy_range(1)&E<=cfg.preparation.energy_range(2);E=E(keep);F=F(keep,:,:);
frames=zeros(numel(E),T,numel(units),'single');max_weight_error=0;
for b=1:numel(units)
 width=units(b).q_width;center=units(b).q_Ainv;
 for t=1:T
  left=center+drift(t)-width/2;right=center+drift(t)+width/2;
  w=max(0,min(Q+s.dq/2,right)-max(Q-s.dq/2,left))/s.dq;
  max_weight_error=max(max_weight_error,abs(sum(w)-units(b).source_q_count));
  used=find(w>1e-12);frames(:,t,b)=single(reshape(F(:,t,used),numel(E),[])*w(used).');
 end
 units(b).N=units(b).source_q_count;
end
assert(max_weight_error<1e-8,'q bin transport lost counts');
D=struct('signature',signature,'session',s.id,'E',E,'q_native',Q,'units',units,'frames',frames, ...
 'beam',beam,'dq',s.dq,'h',s.h,'sigma_probe',probe,'drift_Ainv',drift, ...
 'alignment',al,'energy_origin_meV',Eraw(al.reference_pixel),'train_frames',cfg.preparation.train_frames, ...
 'max_weight_error',max_weight_error,'train_sigma',zeros(numel(E),numel(units)), ...
 'test_sigma',zeros(numel(E),numel(units)));
for b=1:numel(units)
 nt=D.train_frames;
 D.train_sigma(:,b)=b1_p1b_noise(frames(:,1:nt,b),E,beam(1:nt),zeros(0,2),'first_diff',21);
 D.test_sigma(:,b)=b1_p1b_noise(frames(:,nt+1:end,b),E,beam(nt+1:end),zeros(0,2),'first_diff',21);
end
save(file,'D','-v7.3');
fprintf('MATLAB PREPARED %s: %d energies, %d frames, %d bins, origin %.1f, probe %.5g\n',s.id,numel(E),T,numel(units),D.energy_origin_meV,probe);
end

function [E,Y,sg,args,K]=problem(D,b,cfg)
E=D.E;Y=sum(double(D.frames(:,1:D.train_frames,b)),2);sg=D.train_sigma(:,b);o=cfg.fit;
u=D.units(b);members=sscanf(u.source_q_Ainv,'%f,');
K=qe_kinematic_prefactor(members,dq=D.dq,h_perp=D.h,sigma_probe=D.sigma_probe, ...
 beam_kV=o.beam_kV,form='2d',n_sub=o.n_sub,n_probe=o.n_probe);
args={'fit_window',o.fit_window,'signal_window',o.signal_window,'peak_model',o.peak_model, ...
 'n_zlp',o.n_zlp,'aux_windows',o.aux_windows,'n_starts',o.n_starts, ...
 'seed',o.seed+b,'prefactor_floor_meV',o.prefactor_floor,'prefactor_on_aux',false, ...
 'peak_E0_bounds',o.peak_E0_bounds,'peak_width_max',o.peak_width_max,'max_iterations',o.max_iterations};
end

function row=compact(f,D,b,n)
row=struct('session',D.session,'bin',b,'q_Ainv',D.units(b).q_Ainv,'N',D.units(b).N, ...
 'n_peaks',n,'success',f.success,'status',f.numerical_status,'low_E0_meV',NaN,'high_E0_meV',NaN, ...
 'low_count_apex_meV',NaN,'low_loss_apex_meV',NaN,'low_width_meV',NaN,'signal_msr',NaN, ...
 'test_msr',NaN,'low_boundary',true,'low_area_fraction',NaN);
if ~f.success,return;end
row.low_E0_meV=f.parameters(1,1);row.low_width_meV=f.parameters(1,2);
if n==2,row.high_E0_meV=f.parameters(2,1);end
m=f.energy_meV>=300&f.energy_meV<=1800;
e=f.energy_meV(m);[~,i]=max(f.peaks(m,1));row.low_count_apex_meV=e(i);
[~,i]=max(f.loss_function(m,1));row.low_loss_apex_meV=e(i);
row.low_boundary=any(f.peak_boundary(1,:))||f.zero_amplitude(1);
row.low_area_fraction=sum(f.peaks(m,1))/max(sum(f.peaks(m,:),'all'),eps);
row.signal_msr=mean(((f.observed(m)-f.prediction(m))./f.sigma(m)).^2);
[tf,ix]=ismember(f.energy_meV,D.E);assert(all(tf),'Fit energy axis mismatch');
scale=sum(D.beam(D.train_frames+1:end))/sum(D.beam(1:D.train_frames));
y=sum(double(D.frames(ix,D.train_frames+1:end,b)),2);sg=D.test_sigma(ix,b);
row.test_msr=mean(((y(m)-scale*f.prediction(m))./sg(m)).^2);
end

function bootstrap(D,b,cfg,folder,target)
assert(isfile(target),'Run main fit first');L=load(target,'fits');base=L.fits{2};if ~base.success,return;end
dest=fullfile(folder,'bootstrap');if ~isfolder(dest),mkdir(dest);end
signature=struct('main_fit_sha256',b1_sha256({'file',target}), ...
 'options',cfg.bootstrap,'seed',cfg.seed+sum(double(D.session)));
sigfile=fullfile(dest,sprintf('bin_%03d_signature.mat',b));
if isfile(sigfile)
 L=load(sigfile,'signature');assert(isequaln(L.signature,signature),'Bootstrap configuration changed');
else
 save(sigfile,'signature');
end
[E,~,sg,args,K]=problem(D,b,cfg);j=find(strcmp(args(1:2:end),'n_starts'));args{2*j}=cfg.bootstrap.n_starts;
T=D.train_frames;rs=RandStream('mt19937ar','Seed',signature.seed);F=double(D.frames(:,1:T,b));rows=struct([]);
draws=zeros(T,cfg.bootstrap.repeats);
for r=1:cfg.bootstrap.repeats
 starts=randi(rs,T,ceil(T/cfg.bootstrap.block_frames),1)-1;
 ix=mod(starts+(0:cfg.bootstrap.block_frames-1),T)+1;ix=reshape(ix.',[],1);draws(:,r)=ix(1:T);
end
for r=1:cfg.bootstrap.repeats
 file=fullfile(dest,sprintf('bin_%03d_rep_%03d.mat',b,r));
 if isfile(file),R=load(file,'row');row=R.row;
 else
  draw=draws(:,r);
  Y=sum(F(:,draw),2);sboot=sg.*sqrt(max(Y,1)./max(sum(F,2),1));
  f=qe_zlp_joint_fit(E,Y,args{:},n_peaks=2,noise_sigma=sboot,prefactor=K,warm_u=base.u);
  dose=D;dose.beam(1:T)=D.beam(draw);
  row=compact(f,dose,b,2);row.replicate=r;save(file,'row','draw','-v7');
 end
 if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
writetable(struct2table(rows),fullfile(dest,sprintf('bin_%03d.csv',b)));
fprintf('MATLAB BOOTSTRAP %s q=%+.5f %d/%d successful\n',D.session,D.units(b).q_Ainv,nnz([rows.success]),numel(rows));
end

function report(out_root,cfg)
allrows=struct([]);alltables={};rises=struct([]);
for s=cfg.session
 folder=fullfile(out_root,s.id);files=dir(fullfile(folder,'fit','bin_*.mat'));if isempty(files),continue;end
 rows=struct([]);
 for k=1:numel(files),L=load(fullfile(files(k).folder,files(k).name),'rows');rows=[rows;L.rows];end %#ok<AGROW>
 allrows=[allrows;rows]; %#ok<AGROW>
 L=load(fullfile(folder,'prepared.mat'),'D');D=L.D;
 t=struct2table(rows);
 t.low_count_censored=t.low_count_apex_meV<=min(D.E(D.E>=300))+1e-8| ...
  t.low_count_apex_meV>=max(D.E(D.E<=1800))-1e-8;
 t.low_loss_censored=t.low_loss_apex_meV<=min(D.E(D.E>=300))+1e-8| ...
  t.low_loss_apex_meV>=max(D.E(D.E<=1800))-1e-8;
 t.components_merge=t.n_peaks==2 & t.high_E0_meV-t.low_E0_meV<cfg.review.merge_gap_meV;
 t.low_component_weak=t.n_peaks==2 & t.low_area_fraction<cfg.review.weak_area_fraction;
 t.test_gain_fraction=nan(height(t),1);
 t.bootstrap_successes=zeros(height(t),1);t.low_count_lo_meV=nan(height(t),1);
 t.low_count_hi_meV=nan(height(t),1);t.low_E0_lo_meV=nan(height(t),1);t.low_E0_hi_meV=nan(height(t),1);
 for j=find(t.n_peaks==2).'
  h0=t(t.bin==t.bin(j)&t.n_peaks==1,:);
  if ~isempty(h0),t.test_gain_fraction(j)=1-t.test_msr(j)/h0.test_msr;end
  bootfile=fullfile(folder,'bootstrap',sprintf('bin_%03d.csv',t.bin(j)));
  if ~isfile(bootfile),continue;end
  B=readtable(bootfile);ok=logical(B.success);
  t.bootstrap_successes(j)=nnz(ok);
  if nnz(ok)<2,continue;end
  c=prctile(B.low_count_apex_meV(ok),[2.5 97.5]);p=prctile(B.low_E0_meV(ok),[2.5 97.5]);
  t.low_count_lo_meV(j)=c(1);t.low_count_hi_meV(j)=c(2);t.low_E0_lo_meV(j)=p(1);t.low_E0_hi_meV(j)=p(2);
 end
 writetable(t,fullfile(folder,'fit_summary.csv'));
 alltables{end+1}=t; %#ok<AGROW>
 fig=figure('Visible','off','Position',[50 50 1200 450]);tl=tiledlayout(fig,1,2);
 for side=[1 -1]
  ax=nexttile(tl);hold(ax,'on');
  colors=[0 .447 .741;.850 .325 .098;.466 .674 .188];
  for N=[1 3 5]
   m=t.n_peaks==2&sign(t.q_Ainv)==side&t.N==N&t.success;
   scatter(ax,abs(t.q_Ainv(m)),t.low_count_apex_meV(m)/1000,30,colors([1 3 5]==N,:), ...
    'DisplayName',sprintf('N=%d',N));
   flagged=m&(t.low_boundary|t.low_count_censored|t.components_merge|t.low_component_weak);
   scatter(ax,abs(t.q_Ainv(flagged)),t.low_count_apex_meV(flagged)/1000,60,'k','x','HandleVisibility','off');
   ids=find(m&isfinite(t.low_count_lo_meV));
   for j=ids.'
    plot(ax,repmat(abs(t.q_Ainv(j)),1,2),[t.low_count_lo_meV(j) t.low_count_hi_meV(j)]/1000, ...
     'Color',[.2 .2 .2],'HandleVisibility','off');
   end
  end
  xlabel(ax,'|q| (A^{-1})');ylabel(ax,'Low-component count apex (eV)');grid(ax,'on');
  ylim(ax,[.3 1.8]);xlim(ax,[0 cfg.preparation.q_max]);
  title(ax,sprintf('%s signed q %+d',s.id,side));legend(ax,'Location','best');
 end
 title(tl,sprintf('MATLAB two-B1 best fits; x: boundary/edge/merged/weak (<%.0f%%); bars: selected bootstrap intervals', ...
  100*cfg.review.weak_area_fraction),'FontSize',13);
 exportgraphics(fig,fullfile(folder,'low_component_counts.png'),'Resolution',180);close(fig);
 % Actual held-out spectra and the two fitted count components.
 fig=figure('Visible','off','Position',[50 50 1600 420]);tl=tiledlayout(fig,1,5);
 for q=[.001 .0025 .005 .010 .015]
  tw=t(t.n_peaks==2&t.N~=3,:);[~,j]=min(abs(tw.q_Ainv-q));b=tw.bin(j);
  F=load(fullfile(folder,'fit',sprintf('bin_%03d.mat',b)),'fits');f=F.fits{2};ax=nexttile(tl);hold(ax,'on');
  if ~f.success,title(ax,'Fit failed');continue;end
  m=f.energy_meV>=300&f.energy_meV<=1800;e=f.energy_meV(m)/1000;
  scale=sum(D.beam(D.train_frames+1:end))/sum(D.beam(1:D.train_frames));
  [~,ix]=ismember(f.energy_meV,D.E);y=sum(double(D.frames(ix,D.train_frames+1:end,b)),2)/scale;
  z=max(f.prediction(m));plot(ax,e,y(m)/z,'.','Color',[.6 .6 .6],'DisplayName','Held-out counts');
  plot(ax,e,f.prediction(m)/z,'k','LineWidth',1.2,'DisplayName','Total');
  plot(ax,e,f.peaks(m,1)/z,'Color',[0 .447 .741],'LineWidth',1.2,'DisplayName','Low component');
  plot(ax,e,f.peaks(m,2)/z,'Color',[.850 .325 .098],'LineWidth',1.2,'DisplayName','High component');
  plot(ax,e,f.background(m)/z,':','Color',[.466 .674 .188],'DisplayName','ZLP + aux');
  title(ax,sprintf('q=%+.4f N=%d',D.units(b).q_Ainv,D.units(b).N));xlabel(ax,'Loss (eV)');
  ylim(ax,[0 1.2]);xlim(ax,[.3 1.8]);
 end
 lg=legend(ax,'Orientation','horizontal','NumColumns',5);lg.Layout.Tile='south';
 ylabel(tl,'Counts / fitted maximum');title(tl,[s.id ' MATLAB held-out raw counts and component fits']);
 exportgraphics(fig,fullfile(folder,'component_spectra.png'),'Resolution',180);close(fig);
 for side=[1 -1]
  a=find(t.n_peaks==2&t.N==1&abs(t.q_Ainv-side*.001)<1e-10);
  b=find(t.n_peaks==2&t.N==1&abs(t.q_Ainv-side*.005)<1e-10);
  if isempty(a)||isempty(b),continue;end
  row=struct('session',s.id,'side',side,'count_rise_meV',t.low_count_apex_meV(b)-t.low_count_apex_meV(a), ...
   'E0_rise_meV',t.low_E0_meV(b)-t.low_E0_meV(a),'bootstrap_pairs',0, ...
   'count_rise_lo_meV',NaN,'count_rise_hi_meV',NaN,'E0_rise_lo_meV',NaN,'E0_rise_hi_meV',NaN, ...
   'endpoint_review_required',any(~t.success([a b])|t.low_component_weak([a b])| ...
   t.low_boundary([a b])|t.low_count_censored([a b])|t.components_merge([a b])), ...
   'min_endpoint_area_fraction',min(t.low_area_fraction([a b])));
  fa=fullfile(folder,'bootstrap',sprintf('bin_%03d.csv',t.bin(a)));
  fb=fullfile(folder,'bootstrap',sprintf('bin_%03d.csv',t.bin(b)));
  if isfile(fa)&&isfile(fb)
   A=readtable(fa);B=readtable(fb);assert(isequal(A.replicate,B.replicate),'Bootstrap pairing changed');
   ok=logical(A.success)&logical(B.success);row.bootstrap_pairs=nnz(ok);
   if nnz(ok)>1
    c=prctile(B.low_count_apex_meV(ok)-A.low_count_apex_meV(ok),[2.5 97.5]);
    p=prctile(B.low_E0_meV(ok)-A.low_E0_meV(ok),[2.5 97.5]);
    row.count_rise_lo_meV=c(1);row.count_rise_hi_meV=c(2);row.E0_rise_lo_meV=p(1);row.E0_rise_hi_meV=p(2);
   end
  end
  if isempty(rises),rises=row;else,rises(end+1)=row;end %#ok<AGROW>
 end
end
if ~isempty(allrows)
 writetable(struct2table(allrows),fullfile(out_root,'best_fit_summary.csv'));
 writetable(vertcat(alltables{:}),fullfile(out_root,'fit_summary.csv'));
 if ~isempty(rises),writetable(struct2table(rises),fullfile(out_root,'low_q_rise.csv'));end
end
end

function verify_data(D,s,cfg,folder)
source=fullfile(cfg.data_root,'20260120 BiSb',s.folder,s.file);
assert(strcmp(b1_sha256({'file',source}),s.sha256),'Raw data changed');
assert(all(isfinite(D.frames(:)))&&all(D.train_sigma(:)>0)&&all(D.test_sigma(:)>0),'Invalid prepared data');
assert(D.max_weight_error<1e-8,'Bin weight mismatch');
files=dir(fullfile(folder,'fit','bin_*.mat'));n=0;
assert(numel(files)==numel(D.units),'Main fit is incomplete');
expected=struct('prepared_sha256',b1_sha256({'file',fullfile(folder,'prepared.mat')}), ...
 'options',cfg.fit,'fitter_sha256',b1_sha256({'file',which('qe_zlp_joint_fit')}), ...
 'kinematic_sha256',b1_sha256({'file',which('qe_kinematic_prefactor')}));
for k=1:numel(files)
 L=load(fullfile(files(k).folder,files(k).name),'fits','rows','stamp');
 assert(isequaln(L.stamp,expected),'Saved fit signature does not match current configuration');
 for j=1:2
  f=L.fits{j};if ~f.success,continue;end
  assert(max(abs(f.prediction-f.background-sum(f.peaks,2)))<1e-8*max(f.prediction),'Component readback failed');
  assert(all(diff(f.parameters(:,1))>=0),'Component energy ordering failed');n=n+1;
 end
end
verification=struct('language','MATLAB','raw_hash_unchanged',true,'prepared_bins',numel(D.units), ...
 'models_attempted',2*numel(files),'models_read_back',n,'q_weight_error',D.max_weight_error,'gui_verified',false);
save(fullfile(folder,'verification.mat'),'verification');
fprintf('MATLAB VERIFIED %s: %d models, raw hash unchanged\n',s.id,n);
end
