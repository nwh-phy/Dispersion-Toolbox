function file = b1_fullq_build(session, run_root)
% B1 full-q task, stage 1: per-frame ZLP-aligned N=3 bins (A1 scheme of v7:
% common integer shift per frame from the |q| <= 0.001 reference in
% [-100, 100] meV) on the grid |q| = 0.0025 + 0.0015 k, both signs, plus the
% measured noise model. Writes <run_root>/stage1/frames_<session>.mat.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
switch session
 case '590_PL2_10w', raw_path=fullfile(root,'20260120 BiSb','590 PL2 10w 0.004 10sx300','Sequence EELS Image 274.npy');
 case 'n0_PL2_10w_repeat', raw_path=fullfile(root,'20260120 BiSb','n0 pl2 10w 0.004 10s x300','Sequence EELS Image 272.npy');
 otherwise, error('b1_fullq_build:Session','Unknown session %s',session);
end
cfg=b1_fullq_config(); out=fullfile(run_root,'stage1'); if ~isfolder(out), mkdir(out); end
l=load(fullfile(p2,session,'L1_minimal.mat')); qe=l.d.qe; E=qe.energy_meV(:);
grid=cfg.q_grid; targets=[-fliplr(grid) grid];
raw=read_npy(raw_path); T=size(raw,1);
assert(size(raw,2)==numel(qe.q_Ainv)&&size(raw,3)==numel(E),'Raw axes mismatch');
invalid=false(1,size(raw,2));
for q=1:size(raw,2), x=raw(:,q,:); invalid(q)=any(~isfinite(x(:))|double(x(:))>=double(intmax('uint32'))); end
invalid=invalid|any(~isfinite(qe.intensity),1);
bins=centered_bins(qe,targets,invalid); ok=[bins.valid];
fprintf('%s: %d/%d bins valid; invalid targets %s\n',session,nnz(ok),numel(ok),mat2str(targets(~ok)));
bins=bins(ok);
refq=find(abs(qe.q_Ainv)<=cfg.reference_abs_q & ~invalid);
reference=squeeze(sum(double(raw(:,refq,:)),2)).';
X=zeros(numel(E),T,numel(bins));
for b=1:numel(bins), X(:,:,b)=squeeze(sum(double(raw(:,bins(b).source_channel,:)),2)).'; end
elastic=squeeze(sum(sum(double(raw(:,:,abs(E)<=12)),1),3)); elastic(invalid)=0;
clear raw
al=qe_zlp_integer_align(E,reference,X,cfg.zlp_window); clear X
assert(all(al.valid),'Frames without a valid ZLP reference');
Ea=al.E(:); F=al.aligned; Ysum=squeeze(sum(F,2));

% Noise: per-frame variance / mean after removing the common-mode intensity
% (each frame scaled to the mean total), lag correlations along E.
nb=numel(bins); ratio=nan(numel(Ea),nb); ratio_raw=ratio; rho=nan(3,nb,2); drift=nan(T,nb);
lossr=Ea>=300&Ea<=1800; gainr=Ea<=-60;
for b=1:nb
 Fb=F(:,:,b); s=sum(Fb,1); drift(:,b)=s/mean(s);
 m=mean(Fb,2); ratio_raw(:,b)=var(Fb,0,2)./max(m,eps);
 Fc=Fb.*(mean(s)./s); vc=var(Fc,0,2); ratio(:,b)=vc./max(m,eps);
 R=(Fc-mean(Fc,2))./sqrt(max(vc,eps));
 for k=1:3
  c=mean(R(1:end-k,:).*R(1+k:end,:),2);
  rho(k,b,1)=median(c(lossr(1:end-k))); rho(k,b,2)=median(c(gainr(1:end-k)));
 end
end
g=movmedian(ratio,21,1,'omitnan'); g(~isfinite(g)|g<=0)=1;
sigma=sqrt(g.*max(Ysum,1));

% Probe angular spread from the elastic q-profile (central five channels).
[~,i0]=max(elastic); k=-2:2; w=elastic(i0+k); w=w/sum(w); mu=sum(w.*k);
sigma_probe=sqrt(max(sum(w.*(k-mu).^2)-1/12,0))*qe.dq_Ainv;

D=struct('session',session,'E',Ea,'bins',bins,'q_center',[bins.q_Ainv],'Ysum',Ysum,'sigma',sigma, ...
 'noise_ratio',ratio,'noise_ratio_raw',ratio_raw,'noise_rho_loss',squeeze(rho(:,:,1)),'noise_rho_gain',squeeze(rho(:,:,2)), ...
 'frame_drift',drift,'offsets',al.measured_offset_pixels,'support',al.support,'sigma_probe',sigma_probe, ...
 'elastic_q_profile',elastic,'q_axis',qe.q_Ainv,'T',T,'raw_path',raw_path);
file=fullfile(out,sprintf('frames_%s.mat',session));
frames=single(F); %#ok<NASGU>
save(file,'D','frames','-v7.3');
fprintf('saved %s: %d bins, E %g..%g meV, sigma_probe %.2e A^-1\n',file,nb,Ea(1),Ea(end),sigma_probe);
fprintf('noise ratio (loss 300-1800, median over bins) %.3f, raw %.3f; lag-1/2/3 rho loss %s gain %s\n', ...
 median(ratio(lossr,:),'all','omitnan'),median(ratio_raw(lossr,:),'all','omitnan'), ...
 mat2str(round(median(rho(:,:,1),2,'omitnan')',3)),mat2str(round(median(rho(:,:,2),2,'omitnan')',3)));
end


function bins=centered_bins(qe,targets,invalid)
% N=3 bins centred on a native channel: members c-1..c+1, same sign, q ~= 0,
% no invalid channel. (qe_prepare_count_bins is limited to |q| <= 0.015.)
q=qe.q_Ainv(:).'; src=1:numel(q); if isfield(qe,'source_channel'), src=qe.source_channel(:).'; end
bins=struct('target_q',{},'valid',{},'q_Ainv',{},'source_channel',{},'q_members',{});
for t=targets
 [d,c]=min(abs(q-t)); ii=c-1:c+1; okb=d<qe.dq_Ainv*1e-6 && ii(1)>=1 && ii(end)<=numel(q);
 if okb, okb=all(sign(q(ii))==sign(t)) && all(q(ii)~=0) && ~any(invalid(src(ii))); end
 b=struct('target_q',t,'valid',okb,'q_Ainv',NaN,'source_channel',[],'q_members',[]);
 if okb, b.q_Ainv=mean(q(ii)); b.source_channel=src(ii); b.q_members=q(ii); end
 bins(end+1)=b; %#ok<AGROW>
end
end
