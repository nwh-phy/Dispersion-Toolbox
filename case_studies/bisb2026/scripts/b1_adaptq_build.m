function file = b1_adaptq_build(session, run_root, qmax)
% B1 adaptive-binning task, stage 1: every native q channel with
% 0 < |q| <= qmax (both signs), per frame, ZLP-aligned with the A1 scheme of
% v7 / b1_fullq (common integer shift per frame from the |q| <= 0.001
% reference in [-100, 100] meV). The shift is common to all channels of a
% frame, so any later binning is a plain sum over aligned channels.
% Writes <run_root>/stage1/channels_<session>.mat.
if nargin < 3, qmax = 0.035; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
switch session
 case '590_PL2_10w', raw_path=fullfile(root,'20260120 BiSb','590 PL2 10w 0.004 10sx300','Sequence EELS Image 274.npy');
 case 'n0_PL2_10w_repeat', raw_path=fullfile(root,'20260120 BiSb','n0 pl2 10w 0.004 10s x300','Sequence EELS Image 272.npy');
 case 'no_PL2_20w_2film', raw_path=fullfile(root,'20260120 BiSb','no pl2 20w 0.004 10sx300 2film','Sequence EELS Image 273.npy');
 otherwise, error('b1_adaptq_build:Session','Unknown session %s',session);
end
cfg=b1_fullq_config(); out=fullfile(run_root,'stage1'); if ~isfolder(out), mkdir(out); end
l=load(fullfile(p2,session,'L1_minimal.mat')); qe=l.d.qe; E=qe.energy_meV(:);
raw=read_npy(raw_path); T=size(raw,1);
assert(size(raw,2)==numel(qe.q_Ainv)&&size(raw,3)==numel(E),'Raw axes mismatch');
src=1:size(raw,2); if isfield(qe,'source_channel'), src=qe.source_channel(:).'; end
q=qe.q_Ainv(:).';
invalid=false(1,size(raw,2));
for c=1:size(raw,2), x=raw(:,c,:); invalid(c)=any(~isfinite(x(:))|double(x(:))>=double(intmax('uint32'))); end
invalid=invalid|any(~isfinite(qe.intensity),1);
ch=find(abs(q)>0 & abs(q)<=qmax+1e-9 & ~invalid(src));
refq=find(abs(q)<=cfg.reference_abs_q & ~invalid(src));
reference=squeeze(sum(double(raw(:,src(refq),:)),2)).';
X=zeros(numel(E),T,numel(ch));
for k=1:numel(ch), X(:,:,k)=squeeze(double(raw(:,src(ch(k)),:))).'; end
elastic=squeeze(sum(sum(double(raw(:,:,abs(E)<=12)),1),3)); elastic(invalid)=0;
clear raw
al=qe_zlp_integer_align(E,reference,X,cfg.zlp_window); clear X
assert(all(al.valid),'Frames without a valid ZLP reference');

% Probe angular spread from the elastic q-profile (central five channels), as b1_fullq_build.
[~,i0]=max(elastic); kk=-2:2; w=elastic(i0+kk); w=w/sum(w); mu=sum(w.*kk);
sigma_probe=sqrt(max(sum(w.*(kk-mu).^2)-1/12,0))*qe.dq_Ainv;

C=struct('session',session,'E',al.E(:),'q',q(ch),'source_channel',src(ch),'dq',qe.dq_Ainv, ...
 'offsets',al.measured_offset_pixels,'support',al.support,'sigma_probe',sigma_probe, ...
 'elastic_q_profile',elastic,'q_axis',q,'T',T,'raw_path',raw_path,'qmax',qmax);
frames=single(al.aligned); %#ok<NASGU>
file=fullfile(out,sprintf('channels_%s.mat',session));
save(file,'C','frames','-v7.3');
fprintf('saved %s: %d channels (%d per side), E %g..%g meV, sigma_probe %.2e A^-1\n', ...
 file,numel(ch),nnz(q(ch)>0),C.E(1),C.E(end),sigma_probe);
end
