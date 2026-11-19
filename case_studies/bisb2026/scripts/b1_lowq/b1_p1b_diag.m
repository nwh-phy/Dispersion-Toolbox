function out = b1_p1b_diag(run_root, session, cfg)
% P1b diagnostics from the frame data (no fits):
%  (T3) noise coefficient g(E) = var/mean of the frames for four estimators:
%       raw, common-beam corrected (whole-frame reference), per-channel corrected
%       (the model used up to P1) and first differences of the common-beam
%       corrected frames (removes slow drifts);
%  (T7) near-zero-loss 2D response including the central channel: ZLP centre
%       c(q), FWHM and tails, and separability I_el(q,E) ~ P(q) R(E) after the
%       c(q) bending is removed;
%  (T8) frame-to-frame relations of the elastic core, low-energy bands and the
%       B1 band (common-beam normalised) with ZLP position and width.
% The common-beam reference (all valid channels, all energies except the
% 150-250 meV holdout) is saved for the P1b fits.
% Output: <run_root>/p1b/diag_<session>/.
if nargin<1, run_root=[]; end
if nargin<2, session=[]; end
if nargin<3, cfg=struct('session',session); end
paths=b1_forward_paths(cfg,run_root); session=paths.session;
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
out=paths.diag_dir; if ~isfolder(out), mkdir(out); end
ho=[150 250]; if isfield(cfg,'holdout'), ho=cfg.holdout; end
cf=paths.stage1_file; S=load(cf,'C'); C=S.C; mf=matfile(cf);
assert(strcmp(C.session,session),'Stage1 has another session.');
E=C.E(:); T=C.T; off=C.offsets(:); sup=C.support(:);
prepared=isfield(C,'native_prepared') && C.native_prepared;
hidx=sup(E>=ho(1)&E<=ho(2));
b1idx=sup(E>=300&E<=1800);
ref=zeros(T,1); refB1=zeros(T,1);
if prepared
 D=load(cf,'detector_counts'); validq=C.valid_channels;
 for t=1:T
  y=D.detector_counts(:,t);
  ref(t)=sum(y)-sum(y(hidx+off(t)));
  refB1(t)=sum(y(b1idx+off(t)));
 end
else
 d=bisb_datasets(session);
 raw=read_npy(d.path);
 l=load(paths.l1_file); qe=l.d.qe; q=qe.q_Ainv(:).';
 src=1:size(raw,2); if isfield(qe,'source_channel'), src=qe.source_channel(:).'; end
 invalid=false(1,size(raw,2));
 for c=1:size(raw,2), x=raw(:,c,:); invalid(c)=any(~isfinite(x(:))|double(x(:))>=double(intmax('uint32'))); end
 invalid=invalid|any(~isfinite(qe.intensity),1);
 validq=setdiff(src(~invalid(src)),[]);
 for t=1:T
  fr=double(squeeze(raw(t,validq,:)));
  ref(t)=sum(fr(:))-sum(sum(fr(:,hidx+off(t))));
  refB1(t)=sum(sum(fr(:,b1idx+off(t))));
 end
end
assert(all(isfinite(ref) & ref>0),'Invalid common-beam reference.');
save(fullfile(out,'common_ref.mat'),'ref','refB1','ho','validq');

% ---- (T7) near-zero-loss map, |q| <= 0.006 including the central channel
ew=find(E>=-200&E<=300); Ez=E(ew);
if prepared
 [A,qs]=b1_forward_stage_frames(cf,C,0.006);
 I=sum(A(ew,:,:),3);
 Z0=reshape(A(ew,abs(qs)<1e-9,:),numel(ew),T);
 clear A
else
 sel=find(abs(q)<=0.006+1e-9 & ~invalid(src)); [~,o]=sort(q(sel)); sel=sel(o); qs=q(sel);
 I=zeros(numel(ew),numel(qs));
 for t=1:T
  fr=double(squeeze(raw(t,src(sel),:)));
  I=I+fr(:,sup(ew)+off(t)).';
 end
 c0=find(abs(q)<1e-9 & ~invalid(src)); assert(~isempty(c0),'no central channel');
 Z0=zeros(numel(ew),T);
 for t=1:T
  y=double(squeeze(raw(t,src(c0),:))); Z0(:,t)=y(sup(ew)+off(t));
 end
 clear raw
end
cq=nan(1,numel(qs)); fw=cq; tg=cq; tl=cq;
for k=1:numel(qs)
 y=I(:,k); [~,i0]=max(y.*(abs(Ez)<40));
 if i0>1 && i0<numel(y)
  ly=log(y(i0-1:i0+1)); d=(ly(3)-ly(1))/(2*(2*ly(2)-ly(1)-ly(3))); cq(k)=Ez(i0)+d*(Ez(2)-Ez(1));
 end
 h=y(i0)/2; l1=find(y(1:i0)<h,1,'last'); r1=i0-1+find(y(i0:end)<h,1,'first');
 if ~isempty(l1)&&~isempty(r1), fw(k)=interp1(y([l1 l1+1]),Ez([l1 l1+1]),h)-0; fw(k)=interp1(y([r1-1 r1]),Ez([r1-1 r1]),h)-fw(k); end
 tg(k)=mean(y(Ez>=-150&Ez<=-60))/y(i0); tl(k)=mean(y(Ez>=60&Ez<=150))/y(i0);
end
% separability after removing c(q): rows = q, columns = x = E - c(q)
xg=(-180:2:280).'; R=nan(numel(xg),numel(qs));
for k=1:numel(qs)
 if isnan(cq(k)), continue; end
 y=interp1(Ez-cq(k),I(:,k),xg,'linear',NaN); R(:,k)=y/max(y(abs(xg)<=4));
end
win=struct('name',{'core |x|<=30','gain -180..-40','loss 40..150','loss 150..280'},'r',{[-30 30],[-180 -40],[40 150],[150 280]});
sepTab=table(); kc=find(abs(qs)<1e-9);
for w=1:numel(win)
 m=xg>=win(w).r(1)&xg<=win(w).r(2); A=R(m,:); A=A(:,all(isfinite(A),1));
 s=svd(A); frac=s(1)^2/sum(s.^2);
 dev=max(abs(R(m,:)./R(m,kc)-1),[],1,'omitnan');
 sepTab=[sepTab; table(string(win(w).name),frac,median(dev,'omitnan'),max(dev),'VariableNames',{'window','rank1_fraction','median_maxdev_vs_q0','worst_maxdev_vs_q0'})]; %#ok<AGROW>
end
writetable(sepTab,fullfile(out,'nearzlp_separability.csv'));
writetable(table(qs(:),cq(:),fw(:),tg(:),tl(:),'VariableNames',{'q','zlp_centre_meV','zlp_fwhm_meV','gain_tail_-150_-60_over_peak','loss_60_150_over_peak'}),fullfile(out,'nearzlp_per_channel.csv'));
f=figure('Visible','off','Position',[0 0 1700 950]); tlo=tiledlayout(f,2,2,'TileSpacing','compact');
ax=nexttile(tlo); imagesc(ax,qs,Ez,log10(max(I,1))); axis(ax,'xy'); colorbar(ax); hold(ax,'on'); plot(ax,qs,cq,'w.-');
xlabel(ax,'q (A^{-1})'); ylabel(ax,'E (meV)'); title(ax,'log_{10} counts near zero loss (all frames, aligned); white: ZLP centre c(q)');
ax=nexttile(tlo); yyaxis(ax,'left'); plot(ax,qs,cq,'o-'); ylabel(ax,'c(q) (meV)'); yyaxis(ax,'right'); plot(ax,qs,fw,'s-'); ylabel(ax,'FWHM (meV)');
xlabel(ax,'q (A^{-1})'); grid(ax,'on'); title(ax,'ZLP centre and width per channel');
ax=nexttile(tlo); hold(ax,'on'); cm=turbo(numel(qs));
for k=1:numel(qs), semilogy(ax,xg,R(:,k),'Color',cm(k,:)); end
set(ax,'YScale','log'); xlim(ax,[-180 280]); ylim(ax,[1e-6 2]); grid(ax,'on'); xlabel(ax,'E - c(q) (meV)'); ylabel(ax,'normalised to peak');
title(ax,'line shapes after removing c(q) (colour: q from - to +)'); colormap(ax,turbo); cb=colorbar(ax); clim(ax,[min(qs) max(qs)]); cb.Label.String='q (A^{-1})';
ax=nexttile(tlo); hold(ax,'on'); xs=[-120 -60 -30 -10 0 10 30 60 120]; cm2=lines(numel(xs));
for j=1:numel(xs), [~,ix]=min(abs(xg-xs(j))); pr=R(ix,:)./max(R(ix,:)); plot(ax,qs,pr,'o-','Color',cm2(j,:),'DisplayName',sprintf('E-c = %+d meV',xs(j))); end
set(ax,'YScale','log'); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'q profile, normalised'); legend(ax,'Location','southoutside','NumColumns',5,'FontSize',7);
title(ax,'q profiles of the line-shape-normalised map at fixed E - c (equal if separable)');
exportgraphics(f,fullfile(out,'nearzlp_2d.png'),'Resolution',120); close(f);

% ---- per-frame ZLP descriptors from the central channel
zc=nan(T,1); zw=zc;
for t=1:T
 y=Z0(:,t); [~,i0]=max(y.*(abs(Ez)<40));
 ly=log(max(y(i0-1:i0+1),1)); d=(ly(3)-ly(1))/(2*(2*ly(2)-ly(1)-ly(3))); zc(t)=Ez(i0)+d*(Ez(2)-Ez(1));
 h=y(i0)/2; l1=find(y(1:i0)<h,1,'last'); r1=i0-1+find(y(i0:end)<h,1,'first');
 if ~isempty(l1)&&~isempty(r1), zw(t)=interp1(y([r1-1 r1]),Ez([r1-1 r1]),h)-interp1(y([l1 l1+1]),Ez([l1 l1+1]),h); end
end

% ---- (T8) band intensities per frame, 8 channels, and (T3) noise coefficients
qlist=[0.0005 0.001 0.002 0.005]; if isfield(cfg,'noise_q'), qlist=cfg.noise_q; end
bands=struct('name',{'core |E|<=20','gain -150..-40','low 40..150','holdout 150..250','250..400','B1 400..1800'},'r',{[-20 20],[-150 -40],[40 150],[150 250],[250 400],[400 1800]});
corr_rows={}; gtab={}; Eq=[-100 45 100 200 392 800 1500];
fT=figure('Visible','off','Position',[0 0 1800 1000]); tlT=tiledlayout(fT,2,numel(qlist),'TileSpacing','compact');
fG=figure('Visible','off','Position',[0 0 1800 1000]); tlG=tiledlayout(fG,2,numel(qlist),'TileSpacing','compact');
if prepared
 q_noise=reshape([-qlist; qlist],1,[]);
 sigma=zeros(numel(E),numel(q_noise)); g=sigma;
end
for qa=qlist
 for sg=[-1 1]
  k=find(abs(C.q-sg*qa)<1e-7);
  assert(numel(k)==1,'Missing P1b diagnostic channel %+.6f A^-1.',sg*qa);
  Fb=double(mf.frames(:,:,k));          % E x T
  if prepared
   col=find(q_noise==sg*qa);
   [sigma(:,col),g(:,col)]=b1_p1b_noise(Fb,E,ref,ho,'first_diff',21);
  end
  % noise coefficients
  m=mean(Fb,2); s_ch=sum(Fb,1);
  Fcom=Fb.*(mean(ref)./ref.'); Fch=Fb.*(mean(s_ch)./s_ch);
  g_raw=var(Fb,0,2)./max(m,eps); g_com=var(Fcom,0,2)./max(m,eps); g_ch=var(Fch,0,2)./max(m,eps);
  g_dif=var(diff(Fcom,1,2),0,2)/2./max(m,eps);
  sm=@(g) movmedian(g,21,1,'omitnan');
  G=[sm(g_raw) sm(g_com) sm(g_ch) sm(g_dif)];
  for e=Eq, [~,ie]=min(abs(E-e)); gtab(end+1,:)={sg*qa,e,G(ie,1),G(ie,2),G(ie,3),G(ie,4)}; end %#ok<AGROW>
  ax=nexttile(tlG); semilogy(ax,E,G,'LineWidth',1); grid(ax,'on'); xlim(ax,[-180 1800]); ylim(ax,[0.5 50]);
  xregion(ax,ho(1),ho(2)); title(ax,sprintf('q=%+.4f: g(E)=var/mean (21-pt median)',sg*qa),'FontSize',9);
  if qa==qlist(1)&&sg<0, legend(ax,{'raw','common beam','per channel (used to P1)','first differences, common beam'},'FontSize',7,'Location','north'); end
  % band intensities, common-beam normalised
  Bv=zeros(T,numel(bands));
  for b=1:numel(bands), w=E>=bands(b).r(1)&E<=bands(b).r(2); Bv(:,b)=sum(Fcom(w,:),1).'; end
  Z=(Bv-mean(Bv))./std(Bv);
  X=[off zc zw]; Xn={'ZLP offset (px)','ZLP centre in frame (meV)','ZLP FWHM (meV)'};
  for b=1:numel(bands)
   rr=[pcorr(Bv(:,b),X(:,1)) pcorr(Bv(:,b),X(:,2)) pcorr(Bv(:,b),X(:,3)) pcorr(Bv(:,b),Bv(:,1)) pcorr(Bv(:,b),Bv(:,6))];
   corr_rows(end+1,:)={sg*qa,bands(b).name,std(Bv(:,b))/mean(Bv(:,b)),sqrt(mean(Bv(:,b)))/mean(Bv(:,b)),rr(1),rr(2),rr(3),rr(4),rr(5)}; %#ok<AGROW>
  end
  ax=nexttile(tlT); hold(ax,'on'); nb=10; tb=reshape(1:floor(T/nb)*nb,nb,[]);
  for b=[1 3 4 5 6], zb=mean(reshape(Z(tb(:),b),nb,[]),1); plot(ax,mean(tb,1),zb,'DisplayName',bands(b).name); end
  plot(ax,mean(tb,1),mean(reshape((zw(tb(:))-mean(zw,'omitnan'))/std(zw,'omitnan'),nb,[]),1),'k--','DisplayName','ZLP FWHM (q=0)');
  grid(ax,'on'); xlabel(ax,'frame'); ylabel(ax,'z-score (10-frame means)'); title(ax,sprintf('q=%+.4f, common-beam normalised',sg*qa),'FontSize',9);
  if qa==qlist(1)&&sg<0, legend(ax,'FontSize',7,'Location','southoutside','NumColumns',2); end
 end
end
exportgraphics(fG,fullfile(out,'noise_coefficients.png'),'Resolution',110); close(fG);
exportgraphics(fT,fullfile(out,'frame_bands.png'),'Resolution',110); close(fT);
writetable(cell2table(gtab,'VariableNames',{'q','E_meV','g_raw','g_common_beam','g_per_channel','g_first_diff'}),fullfile(out,'noise_coefficients.csv'));
writetable(cell2table(corr_rows,'VariableNames',{'q','band','rel_std','poisson_rel_std','r_offset','r_zlp_centre','r_zlp_fwhm','r_core_same_channel','r_B1_same_channel'}),fullfile(out,'frame_band_correlations.csv'));
writetable(table((1:T).',ref,refB1,off,zc,zw,'VariableNames',{'frame','common_ref','B1_ref_300_1800','offset_px','zlp_centre_meV','zlp_fwhm_meV'}),fullfile(out,'frame_descriptors.csv'));
if prepared, save(fullfile(out,'noise_model.mat'),'E','q_noise','g','sigma','ho'); end
fprintf('diagnostics written to %s; common ref rel std %.4f, r(ref, B1 ref) %.3f\n',out,std(ref)/mean(ref),pcorr(ref,refB1));
end


function r = pcorr(a, b)
% Pearson correlation over the finite pairs (no Statistics Toolbox needed).
m = isfinite(a) & isfinite(b); a = a(m) - mean(a(m)); b = b(m) - mean(b(m));
r = sum(a .* b) / sqrt(sum(a .^ 2) * sum(b .^ 2));
end
