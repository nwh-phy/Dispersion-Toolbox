function out = b1_p2a_inputs(run_root)
% P2a step 1: input checks from the frame data of 590 (no B1 fits).
%  (1) Angular response of the elastic core per frame: A1-aligned frames,
%      window |E - c(q) - dE_t| <= w around the ZLP centre c(q) of the summed
%      map (w = 8 and 16 meV), channels |q| <= 0.0065 including q = 0;
%      centroid and width by moments (fixed central five channels, minus the
%      channel width^2/12) and by a channel-integrated Gaussian (central five
%      channels, relative-error weights); 6 blocks of 50 frames. On the summed
%      and block profiles (|q| <= 0.0065): one split Pearson VII (reference) and
%      two split Pearson VII components (narrow + broad, separate positions; the
%      profile model for W); a joint block fit with a common two-component shape
%      and per-block positions (proposed x-kernel of the blocked operator); the
%      per-frame beam position from the fixed-shape fit (primary dq_t for (2b),
%      (4) and (5)).
%  (2) Energy response per frame: ZLP centre per channel (log-parabola), frame
%      shift dE_t (count-weighted over |q| <= 0.001; mean of the wing channels
%      as an alternative that does not see the bending), FWHM and RMS width of
%      the central channel, centre and FWHM of the line summed over |q| <= 0.0065;
%      whether the bending c(q) moves with the beam (slope of c_t(q) on dq_t
%      against -c'(q)).
%  (3) sigma_probe: old definition (b1_adaptq_build: unaligned, |E| <= 12 meV,
%      five channels around the maximum, moments), the aligned / c(q)-centred
%      profiles (moments, Gaussian, split Pearson VII), block synthesis
%      sqrt(<s_b^2> + var(mu_b)).
%  (4) Frame-to-frame changes of band intensities (common-beam normalised),
%      |q| <= 0.0055: prediction -dq dB/dq + (dq^2/2 + d(s_x^2)/2) d2B/dq2
%      - dE dB/dE + d(s_E^2)/2 d2B/dE2 from the mean map, regression on the
%      same descriptors, Poisson ceiling.
%  (5) Block-to-block spectral shape of the key channels: leave-one-out
%      ratios, common-beam and shape-only chi^2 per energy region, before and
%      after the ratio predicted by the block centroid shifts (mean map moved
%      rigidly in q, interpolated across channels at each energy bin).
% Output: <run_root>/p2a/inputs_<session>_<config hash>/.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_p2a_inputs_config(); session=cfg.session; tag=b1_sha256(cfg); t_start=tic;
out=fullfile(run_root,'p2a',sprintf('inputs_%s_%s',session,tag(1:8))); if ~isfolder(out), mkdir(out); end
S=load(fullfile(run_root,'stage1',sprintf('channels_%s.mat',session)),'C'); C=S.C;
E=C.E(:); T=C.T; off=C.offsets(:); sup=C.support(:); dE=E(2)-E(1);
L=load(fullfile(run_root,'p1b',['diag_' session],'common_ref.mat'),'ref'); ref=L.ref(:); nrm=mean(ref)./ref;
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
l=load(fullfile(p2,session,'L1_minimal.mat')); qe=l.d.qe; q=qe.q_Ainv(:).'; dq=qe.dq_Ainv; Er=qe.energy_meV(:);
assert(abs(dq-C.dq)<1e-12,'b1_p2a_inputs:dq','dq of L1 and stage 1 differ');
raw=read_npy(C.raw_path);
src=1:size(raw,2); if isfield(qe,'source_channel'), src=qe.source_channel(:).'; end
assert(isequal(src,1:size(raw,2)),'b1_p2a_inputs:src','non-identity source channels are not handled');
invalid=false(1,size(raw,2));
for c=1:size(raw,2), x=raw(:,c,:); invalid(c)=any(~isfinite(x(:))|double(x(:))>=double(intmax('uint32'))); end
invalid=invalid|any(~isfinite(qe.intensity),1);
opts=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-14,'StepTolerance',1e-12, ...
 'OptimalityTolerance',1e-12,'MaxIterations',500,'MaxFunctionEvaluations',5000);

% ---- (3i) old sigma_probe definition: unaligned frames, raw |E| <= 12 meV, argmax +/- 2 channels
el_old=squeeze(sum(sum(double(raw(:,:,abs(Er)<=cfg.old_window)),1),3)); el_old(invalid)=0;
[~,i0]=max(el_old); kk=-2:2; w=el_old(i0+kk); w=w/sum(w); m1=sum(w.*kk);
sig_old=sqrt(max(sum(w.*(kk-m1).^2)-1/12,0))*dq; mu_old=q(i0)+m1*dq;
fprintf('old sigma_probe %.4e (stage 1 value %.4e), centroid %+.3e\n',sig_old,C.sigma_probe,mu_old);

% ---- aligned frames of the channels |q| <= q_max (incl. q = 0): E x channel x frame
sel=find(abs(q)<=cfg.q_max+1e-9 & ~invalid); [~,o]=sort(q(sel)); sel=sel(o); qs=q(sel); nk=numel(sel);
k0=find(abs(qs)<1e-12); prof_old=el_old(sel);
A=zeros(numel(E),nk,T);
for t=1:T, fr=double(squeeze(raw(t,sel,:))); A(:,:,t)=fr(:,sup+off(t)).'; end
clear raw
prof_al12=squeeze(sum(sum(A(abs(E)<=cfg.old_window,:,:),1),3));

% ---- (2) energy centres: summed map c(q), per frame c_t(q), frame shifts
Isum=sum(A,3); cq=zeros(1,nk); fwq=cq;
for k=1:nk, [cq(k),fwq(k)]=local_centre(Isum(:,k),E); end
ct=zeros(T,nk); fwt=ct;
for t=1:T, for k=1:nk, [ct(t,k),fwt(t,k)]=local_centre(A(:,k,t),E); end, end
kc=find(abs(qs)<=cfg.centre_q_max+1e-9); wc=sum(Isum(abs(E)<=16,kc),1); wc=wc(:)/sum(wc);
dE_core=(ct(:,kc)-cq(kc))*wc;
kwg=find(abs(qs)>=cfg.wing_centre_q(1)-1e-9 & abs(qs)<=cfg.wing_centre_q(2)+1e-9);
dE_wing=mean(ct(:,kwg)-cq(kwg),2); dE_0=ct(:,k0)-cq(k0);
rmsE=zeros(T,1);
for t=1:T
 y=A(:,k0,t); wv=local_window(E,dE,ct(t,k0),16); rmsE(t)=sqrt(sum((E-ct(t,k0)).^2.*y.*wv)/sum(y.*wv));
end
call=zeros(T,1); fwall=call;
for t=1:T, [call(t),fwall(t)]=local_centre(sum(A(:,:,t),2),E); end
[c_all,fw_all]=local_centre(sum(Isum,2),E); dE_all=call-c_all;
P1=readtable(fullfile(run_root,'p1b',['diag_' session],'frame_descriptors.csv'));
fprintf('check vs P1b frame descriptors: max |centre diff| %.2e meV, max |FWHM diff| %.2e meV\n', ...
 max(abs(ct(:,k0)-P1.zlp_centre_meV)),max(abs(fwt(:,k0)-P1.zlp_fwhm_meV)));

% ---- (1) core-window angular profiles per frame (window follows c(q) + dE_t; 'fix' = c(q) only)
wl=cfg.window_list; nwl=numel(wl); P=zeros(T,nk,nwl); Pfix=P;
for t=1:T
 At=A(:,:,t);
 for j=1:nwl
  P(t,:,j)=sum(At.*local_window(E,dE,cq+dE_core(t),wl(j)),1);
  Pfix(t,:,j)=sum(At.*local_window(E,dE,cq,wl(j)),1);
 end
end
km=find(abs(qs)<=cfg.moment_q_max+1e-9); kg=find(abs(qs)<=cfg.gauss_q_max+1e-9);
jw=arrayfun(@(v) find(wl==v),cfg.core_w); nw=numel(jw);
mom_mu=zeros(T,nw); mom_s=mom_mu; g_mu=mom_mu; g_s=mom_mu; g_rms=mom_mu; momfix_mu=mom_mu; momfix_s=mom_mu;
for j=1:nw
 [mom_mu(:,j),mom_s(:,j)]=local_moments(P(:,km,jw(j)),qs(km),dq);
 [momfix_mu(:,j),momfix_s(:,j)]=local_moments(Pfix(:,km,jw(j)),qs(km),dq);
 for t=1:T
  [g_mu(t,j),g_s(t,j),r]=local_gauss(P(t,kg,jw(j)),qs(kg),dq,mom_mu(t,j),mom_s(t,j),opts); g_rms(t,j)=sqrt(mean(r.^2));
 end
end
nb=cfg.n_blocks; blk=ceil((1:T).'/(T/nb));

% block and total profiles: moments, Gaussian (central five), one split Pearson VII
% (reference) and the profile model for W: two split Pearson VII components
% (narrow core + broad, separate positions), |q| <= q_max, log residuals; block
% test with the two-component shape of the all-frame profile fixed.
bt=struct([]); qmr=max(abs(qs))+dq/2;
for j=1:nw
 for b=0:nb
  if b==0, m=true(T,1); else, m=blk==b; end
  Pb=sum(P(m,:,jw(j)),1); [mb,sb]=local_moments(Pb(km),qs(km),dq);
  [gm,gs,rg]=local_gauss(Pb(kg),qs(kg),dq,mb,sb,opts);
  [ps,~,rs]=local_fit_spv(Pb,qs,dq,gm,opts,[]);
  if b==0
   [pd,~,rd]=local_fit_dspv(Pb,qs,dq,local_dspv_start(Pb,gm,dq),opts,[]); dtot=pd; pf=pd; rf=rd; Pmax0=max(Pb); g0=gm;
  else
   [pd,~,rd]=local_fit_dspv(Pb,qs,dq,dtot,opts,[]);
   [pf,~,rf]=local_fit_dspv(Pb,qs,dq,[dtot(1)+log(max(Pb)/Pmax0) dtot(2)+(gm-g0)/dq],opts,dtot(3:10));
  end
  wt=sum(P(m,kg,jw(j)),2); wt=wt/sum(wt);
  bt(end+1).w=cfg.core_w(j); bt(end).block=b; bt(end).profile=Pb; %#ok<AGROW>
  bt(end).mom=[mb sb]; bt(end).gauss=[gm gs]; bt(end).spv=ps; bt(end).dspv=pd; bt(end).dspv_fixed=pf;
  bt(end).dspv_stats=local_dspv_stats(pd,dq,qmr);
  bt(end).logrms=[sqrt(mean(rg.^2)) sqrt(mean(rs.^2)) sqrt(mean(rd.^2)) sqrt(mean(rf.^2))];
  bt(end).tail=sum(Pb(abs(qs)>=cfg.tail_q_min-1e-9))/sum(Pb);
  bt(end).within_std_g=sqrt(sum(wt.*(g_mu(m,j)-sum(wt.*g_mu(m,j))).^2));
 end
end

% per-frame beam position: two-component shape of the all-frame profile fixed,
% amplitude and position free (position = mode of the narrow component)
pf_mu=zeros(T,nw);
for j=1:nw
 tot=bt([bt.w]==cfg.core_w(j) & [bt.block]==0); d0=tot.dspv;
 for t=1:T
  pp=local_fit_dspv(P(t,:,jw(j)),qs,dq,[d0(1)+log(max(P(t,:,jw(j)))/max(tot.profile)) d0(2)+(g_mu(t,j)-tot.gauss(1))/dq],opts,d0(3:10));
  pf_mu(t,j)=pp(2)*dq;
 end
end

% joint fit of the six block profiles: common two-component shape, amplitude and
% position per block. This is the x-kernel proposed for the blocked operator
% (within-block jitter inside the shape, between-block drift only as the shifts).
jb=struct([]); wrow={}; [lbd,ubd]=local_dspv_bounds(); n_b=accumarray(blk,ref,[nb 1]).';
for j=1:nw
 Bj=bt([bt.w]==cfg.core_w(j) & [bt.block]>0); Pj=cat(1,Bj.profile); tot=bt([bt.w]==cfg.core_w(j) & [bt.block]==0);
 u0=tot.dspv(3:10); for b=1:nb, u0=[u0 Bj(b).dspv_fixed(1:2)]; end %#ok<AGROW>
 u=lsqnonlin(@(u) local_joint_res(u,Pj,qs,dq,nb),u0,[lbd(3:10) repmat([-Inf -3],1,nb)],[ubd(3:10) repmat([Inf 3],1,nb)],opts);
 r=local_joint_res(u,Pj,qs,dq,nb);
 jb(j).w=cfg.core_w(j); jb(j).shape=u(1:8); jb(j).logA=u(9:2:end); jb(j).pos=u(10:2:end)*dq; jb(j).logrms=sqrt(mean(r.^2,1));
 jb(j).pos_mean=sum(n_b.*jb(j).pos)/sum(n_b);
 wrow(end+1,:)=[{cfg.core_w(j),'all frames (time-integrated operator)',0} num2cell(local_dspv_row(tot.dspv,dq,qmr)) {tot.logrms(3)}]; %#ok<AGROW>
 for b=1:nb
  wrow(end+1,:)=[{cfg.core_w(j),'blocks, common shape',b} num2cell(local_dspv_row([jb(j).logA(b) jb(j).pos(b)/dq jb(j).shape],dq,qmr)) {jb(j).logrms(b)}]; %#ok<AGROW>
 end
end
writetable(cell2table(wrow,'VariableNames',{'w_meV','version','block','narrow_mode_Ainv','mean_position_Ainv', ...
 'narrow_gamma_minus_Ainv','narrow_gamma_plus_Ainv','narrow_2m','narrow_fwhm_Ainv','broad_offset_Ainv','broad_gamma_minus_Ainv', ...
 'broad_gamma_plus_Ainv','broad_2m','broad_fraction_total','fraction_outside_qmax','logrms'}),fullfile(out,'w_kernel_x.csv'));

% ---- (3) sigma_probe comparison on summed profiles
prof=struct('name',{},'w',{},'P',{});
prof(end+1)=struct('name',"old: unaligned, raw |E|<=12",'w',12,'P',prof_old);
prof(end+1)=struct('name',"aligned, |E|<=12 (no c(q))",'w',12,'P',prof_al12);
for j=1:nwl
 prof(end+1)=struct('name',"aligned, |E-c(q)-dE_t|<=w",'w',wl(j),'P',sum(P(:,:,j),1)); %#ok<AGROW>
 prof(end+1)=struct('name',"aligned, |E-c(q)|<=w (no dE_t)",'w',wl(j),'P',sum(Pfix(:,:,j),1)); %#ok<AGROW>
end
srow={};
for i=1:numel(prof)
 Pv=prof(i).P(:).'; [ma,sa]=local_moments_argmax(Pv,qs,dq); [mf,sf]=local_moments(Pv(km),qs(km),dq);
 [gm,gs]=local_gauss(Pv(kg),qs(kg),dq,mf,sf,opts); [ps,~,rs]=local_fit_spv(Pv,qs,dq,gm,opts,[]);
 [pd,~,rd]=local_fit_dspv(Pv,qs,dq,local_dspv_start(Pv,gm,dq),opts,[]); st=local_dspv_stats(pd,dq,qmr);
 srow(end+1,:)={prof(i).name,prof(i).w,'moments, argmax +/-2 (old)',ma,sa,NaN,NaN,NaN}; %#ok<AGROW>
 srow(end+1,:)={prof(i).name,prof(i).w,'moments, fixed central 5',mf,sf,NaN,NaN,NaN}; %#ok<AGROW>
 srow(end+1,:)={prof(i).name,prof(i).w,'Gaussian, central 5',gm,gs,NaN,NaN,NaN}; %#ok<AGROW>
 srow(end+1,:)={prof(i).name,prof(i).w,'one split Pearson VII, |q|<=0.0065: FWHM/2.355; tail power',ps(2)*dq,local_spv_fwhm(ps,dq)/(2*sqrt(2*log(2))),NaN,2*exp(ps(5)),sqrt(mean(rs.^2))}; %#ok<AGROW>
 srow(end+1,:)={prof(i).name,prof(i).w,'two split Pearson VII, |q|<=0.0065: narrow FWHM/2.355; broad fraction and power',pd(2)*dq,st(1)/(2*sqrt(2*log(2))),st(2),2*exp(pd(10)),sqrt(mean(rd.^2))}; %#ok<AGROW>
end
for j=1:nw
 B0=bt([bt.w]==cfg.core_w(j)); Bb=B0([B0.block]>0); wb=arrayfun(@(x) sum(x.profile(kg)),Bb); wb=wb/sum(wb);
 for est=1:2
  if est==1, mu_b=arrayfun(@(x) x.mom(1),Bb); s_b=arrayfun(@(x) x.mom(2),Bb); nm='moments'; mt=mom_mu(:,j); st=mom_s(:,j);
  else, mu_b=arrayfun(@(x) x.gauss(1),Bb); s_b=arrayfun(@(x) x.gauss(2),Bb); nm='Gaussian'; mt=g_mu(:,j); st=g_s(:,j); end
  mub=sum(wb.*mu_b); syn=sqrt(sum(wb.*s_b.^2)+sum(wb.*(mu_b-mub).^2));
  wt=sum(P(:,kg,jw(j)),2); wt=wt/sum(wt); mtt=sum(wt.*mt); synf=sqrt(sum(wt.*st.^2)+sum(wt.*(mt-mtt).^2));
  srow(end+1,:)={"synthesis: blocks sqrt(<s_b^2>+var mu_b)",cfg.core_w(j),nm,mub,syn,NaN,NaN,NaN}; %#ok<AGROW>
  srow(end+1,:)={"blocks: mean within-block width sqrt(<s_b^2>)",cfg.core_w(j),nm,mub,sqrt(sum(wb.*s_b.^2)),NaN,NaN,NaN}; %#ok<AGROW>
  srow(end+1,:)={"blocks: between-block centroid std",cfg.core_w(j),nm,mub,sqrt(sum(wb.*(mu_b-mub).^2)),NaN,NaN,NaN}; %#ok<AGROW>
  srow(end+1,:)={"synthesis: frames sqrt(<s_t^2>+var mu_t)",cfg.core_w(j),nm,mtt,synf,NaN,NaN,NaN}; %#ok<AGROW>
  srow(end+1,:)={"frames: centroid std (all 300)",cfg.core_w(j),nm,mtt,sqrt(sum(wt.*(mt-mtt).^2)),NaN,NaN,NaN}; %#ok<AGROW>
 end
 pj=jb(j).pos; srow(end+1,:)={"blocks: two-Pearson joint fit, between-block position std",cfg.core_w(j),'two split Pearson VII',jb(j).pos_mean, ...
  sqrt(sum(n_b.*(pj-jb(j).pos_mean).^2)/sum(n_b)),NaN,NaN,NaN}; %#ok<AGROW>
end
ST=cell2table(srow,'VariableNames',{'profile','w_meV','estimator','position_Ainv','sigma_Ainv','broad_fraction_total','tail_power_2m','logrms_27ch'});
ST.sigma_channels=ST.sigma_Ainv/dq; ST.position_channels=ST.position_Ainv/dq;
writetable(ST,fullfile(out,'sigma_probe_comparison.csv'));

% ---- (2b) does the bending move with the beam? c_t(q) = c(q - dq_t) + eps_t predicts slope -c'(q)
dqv=pf_mu(:,1)-mean(pf_mu(:,1)); cpr=gradient(cq,dq); eps_t=dE_wing-mean(dE_wing);
kt=find(abs(qs)<=cfg.tie_q_max+1e-9); tie=zeros(numel(kt),8);
for i=1:numel(kt)
 k=kt(i); X=[ones(T,1) dqv eps_t]; bb=X\ct(:,k); r=ct(:,k)-X*bb; R2=1-sum(r.^2)/sum((ct(:,k)-mean(ct(:,k))).^2);
 tie(i,:)=[qs(k) cq(k) cpr(k)*dq -cpr(k)*dq bb(2)*dq bb(3) R2 pcorr(ct(:,k),dqv)];
end
TT=array2table(tie,'VariableNames',{'q','c_q_meV','dc_dq_meV_per_channel','beam_tied_prediction_meV_per_channel', ...
 'measured_slope_on_dq_meV_per_channel','slope_on_wing_shift','R2','r_c_t_vs_dq'});
writetable(TT,fullfile(out,'bending_vs_beam_position.csv'));
cb=zeros(nb,nk); fwb=cb;
for b=1:nb, Ib=sum(A(:,:,blk==b),3); for k=1:nk, [cb(b,k),fwb(b,k)]=local_centre(Ib(:,k),E); end, end
writetable(array2table([qs(:) cq(:) fwq(:) cb.' fwb.'],'VariableNames',[{'q','c_all_meV','fwhm_all_meV'}, ...
 arrayfun(@(b) sprintf('c_block%d',b),1:nb,'UniformOutput',false),arrayfun(@(b) sprintf('fwhm_block%d',b),1:nb,'UniformOutput',false)]), ...
 fullfile(out,'block_channel_centres.csv'));

% ---- frame and block descriptor tables
FT=table((1:T).',blk,ref,off,dE_core,dE_0,dE_wing,dE_all,fwt(:,k0),rmsE,fwall,'VariableNames', ...
 {'frame','block','common_ref','offset_px','dE_core_meV','dE_central_meV','dE_wing_meV','dE_allq_meV','fwhm_central_meV','rmsE_central_16_meV','fwhm_allq_meV'});
for j=1:nw
 s=sprintf('_w%d',cfg.core_w(j));
 FT.(['mu_mom' s])=mom_mu(:,j); FT.(['s_mom' s])=mom_s(:,j); FT.(['mu_gauss' s])=g_mu(:,j); FT.(['s_gauss' s])=g_s(:,j);
 FT.(['gauss_logrms' s])=g_rms(:,j); FT.(['pos_two_pearson_fixed' s])=pf_mu(:,j); FT.(['mu_mom_nodE' s])=momfix_mu(:,j); FT.(['core_counts_central5' s])=sum(P(:,km,jw(j)),2);
end
writetable(FT,fullfile(out,'frame_descriptors_p2a.csv'));
brow={};
for i=1:numel(bt)
 x=bt(i); if x.block==0, m=true(T,1); else, m=blk==x.block; end
 brow(end+1,:)={x.w,x.block,find(m,1),find(m,1,'last'),x.mom(1),x.mom(2),x.gauss(1),x.gauss(2),x.within_std_g, ...
  x.spv(2)*dq,2*exp(x.spv(5)),x.dspv(2)*dq,x.dspv_stats(4),x.dspv_stats(1),x.dspv_stats(2),x.dspv(4)*dq,2*exp(x.dspv(10)), ...
  x.dspv_fixed(2)*dq,x.tail,x.logrms(1),x.logrms(2),x.logrms(3),x.logrms(4), ...
  mean(dE_core(m)),mean(dE_wing(m)),mean(dE_all(m)),mean(fwt(m,k0)),mean(fwall(m))}; %#ok<AGROW>
end
BT=cell2table(brow,'VariableNames',{'w_meV','block','first_frame','last_frame','mu_mom','s_mom','mu_gauss','s_gauss', ...
 'within_block_std_mu_gauss','mu_one_pearson','one_pearson_2m','narrow_mode_two_pearson','mean_position_two_pearson', ...
 'narrow_fwhm_two_pearson','broad_fraction_two_pearson','broad_offset_two_pearson','broad_2m_two_pearson', ...
 'narrow_mode_shape_fixed','frac_counts_absq_ge_tail_q_min','logrms_gauss5','logrms_one_pearson','logrms_two_pearson', ...
 'logrms_two_pearson_shape_fixed','dE_core_mean','dE_wing_mean','dE_allq_mean','fwhm_central_mean','fwhm_allq_mean'});
writetable(BT,fullfile(out,'block_descriptors.csv'));
PT=table(qs(:),'VariableNames',{'q'});
for i=1:numel(bt), PT.(sprintf('w%d_block%d',bt(i).w,bt(i).block))=bt(i).profile(:); end
PT.old_unaligned_E12=prof_old(:); PT.aligned_E12=prof_al12(:);
writetable(PT,fullfile(out,'angular_profiles.csv'));

% ---- (4) attribution of frame-to-frame band changes
ka=find(abs(qs)<=cfg.attrib_q_max+1e-9); bands=cfg.bands; nbd=size(bands,1);
B=zeros(T,nk,nbd);
for b=1:nbd, m=E>=bands(b,1)&E<=bands(b,2); B(:,:,b)=squeeze(sum(A(m,:,:),1)).'; end
Bn=B.*nrm; Bbar=squeeze(mean(Bn,1));
Ibar=mean(A.*reshape(nrm,1,1,T),3); [~,dIdE]=gradient(Ibar,1,dE); [~,d2IdE2]=gradient(dIdE,1,dE);
dq2=dqv.^2-mean(dqv.^2); dsx2=g_s(:,1).^2-mean(g_s(:,1).^2); dEv=dE_all-mean(dE_all);
sE2=(fwall/(2*sqrt(2*log(2)))).^2; dsE2=sE2-mean(sE2);
dqg=g_mu(:,1)-mean(g_mu(:,1)); dqg2=dqg.^2-mean(dqg.^2);
j16=jw(cfg.core_w==16); arow={};
for k=ka(:).'
 for b=1:nbd
  dB=Bn(:,k,b)-Bbar(k,b); V=sum(dB.^2);
  d1=(Bbar(k+1,b)-Bbar(k-1,b))/(2*dq); d2=(Bbar(k+1,b)-2*Bbar(k,b)+Bbar(k-1,b))/dq^2;
  m=E>=bands(b,1)&E<=bands(b,2); e1=sum(dIdE(m,k)); e2=sum(d2IdE2(m,k));
  % prediction from the mean map: beam shift (two orders) and energy terms; the
  % apparent Gaussian width is a sub-pixel artefact and stays only in the regression
  predq=-dqv*d1+0.5*dq2*d2; pred=predq-dEv*e1+0.5*dsE2*e2; predqg=-dqg*d1+0.5*dqg2*d2;
  predg=predqg-dEv*e1+0.5*dsE2*e2;
  X=[dqv dq2 dsx2 dEv dsE2]; beta=X\dB; Xg=[dqg dqg2 dsx2 dEv dsE2]; betag=Xg\dB;
  vp=mean(B(:,k,b).*nrm.^2); vt=V/T;
  core=P(:,k,j16).*nrm; rc=pcorr(dB,core);
  arow(end+1,:)={qs(k),sprintf('%d-%d',bands(b,1),bands(b,2)),Bbar(k,b),sqrt(vt)/Bbar(k,b),sqrt(vp)/Bbar(k,b), ...
   1-vp/vt,1-sum((dB-predq).^2)/V,1-sum((dB-pred).^2)/V,1-sum((dB-X*beta).^2)/V,rc^2,sum(dB.*predq)/sum(predq.^2), ...
   1-sum((dB-predg).^2)/V,1-sum((dB-Xg*betag).^2)/V,sum(dB.*predqg)/sum(predqg.^2)}; %#ok<AGROW>
 end
end
AT=cell2table(arow,'VariableNames',{'q','band_meV','mean_counts','rel_std','poisson_rel_std','R2_ceiling', ...
 'R2_pred_dq_only','R2_pred_shift_energy','R2_regression','R2_core_proxy','amp_ratio_meas_over_shift_pred', ...
 'R2_pred_shift_energy_gauss','R2_regression_gauss_centroid','amp_ratio_gauss_centroid'});
writetable(AT,fullfile(out,'deltaI_explained.csv'));

% ---- (5) block-to-block shape of the key channels
edges=cfg.shape_edges; ec=(edges(1:end-1)+edges(2:end))/2; nbin=numel(ec); reg=cfg.shape_regions;
ib=discretize(E,edges); kq=arrayfun(@(v) find(abs(qs-v)<1e-9),[-fliplr(cfg.key_q) cfg.key_q]);
okb=~isnan(ib); Mb=zeros(nbin,nk); for k=1:nk, Mb(:,k)=accumarray(ib(okb),Ibar(okb,k),[nbin 1]); end
mu_blk=jb(1).pos; mu_all=jb(1).pos_mean;
n_blk=accumarray(blk,ref,[nb 1]).';
% leave-one-block-out templates: mean map of the other blocks (common-beam normalised), binned
An=A.*reshape(nrm,1,1,T); Sall=sum(An,3); Mo=zeros(nbin,nk,nb);
for b=1:nb
 Io=(Sall-sum(An(:,:,blk==b),3))/nnz(blk~=b);
 for k=1:nk, Mo(:,k,b)=accumarray(ib(okb),Io(okb,k),[nbin 1]); end
end
clear An
shp=struct('q',{},'r',{},'sr',{},'rp',{}); crow={};
for k=kq
 R=zeros(nb,nbin); SR=R; RP=R;
 for b=1:nb
  sb=accumarray(ib(okb),sum(A(okb,k,blk==b),3),[nbin 1]).';
  so=accumarray(ib(okb),sum(A(okb,k,blk~=b),3),[nbin 1]).';
  nb_=sum(ref(blk==b)); no_=sum(ref(blk~=b)); r=(sb/nb_)./(so/no_); sr=r.*sqrt(1./sb+1./so); R(b,:)=r; SR(b,:)=sr;
  % shift prediction from the other blocks only: their template moved by mu_b - mu_o
  oth=setdiff(1:nb,b); mu_o=sum(n_blk(oth).*mu_blk(oth))/sum(n_blk(oth));
  lb_=arrayfun(@(i) pchip(qs,log(Mo(i,:,b)),qs(k)-(mu_blk(b)-mu_o)),1:nbin);
  rp=exp(lb_-log(Mo(:,k,b)).'); RP(b,:)=rp;
  for g=1:size(reg,1)
   mm=ec>=reg(g,1)&ec<=reg(g,2); a=sum(sb(mm))/sum(so(mm)*nb_/no_); rs=r(mm)/a; srs=sr(mm)/a;
   ap=sum(rp(mm).*so(mm))/sum(so(mm)); rps=rp(mm)/ap;
   crow(end+1,:)={qs(k),b,sprintf('%d-%d',reg(g,1),reg(g,2)),nnz(mm),a,ap,sum(((r(mm)-1)./sr(mm)).^2)/nnz(mm), ...
    sum(((r(mm)-rp(mm))./sr(mm)).^2)/nnz(mm),sum(((rs-1)./srs).^2)/(nnz(mm)-1),sum(((rs-rps)./srs).^2)/(nnz(mm)-1), ...
    sqrt(mean((rs-1).^2)),sqrt(mean((rs-rps).^2)),sqrt(mean(srs.^2))}; %#ok<AGROW>
  end
 end
 shp(end+1)=struct('q',qs(k),'r',R,'sr',SR,'rp',RP); %#ok<AGROW>
end
CT=cell2table(crow,'VariableNames',{'q','block','region_meV','n_bins','amplitude_vs_others','amplitude_predicted', ...
 'chi2dof_common_beam','chi2dof_common_beam_after_shift','chi2dof_shape_only','chi2dof_shape_only_after_shift', ...
 'rms_shape_dev','rms_shape_dev_after_shift','rms_poisson_dev'});
writetable(CT,fullfile(out,'block_shape_consistency.csv'));

% ---- figures
blkmean=@(x) accumarray(blk,x,[nb 1],@mean);
f=figure('Visible','off','Position',[0 0 1700 1050]); tl=tiledlayout(f,2,2,'TileSpacing','compact');
ax=nexttile(tl); hold(ax,'on'); cl=lines(4); lab={};
for j=1:nw
 plot(ax,1:T,g_mu(:,j)/dq,'-','Color',cl(j,:)); plot(ax,1:T,mom_mu(:,j)/dq,':','Color',cl(j,:));
 lab=[lab {sprintf('Gaussian, w=%d',cfg.core_w(j)),sprintf('moments, w=%d',cfg.core_w(j))}]; %#ok<AGROW>
end
plot(ax,1:T,pf_mu(:,1)/dq,'-','Color',[0.2 0.6 0.2]); lab=[lab {'two-Pearson position (fixed shape), w=8'}];
bm=blkmean(pf_mu(:,1)); stairs(ax,0:T/nb:T,[bm; bm(end)]/dq,'k-','LineWidth',1.5); lab=[lab {'block mean (two-Pearson, w=8)'}];
grid(ax,'on'); xlabel(ax,'frame'); ylabel(ax,'angular centroid (channels, 1 = 0.0005 A^{-1})'); legend(ax,lab,'FontSize',7,'Location','best');
title(ax,'(a) elastic-core angular centroid \mu_t');
ax=nexttile(tl); hold(ax,'on');
for j=1:nw, plot(ax,1:T,g_s(:,j)/dq,'-','Color',cl(j,:)); plot(ax,1:T,mom_s(:,j)/dq,':','Color',cl(j,:)); end
bm=blkmean(g_s(:,1)); stairs(ax,0:T/nb:T,[bm; bm(end)]/dq,'k-','LineWidth',1.5);
grid(ax,'on'); xlabel(ax,'frame'); ylabel(ax,'angular width s_x (channels)'); title(ax,'(b) elastic-core angular width (channel width removed)');
ax=nexttile(tl); hold(ax,'on');
plot(ax,1:T,dE_core,'-',1:T,dE_0,':',1:T,dE_wing,'-');
bm=blkmean(dE_core); stairs(ax,0:T/nb:T,[bm; bm(end)],'k-','LineWidth',1.5);
grid(ax,'on'); xlabel(ax,'frame'); ylabel(ax,'energy shift after A1 alignment (meV)');
legend(ax,{'dE_t: core channels |q|<=0.001 (count weighted)','central channel','wing channels 0.002<=|q|<=0.004','block mean'},'FontSize',7,'Location','best');
title(ax,'(c) residual energy shift per frame');
ax=nexttile(tl); yyaxis(ax,'left'); plot(ax,1:T,fwt(:,k0),'-',1:T,fwall,'--'); ylabel(ax,'FWHM (meV): central channel (solid), |q|<=0.0065 summed (dashed)');
yyaxis(ax,'right'); plot(ax,1:T,ref/mean(ref),'-'); ylabel(ax,'common-beam reference / mean'); grid(ax,'on'); xlabel(ax,'frame');
title(ax,sprintf('(d) energy width, beam current; r(\\mu_t, FWHM): central %.2f, summed %.2f',pcorr(fwt(:,k0),g_mu(:,1)),pcorr(fwall,g_mu(:,1))));
exportgraphics(f,fullfile(out,'fig1_time_series.png'),'Resolution',110); close(f);

f=figure('Visible','off','Position',[0 0 1700 1050]); tl=tiledlayout(f,2,2,'TileSpacing','compact'); cb6=turbo(nb+2);
for j=1:nw
 B0=bt([bt.w]==cfg.core_w(j)); tot=B0([B0.block]==0);
 ax=nexttile(tl); hold(ax,'on');
 for b=1:nb, x=B0([B0.block]==b); semilogy(ax,qs,x.profile/sum(x.profile),'o-','Color',cb6(b+1,:),'MarkerSize',3,'DisplayName',sprintf('block %d',b)); end
 d=tot.dspv; S0=sum(tot.profile);
 semilogy(ax,qs,tot.profile/S0,'k.','MarkerSize',14,'DisplayName','all frames');
 semilogy(ax,qs,local_dspv(qs(:),dq,d)/S0,'k-','DisplayName',sprintf('two split Pearson VII (log rms %.3f)',tot.logrms(3)));
 semilogy(ax,qs,local_spv(qs(:),dq,[d(1) d(2) d(5:7)])/S0,'b--','DisplayName','narrow component');
 semilogy(ax,qs,local_spv(qs(:),dq,[d(1)+d(3) d(2)+d(4) d(8:10)])/S0,'b:','DisplayName',sprintf('broad component (%.3f of total)',tot.dspv_stats(2)));
 semilogy(ax,qs,local_spv(qs(:),dq,tot.spv)/S0,'-','Color',[0.6 0.6 0.6],'DisplayName',sprintf('one split Pearson (log rms %.3f)',tot.logrms(2)));
 set(ax,'YScale','log'); ylim(ax,[1e-5 1]); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'core counts / sum');
 legend(ax,'FontSize',7,'Location','southoutside','NumColumns',5);
 title(ax,sprintf('(%s) elastic-core angular profile, w = %d meV; counts at |q| >= %.4f: %.3f',char('a'+2*(j-1)),cfg.core_w(j),cfg.tail_q_min,tot.tail));
 ax=nexttile(tl); hold(ax,'on');
 for b=1:nb, x=B0([B0.block]==b); plot(ax,qs,(x.profile/sum(x.profile))./(tot.profile/sum(tot.profile))-1,'o-','Color',cb6(b+1,:),'MarkerSize',3); end
 grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'block / all - 1 (both normalised to sum)'); ylim(ax,[-0.6 0.6]);
 title(ax,sprintf('(%s) block profile differences, w = %d meV',char('b'+2*(j-1)),cfg.core_w(j)));
end
exportgraphics(f,fullfile(out,'fig2_block_profiles.png'),'Resolution',110); close(f);

f=figure('Visible','off','Position',[0 0 1700 650]); tl=tiledlayout(f,1,2,'TileSpacing','compact');
ax=nexttile(tl); hold(ax,'on'); sh=[1 2 3 5 7 9]; cl6=lines(numel(sh));
for i=1:numel(sh), Pv=prof(sh(i)).P; semilogy(ax,qs,Pv/sum(Pv),'o-','Color',cl6(i,:),'MarkerSize',3,'DisplayName',sprintf('%s, w=%d',prof(sh(i)).name,prof(sh(i)).w)); end
set(ax,'YScale','log'); xlim(ax,[-0.0035 0.0035]); ylim(ax,[1e-4 1]); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'normalised to sum');
legend(ax,'FontSize',7,'Location','southoutside'); title(ax,'(a) summed elastic-core profiles used for \sigma_{probe}');
ax=nexttile(tl); hold(ax,'on'); ref3=prof(3).P/sum(prof(3).P);
for i=1:numel(sh), Pv=prof(sh(i)).P; plot(ax,qs,(Pv/sum(Pv))./ref3-1,'o-','Color',cl6(i,:),'MarkerSize',3); end
xlim(ax,[-0.0035 0.0035]); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'relative to aligned, c(q)+dE_t, w = 8');
title(ax,'(b) profile differences between definitions');
exportgraphics(f,fullfile(out,'fig3_sigma_probe_profiles.png'),'Resolution',110); close(f);

f=figure('Visible','off','Position',[0 0 1700 1050]); tl=tiledlayout(f,2,2,'TileSpacing','compact');
ax=nexttile(tl); hold(ax,'on'); for b=1:nb, plot(ax,qs,cb(b,:),'o-','Color',cb6(b+1,:),'MarkerSize',3); end
plot(ax,qs,cq,'k.-','MarkerSize',12); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'ZLP centre (meV)');
legend(ax,[arrayfun(@(b) sprintf('block %d (position %+.2f ch)',b,jb(1).pos(b)/dq),1:nb,'UniformOutput',false) {'all frames'}],'FontSize',7,'Location','best');
title(ax,'(a) energy-centre bending c(q) per block');
ax=nexttile(tl); hold(ax,'on'); plot(ax,tie(:,1),tie(:,4),'ks-','DisplayName','beam-tied prediction -c''(q)'); plot(ax,tie(:,1),tie(:,5),'ro-','DisplayName','measured slope of c_t(q) on \mu_t');
grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'meV per channel of beam shift'); legend(ax,'FontSize',8,'Location','best');
title(ax,'(b) does the bending move with the beam?');
ax=nexttile(tl); hold(ax,'on'); kpm=[find(abs(qs+dq)<1e-9) find(abs(qs-dq)<1e-9)];
for i=1:2, scatter(ax,dqv/dq,ct(:,kpm(i))-eps_t,8,'filled','DisplayName',sprintf('q = %+.4f',qs(kpm(i)))); end
grid(ax,'on'); xlabel(ax,'\mu_t - <\mu> (channels)'); ylabel(ax,'c_t(q) - wing shift (meV)'); legend(ax,'FontSize',8,'Location','best');
title(ax,'(c) per-frame energy centre against beam position');
ax=nexttile(tl); hold(ax,'on'); for b=1:nb, plot(ax,qs,fwb(b,:),'o-','Color',cb6(b+1,:),'MarkerSize',3); end
plot(ax,qs,fwq,'k.-','MarkerSize',12); yline(ax,fw_all,'k--',sprintf('summed |q|<=0.0065: %.1f meV',fw_all)); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'FWHM (meV)');
title(ax,sprintf('(d) ZLP FWHM per channel and block; summed line: r(centre, \\mu_t) = %.2f, r(FWHM, \\mu_t) = %.2f',pcorr(call,g_mu(:,1)),pcorr(fwall,g_mu(:,1))));
exportgraphics(f,fullfile(out,'fig4_energy_response.png'),'Resolution',110); close(f);

f=figure('Visible','off','Position',[0 0 1800 900]); tl=tiledlayout(f,1,nbd,'TileSpacing','compact');
for b=1:nbd
 ax=nexttile(tl); hold(ax,'on'); m=strcmp(AT.band_meV,sprintf('%d-%d',bands(b,1),bands(b,2)));
 plot(ax,AT.q(m),AT.R2_ceiling(m),'k-','LineWidth',1.5); plot(ax,AT.q(m),AT.R2_regression(m),'bo-','MarkerSize',3);
 plot(ax,AT.q(m),AT.R2_pred_shift_energy(m),'rs-','MarkerSize',3); plot(ax,AT.q(m),AT.R2_pred_dq_only(m),'m^:','MarkerSize',3);
 plot(ax,AT.q(m),AT.R2_core_proxy(m),'g--');
 ylim(ax,[-0.5 1]); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); title(ax,sprintf('%d-%d meV',bands(b,1),bands(b,2)));
 if b==1, ylabel(ax,'fraction of frame-to-frame variance'); legend(ax,{'ceiling (1 - Poisson/total)','regression on descriptors', ...
   'prediction from mean map (shift + energy terms)','prediction, beam shift only','same-channel core (r^2)'},'FontSize',7,'Location','southoutside'); end
end
exportgraphics(f,fullfile(out,'fig5_deltaI_explained.png'),'Resolution',110); close(f);

f=figure('Visible','off','Position',[0 0 1900 950]); tl=tiledlayout(f,2,4,'TileSpacing','compact');
for i=1:numel(shp)
 ax=nexttile(tl); hold(ax,'on'); x=shp(i);
 fill(ax,[ec fliplr(ec)],[1+x.sr(1,:) fliplr(1-x.sr(1,:))],[0.85 0.85 0.85],'EdgeColor','none');
 for b=1:nb, plot(ax,ec,x.r(b,:),'-','Color',cb6(b+1,:)); end
 yline(ax,1,'k-'); xlim(ax,[40 1800]); ylim(ax,[0.8 1.2]); grid(ax,'on'); set(ax,'XScale','log'); xticks(ax,[40 100 150 250 400 800 1800]);
 title(ax,sprintf('q = %+.4f',x.q)); if i>4, xlabel(ax,'E (meV)'); end
 if mod(i,4)==1, ylabel(ax,'block / other blocks (common beam)'); end
 if i==1, legend(ax,[{'\pm1\sigma Poisson (block 1)'} arrayfun(@(b) sprintf('block %d',b),1:nb,'UniformOutput',false)],'FontSize',7,'Location','north','NumColumns',2); end
end
exportgraphics(f,fullfile(out,'fig6_block_shapes.png'),'Resolution',110); close(f);
f=figure('Visible','off','Position',[0 0 1900 950]); tl=tiledlayout(f,2,4,'TileSpacing','compact');
for i=1:numel(shp)
 ax=nexttile(tl); hold(ax,'on'); x=shp(i);
 for b=1:nb, plot(ax,ec,x.r(b,:),'-','Color',[cb6(b+1,:) 0.5]); plot(ax,ec,x.rp(b,:),'--','Color',cb6(b+1,:),'LineWidth',1.8); end
 yline(ax,1,'k-'); xlim(ax,[40 1800]); ylim(ax,[0.75 1.25]); grid(ax,'on'); set(ax,'XScale','log'); xticks(ax,[40 100 150 250 400 800 1800]);
 title(ax,sprintf('q = %+.4f',x.q)); if i>4, xlabel(ax,'E (meV)'); end
 if mod(i,4)==1, ylabel(ax,'block / other blocks: measured (thin), shift prediction (dashed)'); end
end
exportgraphics(f,fullfile(out,'fig7_block_shapes_vs_shift.png'),'Resolution',110); close(f);

% ---- count cube for the step-2 operator (all channels |q| <= q_max, full energy axis)
cube=struct('E',E,'q',qs,'dq',dq,'note','J(E,channel,block): A1-aligned counts summed over the block frames; d: common-beam reference summed over the block');
for nbk=[nb cfg.n_blocks_check]
 bk=ceil((1:T).'/(T/nbk)); J=zeros(numel(E),nk,nbk); d=zeros(1,nbk); cE=d; fE=d;
 for b=1:nbk, J(:,:,b)=sum(A(:,:,bk==b),3); d(b)=sum(ref(bk==b)); [cE(b),fE(b)]=local_centre(sum(J(:,:,b),2),E); end
 Pk=zeros(nbk,nk); pos=zeros(1,nbk); d0=jb(1).shape;
 for b=1:nbk
  Pk(b,:)=sum(P(bk==b,:,jw(1)),1);
  pp=local_fit_dspv(Pk(b,:),qs,dq,[log(max(Pk(b,:))) mean(g_mu(bk==b,1))/dq-0.35],opts,d0); pos(b)=pp(2)*dq;
 end
 cube.(sprintf('b%d',nbk))=struct('blk',bk,'J',J,'d',d,'c_allq',cE,'fwhm_allq',fE,'core_profile_w8',Pk,'pos_narrow_mode',pos, ...
  'shape_w8',d0,'note','pos from the six-block common two-Pearson shape (w = 8), amplitude and position fitted per block');
end
cube.jb=jb; cube.all_frame_two_pearson_w8=bt([bt.w]==cfg.core_w(1) & [bt.block]==0).dspv; cube.c_q=cq; cube.ct=ct;
fi=dir(C.raw_path); [~,hs]=system(['shasum -a 256 "' C.raw_path '"']); hs=strtok(strtrim(hs));
cube.raw=struct('path',C.raw_path,'bytes',fi.bytes,'datenum',fi.datenum,'sha256',hs);
save(fullfile(out,'p2a_cube.mat'),'cube','-v7.3');

% ---- save
save(fullfile(out,'inputs_p2a.mat'),'cfg','qs','cq','fwq','ct','fwt','dE_core','dE_0','dE_wing','dE_all','call','fwall','c_all','fw_all','rmsE','P','Pfix','wl','Mb','ec','pf_mu','mu_blk','mu_all','jb', ...
 'mom_mu','mom_s','g_mu','g_s','g_rms','momfix_mu','momfix_s','bt','blk','sig_old','mu_old','prof_old','prof_al12','cb','fwb','tie','shp','dq','E');
meta=struct('config',cfg,'config_sha256',tag,'code_sha256',b1_sha256({'file',[mfilename('fullpath') '.m']}), ...
 'config_code_sha256',b1_sha256({'file',which('b1_p2a_inputs_config')}),'raw_path',C.raw_path,'run_root',run_root, ...
 'created',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')),'runtime_s',toc(t_start));
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
fprintf('P2a step 1 written to %s (%.0f s)\n',out,toc(t_start));
end


function [c, fw] = local_centre(y, E)
% ZLP centre (log-parabola through the maximum within |E| < 40 meV) and FWHM
% (linear interpolation), as b1_p1b_diag.
y = y(:); dE = E(2) - E(1); [~, i0] = max(y .* (abs(E) < 40));
ly = log(max(y(i0-1:i0+1), 1)); d = (ly(3) - ly(1)) / (2 * (2 * ly(2) - ly(1) - ly(3))); c = E(i0) + d * dE;
h = y(i0) / 2; l1 = find(y(1:i0) < h, 1, 'last'); r1 = i0 - 1 + find(y(i0:end) < h, 1, 'first'); fw = NaN;
if ~isempty(l1) && ~isempty(r1)
    fw = interp1(y([r1-1 r1]), E([r1-1 r1]), h) - interp1(y([l1 l1+1]), E([l1 l1+1]), h);
end
end

function wv = local_window(E, dE, c, w)
% Fraction of each energy channel [E - dE/2, E + dE/2] inside [c - w, c + w];
% E is a column, c a row (one centre per q channel).
wv = max(min(E + dE / 2, c + w) - max(E - dE / 2, c - w), 0) / dE;
end

function [mu, s] = local_moments(P, qv, dq)
% Centroid and width of each row of P over the channels qv (channel width^2/12 removed).
W = P ./ sum(P, 2); qv = qv(:); mu = W * qv; s = sqrt(max(W * (qv .^ 2) - mu .^ 2 - dq ^ 2 / 12, 0));
end

function [mu, s] = local_moments_argmax(P, qv, dq)
% Old definition (b1_adaptq_build): five channels around the maximum.
[~, i0] = max(P); kk = -2:2; [mu, s] = local_moments(P(i0 + kk), qv(i0 + kk), dq);
end

function y = local_gint(qv, dq, mu, s, A)
% Counts of a Gaussian (total A) in channels [qv - dq/2, qv + dq/2] (erfc for the far tails).
a = (qv - dq / 2 - mu) / (sqrt(2) * s); b = (qv + dq / 2 - mu) / (sqrt(2) * s);
y = A / 2 * (erf(b) - erf(a));
hi = a > 0; y(hi) = A / 2 * (erfc(a(hi)) - erfc(b(hi)));
lo = b < 0; y(lo) = A / 2 * (erfc(-b(lo)) - erfc(-a(lo)));
end

function [mu, s, r] = local_gauss(P, qv, dq, mu0, s0, opts)
% Channel-integrated Gaussian fitted to the profile with relative-error weights (log residuals).
P = P(:); qv = qv(:);
f = @(p) log(local_gint(qv, dq, p(2) * dq, exp(p(3)) * dq, exp(p(1)))) - log(P);
p = lsqnonlin(f, [log(sum(P)) mu0 / dq log(max(s0, 0.3 * dq) / dq)], [-Inf -3 log(0.1)], [Inf 3 log(5)], opts);
mu = p(2) * dq; s = exp(p(3)) * dq; r = f(p);
end

function d = local_spv_density(Q, dq, p)
% Split Pearson VII in counts per channel at the points Q; p = [log A, mu/dq,
% log(gamma_minus/dq), log(gamma_plus/dq), log m], A = counts per channel at mu.
mu = p(2) * dq; g = exp(p(3)) * dq * ones(size(Q)); g(Q > mu) = exp(p(4)) * dq;
d = exp(p(1)) * (1 + ((Q - mu) ./ g) .^ 2) .^ (-exp(p(5)));
end

function m = local_spv(qv, dq, p)
% Channel-integrated split Pearson VII (21 samples per channel).
ns = 21; xs = ((1:ns) - 0.5) / ns - 0.5; m = mean(local_spv_density(qv(:) + dq * xs, dq, p), 2);
end

function fw = local_spv_fwhm(p, dq)
% FWHM of the split Pearson VII (before channel integration).
fw = (exp(p(3)) + exp(p(4))) * dq * sqrt(2 ^ (1 / exp(p(5))) - 1);
end

function [p, mdl, r] = local_fit_spv(P, qv, dq, mu0, opts, shape)
% Split Pearson VII fitted with log residuals; shape = [log(gm/dq) log(gp/dq) log m]
% fixes the shape (amplitude and centroid free).
P = P(:); qv = qv(:);
if isempty(shape)
    p0 = [log(max(P)) mu0 / dq log(0.8) log(0.8) log(1.35)];
    lb = [-Inf -3 log(0.1) log(0.1) log(0.55)]; ub = [Inf 3 log(20) log(20) log(10)];
    p = lsqnonlin(@(u) log(local_spv(qv, dq, u)) - log(P), p0, lb, ub, opts);
else
    u = lsqnonlin(@(u) log(local_spv(qv, dq, [u shape])) - log(P), [log(max(P)) mu0 / dq], [-Inf -3], [Inf 3], opts);
    p = [u shape];
end
mdl = local_spv(qv, dq, p); r = log(mdl) - log(P);
end

function m = local_dspv(qv, dq, p)
% Two split Pearson VII components, channel-integrated (21 samples per channel);
% p = [log A1, mu1/dq, log(A2/A1), (mu2 - mu1)/dq, log(g1m/dq), log(g1p/dq), log m1,
% log(g2m/dq), log(g2p/dq), log m2]; A in counts per channel at the mode.
ns = 21; xs = ((1:ns) - 0.5) / ns - 0.5; m = mean(local_dspv_density(qv(:) + dq * xs, dq, p), 2);
end

function d = local_dspv_density(Q, dq, p)
d = local_spv_density(Q, dq, [p(1) p(2) p(5:7)]) + local_spv_density(Q, dq, [p(1) + p(3) p(2) + p(4) p(8:10)]);
end

function [lb, ub] = local_dspv_bounds()
lb = [-Inf -3 log(1e-4) -8 log(0.05) log(0.05) log(0.55) log(1) log(1) log(0.55)];
ub = [Inf 3 log(1) 8 log(5) log(5) log(20) log(40) log(40) log(10)];
end

function p0 = local_dspv_start(P, gm, dq)
p0 = [log(max(P)) gm / dq - 0.35 log(0.01) -1.85 log(1.1) log(1.8) log(2.5) log(2.4) log(2.7) log(1.2)];
end

function [p, mdl, r] = local_fit_dspv(P, qv, dq, p0, opts, shape)
% Two-component profile fitted with log residuals. Free fit: p0 is the full
% start vector. With shape = p(3:10) fixed, p0 = [log A1, mu1/dq] and only the
% amplitude and position are fitted.
P = P(:); qv = qv(:);
if isempty(shape)
    [lb, ub] = local_dspv_bounds(); p0 = min(max(p0, lb), ub);
    p = lsqnonlin(@(u) log(local_dspv(qv, dq, u)) - log(P), p0, lb, ub, opts);
else
    u = lsqnonlin(@(u) log(local_dspv(qv, dq, [u shape])) - log(P), p0, [-Inf -3], [Inf 3], opts); p = [u shape];
end
mdl = local_dspv(qv, dq, p); r = log(mdl) - log(P);
end

function st = local_dspv_stats(p, dq, qmax)
% [narrow FWHM, broad fraction of the total, fraction outside |q| <= qmax,
% mean position, narrow mode], numerically to |q| = 0.5.
qf = (-0.5:dq / 40:0.5).'; d1 = local_spv_density(qf, dq, [p(1) p(2) p(5:7)]); d2 = local_spv_density(qf, dq, [p(1) + p(3) p(2) + p(4) p(8:10)]);
d = d1 + d2; h = max(d1) / 2; i1 = find(d1 >= h, 1, 'first'); i2 = find(d1 >= h, 1, 'last');
st = [qf(i2) - qf(i1), sum(d2) / sum(d), sum(d(abs(qf) > qmax)) / sum(d), sum(qf .* d) / sum(d), p(2) * dq];
end

function row = local_dspv_row(p, dq, qmax)
% Values for w_kernel_x.csv.
st = local_dspv_stats(p, dq, qmax);
row = [p(2) * dq, st(4), exp(p(5)) * dq, exp(p(6)) * dq, 2 * exp(p(7)), st(1), p(4) * dq, exp(p(8)) * dq, exp(p(9)) * dq, 2 * exp(p(10)), st(2), st(3)];
end

function r = local_joint_res(u, Pj, qs, dq, nb)
% Log residuals of the block profiles (rows of Pj) for a common two-component shape u(1:8).
r = zeros(size(Pj, 2), nb);
for b = 1:nb, r(:, b) = log(local_dspv(qs, dq, [u(7 + 2 * b) u(8 + 2 * b) u(1:8)])) - log(Pj(b, :).'); end
end

function r = pcorr(a, b)
% Pearson correlation over the finite pairs (no Statistics Toolbox needed).
m = isfinite(a) & isfinite(b); a = a(m) - mean(a(m)); b = b(m) - mean(b(m));
r = sum(a .* b) / sqrt(sum(a .^ 2) * sum(b .^ 2));
end
