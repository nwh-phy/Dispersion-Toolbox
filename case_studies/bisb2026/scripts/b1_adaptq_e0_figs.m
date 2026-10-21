function T = b1_adaptq_e0_figs(run_root, session, stage, key_ref)
% B1 adaptive-binning task: figures to inspect the fitted Drude-Lorentz E0
% (= omega_p in Do et al. 2025) against the data, for one fit stage.
% Measured loss function per bin: (observed - ZLP - aux - constant) / I_kin,
% with I_kin the 2D kinematic prefactor used in that fit; the fitted loss
% function is the sum of the two DL components. Writes
% <run_root>/stage5/figures_<stage>/ and e0_<session>.csv.
if nargin<2, session='590_PL2_10w'; end
if nargin<3, stage='stage2_2d_h5'; end
if nargin<4, key_ref='lorentz_n2_2d_h20'; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
figd=fullfile(run_root,'stage5',['figures_' stage]); if ~isfolder(figd), mkdir(figd); end
d=dir(fullfile(run_root,stage,session,'bin_*.mat'));
cf=fullfile(run_root,'stage1',sprintf('channels_%s.mat',session)); S=load(cf,'C'); C=S.C; mf=matfile(cf);
Y=zeros(numel(C.E),numel(C.q)); for k=1:numel(C.q), Y(:,k)=sum(double(mf.frames(:,:,k)),2); end
Ex=C.E(:); rx=Ex>=250&Ex<=1800; Ex=Ex(rx); dl=peak_models('lorentz');
B=struct('q',{},'N',{},'members',{},'E',{},'Lm',{},'Lmod',{},'Lc',{},'E0',{},'G',{},'E0sd',{},'apex_c',{},'apex_m',{},'apex_mod',{},'status',{},'E0ref',{},'Ex',{},'Lmx',{},'Lcx',{});
for k=1:numel(d)
 L=load(fullfile(d(k).folder,d(k).name),'R'); R=L.R; m=R.fits(1); c=m.curves;
 if ~m.success || isempty(fieldnames(c)), continue; end
 E=c.E(:); bg=c.background(:); Kp=c.prefactor(:);
 Lm=(c.observed(:)-bg)./Kp; Lc=c.loss_function; Lmod=sum(Lc,2);
 E0=m.parameters(:,1)'; G=m.parameters(:,2)';
 ac=arrayfun(@(i) apex_e(E0(i),G(i)),1:2);
 % Cross-check on the fit window: measured counts of the bin over the 2D prefactor
 % recomputed with the fit's own settings (background < 1 % above 1500 meV,
 % neglected beyond the window), and the two DL evaluated analytically.
 ch=arrayfun(@(x) find(abs(C.q-x)<1e-7,1),R.members); Yb=sum(Y(:,ch),2);
 Kf=qe_kinematic_prefactor(R.members,dq=C.dq,form=m.form,h_perp=m.h,sigma_probe=C.sigma_probe,beam_kV=30);
 Kx=reshape(Kf(max(Ex,250)),[],1); w0=E>=300&E<=1800;
 kchk=median(Kp(w0)./reshape(Kf(max(E(w0),250)),[],1)); assert(abs(kchk-1)<1e-6,'prefactor mismatch %g',kchk);
 bgx=interp1(E,bg,Ex,'linear',0); Lmx=(Yb(rx)-bgx)./Kx;
 Lcx=zeros(numel(Ex),2);
 for i=1:2, u1=dl.model_fn(E0(i),G(i),1,E(w0)); sc_i=median(Lc(w0,i)./u1); Lcx(:,i)=sc_i*dl.model_fn(E0(i),G(i),1,Ex); end
 w=E>=300&E<=1800; s=movmean(Lm,51); [~,i]=max(s.*w-1e30*~w); apm=E(i);
 [~,i]=max(Lmod.*w-1e30*~w); apmod=E(i);
 sd=nan(1,2); if ~isempty(R.boot), P=nan(numel(R.boot),2); for j=1:numel(R.boot), if R.boot(j).success, P(j,:)=R.boot(j).parameters(:,1)'; end, end, sd=std(P,0,1,'omitnan'); end
 e0r=nan(1,2); f2=fullfile(run_root,'stage2',session,d(k).name);
 if isfile(f2), L2=load(f2,'R'); F2=L2.R.fits; j=find(strcmp({F2.key},key_ref),1); if ~isempty(j)&&F2(j).success, e0r=F2(j).parameters(:,1)'; end, end
 B(end+1)=struct('q',R.q,'N',R.N,'members',R.members,'E',E,'Lm',Lm,'Lmod',Lmod,'Lc',Lc,'E0',E0,'G',G,'E0sd',sd, ...
  'apex_c',ac,'apex_m',apm,'apex_mod',apmod,'status',m.status,'E0ref',e0r,'Ex',Ex,'Lmx',Lmx,'Lcx',Lcx); %#ok<AGROW>
end
[~,o]=sort([B.q]); B=B(o); q=[B.q]; aq=abs(q);
E0=reshape([B.E0],2,[])'; G=reshape([B.G],2,[])'; sd=reshape([B.E0sd],2,[])'; ac=reshape([B.apex_c],2,[])'; e0r=reshape([B.E0ref],2,[])';
T=table(q(:),[B.N]',{B.status}',E0(:,1),sd(:,1),G(:,1),ac(:,1),E0(:,2),sd(:,2),G(:,2),ac(:,2),[B.apex_m]',[B.apex_mod]',e0r(:,1),e0r(:,2), ...
 'VariableNames',{'q','N','status','E0_low','E0_low_sd','G_low','apex_low','E0_high','E0_high_sd','G_high','apex_high','apex_measured_LF','apex_model_LF','E0_low_h20','E0_high_h20'});
writetable(T,fullfile(run_root,'stage5',sprintf('e0_%s_%s.csv',stage,session)));
cl=[0.10 0.35 0.75; 0.80 0.25 0.10]; sel=aq<=0.0241;

% ---- Figure 1: raw counts map and measured loss function map, E0 on top
f=figure('Visible','off','Position',[0 0 1700 700]); tl=tiledlayout(f,1,2,'TileSpacing','compact');
ax=nexttile(tl); [qs,o2]=sort(C.q); u=abs(qs)<=0.0245; M=movmean(Y(:,o2(u)),9,1); r=C.E>=250&C.E<=1800; M=M(r,:);
M=M./max(M(C.E(r)>=400&C.E(r)<=1800,:),[],1); imagesc(ax,qs(u),C.E(r),M); axis(ax,'xy'); colormap(ax,turbo); clim(ax,[0.4 1]); hold(ax,'on');
for b=1:2, plot(ax,q(sel),E0(sel,b),'o','MarkerSize',5,'MarkerFaceColor',cl(b,:),'MarkerEdgeColor','w'); end
title(ax,'(a) raw counts, native channels, column-normalised; dots = fitted E0 (2D, h=0.0005)'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'energy loss (meV)'); ylim(ax,[250 1800]);
ax=nexttile(tl); hold(ax,'on'); Eg=B(1).E; r=Eg>=250&Eg<=1800; Er=Eg(r);
ed=zeros(1,numel(B)+1); % bin edges along q
for b=1:numel(B), w=0.0005*B(b).N; ed(b)=B(b).q-w/2; ed(b+1)=B(b).q+w/2; end
Mq=nan(nnz(r),numel(B));
for b=1:numel(B), s=movmean(B(b).Lm,25); s=s(r); Mq(:,b)=s/max(s(Er>=400&Er<=1800)); end
us=find(sel); xe=[ed(us) ed(us(end)+1)]; ye=[Er; Er(end)+4];
Cm=Mq(:,us); Cm=[Cm nan(size(Cm,1),1); nan(1,size(Cm,2)+1)];
pcolor(ax,xe,ye,Cm); shading(ax,'flat'); colormap(ax,turbo); clim(ax,[0.4 1]); colorbar(ax);
for b=1:2
 plot(ax,q(sel),E0(sel,b),'o','MarkerSize',5,'MarkerFaceColor',cl(b,:),'MarkerEdgeColor','w');
 plot(ax,q(sel),ac(sel,b),'^','MarkerSize',5,'MarkerEdgeColor',cl(b,:),'MarkerFaceColor','none','LineWidth',1);
end
title(ax,'(b) measured loss function (obs - ZLP - aux)/I_{kin}^{2D}, column-normalised; dots E0, triangles component apex');
xlabel(ax,'q (A^{-1}), bin edges'); ylim(ax,[250 1800]); xlim(ax,[-0.0245 0.0245]);
exportgraphics(f,fullfile(figd,sprintf('maps_E0_%s.png',session)),'Resolution',130); close(f);

% ---- Figure 2: single spectra, measured vs fitted loss function, E0 marked
pick=[0.0005 0.0010 0.0015 0.0020 0.0030 0.0040 0.0050 0.0060 0.0085 0.0115 0.0155 0.0230];
for sg=[1 -1]
 f=figure('Visible','off','Position',[0 0 1800 1150]); tl=tiledlayout(f,3,4,'TileSpacing','compact','Padding','compact');
 for p=pick
  b=find(abs(q-sg*p)<1e-6,1); ax=nexttile(tl); hold(ax,'on'); box(ax,'on');
  if isempty(b), continue; end
  X=B(b); in=X.E>=250&X.E<=1800; E=X.E(in); sc=max(X.Lmod(X.E>=300&X.E<=1800));
  plot(ax,E,X.Lm(in)/sc,'Color',[0.75 0.75 0.75]);
  plot(ax,E,movmean(X.Lm(in),25)/sc,'k','LineWidth',1);
  plot(ax,E,X.Lmod(in)/sc,'m','LineWidth',1.4);
  for c=1:2
   plot(ax,E,X.Lc(in,c)/sc,'Color',cl(c,:),'LineWidth',1.2);
   xline(ax,X.E0(c),'-','Color',cl(c,:),'LineWidth',1.6);
   xline(ax,X.apex_c(c),':','Color',cl(c,:),'LineWidth',1.4);
  end
  ylim(ax,[-0.1 1.3]); xlim(ax,[250 1800]); grid(ax,'on');
  title(ax,sprintf('q = %+.4f (N=%d)  E0 = %.0f / %.0f,  \\Gamma = %.0f / %.0f',X.q,X.N,X.E0(1),X.E0(2),X.G(1),X.G(2)),'FontSize',9);
  if p==pick(1)
   lg=legend(ax,{'measured LF','measured LF, 100-meV smooth','fitted LF (sum of 2 DL)','lower DL','E0 lower','apex lower DL','upper DL','E0 upper','apex upper DL'},'FontSize',8);
   lg.Layout.Tile='east';
  end
 end
 xlabel(tl,'energy loss (meV)'); ylabel(tl,'loss function (normalised to fitted max in 300-1800 meV)');
 title(tl,sprintf('%s, q %s 0: measured vs fitted loss function in the fit window; solid lines = E0, dotted = apex of each DL (2D I_{kin}, h = 0.0005)',strrep(session,'_','\_'),char('<'*(sg<0)+'>'*(sg>0))));
 exportgraphics(f,fullfile(figd,sprintf('spectra_E0_%s_%s.png',session,char('m'*(sg<0)+'p'*(sg>0)))),'Resolution',120); close(f);
end

% ---- Figure 3: E0(q)
f=figure('Visible','off','Position',[0 0 1000 650]); ax=axes(f); hold(ax,'on'); grid(ax,'on');
patch(ax,[0 0.0025 0.0025 0],[0 0 2000 2000],[0.93 0.93 0.93],'EdgeColor','none','DisplayName','q < q_E(1.6 eV)');
mk={'o','s'}; nm={'lower','upper'}; hv={'off','on'};
for sg=[-1 1]
 u=sel&sign(q)==sg; ii=(sg>0)+1;
 for b=1:2
  errorbar(ax,aq(u),E0(u,b),sd(u,b),mk{ii},'Color',cl(b,:),'MarkerFaceColor',cl(b,:),'MarkerSize',5,'CapSize',3, ...
   'DisplayName',sprintf('E0 %s, q%s0 (2D, h=0.0005; bars: frame bootstrap)',nm{b},char('<'*(sg<0)+'>'*(sg>0))));
  plot(ax,aq(u),e0r(u,b),mk{ii},'Color',cl(b,:),'MarkerSize',4,'HandleVisibility','off');
 end
 am=[B(u).apex_m]; au=aq(u); e=am>=1790;
 plot(ax,au(~e),am(~e),'kx','MarkerSize',7,'LineWidth',1.2,'HandleVisibility',hv{(sg>0)+1},'DisplayName','apex of measured LF in 300-1800 meV (100-meV smooth)');
 plot(ax,au(e),am(e),'k^','MarkerSize',6,'HandleVisibility',hv{(sg>0)+1},'DisplayName','measured LF still rising at 1800 meV (no apex in window)');
end
plot(ax,nan,nan,'o','Color',[0.4 0.4 0.4],'DisplayName','open: E0 with h = 0.002');
xlabel(ax,'|q| (A^{-1})'); ylabel(ax,'energy (meV)'); xlim(ax,[0 0.0245]); ylim(ax,[300 2000]);
legend(ax,'Location','eastoutside','FontSize',7);
title(ax,sprintf('%s: fitted DL E0 (= \\omega_p in Do et al. 2025) vs |q|',strrep(session,'_','\_')));
exportgraphics(f,fullfile(figd,sprintf('E0_vs_q_%s.png',session)),'Resolution',140); close(f);
end

function a=apex_e(E0,G)
c=2*E0.^2-G.^2; a=sqrt((c+sqrt(c.^2+12*E0.^4))/6);
end
