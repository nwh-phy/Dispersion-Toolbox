function out = b1_fullq_fig_twopeak(run_root, session, abs_q)
% Proof figure for the two-component B1: single spectra at +/-|q| (590 by
% default), each with the two DL components of the main fit (3D prefactor,
% aux2f, h = 0.002) on top of the ZLP tail + low-energy losses, and a
% residual strip comparing the one- and two-peak fits (residuals in units of
% sigma, averaged over 40 meV and scaled by sqrt(N) so that pure noise has
% unit spread). Both fits are recomputed from the stored optima and checked
% against the stored cost. Writes PNG and PDF into <run_root>/stage6_figures.
arguments
 run_root char
 session char = '590_PL2_10w'
 abs_q (1,:) double = [0.0025 0.0055 0.0085 0.0115]
end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_fullq_config(); outd=fullfile(run_root,'stage6_figures'); if ~isfolder(outd), mkdir(outd); end
cU=[0.165 0.471 0.839]; cL=[0.922 0.408 0.204]; cB=[0.80 0.80 0.80]; c1=[0.80 0.20 0.20];
nc=numel(abs_q); signs=[-1 1];
f=figure('Visible','off','Units','centimeters','Position',[2 2 32 15.5]); if isprop(f,'Theme'), f.Theme='light'; end
set(f,'Color','w');
T=tiledlayout(f,2,nc,'TileSpacing','compact','Padding','compact');
for si=1:2
 for k=1:nc
  q=signs(si)*abs_q(k);
  L=load(fullfile(run_root,'stage2',session,sprintf('bin_%+.4f.mat',q)),'R'); R=L.R;
  m=R.fits([R.fits.main]); key1=strrep(m.key,'_n2_','_n1_'); n1=R.fits(strcmp({R.fits.key},key1));
  c=m.curves; E=c.E; Y=c.observed; sg=c.sigma;
  Kfun=qe_kinematic_prefactor(R.members,form=m.form,h_perp=m.h,sigma_probe=R.sigma_probe,beam_kV=cfg.beam_kV);
  fit=@(n,u) qe_zlp_joint_fit(E,Y,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
   'peak_model',m.model,'n_peaks',n,'n_starts',1,'seed',cfg.seed,'aux_windows',cfg.aux{1,2},'n_zlp',cfg.n_zlp, ...
   'prefactor',Kfun,'prefactor_floor_meV',cfg.prefactor_floor,'prefactor_on_aux',cfg.prefactor_on_aux, ...
   'noise_sigma',sg,'warm_u',u);
  f2=fit(2,m.u); f1=fit(1,n1.u);
  fprintf('%s q=%+.4f: cost n2 %.1f (stored %.1f), n1 %.1f (stored %.1f), dchi2 %.0f; E0 %s\n',session,q, ...
   f2.cost,m.cost,f1.cost,n1.cost,f1.cost-f2.cost,mat2str(round(f2.parameters(:,1)')));
  w=E>=200&E<=cfg.window(2); x=E(w); s=1e-3;
  bg=f2.zlp(w)+sum(f2.aux_peaks(w,:),2); lo=f2.peaks(w,1); up=f2.peaks(w,2);
  inner=tiledlayout(T,4,1,'TileSpacing','none','Padding','none'); inner.Layout.Tile=(si-1)*nc+k;
  ax=nexttile(inner,[3 1]); hold(ax,'on');
  fill(ax,[x;flipud(x)],s*[zeros(size(bg));flipud(bg)],cB,'EdgeColor','none','FaceAlpha',0.6);
  fill(ax,[x;flipud(x)],s*[bg;flipud(bg+up)],cU,'EdgeColor','none','FaceAlpha',0.18);
  fill(ax,[x;flipud(x)],s*[bg;flipud(bg+lo)],cL,'EdgeColor','none','FaceAlpha',0.30);
  plot(ax,x,s*Y(w),'.','Color',[0.35 0.35 0.35],'MarkerSize',4);
  plot(ax,x,s*(bg+lo),'-','Color',cL,'LineWidth',1.3);
  plot(ax,x,s*(bg+up),'-','Color',cU,'LineWidth',1.3);
  plot(ax,x,s*f2.prediction(w),'k-','LineWidth',1.2);
  xlim(ax,[200 cfg.window(2)]); ylim(ax,[0 1.42*s*max(Y(w))]);
  set(ax,'XTickLabel',[],'FontName','Helvetica','FontSize',9,'Box','on','TickDir','in','LineWidth',0.6);
  yt=get(ax,'YTick'); set(ax,'YTick',yt(yt>0));
  text(ax,0.04,0.93,sprintf('\\itq\\rm = %+.4f Å^{-1}',q),'Units','normalized','FontSize',9.5,'FontName','Helvetica');
  text(ax,0.96,0.93,sprintf('\\Delta\\chi^2 = %.0f',f1.cost-f2.cost),'Units','normalized','FontSize',9.5, ...
   'HorizontalAlignment','right','FontName','Helvetica');
  text(ax,0.96,0.81,sprintf('\\chi^2_\\nu: %.2f \\rightarrow %.2f',f1.chi2_red_signal,f2.chi2_red_signal),'Units','normalized', ...
   'FontSize',8.5,'HorizontalAlignment','right','FontName','Helvetica','Color',[0.3 0.3 0.3]);
  text(ax,0.04,0.81,sprintf('\\itE\\rm_0: %.0f / %.0f meV',f2.parameters(1,1),f2.parameters(2,1)),'Units','normalized', ...
   'FontSize',8.5,'FontName','Helvetica','Color',[0.3 0.3 0.3]);
  if k==1, ylabel(ax,'Counts (10^3)'); end
  % residual strip: 40 meV bins (10 points), mean * sqrt(N)
  axr=nexttile(inner); hold(axr,'on');
  sel=E>=cfg.window(1)&E<=cfg.window(2); Es=E(sel); r2=(Y(sel)-f2.prediction(sel))./sg(sel); r1=(Y(sel)-f1.prediction(sel))./sg(sel);
  nb=floor(numel(Es)/10); idx=reshape(1:10*nb,10,nb);
  eb=mean(Es(idx),1); b2=mean(r2(idx),1)*sqrt(10); b1=mean(r1(idx),1)*sqrt(10);
  yline(axr,0,'-','Color',[0.6 0.6 0.6]);
  plot(axr,eb,b1,'-o','Color',c1,'MarkerSize',3,'MarkerFaceColor',c1,'LineWidth',1);
  plot(axr,eb,b2,'-o','Color','k','MarkerSize',3,'MarkerFaceColor','k','LineWidth',1);
  xlim(axr,[200 cfg.window(2)]); lim=max(4,ceil(max(abs([b1 b2]))*1.2)); ylim(axr,[-lim lim]);
  nice=[2 3 5 10 20 50]; t=nice(find(nice<=0.7*lim,1,'last')); if isempty(t), t=2; end; set(axr,'YTick',[-t 0 t]);
  set(axr,'FontName','Helvetica','FontSize',9,'Box','on','TickDir','in','LineWidth',0.6);
  if si==2, xlabel(axr,'Energy loss \itE\rm (meV)'); else, set(axr,'XTickLabel',[]); end
  if k==1, ylabel(axr,'res. (\sigma)'); end
 end
end
% shared legend
lg=axes(f,'Position',[0 0 1 1],'Visible','off'); hold(lg,'on');
h=[plot(lg,nan,nan,'.','Color',[0.35 0.35 0.35],'MarkerSize',10), plot(lg,nan,nan,'k-','LineWidth',1.2), ...
   fill(lg,nan,nan,cL,'FaceAlpha',0.4,'EdgeColor',cL), fill(lg,nan,nan,cU,'FaceAlpha',0.4,'EdgeColor',cU), ...
   fill(lg,nan,nan,cB,'EdgeColor','none'), plot(lg,nan,nan,'-o','Color',c1,'MarkerFaceColor',c1,'MarkerSize',3), ...
   plot(lg,nan,nan,'-o','Color','k','MarkerFaceColor','k','MarkerSize',3)];
legend(lg,h,{'data','two-peak fit','lower branch','upper branch','ZLP tail + low-energy losses','residual, one peak','residual, two peaks'}, ...
 'Orientation','horizontal','Location','north','Box','off','FontName','Helvetica','FontSize',9);
T.Padding='loose';
base=fullfile(outd,sprintf('twopeak_proof_%s',session));
exportgraphics(f,[base '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(f,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
close(f); out=base; fprintf('saved %s.png/.pdf\n',base);
end
