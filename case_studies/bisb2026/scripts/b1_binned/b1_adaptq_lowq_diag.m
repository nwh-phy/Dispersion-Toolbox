function b1_adaptq_lowq_diag(run_root, session, stage, qlist)
% Low-q diagnostic of the joint fit in count space: observed spectrum, ZLP,
% auxiliary DL peaks (30-80 and 80-300 meV, no kinematic factor), the two B1
% DL components (E0 bounded to the signal window 300-1800 meV, times I_kin),
% total model and normalised residual. Linear (0-1000 meV) and log (-200-1800
% meV) rows. Writes <run_root>/stage5/figures_<stage>/lowq_diag_<session>.png.
if nargin<2, session='590_PL2_10w'; end
if nargin<3, stage='stage2_2d_h5'; end
if nargin<4, qlist=[0.0005 0.0010 0.0020 0.0050]; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
figd=fullfile(run_root,'stage5',['figures_' stage]); if ~isfolder(figd), mkdir(figd); end
f=figure('Visible','off','Position',[0 0 1900 1050]); tl=tiledlayout(f,3,numel(qlist),'TileSpacing','compact','Padding','compact');
cA=[0.20 0.60 0.20; 0.55 0.80 0.30]; cB=[0.10 0.35 0.75; 0.80 0.25 0.10];
for row=1:3
 for qa=qlist
  d=dir(fullfile(run_root,stage,session,sprintf('bin_%+.4f_N*.mat',qa))); ax=nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
  L=load(fullfile(d(1).folder,d(1).name),'R'); m=L.R.fits(1); c=m.curves; E=c.E;
  B1=c.peaks; A=c.aux_peaks; z=c.zlp;
  if row<3
   plot(ax,E,c.observed,'Color',[0.6 0.6 0.6],'DisplayName','observed');
   plot(ax,E,z,'k--','DisplayName','ZLP (2 Pearson VII)');
   for j=1:size(A,2), plot(ax,E,A(:,j),'Color',cA(j,:),'LineWidth',1.4,'DisplayName',sprintf('aux DL %d: E0 %.0f, \\Gamma %.0f',j,m.aux_parameters(j,1),m.aux_parameters(j,2))); end
   for j=1:2, plot(ax,E,B1(:,j),'Color',cB(j,:),'LineWidth',1.4,'DisplayName',sprintf('B1 DL %d x I_{kin}: E0 %.0f, \\Gamma %.0f',j,m.parameters(j,1),m.parameters(j,2))); end
   plot(ax,E,c.prediction,'m','LineWidth',1.1,'DisplayName','total model');
   xline(ax,300,':','B1 E0 lower bound 300','HandleVisibility','off'); xline(ax,250,':','I_{kin} floor 250','LabelVerticalAlignment','bottom','HandleVisibility','off');
  end
  if row==1
   xlim(ax,[0 1000]); w=E>=60&E<=1000; ylim(ax,[0 1.3*max(c.observed(w))]);
   title(ax,sprintf('q = %+.4f (N=%d), linear',L.R.q,L.R.N));
   if qa==qlist(1), legend(ax,'Location','northeast','FontSize',7); end
   % share of 100-400 meV counts assigned to each part
   w=E>=100&E<=400; tot=sum(c.observed(w));
   text(ax,520,0.55*1.3*max(c.observed(E>=60&E<=1000)),sprintf('100-400 meV counts:\n ZLP %.0f%%, aux %.0f%%, B1 %.0f%%', ...
    100*sum(z(w))/tot,100*sum(sum(A(w,:),2))/tot,100*sum(sum(B1(w,:),2))/tot),'FontSize',8);
  elseif row==2
   set(ax,'YScale','log'); xlim(ax,[-200 1800]); ylim(ax,[max(1,0.3*min(c.observed(E<=1800))) 2*max(c.observed)]);
   title(ax,'log scale');
  else
   r=(c.observed-c.prediction)./c.sigma; plot(ax,E,r,'Color',[0.4 0.4 0.4]); plot(ax,E,movmean(r,25),'k','LineWidth',1.2);
   xlim(ax,[-200 1800]); ylim(ax,[-6 6]); yline(ax,0,'r'); title(ax,'(obs - model)/\sigma, 100-meV running mean in black');
   xline(ax,300,':'); xline(ax,250,':');
  end
  xlabel(ax,'energy loss (meV)');
 end
end
title(tl,sprintf('%s low-q joint fits (%s: 2D I_{kin}, h = 0.0005); counts space',strrep(session,'_','\_'),strrep(stage,'_','\_')));
exportgraphics(f,fullfile(figd,sprintf('lowq_diag_%s.png',session)),'Resolution',120); close(f);
end
