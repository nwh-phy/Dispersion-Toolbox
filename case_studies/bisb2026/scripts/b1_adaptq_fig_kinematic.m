function out = b1_adaptq_fig_kinematic(run_root, session, stage)
% Why the kinematic factor matters at 30 kV and small q.
% (a) 2D kinematic factor I_kin(E) of single native channels (dq, h, probe
%     spread as in the fits), normalised at 1.5 eV, for several |q|; dashed:
%     80 kV, point q = 0.008 (conditions of Do et al. 2025) for comparison.
% (b) B1 signal in counts: observed minus the fitted ZLP tail and low-energy
%     losses, per bin of the 2D, h = 0.0005 fits (stored curves).
% (c) the same divided by the kinematic factor of each fit (measured loss
%     function). (b) and (c) differ only by that factor: both get the same
%     100 meV moving mean along E, and each column is normalised to its
%     maximum in 0.3-1.8 eV.
% Prints I_kin(0.7 eV)/I_kin(1.4 eV) per |q|. Writes
% stage5/figures_<stage>/kinematic_effect_<session>.png/.pdf.
arguments
 run_root char
 session char = '590_PL2_10w'
 stage char = 'stage2_2d_h5'
end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
S=load(fullfile(run_root,'stage1',sprintf('channels_%s.mat',session)),'C'); C=S.C;
d=dir(fullfile(run_root,stage,session,'bin_*.mat'));
B=struct('q',{},'N',{},'sig',{},'lf',{},'h',{}); E=[];
for k=1:numel(d)
 L=load(fullfile(d(k).folder,d(k).name),'R'); R=L.R; m=R.fits([R.fits.main]);
 if ~m.success || isempty(fieldnames(m.curves)), continue; end
 c=m.curves;
 if isempty(E), E=c.E(:); else, assert(isequal(E,c.E(:)),'b1_adaptq_fig_kinematic:Grid','bins use different energy grids'); end
 s=c.observed(:)-c.background(:);
 B(end+1)=struct('q',R.q,'N',R.N,'sig',s,'lf',s./c.prefactor(:),'h',m.h); %#ok<AGROW>
end
[~,o]=sort([B.q]); B=B(o); h=unique([B.h]); assert(isscalar(h),'b1_adaptq_fig_kinematic:H','mixed h in %s',stage);
% Display grid offset by dq/4 so that every bin (q +/- N*dq/2) covers whole
% grid cells; both maps get the same 100 meV moving mean along E.
qmax=0.0245; qs=(-qmax+C.dq/4):C.dq/2:qmax; r=E>=300&E<=1800; Er=E(r)/1e3;
nsm=max(1,round(100/median(diff(E))));
Ms=nan(nnz(r),numel(qs)); Ml=Ms;
for k=1:numel(B)
 half=B(k).N*C.dq/2; col=qs>B(k).q-half & qs<B(k).q+half;
 s=movmean(B(k).sig,nsm); l=movmean(B(k).lf,nsm); s=s(r); l=l(r);
 Ms(:,col)=repmat(s/max(s),1,nnz(col)); Ml(:,col)=repmat(l/max(l),1,nnz(col));
end

f=figure('Visible','off','Units','centimeters','Position',[2 2 33 10.5],'Color','w'); if isprop(f,'Theme'), f.Theme='light'; end
T=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
% (a) kinematic factor
ax=nexttile(T); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
Eg=linspace(300,1800,301); qa=[0.001 0.002 0.003 0.005 0.010 0.020]; cm=parula(numel(qa)+1);
fprintf('I_kin(0.7 eV)/I_kin(1.4 eV), 30 kV, 2D, single channel, h = %g, sigma_probe = %.3g:\n',h,C.sigma_probe);
for i=1:numel(qa)
 K=qe_kinematic_prefactor(qa(i),dq=C.dq,form='2d',h_perp=h,sigma_probe=C.sigma_probe,beam_kV=30);
 y=K(Eg)/K(1500); plot(ax,Eg/1e3,y,'-','Color',cm(i,:),'LineWidth',1.4,'DisplayName',sprintf('|q| = %.3f',qa(i)));
 fprintf('  |q| = %.3f: %.2f\n',qa(i),K(700)/K(1400));
end
K80=qe_kinematic_prefactor(0.008,dq=1e-7,form='2d',h_perp=0,sigma_probe=0,beam_kV=80);
plot(ax,Eg/1e3,K80(Eg)/K80(1500),'k--','LineWidth',1.2,'DisplayName','80 kV, q = 0.008');
fprintf('  80 kV, q = 0.008 (point): %.2f\n',K80(700)/K80(1400));
set(ax,'YScale','log','FontName','Helvetica','FontSize',9,'TickDir','in','LineWidth',0.6); xlim(ax,[0.3 1.8]); ylim(ax,[0.8 40]);
set(ax,'YTick',[1 2 5 10 20 40]);
xlabel(ax,'Energy loss (eV)'); ylabel(ax,'Kinematic factor / value at 1.5 eV');
title(ax,'(a) Kinematic factor, 30 kV (q in Å^{-1})','FontWeight','normal','FontSize',10);
legend(ax,'Location','northeast','FontSize',8,'Box','off');
% (b), (c) maps
titles={'(b) Counts, ZLP tail removed','(c) Divided by the kinematic factor'};
maps={Ms,Ml};
for p=1:2
 ax=nexttile(T); M=maps{p};
 imagesc(ax,qs,Er,M,'AlphaData',~isnan(M)); axis(ax,'xy'); colormap(ax,turbo); clim(ax,[0.3 1]);
 set(ax,'Color','w','FontName','Helvetica','FontSize',9,'TickDir','out','LineWidth',0.6,'Layer','top'); box(ax,'on');
 xlim(ax,[-qmax qmax]); ylim(ax,[0.3 1.8]); xlabel(ax,'q (Å^{-1})'); if p==1, ylabel(ax,'Energy loss (eV)'); end
 title(ax,titles{p},'FontWeight','normal','FontSize',10);
end
cb=colorbar(ax); cb.Label.String='normalised in each column'; cb.FontSize=8;
outd=fullfile(run_root,'stage5',['figures_' stage]); if ~isfolder(outd), mkdir(outd); end
out=fullfile(outd,sprintf('kinematic_effect_%s',session));
exportgraphics(f,[out '.png'],'Resolution',250,'BackgroundColor','white');
exportgraphics(f,[out '.pdf'],'ContentType','image','Resolution',300,'BackgroundColor','white');
close(f); fprintf('saved %s.png/.pdf\n',out);
end
