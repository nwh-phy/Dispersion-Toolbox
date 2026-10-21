function T = b1_adaptq_aggregate(run_root, fullq_root, sessions)
% B1 adaptive-binning task, stage 5: tables and figures. Besides the DL
% parameter E0, every component is reported by its mode energy:
%   pole  = Re of the pole of A*E*G/((E^2-E0^2)^2+E^2*G^2) = sqrt(E0^2-G^2/4)
%   apex  = maximum of that line shape
% With widths comparable to E0 the two lie well below E0, so a flat E0 with
% a narrowing width is a dispersing mode.
if nargin<3, sessions={'590_PL2_10w'}; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
out=fullfile(run_root,'stage5'); figd=fullfile(out,'figures'); if ~isfolder(figd), mkdir(figd); end
rows={};
for s=1:numel(sessions)
 ses=sessions{s}; d=dir(fullfile(run_root,'stage2',ses,'bin_*.mat'));
 for k=1:numel(d)
  L=load(fullfile(d(k).folder,d(k).name),'R'); R=L.R; F=R.fits;
  m=F(1); m2=pick(F,'lorentz_n2_2d_h20'); n1=pick(F,'lorentz_n1_3d_h20'); n12=pick(F,'lorentz_n1_2d_h20');
  [E0,G,ap,po]=modes(m); [E0b,Gb,apb,pob]=modes(m2); [~,~,ap1,po1]=modes(n1); [~,~,ap12,~]=modes(n12);
  wl=NaN; if m.success, a=m.loss_area; wl=a(1)/sum(a); end
  bp=nan(numel(R.boot),2); ba=bp; bw=nan(numel(R.boot),1);
  for j=1:numel(R.boot)
   [~,~,aa,pp]=modes(R.boot(j)); bp(j,:)=pp; ba(j,:)=aa;
   if R.boot(j).success, a=R.boot(j).loss_area; bw(j)=a(1)/sum(a); end
  end
  hv=F([F.h_variant]); hp=nan(1,2);
  if ~isempty(hv), P=nan(numel(hv),2); for j=1:numel(hv), [~,~,~,P(j,:)]=modes(hv(j)); end, hp=max(abs(P-po),[],1,'omitnan'); end
  dchi=NaN; if n1.success && m.success, dchi=n1.cost-m.cost; end
  dchi2d=NaN; if n12.success && m2.success, dchi2d=n12.cost-m2.cost; end
  rows(end+1,:)={ses,R.q,abs(R.q),R.N,m.q_mean,m.status,any(m.peak_boundary(:)), ...
   E0(1),G(1),ap(1),po(1),E0(2),G(2),ap(2),po(2),wl, ...
   std(bp(:,1),'omitnan'),std(bp(:,2),'omitnan'),std(ba(:,1),'omitnan'),std(ba(:,2),'omitnan'),std(bw,'omitnan'), ...
   E0b(1),Gb(1),pob(1),E0b(2),Gb(2),pob(2),hp(1),hp(2),ap1(1),po1(1),ap12(1),dchi,dchi2d,m.chi2_signal,m.chi2_gain}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'session','q','abs_q','N','q_mean','status','any_boundary', ...
 'E0_low','G_low','apex_low','pole_low','E0_high','G_high','apex_high','pole_high','weight_low', ...
 'boot_pole_low','boot_pole_high','boot_apex_low','boot_apex_high','boot_weight_low', ...
 'E0_low_2d','G_low_2d','pole_low_2d','E0_high_2d','G_high_2d','pole_high_2d','h_sys_pole_low','h_sys_pole_high', ...
 'n1_apex','n1_pole','n1_apex_2d','dchi2_n1_n2','dchi2_n1_n2_2d','chi2_signal','chi2_gain'});
T=sortrows(T,{'session','q'});
writetable(T,fullfile(out,'adaptq_main.csv'));

% fullq main configuration converted the same way, for the overlay.
Fq=readtable(fullfile(fullq_root,'stage5','b1_dispersion_main.csv'));
Fq.pole_high=pole_e(Fq.E0_high,Fq.W_high); Fq.pole_low=pole_e(Fq.E0_low,Fq.W_low);

for s=1:numel(sessions)
 ses=sessions{s}; t=T(strcmp(T.session,ses)&T.abs_q<=0.0300,:); fq=Fq(strcmp(Fq.session,ses)&Fq.abs_q<=0.0300,:);
 ok=strcmp(t.status,'converged')&~t.any_boundary;
 sig=ok&t.dchi2_n1_n2>=11.3;
 f=figure('Visible','off','Position',[0 0 1500 560]); tl=tiledlayout(f,1,3,'TileSpacing','compact');
 % (a) mode energies (pole), main 3D, both signs
 ax=nexttile(tl); hold(ax,'on'); grid(ax,'on');
 for sg=[-1 1]
  u=sig&sign(t.q)==sg; mk='o'; if sg>0, mk='s'; end
  errorbar(ax,t.abs_q(u),t.pole_high(u),t.boot_pole_high(u),mk,'Color',[0.80 0.25 0.10],'MarkerFaceColor',[0.80 0.25 0.10],'MarkerSize',5,'DisplayName',sprintf('upper, q%s0',char('<'*(sg<0)+'>'*(sg>0))));
  errorbar(ax,t.abs_q(u),t.pole_low(u),t.boot_pole_low(u),mk,'Color',[0.10 0.35 0.75],'MarkerFaceColor',[0.10 0.35 0.75],'MarkerSize',5,'DisplayName',sprintf('lower, q%s0',char('<'*(sg<0)+'>'*(sg>0))));
 end
 plot(ax,fq.abs_q,fq.pole_high,'x','Color',[0.5 0.5 0.5],'DisplayName','b1\_fullq (N=3), same conversion');
 plot(ax,fq.abs_q,fq.pole_low,'x','Color',[0.5 0.5 0.5],'HandleVisibility','off');
 xline(ax,0.00625,':','N=1 | N=3','HandleVisibility','off'); xline(ax,0.01425,':','N=5','HandleVisibility','off'); xline(ax,0.02425,':','N=7','HandleVisibility','off');
 xlabel(ax,'|q| (A^{-1}), bin centre'); ylabel(ax,'mode energy Re(pole) (meV)'); title(ax,'(a) two-component fit, 3D kinematic factor');
 legend(ax,'Location','southeast','FontSize',7); xlim(ax,[0 0.031]); ylim(ax,[300 1900]);
 % (b) kinematic-form sensitivity
 ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); u=ok&t.dchi2_n1_n2>=11.3;
 plot(ax,t.abs_q(u),t.pole_high(u),'o','Color',[0.80 0.25 0.10],'MarkerFaceColor',[0.80 0.25 0.10],'DisplayName','upper, 3D');
 plot(ax,t.abs_q(u),t.pole_high_2d(u),'o','Color',[0.80 0.25 0.10],'DisplayName','upper, 2D');
 plot(ax,t.abs_q(u),t.pole_low(u),'o','Color',[0.10 0.35 0.75],'MarkerFaceColor',[0.10 0.35 0.75],'DisplayName','lower, 3D');
 plot(ax,t.abs_q(u),t.pole_low_2d(u),'o','Color',[0.10 0.35 0.75],'DisplayName','lower, 2D');
 xlabel(ax,'|q| (A^{-1})'); ylabel(ax,'Re(pole) (meV)'); title(ax,'(b) 3D vs 2D kinematic factor'); legend(ax,'Location','southeast','FontSize',7);
 xlim(ax,[0 0.031]); ylim(ax,[300 1900]);
 % (c) single-component B1 (n=1) apex and low-branch weight
 ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); yyaxis(ax,'left');
 plot(ax,t.abs_q(ok),t.n1_apex(ok),'k.','MarkerSize',12,'DisplayName','single B1 apex, 3D');
 plot(ax,t.abs_q(ok),t.n1_apex_2d(ok),'ko','MarkerSize',4,'DisplayName','single B1 apex, 2D'); ylabel(ax,'apex (meV)'); ylim(ax,[300 1900]);
 yyaxis(ax,'right'); plot(ax,t.abs_q(sig),t.weight_low(sig),'^','DisplayName','lower-branch weight'); ylabel(ax,'lower-branch share of B1 loss'); ylim(ax,[0 1]);
 xlabel(ax,'|q| (A^{-1})'); title(ax,'(c) one-component B1 and weight transfer'); legend(ax,'Location','north','FontSize',7); xlim(ax,[0 0.031]);
 exportgraphics(f,fullfile(figd,sprintf('dispersion_modes_%s.png',ses)),'Resolution',140); close(f);
 % Fits of selected bins (main configuration): data, model, ZLP, components.
 sel=[0.0010 0.0020 0.0030 0.0045 0.0060 0.0100]; f=figure('Visible','off','Position',[0 0 1500 800]); tl=tiledlayout(f,2,3,'TileSpacing','compact');
 for qa=sel
  d=dir(fullfile(run_root,'stage2',ses,sprintf('bin_%+.4f_N*.mat',qa))); ax=nexttile(tl); hold(ax,'on'); grid(ax,'on');
  if isempty(d), continue; end
  L=load(fullfile(d(1).folder,d(1).name),'R'); c=L.R.fits(1).curves; if isempty(fieldnames(c)), continue; end
  r=c.E>=150&c.E<=2000; plot(ax,c.E(r),c.observed(r),'Color',[0.6 0.6 0.6]); plot(ax,c.E(r),c.prediction(r),'k','LineWidth',1.2);
  plot(ax,c.E(r),c.zlp(r),'--','Color',[0.4 0.4 0.4]); cl=[0.10 0.35 0.75; 0.80 0.25 0.10];
  for j=1:size(c.peaks,2), plot(ax,c.E(r),c.peaks(r,j),'Color',cl(min(j,2),:),'LineWidth',1.2); end
  p=L.R.fits(1).parameters; title(ax,sprintf('q=%+.4f (N=%d): E0 %s, \\Gamma %s',L.R.q,L.R.N,mat2str(round(p(:,1)')),mat2str(round(p(:,2)'))),'FontSize',9);
  xlabel(ax,'meV'); ylabel(ax,'counts'); yl=max(c.observed(c.E>=400&c.E<=2000)); ylim(ax,[0 1.6*yl]); xlim(ax,[150 2000]);
 end
 exportgraphics(f,fullfile(figd,sprintf('fits_selected_%s.png',ses)),'Resolution',130); close(f);
end
end

function f=pick(F,key)
i=find(strcmp({F.key},key),1); if isempty(i), f=struct('success',false,'parameters',nan(1,3),'cost',NaN); else, f=F(i); end
end

function [E0,G,ap,po]=modes(f)
n=2; E0=nan(1,n); G=E0; ap=E0; po=E0;
if ~isfield(f,'success') || ~f.success, return; end
P=f.parameters; k=size(P,1); E0(1:k)=P(:,1)'; G(1:k)=P(:,2)';
for i=1:k, ap(i)=apex_e(E0(i),G(i)); po(i)=pole_e(E0(i),G(i)); end
end

function a=apex_e(E0,G)
c=2*E0.^2-G.^2; a=sqrt((c+sqrt(c.^2+12*E0.^4))/6);
end

function p=pole_e(E0,G)
d=E0.^2-G.^2/4; p=sqrt(max(d,0)); p(d<=0)=NaN;
end
