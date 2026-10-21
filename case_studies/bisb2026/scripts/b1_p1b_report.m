function b1_p1b_report(outdir, diagdir)
% Tables and figures for a P1b run folder (b1_p1b_lowq). All models (baseline,
% P1, P1b KA/KB with floor 250 and floor 0) are scored with the P1b noise model
% (common-beam reference, first differences) on the same windows; component
% shares are reported per component; quadrature convergence of the refined
% quadrature is checked against a still finer one.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_p1b_config(); sw=cfg.score_windows;
d=dir(fullfile(outdir,'q_*.mat')); Rs={}; for k=1:numel(d), L=load(fullfile(outdir,d(k).name),'R'); Rs{end+1}=L.R; end %#ok<AGROW>
q=cellfun(@(r) r.q,Rs); [~,o]=sort(abs(q)+1e-9*sign(q)); Rs=Rs(o); q=q(o);
hi2=struct('n_sub',45,'n_probe',61,'probe_extent',5.5,'q_round',1e-8);
Sc={}; Pa={}; Sh={}; Kt={}; Qc={};
for i=1:numel(Rs)
 R=Rs{i}; E=R.E; Y=R.Y;
 M=models(R);
 for m=1:numel(M)
  for sg={'sigma','sigma_channel'}
   s=R.(sg{1}); r=(Y-M(m).pred)./s;
   for w=1:numel(sw)
    k=E>=sw(w).range(1)&E<=sw(w).range(2);
    Sc(end+1,:)={R.q,M(m).name,sg{1},sw(w).name,nnz(k),mean(r(k).^2),sqrt(mean(r(k).^2)),mean(r(k)),M(m).in_data(w)}; %#ok<AGROW>
   end
  end
  for e=cfg.share_energies
   [~,ie]=min(abs(E-e)); ii=max(1,ie-2):min(numel(E),ie+2);
   Sh(end+1,:)={R.q,M(m).name,e,sum(M(m).zlp(ii))/sum(Y(ii)),sum(sum(M(m).aux(ii,:),2))/sum(Y(ii)),sum(sum(M(m).b1(ii,:),2))/sum(Y(ii))}; %#ok<AGROW>
  end
  p=M(m).par; ap=arrayfun(@(j) apex_e(p(j,1),p(j,2)),1:2);
  Pa(end+1,:)={R.q,M(m).name,M(m).status,p(1,1),p(1,2),ap(1),p(2,1),p(2,2),ap(2),M(m).cA,M(m).cB,M(m).auxp(1,1),M(m).auxp(1,2),M(m).auxp(2,1),M(m).auxp(2,2),M(m).auxbd}; %#ok<AGROW>
 end
 for kv=1:numel(R.fits)
  F=R.fits(kv); a=F.f250;
  Kn=a.kernel.K/sum(a.kernel.K); lg=a.kernel.lag_meV; pk=max(Kn);
  hm=lg(find(Kn>=pk/2,1,'first')); hp=lg(find(Kn>=pk/2,1,'last'));
  Kt(end+1,:)={R.q,F.kernel,hp-hm,sum(Kn(abs(lg-a.zlp_parameters.center_meV)>60)),a.stageA.zlp_parameters.center_meV,a.zlp_parameters.center_meV,a.stageA.chi2_red_data, ...
   F.forward.floor_max_sigma,F.forward.quad_max_sigma_floor250,F.forward.quad_max_sigma_floor0,F.f0_quadrature,F.forward.rebuild_check}; %#ok<AGROW>
  % finer quadrature check of the one used by the floor-0 fit
  b=F.f0; if ~b.success, continue; end
  ql=cfg.quadrature_lo; qh=cfg.quadrature_hi; if strcmp(F.f0_quadrature,'lo'), qa=ql; qb=qh; else, qa=qh; qb=hi2; end
  Ka=qe_kinematic_prefactor(R.members,dq=R.dq,form=cfg.kinematic_form,h_perp=cfg.h_perp,sigma_probe=R.sigma_probe,beam_kV=cfg.beam_kV,n_sub=qa.n_sub,n_probe=qa.n_probe,probe_extent=qa.probe_extent,q_round=qa.q_round);
  Kb=qe_kinematic_prefactor(R.members,dq=R.dq,form=cfg.kinematic_form,h_perp=cfg.h_perp,sigma_probe=R.sigma_probe,beam_kV=cfg.beam_kV,n_sub=qb.n_sub,n_probe=qb.n_probe,probe_extent=qb.probe_extent,q_round=qb.q_round);
  [ya,yb]=forward_pair(b,Ka,Kb); s=R.sigma; chk=E>=-180&E<=1800;
  Qc(end+1,:)={R.q,F.kernel,F.f0_quadrature,max(abs(yb(chk)-ya(chk))./s(chk))}; %#ok<AGROW>
 end
end
St=cell2table(Sc,'VariableNames',{'q','model','sigma_model','window','n','mean_r2','rms_r','mean_r','in_data_term'});
Pt=cell2table(Pa,'VariableNames',{'q','model','status','E0_low','G_low','apex_low','E0_high','G_high','apex_high','c_stageA','c_stageB','aux1_E0','aux1_G','aux2_E0','aux2_G','aux_boundary'});
Ht=cell2table(Sh,'VariableNames',{'q','model','E_meV','zlp_share','aux_share','b1_share'});
Kk=cell2table(Kt,'VariableNames',{'q','kernel','kernel_fwhm_meV','kernel_weight_beyond_60meV','c_stageA','c_stageB','stageA_chi2_gain','floor_fwd_max_sigma','quad_max_sigma_floor250','quad_max_sigma_floor0','f0_quadrature','rebuild_rel_err'});
Qt=cell2table(Qc,'VariableNames',{'q','kernel','f0_quadrature','max_change_vs_finer_sigma'});
writetable(St,fullfile(outdir,'p1b_scores.csv')); writetable(Pt,fullfile(outdir,'p1b_params.csv')); writetable(Ht,fullfile(outdir,'p1b_shares.csv'));
writetable(Kk,fullfile(outdir,'p1b_kernels_forward.csv')); writetable(Qt,fullfile(outdir,'p1b_quadrature_check.csv'));
pt={}; for i=1:numel(Rs), if ~isempty(Rs{i}.perturbation), p=Rs{i}.perturbation; pt(end+1,:)={Rs{i}.q,p.sigma_outside_identical,p.parameters_identical,p.aux_identical,p.centre_identical,p.stageA_identical,p.kernel_identical,p.max_abs_dparam}; end, end %#ok<AGROW>
writetable(cell2table(pt,'VariableNames',{'q','sigma_outside','parameters','aux','centre','stageA','kernel','max_abs_dparam'}),fullfile(outdir,'p1b_perturbation.csv'));
figs(Rs,q,cfg,outdir,St,diagdir);
fprintf('wrote P1b tables and figures to %s\n',outdir);
end

function M=models(R)
% Prediction and components of every model on R.E; in_data per score window.
nw=6; allin=true(1,nw); p1in=[true true false true true true];
M=struct('name',{},'pred',{},'zlp',{},'aux',{},'b1',{},'par',{},'status',{},'cA',{},'cB',{},'auxp',{},'auxbd',{},'in_data',{});
M(end+1)=struct('name','baseline','pred',R.baseline.prediction,'zlp',R.baseline.zlp,'aux',R.baseline.aux_peaks,'b1',R.baseline.peaks, ...
 'par',R.baseline.parameters,'status','','cA',NaN,'cB',NaN,'auxp',nan(2,2),'auxbd',NaN,'in_data',allin);
M(end+1)=struct('name','P1','pred',R.p1.prediction,'zlp',R.p1.zlp,'aux',R.p1.aux_peaks,'b1',R.p1.peaks, ...
 'par',R.p1.parameters,'status','','cA',NaN,'cB',NaN,'auxp',nan(2,2),'auxbd',NaN,'in_data',p1in);
for kv=1:numel(R.fits)
 for fl={'f250','f0'}
  f=R.fits(kv).(fl{1}); if ~f.success, continue; end
  c=f.curves;
  M(end+1)=struct('name',sprintf('%s_%s',R.fits(kv).kernel,fl{1}),'pred',c.prediction,'zlp',c.zlp,'aux',c.aux_peaks,'b1',c.peaks, ...
   'par',f.parameters,'status',f.status,'cA',f.stageA.zlp_parameters.center_meV,'cB',f.zlp_parameters.center_meV, ...
   'auxp',f.aux_parameters(:,1:2),'auxbd',any(f.aux_boundary(:)),'in_data',p1in); %#ok<AGROW>
 end
end
end

function [ya,yb]=forward_pair(f,Ka,Kb)
% B1 contribution at the fitted parameters (floor 0) with two quadratures.
E=f.curves.E; dE=median(diff(E)); lag=f.kernel.lag_meV; Kn=f.kernel.K/sum(f.kernel.K);
Ex=(E(1)-lag(end))+(0:(numel(E)+numel(lag)-2)).'*dE; Ea=abs(Ex);
m=peak_models('lorentz'); kT=0.08617333262*300; L=zeros(size(Ex));
for j=1:size(f.parameters,1)
 l=f.parameters(j,3)*m.model_fn(f.parameters(j,1),f.parameters(j,2),1,Ea); l(Ex<0)=l(Ex<0).*exp(-Ea(Ex<0)/kT); l(~isfinite(l))=0; L=L+l;
end
ya=conv(L.*Ka(Ea),Kn,'valid'); yb=conv(L.*Kb(Ea),Kn,'valid');
end

function figs(Rs,q,cfg,outdir,St,diagdir) %#ok<INUSD>
ho=cfg.holdout_windows; cB=[0.10 0.35 0.75; 0.80 0.25 0.10]; cA=[0.20 0.60 0.20; 0.55 0.80 0.30];
for sgn=[1 -1]
 ii=find(sign(q)==sgn); tag=char('m'*(sgn<0)+'p'*(sgn>0));
 f=figure('Visible','off','Position',[0 0 1900 1150]); tl=tiledlayout(f,3,numel(ii),'TileSpacing','compact','Padding','compact');
 for row=1:3
  for i=ii
   R=Rs{i}; E=R.E; Y=R.Y; kb=R.fits(2); a=kb.f250.curves; z=kb.f0.curves; ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
   if row<3
    xregion(ax,ho(1),ho(2),'FaceColor',[0.85 0.85 0.85],'HandleVisibility','off');
    plot(ax,E,Y,'Color',[0.6 0.6 0.6],'DisplayName','data');
    plot(ax,E,R.p1.prediction,'--','Color',[0.4 0.4 0.4],'LineWidth',1.1,'DisplayName','P1 total');
    plot(ax,E,a.prediction,'m','LineWidth',1.3,'DisplayName','P1b KB, floor 250: total');
    plot(ax,E,z.prediction,'c--','LineWidth',1.2,'DisplayName','P1b KB, floor 0: total');
    plot(ax,E,a.zlp,'k:','LineWidth',1.1,'DisplayName','KB f250: elastic a*R');
    for j=1:size(a.aux_peaks,2), plot(ax,E,a.aux_peaks(:,j),'Color',cA(j,:),'LineWidth',1.1,'DisplayName',sprintf('KB f250: aux %d',j)); end
    for j=1:2, plot(ax,E,a.peaks(:,j),'Color',cB(j,:),'LineWidth',1.2,'DisplayName',sprintf('KB f250: B1 %d (x I_{kin})',j)); end
   end
   if row==1
    xlim(ax,[0 1000]); w=E>=60&E<=1000; ylim(ax,[0 1.25*max(Y(w))]);
    p=kb.f250.parameters; title(ax,sprintf('q=%+.4f  KB f250: E0 %.0f/%.0f, \\Gamma %.0f/%.0f',R.q,p(1,1),p(2,1),p(1,2),p(2,2)),'FontSize',9);
    if i==ii(1), legend(ax,'Location','northeast','FontSize',7); end
   elseif row==2
    set(ax,'YScale','log'); xlim(ax,[-200 1800]); ylim(ax,[max(1,0.3*min(Y(Y>0))) 2*max(Y)]); title(ax,'log; grey band = held out','FontSize',9);
   else
    xregion(ax,ho(1),ho(2),'FaceColor',[0.85 0.85 0.85]); s=R.sigma;
    plot(ax,E,movmean((Y-R.p1.prediction)./s,25),'Color',[0.55 0.55 0.55],'LineWidth',1.1);
    plot(ax,E,movmean((Y-a.prediction)./s,25),'k','LineWidth',1.3);
    plot(ax,E,movmean((Y-z.prediction)./s,25),'c','LineWidth',1.1);
    h=E>=ho(1)&E<=ho(2); plot(ax,E(h),movmean((Y(h)-a.prediction(h))./s(h),5),'r','LineWidth',1.5);
    yline(ax,0,'b:'); xlim(ax,[-200 1800]); ylim(ax,[-12 12]);
    sq=St(abs(St.q-R.q)<1e-9 & strcmp(St.sigma_model,'sigma'),:); g=@(m,w) sq.rms_r(strcmp(sq.model,m)&strcmp(sq.window,w));
    title(ax,sprintf('25-pt mean (obs-model)/\\sigma (P1b noise); grey P1, black KB f250, cyan KB f0, red KB f250 held out\nRMS 300-1800: P1 %.2f, KB f250 %.2f, f0 %.2f; held out: %.2f, %.2f, %.2f', ...
     g('P1','b1_300_1800'),g('KB_f250','b1_300_1800'),g('KB_f0','b1_300_1800'),g('P1','holdout_150_250'),g('KB_f250','holdout_150_250'),g('KB_f0','holdout_150_250')),'FontSize',8);
   end
   xlabel(ax,'energy loss (meV)');
  end
 end
 title(tl,sprintf('590 single channels, q %s 0: P1 vs P1b (one response for elastic and inelastic, kernel KB), counts space',char('<'*(sgn<0)+'>'*(sgn>0))));
 exportgraphics(f,fullfile(outdir,sprintf('p1b_counts_%s.png',tag)),'Resolution',115); close(f);
end
% floor and quadrature, KB, fixed parameters
f=figure('Visible','off','Position',[0 0 1700 800]); tl=tiledlayout(f,2,4,'TileSpacing','compact');
for i=1:numel(Rs)
 R=Rs{i}; fw=R.fits(2).forward; s=R.sigma; ax=nexttile(tl); hold(ax,'on'); grid(ax,'on');
 plot(ax,fw.E,fw.d_floor./s,'b','DisplayName','floor 0 - floor 250'); plot(ax,fw.E,10*fw.d_quad0./s,'r','DisplayName','10 x (refined - default quadrature), floor 0');
 xregion(ax,ho(1),ho(2),'FaceColor',[0.85 0.85 0.85],'HandleVisibility','off'); xlim(ax,[-200 1800]);
 title(ax,sprintf('q=%+.4f KB: B1 change / \\sigma at fixed parameters',R.q),'FontSize',9); if i==1, legend(ax,'FontSize',7,'Location','northeast'); end
end
exportgraphics(f,fullfile(outdir,'p1b_floor_quadrature.png'),'Resolution',115); close(f);
% kernels KA vs KB, centred
f=figure('Visible','off','Position',[0 0 1700 800]); tl=tiledlayout(f,2,4,'TileSpacing','compact');
for i=1:numel(Rs)
 R=Rs{i}; ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); set(ax,'YScale','log');
 for kv=1:numel(R.fits)
  a=R.fits(kv).f250; Kn=a.kernel.K/max(a.kernel.K); plot(ax,a.kernel.lag_meV-a.zlp_parameters.center_meV,Kn,'LineWidth',1.2,'DisplayName',sprintf('%s (FWHM from table)',R.fits(kv).kernel));
 end
 E=R.E; Y=R.Y; [~,i0]=max(Y.*(abs(E)<40)); plot(ax,E-E(i0),Y/Y(i0),'Color',[0.6 0.6 0.6],'DisplayName','data, normalised');
 xlim(ax,[-200 300]); ylim(ax,[1e-6 2]); title(ax,sprintf('q=%+.4f: kernels KA, KB vs data',R.q),'FontSize',9); if i==1, legend(ax,'FontSize',7,'Location','northeast'); end
end
exportgraphics(f,fullfile(outdir,'p1b_kernels.png'),'Resolution',115); close(f);
% parameters
Pt=readtable(fullfile(outdir,'p1b_params.csv'));
f=figure('Visible','off','Position',[0 0 1500 450]); tl=tiledlayout(f,1,3,'TileSpacing','compact');
mods={'baseline','P1','KA_f250','KB_f250','KA_f0','KB_f0'}; mk={'s','d','o','o','^','^'}; fillm=[0 0 0 1 0 1];
flds={'E0','G','apex'}; ttl={'E0 (meV)','\Gamma (meV)','DL apex (meV)'};
for k=1:3
 ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
 for mm=1:numel(mods)
  u=strcmp(Pt.model,mods{mm}); x=abs(Pt.q(u)).*(1+0.03*(mm-3.5)/3).*(1+0.02*sign(Pt.q(u)));
  for br={'low','high'}
   col=cB(1+strcmp(br{1},'high'),:); y=Pt.([flds{k} '_' br{1}])(u);
   h=plot(ax,x,y,mk{mm},'Color',col,'DisplayName',sprintf('%s %s',mods{mm},br{1}));
   if fillm(mm), h.MarkerFaceColor=col; end
  end
 end
 set(ax,'XScale','log'); xlim(ax,[4e-4 6e-3]); xlabel(ax,'|q| (A^{-1})'); ylabel(ax,ttl{k});
 if k==1, legend(ax,'Location','eastoutside','FontSize',7); end
end
title(tl,'590: parameters of baseline, P1 and P1b (KA/KB, floor 250 / 0); symbols offset slightly in |q| for legibility');
exportgraphics(f,fullfile(outdir,'p1b_params.png'),'Resolution',130); close(f);
end

function a=apex_e(E0,G)
c=2*E0.^2-G.^2; a=sqrt((c+sqrt(c.^2+12*E0.^4))/6);
end
