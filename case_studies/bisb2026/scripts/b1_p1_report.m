function b1_p1_report(outdir)
% Tables and figures for a P1 run folder (b1_p1_lowq): baseline vs P1 on the
% same data and the same scoring windows, Stage A (old: ZLP only on gain side
% + core; P1: ZLP + aux on gain side, core and loss side <= 150 meV),
% parameters with bootstrap spread, retardation flags.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
meta=jsondecode(fileread(fullfile(outdir,'config.json'))); cfg=meta.cfg; sw=cfg.score_windows;
d=dir(fullfile(outdir,'q_*.mat')); Rs={}; for k=1:numel(d), L=load(fullfile(outdir,d(k).name),'R'); Rs{end+1}=L.R; end %#ok<AGROW>
q=cellfun(@(r) r.q,Rs); [~,o]=sort(abs(q)+1e-9*sign(q)); Rs=Rs(o); q=q(o);
hc=1.97327e6;                                  % hbar*c in meV*A: q_gamma = E/hc
P=[]; Sc=[]; SA=[];
for i=1:numel(Rs)
 R=Rs{i}; E=R.E; Y=R.Y; s=R.sigma;
 models={'baseline',R.baseline.curves.prediction; 'P1',R.p1.curves.prediction};
 for m=1:2
  pr=models{m,2}; r=(Y-pr)./s;
  for w=1:numel(sw)
   k=E>=sw(w).range(1)&E<=sw(w).range(2);
   used=true; if m==2, used=any(R.p1.data_mask(k)); end
   Sc=[Sc; {R.q,models{m,1},sw(w).name,nnz(k),mean(r(k).^2),mean(r(k)),used}]; %#ok<AGROW>
  end
 end
 % parameters
 for m=1:2
  if m==1, pa=R.baseline.parameters; ax=R.baseline.aux_parameters; lf=R.baseline.curves.loss_function; Ec=R.baseline.curves.E; st=R.baseline.status; bd=NaN; sd=nan(2,2); asd=nan(1,2);
  else, pa=R.p1.parameters; ax=R.p1.aux_parameters; lf=R.p1.curves.loss_function; Ec=R.p1.curves.E; st=R.p1.status; bd=any(R.p1.peak_boundary(:))+2*any(R.p1.aux_boundary(:));
   sd=nan(2,2); asd=nan(1,2);
   if ~isempty(R.boot), Bp=cat(3,R.boot.parameters); sd=std(Bp(:,1:2,:),0,3); Ba=cat(3,R.boot.aux_parameters); asd=squeeze(std(Ba(:,1,:),0,3)).'; end
  end
  ep=Ec>0&Ec<=1800; ar=trapz(Ec(ep),lf(ep,:)); ap=arrayfun(@(j) apex_e(pa(j,1),pa(j,2)),1:2);
  rt=pa(:,1).'/hc/abs(R.q); rtm=pa(:,1).'/hc/R.q_mean;
  mn={'baseline','P1'};
  P=[P; {R.q,R.q_mean,mn{m},st,bd,pa(1,1),sd(1,1),pa(1,2),sd(1,2),ap(1),pa(2,1),sd(2,1),pa(2,2),sd(2,2),ap(2),ar(1)/sum(ar), ...
   ax(1,1),ax(1,2),ax(2,1),ax(2,2),asd(1),asd(2),rt(1),rt(2),rtm(1),rtm(2),max(rt)>cfg.retardation_flag_ratio}]; %#ok<AGROW>
 end
 % Stage A comparison
 A={R.stageA_old,'old: ZLP only, gain+core'; R.p1.stageA,'P1: ZLP+aux, gain+core+loss<=150'};
 for m=1:2
  a=A{m,1}; z=a.zlp_parameters; ex=@(e0) mean(a.prediction(abs(E-e0)<=10))/mean(Y(abs(E-e0)<=10));
  axp=nan(2,2); if ~isempty(a.aux_parameters), axp=a.aux_parameters(:,1:2); end
  SA=[SA; {R.q,A{m,2},a.chi2_red_data,z.center_meV,z.sL(1),z.sR(1),z.mL(1),z.mR(1),z.sL(2),z.sR(2),z.mL(2),z.mR(2),z.fwhm_meV, ...
   ex(100),ex(200),ex(300),ex(400),axp(1,1),axp(1,2),axp(2,1),axp(2,2)}]; %#ok<AGROW>
 end
end
Pt=cell2table(P,'VariableNames',{'q','q_mean','model','status','boundary_code','E0_low','E0_low_sd','G_low','G_low_sd','apex_low', ...
 'E0_high','E0_high_sd','G_high','G_high_sd','apex_high','low_share_in_window','aux1_E0','aux1_G','aux2_E0','aux2_G','aux1_E0_sd','aux2_E0_sd', ...
 'qgamma_over_q_low','qgamma_over_q_high','qgamma_over_qmean_low','qgamma_over_qmean_high','retardation_to_check'});
St=cell2table(Sc,'VariableNames',{'q','model','window','n','mean_r2','mean_r','in_data_term'});
At=cell2table(SA,'VariableNames',{'q','stageA','chi2_gain','c','sL1','sR1','mL1','mR1','sL2','sR2','mL2','mR2','zlp_fwhm', ...
 'pred_over_data_100','pred_over_data_200','pred_over_data_300','pred_over_data_400','aux1_E0','aux1_G','aux2_E0','aux2_G'});
writetable(Pt,fullfile(outdir,'p1_params.csv')); writetable(St,fullfile(outdir,'p1_scores.csv')); writetable(At,fullfile(outdir,'p1_stageA.csv'));
figs(Rs,q,cfg,outdir,St);
fprintf('wrote tables and figures to %s\n',outdir);
end

function figs(Rs,q,cfg,outdir,St)
ho=cfg.holdout_windows; cB=[0.10 0.35 0.75; 0.80 0.25 0.10]; cA=[0.20 0.60 0.20; 0.55 0.80 0.30];
for sgn=[1 -1]
 ii=find(sign(q)==sgn); tag=char('m'*(sgn<0)+'p'*(sgn>0));
 % ---- F1: counts, baseline vs P1
 f=figure('Visible','off','Position',[0 0 1900 1100]); tl=tiledlayout(f,3,numel(ii),'TileSpacing','compact','Padding','compact');
 for row=1:3
  for i=ii
   R=Rs{i}; E=R.E; Y=R.Y; c=R.p1.curves; b=R.baseline.curves; ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
   if row<3
    xregion(ax,ho(1),ho(2),'FaceColor',[0.85 0.85 0.85],'HandleVisibility','off');
    plot(ax,E,Y,'Color',[0.6 0.6 0.6],'DisplayName','data');
    plot(ax,E,b.prediction,'k--','LineWidth',1.1,'DisplayName','baseline total');
    plot(ax,E,c.prediction,'m','LineWidth',1.3,'DisplayName','P1 total');
    plot(ax,E,c.zlp,'k:','LineWidth',1.1,'DisplayName','P1 ZLP');
    for j=1:size(c.aux_peaks,2), plot(ax,E,c.aux_peaks(:,j),'Color',cA(j,:),'LineWidth',1.2,'DisplayName',sprintf('P1 aux %d',j)); end
    for j=1:2, plot(ax,E,c.peaks(:,j),'Color',cB(j,:),'LineWidth',1.3,'DisplayName',sprintf('P1 B1 %d (x I_{kin})',j)); end
   end
   if row==1
    xlim(ax,[0 1000]); w=E>=60&E<=1000; ylim(ax,[0 1.25*max(Y(w))]);
    p=R.p1.parameters; title(ax,sprintf('q=%+.4f  P1 E0 %.0f/%.0f, \\Gamma %.0f/%.0f',R.q,p(1,1),p(2,1),p(1,2),p(2,2)),'FontSize',9);
    if i==ii(1), legend(ax,'Location','northeast','FontSize',7); end
   elseif row==2
    set(ax,'YScale','log'); xlim(ax,[-200 1800]); ylim(ax,[max(1,0.3*min(Y(Y>0))) 2*max(Y)]); title(ax,'log; grey band = held out (150-250 meV)','FontSize',9);
   else
    xregion(ax,ho(1),ho(2),'FaceColor',[0.85 0.85 0.85]);
    rb=(Y-b.prediction)./R.sigma; rp=(Y-c.prediction)./R.sigma;
    plot(ax,E,movmean(rb,25),'Color',[0.55 0.55 0.55],'LineWidth',1.2); plot(ax,E,movmean(rp,25),'k','LineWidth',1.3);
    h=E>=ho(1)&E<=ho(2); plot(ax,E(h),movmean(rp(h),5),'r','LineWidth',1.5);
    yline(ax,0,'b:'); xlim(ax,[-200 1800]); ylim(ax,[-6 6]);
    s=St(abs(St.q-R.q)<1e-9,:); g=@(m,w) s.mean_r2(strcmp(s.model,m)&strcmp(s.window,w));
    title(ax,sprintf('25-pt mean of (obs-model)/\\sigma; grey baseline, black P1, red P1 in held-out band\n<r^2> 250-1800: %.2f -> %.2f; held-out: %.2f -> %.2f; 30-150: %.2f -> %.2f', ...
     g('baseline','b1_250_1800'),g('P1','b1_250_1800'),g('baseline','holdout_150_250'),g('P1','holdout_150_250'),g('baseline','loss_30_150'),g('P1','loss_30_150')),'FontSize',8);
   end
   xlabel(ax,'energy loss (meV)');
  end
 end
 title(tl,sprintf('590 single channels, q %s 0: baseline (stage2\\_2d\\_h5) vs P1, same data; counts space',char('<'*(sgn<0)+'>'*(sgn>0))));
 exportgraphics(f,fullfile(outdir,sprintf('p1_counts_%s.png',tag)),'Resolution',115); close(f);
 % ---- F2: Stage A, old vs P1
 f=figure('Visible','off','Position',[0 0 1900 800]); tl=tiledlayout(f,2,numel(ii),'TileSpacing','compact','Padding','compact');
 for row=1:2
  for i=ii
   R=Rs{i}; E=R.E; Y=R.Y; a0=R.stageA_old; a1=R.p1.stageA; ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
   xregion(ax,ho(1),ho(2),'FaceColor',[0.85 0.85 0.85],'HandleVisibility','off');
   if row==1
    plot(ax,E,Y,'Color',[0.6 0.6 0.6],'DisplayName','data');
    plot(ax,E,a0.prediction,'b--','LineWidth',1.2,'DisplayName','old Stage A (ZLP only)');
    plot(ax,E,a1.zlp,'k:','LineWidth',1.2,'DisplayName','P1 Stage A ZLP');
    for j=1:size(a1.aux_peaks,2), plot(ax,E,a1.aux_peaks(:,j),'Color',cA(j,:),'LineWidth',1.2,'DisplayName',sprintf('P1 Stage A aux %d (E0 %.0f, \\Gamma %.0f)',j,a1.aux_parameters(j,1),a1.aux_parameters(j,2))); end
    plot(ax,E,a1.prediction,'m','LineWidth',1.2,'DisplayName','P1 Stage A total');
    set(ax,'YScale','log'); xlim(ax,[-190 450]); ylim(ax,[max(1,0.3*min(Y(E<=450&Y>0))) 2*max(Y)]);
    title(ax,sprintf('q=%+.4f Stage A; \\chi^2 gain side: old %.2f, P1 %.2f',R.q,a0.chi2_red_data,a1.chi2_red_data),'FontSize',9);
    legend(ax,'Location','northeast','FontSize',6);
   else
    r0=(Y-a0.prediction)./R.sigma; r1=(Y-a1.prediction)./R.sigma; d0=a0.data_mask; d1=a1.data_mask;
    plot(ax,E(d0),movmean(r0(d0),9),'b','LineWidth',1.1,'DisplayName','old, its data');
    plot(ax,E(d1),movmean(r1(d1),9),'m','LineWidth',1.1,'DisplayName','P1, its data');
    k=E>150&E<=450; plot(ax,E(k),movmean(r1(k),9),'Color',[0.9 0.6 0.9],'DisplayName','P1 Stage A beyond its data (B1 not modelled)');
    yline(ax,0,'k:','HandleVisibility','off'); xlim(ax,[-190 450]); ylim(ax,[-8 8]);
    title(ax,'9-pt mean of (obs - Stage A)/\sigma','FontSize',9); legend(ax,'Location','southwest','FontSize',6);
   end
   xlabel(ax,'energy (meV)');
  end
 end
 title(tl,sprintf('Stage A prefit, q %s 0: old (ZLP only, gain side + core) vs P1 (ZLP + aux, gain side + core + loss side <= 150 meV, 150-250 held out)',char('<'*(sgn<0)+'>'*(sgn>0))));
 exportgraphics(f,fullfile(outdir,sprintf('p1_stageA_%s.png',tag)),'Resolution',115); close(f);
end
% ---- F3: parameters vs |q|
Pt=readtable(fullfile(outdir,'p1_params.csv'));
f=figure('Visible','off','Position',[0 0 1500 450]); tl=tiledlayout(f,1,3,'TileSpacing','compact');
flds={'E0','G','apex'}; ttl={'E0 (meV)','\Gamma (meV)','DL apex (meV)'};
for k=1:3
 ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
 for br={'low','high'}
  b=br{1}; col=cB(1+strcmp(b,'high'),:);
  for m={'baseline','P1'}
   u=strcmp(Pt.model,m{1}); x=abs(Pt.q(u)).*(1+0.04*sign(Pt.q(u))); y=Pt.([flds{k} '_' b])(u);
   if strcmp(m{1},'P1')
    e=zeros(size(y)); if k<3, e=Pt.([flds{k} '_' b '_sd'])(u); end
    errorbar(ax,x,y,e,'o','Color',col,'MarkerFaceColor',col,'DisplayName',sprintf('P1 %s',b));
   else
    plot(ax,x,y,'s','Color',col,'DisplayName',sprintf('baseline %s',b));
   end
  end
 end
 set(ax,'XScale','log'); xlim(ax,[4e-4 6e-3]); xlabel(ax,'|q| (A^{-1}); +q drawn slightly right'); ylabel(ax,ttl{k});
 if k==1, legend(ax,'Location','best','FontSize',8); end
end
title(tl,'590: baseline vs P1 parameters (P1 bars: frame bootstrap SD)');
exportgraphics(f,fullfile(outdir,'p1_params.png'),'Resolution',130); close(f);
end

function a=apex_e(E0,G)
c=2*E0.^2-G.^2; a=sqrt((c+sqrt(c.^2+12*E0.^4))/6);
end
