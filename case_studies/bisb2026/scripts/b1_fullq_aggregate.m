function b1_fullq_aggregate(run_root)
% B1 full-q task, stages 5-6: collect the per-bin fits, build the dispersion
% table (main = DL n=2 with the main aux/kinematic/h settings; statistical
% error = frame-bootstrap SD, systematic = largest shift among the 3D form,
% the other aux layout and the h variants), fit the two branches and draw the
% figures. Writes into <run_root>/stage5.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_fullq_config();
out=fullfile(run_root,'stage5'); figd=fullfile(out,'figures'); if ~isfolder(figd), mkdir(figd); end
sessions={'590_PL2_10w','n0_PL2_10w_repeat'};
ALL={}; DISP={}; S=struct();
for si=1:numel(sessions)
 ses=sessions{si}; d=fullfile(run_root,'stage2',ses); files=dir(fullfile(d,'bin_*.mat'));
 if isempty(files), continue; end
 Rs=cell(numel(files),1); for k=1:numel(files), L=load(fullfile(d,files(k).name),'R'); Rs{k}=L.R; end
 qv=cellfun(@(r)r.q,Rs); [qv,o]=sort(qv); Rs=Rs(o);
 F=load(fullfile(run_root,'stage1',sprintf('frames_%s.mat',ses)),'D'); D=F.D;
 rho=median(D.noise_rho_loss(1,:),'omitnan'); if ~isfinite(rho), rho=0; end
 corr=(1-rho)/(1+rho);
 for k=1:numel(Rs)
  R=Rs{k};
  for j=1:numel(R.fits), ALL(end+1,:)=fit_row(ses,R,R.fits(j)); end %#ok<AGROW>
  [m,n1,vars,ndeg]=pick(R);
  st=boot_sd(R); sy=sys_shift(m,vars);
  sb=synth_summary(run_root,ses,R.q);
  dchi=NaN; if ~isempty(n1) && n1.success && m.success, dchi=n1.cost-m.cost; end
  w=m.loss_area(1)/sum(m.loss_area); wf=m.loss_fsum(1)/sum(m.loss_fsum);
  DISP(end+1,:)={ses,R.q,abs(R.q),sign(R.q),m.q_mean,m.status,any(m.peak_boundary(:)), ...
   m.parameters(1,1),m.parameters(2,1),m.parameters(1,2),m.parameters(2,2),w,wf, ...
   st(1),st(2),st(3),st(4),st(5),sy(1),sy(2),sy(3),sy(4),sy(5), ...
   m.chi2_signal,m.chi2_gain,m.chi2_low,m.background_fraction,dchi,dchi*corr,sb(1),sb(2),sb(3),sb(4),ndeg}; %#ok<AGROW>
 end
 S.(matlab.lang.makeValidName(ses))=struct('D',D,'Rs',{Rs},'name',ses);
end
VA={'session','q','abs_q','sign','q_mean','status','any_boundary','E0_low','E0_high','W_low','W_high','weight_low', ...
 'fsum_low','stat_E0_low','stat_E0_high','stat_W_low','stat_W_high','stat_weight_low','sys_E0_low','sys_E0_high', ...
 'sys_W_low','sys_W_high','sys_weight_low','chi2_signal','chi2_gain','chi2_low','background_fraction', ...
 'dchi2_n1_n2','dchi2_n1_n2_corr','synth_bias_E0_low','synth_bias_E0_high','synth_sd_E0_low','synth_sd_E0_high', ...
 'n_degenerate_variants'};
T=cell2table(DISP,'VariableNames',VA); writetable(T,fullfile(out,'b1_dispersion_main.csv'));
A=cell2table(ALL,'VariableNames',{'session','q','key','model','n','aux','form','h','main','status','cost', ...
 'chi2_signal','chi2_gain','chi2_low','E0_1','W_1','A_1','E0_2','W_2','A_2','area_1','area_2','background_fraction'});
writetable(A,fullfile(out,'b1_all_fits.csv'));

% Dispersion fits on usable bins (main fit converged, no peak at a bound).
use=strcmp(T.status,'converged')&~T.any_boundary;
eT=@(a,b)sqrt(a.^2+b.^2);
fits=struct([]);
for si=0:numel(sessions)
 if si==0, sel=use; name='both'; else, sel=use&strcmp(T.session,sessions{si}); name=sessions{si}; end
 if nnz(sel)<4, continue; end
 q=T.abs_q(sel); sq=sqrt((3*D.q_axis(2)-3*D.q_axis(1))^2/12+D.sigma_probe^2)*ones(size(q));
 hi=fit_branch(q,T.E0_high(sel),eT(T.stat_E0_high(sel),T.sys_E0_high(sel)),sq,'quasi2d');
 lo=fit_branch(q,T.E0_low(sel),eT(T.stat_E0_low(sel),T.sys_E0_low(sel)),sq,'gapped');
 f=struct('name',name,'n_bins',nnz(sel),'high',hi,'low',lo);
 if isempty(fits), fits=f; else, fits(end+1)=f; end %#ok<AGROW>
end
save(fullfile(out,'b1_fullq_summary.mat'),'T','A','fits','cfg');
plots(T,fits,S,cfg,figd);
write_report(T,fits,cfg,out,run_root);
fprintf('aggregate: %d bins, %d usable; figures in %s\n',height(T),nnz(use),figd);
end

function r=fit_row(ses,R,f)
p=nan(2,3); p(1:size(f.parameters,1),:)=f.parameters; a=nan(1,2); a(1:numel(f.loss_area))=f.loss_area;
r={ses,R.q,f.key,f.model,f.n,f.aux,f.form,f.h,f.main,f.status,f.cost,f.chi2_signal,f.chi2_gain,f.chi2_low, ...
 p(1,1),p(1,2),p(1,3),p(2,1),p(2,2),p(2,3),a(1),a(2),f.background_fraction};
end

function [m,n1,vars,ndeg]=pick(R)
% Variants differ from the main fit in exactly one of aux / form / h; those
% that did not converge cleanly (zero amplitude, bound) are counted, not used.
f=R.fits; m=f([f.main]);
same=@(g,model,n)strcmp(g.model,model)&&g.n==n;
n1=[]; vars=struct([]); ndeg=0;
for j=1:numel(f)
 g=f(j);
 if same(g,m.model,1)&&strcmp(g.aux,m.aux)&&strcmp(g.form,m.form)&&g.h==m.h, n1=g; end
 if same(g,m.model,2)&&~g.main&&g.success&&(strcmp(g.aux,m.aux)+strcmp(g.form,m.form)+(g.h==m.h))==2
  if ~strcmp(g.status,'converged'), ndeg=ndeg+1; continue; end
  if isempty(vars), vars=g; else, vars(end+1)=g; end %#ok<AGROW>
 end
end
end

function s=boot_sd(R)
s=nan(1,5); if isempty(R.boot), return; end
ok=[R.boot.success]; if nnz(ok)<5, return; end
P=cat(3,R.boot(ok).parameters); La=cat(1,R.boot(ok).loss_area);
s(1:2)=std(P(:,1,:),0,3,'omitnan'); s(3:4)=std(P(:,2,:),0,3,'omitnan');
s(5)=std(La(:,1)./sum(La,2),'omitnan');
end

function s=sys_shift(m,vars)
s=zeros(1,5); if isempty(vars), return; end
w0=m.loss_area(1)/sum(m.loss_area);
for j=1:numel(vars)
 v=vars(j); x=[v.parameters(1,1)-m.parameters(1,1),v.parameters(2,1)-m.parameters(2,1), ...
  v.parameters(1,2)-m.parameters(1,2),v.parameters(2,2)-m.parameters(2,2),v.loss_area(1)/sum(v.loss_area)-w0];
 s=max(s,abs(x));
end
end

function b=synth_summary(run_root,ses,q)
b=nan(1,4); file=fullfile(run_root,'stage3',ses,sprintf('synth_%+.4f.mat',q));
if ~isfile(file), return; end
L=load(file,'Z'); Z=L.Z; P=cat(3,Z.replicas.parameters);
b(1:2)=squeeze(mean(P(:,1,:),3,'omitnan'))-Z.truth(:,1); b(3:4)=squeeze(std(P(:,1,:),0,3,'omitnan'));
end

function r=fit_branch(q,E,sE,sq,kind)
% Weighted fit with effective variance for the q uncertainty (two passes).
ok=isfinite(q)&isfinite(E)&isfinite(sE); q=q(ok); E=E(ok); sE=max(sE(ok),1); sq=sq(ok);
switch kind
 case 'quasi2d', fun=@(p,x)sqrt(max(p(1)*x./(1+p(2)*x),0)); p0=[max(E)^2/max(q) 10]; lb=[0 0]; ub=[Inf Inf];
 case 'gapped', fun=@(p,x)sqrt(max(p(1)^2+p(2)*x,0)); p0=[min(E) 1e5]; lb=[0 -Inf]; ub=[Inf Inf];
end
o=optimoptions('lsqnonlin','Display','off'); s=sE;
for pass=1:2
 [p,~,res,~,~,~,J]=lsqnonlin(@(p)(fun(p,q)-E)./s,p0,lb,ub,o);
 dq=1e-5; slope=(fun(p,q+dq)-fun(p,q-dq))/(2*dq); s=sqrt(sE.^2+(slope.*sq).^2); p0=p;
end
J=full(J); C=inv(J'*J); dof=max(numel(q)-numel(p),1); chi2=sum(res.^2)/dof;
r=struct('kind',kind,'p',p,'p_sd',sqrt(diag(C))'*sqrt(max(chi2,1)),'chi2_red',chi2,'n',numel(q), ...
 'q',q,'E',E,'sE',s,'fun',fun);
end

function f=newfig(pos)
f=figure('Visible','off','Position',pos); if isprop(f,'Theme'), f.Theme='light'; end
end

function plots(T,fits,S,cfg,figd)
ses=fieldnames(S); col=lines(2); cb=[0 0.447 0.741]; co=[0.85 0.325 0.098];
use=strcmp(T.status,'converged')&~T.any_boundary&T.abs_q<=0.03;
% 1. E(q) with both branches, both sessions and signs.
f=newfig([100 100 900 650]); hold on;
mk={'o','s'};
for si=1:numel(ses)
 t=T(strcmp(matlab.lang.makeValidName(T.session),ses{si})&use,:);
 for sg=[-1 1]
  u=t(t.sign==sg,:); face=col(si,:); if sg<0, face='none'; end
  errorbar(u.abs_q,u.E0_low,hypot(u.stat_E0_low,u.sys_E0_low),mk{si},'Color',col(si,:),'MarkerFaceColor',face,'DisplayName',sprintf('%s low, q%s0',ses{si},char(61+sg)));
  errorbar(u.abs_q,u.E0_high,hypot(u.stat_E0_high,u.sys_E0_high),mk{si},'Color',col(si,:)*0.6,'MarkerFaceColor',face,'DisplayName',sprintf('%s high, q%s0',ses{si},char(61+sg)));
 end
end
for k=1:numel(fits)
 if ~strcmp(fits(k).name,'both'), continue; end
 x=linspace(0,max(T.abs_q),200);
 plot(x,fits(k).high.fun(fits(k).high.p,x),'k-','DisplayName','quasi-2D fit (high)');
 plot(x,fits(k).low.fun(fits(k).low.p,x),'k--','DisplayName','gapped fit (low)');
end
xlim([0 0.031]); xlabel('|q| (Å^{-1})'); ylabel('E_0 (meV)'); grid on; legend('Location','southeast','Interpreter','none');
title('B1 two-branch dispersion (DL, kinematic prefactor, joint ZLP)');
exportgraphics(f,fullfile(figd,'dispersion_E0.png'),'Resolution',140); close(f);
% 2. Widths and lower-branch weight.
f=newfig([100 100 1300 480]); tiledlayout(1,3);
nexttile; hold on; for si=1:numel(ses), t=T(strcmp(matlab.lang.makeValidName(T.session),ses{si})&use,:);
 errorbar(t.abs_q,t.W_low,hypot(t.stat_W_low,t.sys_W_low),'o','Color',col(si,:)); errorbar(t.abs_q,t.W_high,hypot(t.stat_W_high,t.sys_W_high),'s','Color',col(si,:)*0.6); end
xlabel('|q| (Å^{-1})'); ylabel('\Gamma (meV)'); title('widths (o low, s high)'); grid on;
nexttile; hold on; for si=1:numel(ses), t=T(strcmp(matlab.lang.makeValidName(T.session),ses{si})&use,:);
 errorbar(t.abs_q,t.weight_low,hypot(t.stat_weight_low,t.sys_weight_low),'o','Color',col(si,:),'DisplayName',ses{si}); end
xlabel('|q| (Å^{-1})'); ylabel('low-branch share of loss-function area'); grid on; legend('Interpreter','none');
nexttile; hold on; for si=1:numel(ses), t=T(strcmp(matlab.lang.makeValidName(T.session),ses{si}),:);
 semilogy(t.abs_q,max(t.dchi2_n1_n2_corr,1),'o','Color',col(si,:),'DisplayName',ses{si}); end
set(gca,'YScale','log'); yline(11.3,'k--','p=0.01 (3 par.)'); xlabel('|q| (Å^{-1})'); ylabel('\Delta\chi^2 n=1\rightarrow2 (corr.)'); grid on;
exportgraphics(f,fullfile(figd,'widths_weight_dchi2.png'),'Resolution',140); close(f);
% 3. Stacked spectra and 4. loss-function maps per session.
for si=1:numel(ses)
 Rs=S.(ses{si}).Rs; q=cellfun(@(r)r.q,Rs); sel=find(abs(q)<=0.03);
 f=newfig([100 100 1300 900]); tiledlayout(1,2);
 for sg=[-1 1]
  nexttile; hold on; kk=sel(sign(q(sel))==sg); [~,o]=sort(abs(q(kk))); kk=kk(o);
  for j=1:numel(kk)
   m=Rs{kk(j)}.fits([Rs{kk(j)}.fits.main]); c=m.curves; if isempty(fieldnames(c)), continue; end
   w=c.E>=200&c.E<=cfg.window(2); sc=max(c.observed(w)); off=1.1*(j-1);
   plot(c.E(w),c.observed(w)/sc+off,'.','Color',[.4 .4 .4],'MarkerSize',3);
   bg=c.zlp(w)+sum(c.aux_peaks(w,:),2);
   plot(c.E(w),c.prediction(w)/sc+off,'k-'); plot(c.E(w),(c.peaks(w,1)+bg)/sc+off,':','Color',cb,'LineWidth',1.1);
   plot(c.E(w),(c.peaks(w,2)+bg)/sc+off,':','Color',co,'LineWidth',1.1); plot(c.E(w),bg/sc+off,'-','Color',[.6 .6 .6]);
   text(cfg.window(2)+20,off+0.3,sprintf('%+.4f',q(kk(j))),'FontSize',7);
  end
  xlim([200 cfg.window(2)+150]); xlabel('E (meV)'); title(sprintf('%s, q%s0 (each scaled to its max)',S.(ses{si}).name,char(61+sg)),'Interpreter','none');
 end
 exportgraphics(f,fullfile(figd,sprintf('stacked_%s.png',S.(ses{si}).name)),'Resolution',130); close(f);
 f=newfig([100 100 900 600]);
 Eg=200:4:cfg.window(2); M=nan(numel(Eg),numel(Rs));
 for j=1:numel(Rs)
  m=Rs{j}.fits([Rs{j}.fits.main]); c=m.curves; if isempty(fieldnames(c)), continue; end
  L=sum(c.loss_function,2); x=interp1(c.E,L,Eg); M(:,j)=x/max(x);
 end
 imagesc(q,Eg,M); axis xy; colorbar; hold on;
 t=S.(ses{si}); qq=cellfun(@(r)r.q,t.Rs); e1=cellfun(@(r)r.fits(1).parameters(1,1),t.Rs); e2=cellfun(@(r)r.fits(1).parameters(2,1),t.Rs);
 plot(qq,e1,'wo',qq,e2,'ws','MarkerSize',4); xlabel('q (Å^{-1})'); ylabel('E (meV)');
 title(sprintf('%s: fitted loss function (prefactor removed), per-bin max = 1',S.(ses{si}).name),'Interpreter','none');
 exportgraphics(f,fullfile(figd,sprintf('loss_map_%s.png',S.(ses{si}).name)),'Resolution',130); close(f);
end
end

function write_report(T,fits,cfg,out,run_root)
fid=fopen(fullfile(out,'b1_fullq_numbers.md'),'w'); c=onCleanup(@()fclose(fid));
fprintf(fid,'# B1 full-q numbers\n\nrun: %s\n\n',run_root);
fprintf(fid,'main configuration: DL n=2, aux %s, kinematic %s, h_perp = %g A^-1, n_zlp = %d\n\n',cfg.aux{1,1},cfg.kinematic{1},cfg.h_main,cfg.n_zlp);
fprintf(fid,'| session | q | E0 low | E0 high | W low | W high | low share | chi2 sig | dchi2 n1->2 (corr) | status |\n|---|---|---|---|---|---|---|---|---|---|\n');
for k=1:height(T)
 fprintf(fid,'| %s | %+.4f | %.0f ± %.0f ± %.0f | %.0f ± %.0f ± %.0f | %.0f | %.0f | %.3f | %.2f | %.0f | %s |\n',T.session{k},T.q(k), ...
  T.E0_low(k),T.stat_E0_low(k),T.sys_E0_low(k),T.E0_high(k),T.stat_E0_high(k),T.sys_E0_high(k),T.W_low(k),T.W_high(k), ...
  T.weight_low(k),T.chi2_signal(k),T.dchi2_n1_n2_corr(k),T.status{k});
end
fprintf(fid,'\n## dispersion fits\n\n');
for k=1:numel(fits)
 f=fits(k);
 fprintf(fid,'- %s (%d bins): high E^2 = a q/(1+b q): a = %.4g ± %.2g meV^2 A, b = %.4g ± %.2g A (chi2_red %.2f); E_sat = %.0f meV\n', ...
  f.name,f.n_bins,f.high.p(1),f.high.p_sd(1),f.high.p(2),f.high.p_sd(2),f.high.chi2_red,sqrt(f.high.p(1)/max(f.high.p(2),eps)));
 fprintf(fid,'  low E^2 = D^2 + s q: D = %.0f ± %.0f meV, s = %.4g ± %.2g meV^2 A (chi2_red %.2f)\n', ...
  f.low.p(1),f.low.p_sd(1),f.low.p(2),f.low.p_sd(2),f.low.chi2_red);
end
end
