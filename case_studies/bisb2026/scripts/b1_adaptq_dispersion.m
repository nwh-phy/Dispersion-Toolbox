function b1_adaptq_dispersion(run_root, session)
% B1 adaptive-binning task, stage 6: dispersion fits of the mode energies
% Re(pole) from stage 5. Points: converged, no boundary, two components
% required (dchi2 >= 11.3), no degenerate lower component (Gamma >= 300 meV,
% share >= 0.05), |q| <= 0.024. Abscissa: kinematically averaged q_mean.
% Error per point: frame bootstrap (+) kinematic form |3D - 2D| (+) h variants.
if nargin<2, session='590_PL2_10w'; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
out=fullfile(run_root,'stage5'); T=readtable(fullfile(out,'adaptq_main.csv')); T=T(strcmp(T.session,session),:);
ok=strcmp(T.status,'converged')&~T.any_boundary&T.dchi2_n1_n2>=11.3&T.G_low>=300&T.weight_low>=0.05&T.abs_q<=0.0241;
T=T(ok,:); q=T.q_mean;
hs=@(x) fillmissing(x,'constant',0);
err.high=sqrt(T.boot_pole_high.^2+(T.pole_high-T.pole_high_2d).^2+hs(T.h_sys_pole_high).^2);
err.low=sqrt(T.boot_pole_low.^2+(T.pole_low-T.pole_low_2d).^2+hs(T.h_sys_pole_low).^2);
models={ ...
 'quasi2D sqrt(A q/(1+rho0 q))',        @(p,q) sqrt(p(1)*q./(1+p(2)*q)),          [2e9 1000],     [0 0],   [Inf 1e5]; ...
 'gapped quasi2D sqrt(D^2+A q/(1+rho0 q))', @(p,q) sqrt(p(3)^2+p(1)*q./(1+p(2)*q)), [1e9 500 800],  [0 0 0], [Inf 1e5 3000]; ...
 'gapped linear sqrt(D^2+s q)',          @(p,q) sqrt(p(2)^2+p(1)*q),               [5e7 1000],     [-Inf 0],[Inf 3000]; ...
 'constant',                            @(p,q) p(1)+0*q,                          1000,           0,       3000};
opts=optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',2e4,'MaxIterations',2e3);
ranges={'|q|>=0.0035 (3D~2D)',T.abs_q>=0.0035-1e-9; 'all |q|>=0.001',T.abs_q>=0.001-1e-9};
fid=fopen(fullfile(out,sprintf('dispersion_fits_%s.md',session)),'w');
fprintf(fid,'# Dispersion fits of Re(pole), %s\n\nUnits: E in meV, q in 1/A (q_mean). A in meV^2*A, rho0 in A, D in meV, s in meV^2*A.\n\n',session);
for br={'high','low'}
 b=br{1}; E=T.(['pole_' b]); e=err.(b);
 fprintf(fid,'## %s branch\n\n| range | model | n | params | chi2_red |\n|---|---|---|---|---|\n',b);
 for r=1:size(ranges,1)
  u=ranges{r,2}&isfinite(E)&isfinite(e)&e>0; qq=q(u); y=E(u); w=e(u);
  for m=1:size(models,1)
   fn=models{m,2}; res=@(p) (fn(p,qq)-y)./w; best=[]; bc=Inf;
   for st=1:20
    p0=models{m,3}.*exp(0.7*randn(size(models{m,3}))); if m==3, p0(2)=600+400*rand; end
    p0=min(max(p0,models{m,4}),models{m,5});
    try, [p,c]=lsqnonlin(res,p0,models{m,4},models{m,5},opts); catch, continue; end
    if c<bc, bc=c; best=p; end
   end
   dof=numel(y)-numel(best);
   fprintf(fid,'| %s | %s | %d | %s | %.2f |\n',ranges{r,1},models{m,1},numel(y),mat2str(best,4),bc/max(dof,1));
   fits.(b){r,m}=best;
  end
 end
 % weighted linear trend of E vs q, |q| >= 0.0035
 u=ranges{1,2}&isfinite(E)&isfinite(e)&e>0; X=[ones(nnz(u),1) q(u)]; W=diag(1./e(u).^2);
 C=inv(X'*W*X); beta=C*X'*W*E(u);
 fprintf(fid,'\nWeighted linear trend, |q| >= 0.0035: slope %.0f +- %.0f meV per A^-1 (%.1f sigma); change over 0.0035-0.023 = %.0f +- %.0f meV\n\n', ...
  beta(2),sqrt(C(2,2)),beta(2)/sqrt(C(2,2)),beta(2)*(0.023-0.0035),sqrt(C(2,2))*(0.023-0.0035));
end
fclose(fid);
type(fullfile(out,sprintf('dispersion_fits_%s.md',session)));

% Figure: Re(pole) with total error bars and the fits
f=figure('Visible','off','Position',[0 0 900 600]); ax=axes(f); hold(ax,'on'); grid(ax,'on');
col=struct('high',[0.80 0.25 0.10],'low',[0.10 0.35 0.75]); qf=linspace(0.0008,0.025,300);
for br={'high','low'}
 b=br{1}; errorbar(ax,q,T.(['pole_' b]),err.(b),'o','Color',col.(b),'MarkerFaceColor',col.(b),'MarkerSize',4,'CapSize',2,'DisplayName',sprintf('%s branch Re(pole), 3D; bars: stat+kin+h',b));
 p=fits.(b){1,2}; if ~isempty(p), plot(ax,qf,models{2,2}(p,qf),'-','Color',col.(b),'LineWidth',1.2,'DisplayName','gapped quasi-2D fit, |q|>=0.0035'); end
 p=fits.(b){1,1}; if ~isempty(p), plot(ax,qf,models{1,2}(p,qf),'--','Color',col.(b),'LineWidth',1,'DisplayName','quasi-2D fit (no gap), |q|>=0.0035'); end
end
xlabel(ax,'q_{mean} (A^{-1})'); ylabel(ax,'mode energy Re(pole) (meV)'); xlim(ax,[0 0.025]); ylim(ax,[0 1800]);
legend(ax,'Location','southeast','FontSize',8); title(ax,sprintf('%s, adaptive binning (N=1 for |q|<=0.006)',strrep(session,'_','\_')));
exportgraphics(f,fullfile(out,'figures',sprintf('dispersion_fit_%s.png',session)),'Resolution',140); close(f);
end
