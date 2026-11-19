function out = b1_p2b_dispersion(fit_dir, run_cfg)
% Dispersion of the two B1 branches from P2b (shared fit, |q| <= 0.005) and the
% earlier binned fits (0.0055-0.0155): E0, DL apex and pole. Only the nodes that
% lie outside the light cone at the branch energy and inside the B1 window are
% used for the fits (light-cone check per node: Q > E/(hbar c)). Models:
%   gapless 2D: E^2 = A q / (1 + rho q)   (E -> sqrt(A/rho) for q >> 1/rho)
%   gapped:     E^2 = Delta^2 + A q / (1 + rho q)
%   constant:   E = E_c
% For the gapless form the chi2 profile over rho (A refitted) gives the lower
% bound on rho, i.e. the upper bound on the crossover momentum q* = 1/rho,
% allowed by the data (delta chi2 = 4 relative to the constant model).
% Errors: statistical (Jacobian) and the largest variant deviation, in quadrature.
% Output: <fit_dir>/dispersion/.
% run_cfg.node_analysis enables node-only tables with raw chi2, intercept,
% and a profile interval for Delta. run_cfg.variant_dirs supplies matching
% instrument variants stored in separate fit directories.
if nargin<2, run_cfg=struct(); end
paths=b1_forward_paths(run_cfg);
if nargin<1 || isempty(fit_dir), fit_dir=paths.fit_dir; end
S=load(fullfile(fit_dir,'p2b_results.mat'),'res','hq','cfg'); res=S.res; hq=S.hq; cfg=S.cfg;
nodeAnalysis=isfield(run_cfg,'node_analysis') && run_cfg.node_analysis;
if nodeAnalysis
 hq.q=[];
 for f={'E0','sE0','apex','sapex'}, hq.(f{1})=zeros(2,0); end
end
if isfield(run_cfg,'variant_dirs')
 for v=1:numel(run_cfg.variant_dirs)
  V=load(fullfile(run_cfg.variant_dirs{v},'p2b_results.mat'),'res','cfg');
  assert(isequal(V.cfg.nodes,cfg.nodes),'Variant node grids must match.');
  res(end+1)=V.res(1); %#ok<AGROW>
 end
end
out=fullfile(fit_dir,'dispersion');
if isfield(run_cfg,'output_dir'), out=run_cfg.output_dir; end
if nodeAnalysis, assert(~isfolder(out),'Dispersion output already exists: %s',out); end
if ~isfolder(out), mkdir(out); end
window=cfg.variants(1).window;
hbarc=1973.27e3;                                          % meV A
intr=qe_lowq_joint_forward('intrinsic','Qn_E0',cfg.nodes.E0,'Qn_G',cfg.nodes.G,'Qn_A',cfg.nodes.A);
Qn=intr.Qn_E0; meas={'E0','apex','pole'}; brn={'low','high'}; rows={}; prow={}; details={}; D=struct();
for j=1:2
 for k=1:numel(meas)
  [yn,sn,ss]=local_nodes(res,j,meas{k},Qn);
  ok=Qn>yn/hbarc & yn<=window(2) & yn>=window(1);                  % outside the light cone, inside the window
  if nodeAnalysis, ok=ok & isfinite(sn) & sn>0; end
  if isfield(run_cfg,'min_Q'), ok=ok & Qn>=run_cfg.min_Q; end
  x=[Qn(ok) hq.q(:).']; y=[yn(ok) hq.(local_hq(meas{k}))(j,:)]; s=[sn(ok) hq.(['s' local_hq(meas{k})])(j,:)];
  src=[repmat({'P2b'},1,nnz(ok)) repmat({'binned'},1,numel(hq.q))];
  for i=1:numel(Qn), prow(end+1,:)={brn{j},meas{k},Qn(i),yn(i),ss(i),sn(i),Qn(i)>yn(i)/hbarc,ok(i)}; end %#ok<AGROW>
  c=local_const(y,s);
  if nodeAnalysis && isempty(y), c.chi2=NaN; end
  if nodeAnalysis && numel(y)<3
   emptyFit=struct('chi2',NaN,'Delta',NaN,'sDelta',NaN,'f',@(~,q) nan(size(q)),'p',[]);
   g=emptyFit; d=emptyFit; rho_lo=NaN; prof=[];
   if numel(y)==2
    g=local_fit(x,y,s,'gapless');
    [rho_lo,prof]=local_profile(x,y,s,c.chi2,nodeAnalysis);
   end
  else
   g=local_fit(x,y,s,'gapless'); d=local_fit(x,y,s,'gapped');
   [rho_lo,prof]=local_profile(x,y,s,c.chi2,nodeAnalysis);
  end
  if nodeAnalysis
   L=local_linear(x,y,s);
   dp=local_delta_profile(x,y,s,d);
   details(end+1,:)={brn{j},meas{k},numel(y),c.chi2,g.chi2,d.chi2,c.Ec,c.sEc, ...
    d.Delta,d.sDelta,dp.bounds(1),dp.bounds(2),rho_lo,L(1),L(2),L(3),L(4)}; %#ok<AGROW>
   D.(sprintf('%s_delta_profile_%s',brn{j},meas{k}))=dp;
   dof=numel(y)-(1:3);
   dof(dof<=0)=NaN;
   reduced=[c.chi2 g.chi2 d.chi2]./dof;
   slope=L(3)*10;
  else
   reduced=[c.chi2/(numel(y)-1) g.chi2/(numel(y)-2) d.chi2/(numel(y)-3)];
   slope=local_slope(x,y,s);
  end
  D.(sprintf('%s_%s',brn{j},meas{k}))=struct('x',x,'y',y,'s',s,'src',{src},'c',c,'g',g,'d',d,'rho_lo',rho_lo,'prof',prof);
  rows(end+1,:)={brn{j},meas{k},numel(y),c.Ec,c.sEc,reduced(1),reduced(2),reduced(3), ...
   d.Delta,d.sDelta,rho_lo,1/rho_lo,slope}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'branch','measure','n_points','E_const_meV','E_const_err','chi2dof_const','chi2dof_gapless', ...
 'chi2dof_gapped','Delta_meV','Delta_err','rho_lower_bound_A','qstar_upper_bound_Ainv','slope_meV_per_0p01Ainv'});
writetable(T,fullfile(out,'dispersion_summary.csv'));
writetable(cell2table(prow,'VariableNames',{'branch','measure','Q_node','value_meV','sys_meV','total_err_meV','outside_light_cone','used'}), ...
 fullfile(out,'p2b_node_values.csv'));
if nodeAnalysis
 writetable(cell2table(details,'VariableNames',{'branch','measure','n','chi2_constant','chi2_gapless','chi2_gapped', ...
  'constant_meV','constant_err_meV','Delta_meV','Delta_jac_meV','Delta_profile_lo_meV','Delta_profile_hi_meV', ...
  'rho_lower_A','intercept_meV','intercept_err_meV','slope_meV_per_milliQ','slope_err_meV_per_milliQ'}), ...
  fullfile(out,'dispersion_details.csv'));
end
if ~isfield(run_cfg,'make_figure') || run_cfg.make_figure, local_fig(out,D,res,hq,Qn,hbarc); end
save(fullfile(out,'dispersion.mat'),'D');
disp(T);
end


function [y, s, ss] = local_nodes(res, j, meas, Qn)
C = res(1).curves; ip = @(Q, v) interp1(Q, v, min(max(Qn, Q(1)), Q(end)));   % clamp: the curve grid ends at the nodes up to rounding
y = ip(C.Q, C.(meas)(j, :)); st = ip(C.Q, C.(['s' meas])(j, :)); ss = zeros(size(Qn));
for v = 2:numel(res), ss = max(ss, abs(ip(res(v).curves.Q, res(v).curves.(meas)(j, :)) - y)); end
s = hypot(st, ss);
end

function n = local_hq(meas)
n = meas; if strcmp(meas, 'pole'), n = 'apex'; end   % binned fits: pole not stored, apex is within 2% of the pole (Gamma <= E0)
end

function c = local_const(y, s)
w = 1 ./ s .^ 2; Ec = sum(w .* y) / sum(w); c = struct('Ec', Ec, 'sEc', sqrt(1 / sum(w)) * max(1, sqrt(sum(w .* (y - Ec) .^ 2) / max(numel(y) - 1, 1))), ...
    'chi2', sum(w .* (y - Ec) .^ 2));
end

function d = local_fit(x, y, s, model)
o = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 3000, 'MaxFunctionEvaluations', 3e4);
switch model
    case 'gapless'
        f = @(p, q) sqrt(exp(p(1)) * q ./ (1 + exp(p(2)) * q)); best = [];
        for lr = log([30 100 300 1000 3000 1e4 1e5])
            p0 = [log(max(y) ^ 2 * (1 + exp(lr) * 0.005) / 0.005) lr];
            [p, rn] = lsqnonlin(@(p) (f(p, x) - y) ./ s, p0, [-Inf log(1)], [Inf log(1e7)], o);
            if isempty(best) || rn < best.rn, best = struct('p', p, 'rn', rn); end
        end
        d = struct('model', model, 'p', best.p, 'chi2', best.rn, 'f', f, 'Delta', NaN, 'sDelta', NaN);
    case 'gapped'
        f = @(p, q) sqrt(p(3) ^ 2 + exp(p(1)) * q ./ (1 + exp(p(2)) * q)); best = [];
        for lr = log([30 300 3000 3e4])
            p0 = [log(max(y) ^ 2 / 0.01) lr min(y) * 0.9];
            [p, rn, ~, ~, ~, ~, J] = lsqnonlin(@(p) (f(p, x) - y) ./ s, p0, [-Inf log(1) 0], [Inf log(1e7) 2000], o);
            if isempty(best) || rn < best.rn, best = struct('p', p, 'rn', rn, 'J', J); end
        end
        C = pinv(full(best.J.' * best.J)) * max(1, best.rn / max(numel(y) - 3, 1));
        d = struct('model', model, 'p', best.p, 'chi2', best.rn, 'f', f, 'Delta', best.p(3), 'sDelta', sqrt(max(C(3, 3), 0)));
end
end

function [rho_lo, prof] = local_profile(x, y, s, chi2c, analytic)
% gapless form: chi2 as a function of rho (A refitted); lower bound on rho where
% chi2 exceeds the constant-model chi2 by 4
rg = logspace(1, 6, 101);
if analytic
    B=sqrt(x(:)./(1+x(:).*rg))./s(:);
    yw=y(:)./s(:);
    scale=sum(B.*yw,1)./sum(B.^2,1);
    c2=sum((B.*scale-yw).^2,1);
else
    c2 = zeros(size(rg)); o = optimoptions('lsqnonlin', 'Display', 'off');
    for k = 1:numel(rg)
        f = @(a) (sqrt(exp(a) * x ./ (1 + rg(k) * x)) - y) ./ s;
        a0 = log(mean(y .^ 2 .* (1 + rg(k) * x) ./ x)); [~, c2(k)] = lsqnonlin(f, a0, [], [], o);
    end
end
prof = struct('rho', rg, 'chi2', c2); bad = c2 > chi2c + 4; rho_lo = NaN;
if any(~bad), i = find(~bad, 1, 'first'); rho_lo = rg(i); end
end

function b = local_slope(x, y, s)
% weighted linear slope (meV per 0.01 A^-1), a model-free dispersion measure
w = 1 ./ s .^ 2; X = [ones(numel(x), 1) x(:)]; beta = (X.' * (w(:) .* X)) \ (X.' * (w(:) .* y(:))); b = beta(2) * 0.01;
end

function L = local_linear(x,y,s)
L = nan(1,4);
if numel(y)<2, return; end
X = [ones(numel(x),1),1e3*x(:)]./s(:);
u = y(:)./s(:);
b = X\u;
scale = max(1,sum((X*b-u).^2)/max(numel(u)-2,1));
err = sqrt(diag(inv(X.'*X))*scale);
L = [b(1),err(1),b(2),err(2)];
end

function out = local_delta_profile(x,y,s,best)
% For each Delta and rho, optimize A by monotone bisection of the weighted
% least-squares derivative in A. Work with an equivalent scaled A variable.
out=struct('Delta',[],'chi2',[],'bounds',[NaN NaN]);
if numel(y)<4 || ~isfinite(best.chi2), return; end
delta=unique([0:5:2000 best.Delta]);
rho=unique([logspace(0,7,141) exp(best.p(2))]);
q=x(:);
B=q./(1+q.*rho);
B=B./B(end,:);
B=permute(B,[1 3 2]);
y=y(:); w=1./s(:).^2;
lo=zeros(1,numel(delta),numel(rho));
hi=repmat(max(y.^2)./min(B,[],1),1,numel(delta),1);
for k=1:42
 mid=(hi+lo)/2;
 fit=sqrt(delta.^2+B.*mid);
 positive=sum(w.*B.*(1-y./fit),1)>0;
 hi(positive)=mid(positive);
 lo(~positive)=mid(~positive);
end
fit=sqrt(delta.^2+B.*((hi+lo)/2));
chi2=min(sum(w.*(fit-y).^2,1),[],3);
use=chi2<=min(best.chi2,min(chi2))+4;
out.Delta=delta;
out.chi2=chi2;
if any(use), out.bounds=[delta(find(use,1)) delta(find(use,1,'last'))]; end
end

function local_fig(out, D, res, hq, Qn, hbarc)
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1750 1000]); tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
cl = struct('E0', [0 0.3 0.8], 'apex', [0.8 0.2 0.1], 'pole', [0.3 0.6 0.3]); brn = {'low', 'high'};
for sc = 1:2
    for j = [2 1]
        ax = nexttile(tl); hold(ax, 'on'); C = res(1).curves;
        for m = {'E0', 'apex'}
            d = D.(sprintf('%s_%s', brn{j}, m{1})); p = strcmp(d.src, 'P2b');
            errorbar(ax, d.x(p) * 1e3, d.y(p), d.s(p), 's', 'Color', cl.(m{1}), 'MarkerFaceColor', cl.(m{1}), 'DisplayName', sprintf('%s, P2b nodes (used)', m{1}));
            errorbar(ax, d.x(~p) * 1e3, d.y(~p), d.s(~p), 'o', 'Color', cl.(m{1}), 'DisplayName', sprintf('%s, binned fits', m{1}));
            plot(ax, C.Q * 1e3, C.(m{1})(j, :), '-', 'Color', [cl.(m{1}) 0.35], 'LineWidth', 1.5, 'DisplayName', sprintf('%s, P2b curve', m{1}));
            qq = logspace(-4, log10(0.02), 300);
            plot(ax, qq * 1e3, d.d.f(d.d.p, qq), '-.', 'Color', cl.(m{1}), 'DisplayName', sprintf('%s gapped fit', m{1}));
            if isfinite(d.rho_lo)
                a = mean(d.y .^ 2 .* (1 + d.rho_lo * d.x) ./ d.x);
                plot(ax, qq * 1e3, sqrt(a * qq ./ (1 + d.rho_lo * qq)), ':', 'Color', cl.(m{1}), 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%s gapless at the bound rho = %.0f A', m{1}, d.rho_lo));
            end
        end
        qq = logspace(-4, log10(0.02), 200); plot(ax, qq * 1e3, hbarc * qq, 'k--', 'DisplayName', 'light line E = hbar c Q');
        grid(ax, 'on'); ylim(ax, [0 1800]); xlabel(ax, 'Q (10^{-3} A^{-1})'); ylabel(ax, 'energy (meV)');
        if sc == 1, xlim(ax, [0 16]); else, set(ax, 'XScale', 'log'); xlim(ax, [0.2 20]); end
        nm = {'low branch', 'high branch'}; axn = {'linear', 'log'}; title(ax, sprintf('%s (%s Q axis)', nm{j}, axn{sc}));
        if sc == 1, legend(ax, 'FontSize', 6, 'Location', 'southeast'); end
    end
end
exportgraphics(f, fullfile(out, 'fig_dispersion_E0_apex.png'), 'Resolution', 110); close(f);
end
