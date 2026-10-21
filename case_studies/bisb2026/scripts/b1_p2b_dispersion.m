function out = b1_p2b_dispersion(fit_dir)
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
S=load(fullfile(fit_dir,'p2b_results.mat'),'res','hq','cfg'); res=S.res; hq=S.hq; cfg=S.cfg;
out=fullfile(fit_dir,'dispersion'); if ~isfolder(out), mkdir(out); end
hbarc=1973.27e3;                                          % meV A
intr=qe_lowq_joint_forward('intrinsic','Qn_E0',cfg.nodes.E0,'Qn_G',cfg.nodes.G,'Qn_A',cfg.nodes.A);
Qn=intr.Qn_E0; meas={'E0','apex','pole'}; brn={'low','high'}; rows={}; prow={}; D=struct();
for j=1:2
 for k=1:numel(meas)
  [yn,sn,ss]=local_nodes(res,j,meas{k},Qn);
  ok=Qn>yn/hbarc & yn<=1800 & yn>=250;                  % outside the light cone, inside the window
  x=[Qn(ok) hq.q(:).']; y=[yn(ok) hq.(local_hq(meas{k}))(j,:)]; s=[sn(ok) hq.(['s' local_hq(meas{k})])(j,:)];
  src=[repmat({'P2b'},1,nnz(ok)) repmat({'binned'},1,numel(hq.q))];
  for i=1:numel(Qn), prow(end+1,:)={brn{j},meas{k},Qn(i),yn(i),ss(i),sn(i),Qn(i)>yn(i)/hbarc,ok(i)}; end %#ok<AGROW>
  c=local_const(y,s); g=local_fit(x,y,s,'gapless'); d=local_fit(x,y,s,'gapped');
  [rho_lo,prof]=local_profile(x,y,s,c.chi2);
  D.(sprintf('%s_%s',brn{j},meas{k}))=struct('x',x,'y',y,'s',s,'src',{src},'c',c,'g',g,'d',d,'rho_lo',rho_lo,'prof',prof);
  rows(end+1,:)={brn{j},meas{k},numel(y),c.Ec,c.sEc,c.chi2/(numel(y)-1),g.chi2/(numel(y)-2),d.chi2/(numel(y)-3), ...
   d.Delta,d.sDelta,rho_lo,1/rho_lo,local_slope(x,y,s)}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'branch','measure','n_points','E_const_meV','E_const_err','chi2dof_const','chi2dof_gapless', ...
 'chi2dof_gapped','Delta_meV','Delta_err','rho_lower_bound_A','qstar_upper_bound_Ainv','slope_meV_per_0p01Ainv'});
writetable(T,fullfile(out,'dispersion_summary.csv'));
writetable(cell2table(prow,'VariableNames',{'branch','measure','Q_node','value_meV','sys_meV','total_err_meV','outside_light_cone','used'}), ...
 fullfile(out,'p2b_node_values.csv'));
local_fig(out,D,res,hq,Qn,hbarc);
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

function [rho_lo, prof] = local_profile(x, y, s, chi2c)
% gapless form: chi2 as a function of rho (A refitted); lower bound on rho where
% chi2 exceeds the constant-model chi2 by 4
rg = logspace(1, 6, 101); c2 = zeros(size(rg)); o = optimoptions('lsqnonlin', 'Display', 'off');
for k = 1:numel(rg)
    f = @(a) (sqrt(exp(a) * x ./ (1 + rg(k) * x)) - y) ./ s;
    a0 = log(mean(y .^ 2 .* (1 + rg(k) * x) ./ x)); [~, c2(k)] = lsqnonlin(f, a0, [], [], o);
end
prof = struct('rho', rg, 'chi2', c2); bad = c2 > chi2c + 4; rho_lo = NaN;
if any(~bad), i = find(~bad, 1, 'first'); rho_lo = rg(i); end
end

function b = local_slope(x, y, s)
% weighted linear slope (meV per 0.01 A^-1), a model-free dispersion measure
w = 1 ./ s .^ 2; X = [ones(numel(x), 1) x(:)]; beta = (X.' * (w(:) .* X)) \ (X.' * (w(:) .* y(:))); b = beta(2) * 0.01;
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
