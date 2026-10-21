function out = b1_p2b_fit(run_root, cfg, ocfg)
% P2b on 590: one shared loss function L(Q,E) for the 20 native channels
% |q| = 0.0005-0.005 (both signs) and the six time blocks, seen through the
% block operators of P2a step 2 (slice response, measured drift, two-Pearson
% angular weights, true-Q integral of K L). Two DL branches, E0^2 at five and
% the width at three Q nodes, residues proportional to Q x nonnegative hats.
% Data term 250-1800 meV (B1), 150-250 meV held out; sigma^2 = 21-point moving
% mean of the counts; one power-law tail (E/250)^-4 per spectrum, projected out
% (absorbs auxiliary-peak and instrument tails at the low edge). Several starts.
% Variants (one change each against the main model) give the systematic band:
% exit / integrated kernels, broad component off for B1, gapless low-Q
% extension, 2D-mixture y distribution, window 300-1800 without tail, window
% 250-1800 without tail. Block 6 is left out of a second main fit and predicted.
% Reported for every node and on a Q grid: E0, width, DL apex and pole with
% statistical errors (Jacobian, scaled by chi2/dof). The dispersion (E0 and apex)
% is assembled with the earlier binned fits at 0.0055-0.0155 and fitted with
% E^2 = A q / (1 + rho q) (gapless) and E^2 = Delta^2 + A q / (1 + rho q) (gapped).
% cfg and ocfg default to b1_p2b_config and b1_p2a_operator_config; reruns pass
% modified copies (b1_p2b_slot_config). Only cfg names the output folder, so a run
% with another operator config needs its own cfg.version (checked against the
% config.json of an existing folder).
% Output: <run_root>/p2b/fit_590_<config hash>/.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
if nargin<2, cfg=b1_p2b_config(); end
if nargin<3, ocfg=b1_p2a_operator_config(); end
ocfg.exit_q=cfg.exit_q; tag=b1_sha256(cfg); t0=tic;
out=fullfile(run_root,'p2b',sprintf('fit_590_%s',tag(1:8)));
if isfile(fullfile(out,'config.json'))
 old=jsondecode(fileread(fullfile(out,'config.json')));
 assert(isfield(old,'operator_config') && strcmp(jsonencode(jsondecode(jsonencode(ocfg))),jsonencode(old.operator_config)), ...
  'b1_p2b_fit:OperatorConfig','%s holds a run with another operator config; give this one its own cfg.version',out);
end
if ~isfolder(out), mkdir(out); end
ind=fullfile(run_root,'p2a','inputs_590_PL2_10w_66385d96'); rd=fullfile(ind,'response_42ffc295');
L=load(fullfile(ind,'p2a_cube.mat'),'cube'); cube=L.cube; q=cube.q(:).'; dq=cube.dq;
RR=load(fullfile(rd,'response.mat'),'R','J6'); R=RR.R; J6=RR.J6;
ie=arrayfun(@(v) find(abs(q-v)<1e-9),cfg.exit_q); ni=numel(ie);
intr0=qe_lowq_joint_forward('intrinsic','Qn_E0',cfg.nodes.E0,'Qn_G',cfg.nodes.G,'Qn_A',cfg.nodes.A);
base=struct('blocking','b6','variant','slice','ytype','gauss','wide','on','ext','const','n_a',ocfg.n_a, ...
 'n_sx',ocfg.n_sx,'n_sy',ocfg.n_sy,'n_gh',ocfg.n_gh,'Qf',1,'step',ocfg.omega(2),'support',ocfg.support, ...
 'tail_max',ocfg.tail_edges(end),'exits',1:ni,'blocks',[]);
opts=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10,'StepTolerance',1e-8,'OptimalityTolerance',1e-8, ...
 'MaxIterations',300,'MaxFunctionEvaluations',8000);

res=struct('name',{},'fit',{},'curves',{},'chan',{},'perf',{});
for v=1:numel(cfg.variants)
 V=cfg.variants(v); st=base; st.variant=V.kernel; st.wide=V.wide; st.ytype=V.ytype; intr=intr0; intr.ext=V.ext;
 tb=tic; [ops,K,sw,ex]=b1_p2a_build_ops(st,ocfg,cube,R,J6,q,dq,ie);
 cache=qe_lowq_joint_forward('build',ops,K,local_fopt(st,ocfg,sw,ex,V.kernel)); clear ops K; tbuild=toc(tb);
 [Y,S2,Eo]=local_data(cube,cache.E_out,ie,cfg.noise_smooth); nb=size(Y,3);
 mwin=Eo>=V.window(1)&Eo<=V.window(2); mask=repmat(mwin,1,ni,nb)&~isnan(Y);
 lo=abs(abs(cfg.exit_q)-0.0005)<1e-9; mask(Eo>V.lowq_emax,lo,:)=false;     % retardation check: lowest channels only to lowq_emax
 tf=tic; F=local_multistart(cache,Y,S2,mask,Eo,intr,cfg,V.tail,opts); tfit=toc(tf);
 F=local_post(F,cache,Y,S2,Eo,intr,cfg,V,mwin);
 res(v).name=V.name; res(v).fit=F; res(v).curves=local_curves(F,intr,cfg); res(v).chan=F.chan;
 res(v).perf=struct('build_s',tbuild,'fit_s',tfit);
 fprintf('%-30s chi2/dof %.3f  holdout chi2/n %.2f  start %s  (%.0f s)\n',V.name,F.chi2dof,F.ho_chi2n,F.start,tfit);
 if v==1
  % leave block cfg.lobo_block out, predict it with L fixed (tails refitted as nuisance)
  mk=mask; mk(:,:,cfg.lobo_block)=false;
  G=local_multistart(cache,Y,S2,mk,Eo,intr,cfg,V.tail,opts);
  lobo=local_lobo(G,cache,Y,S2,Eo,intr,cfg,V,mwin,mask);
  writetable(struct2table(lobo),fullfile(out,'p2b_lobo.csv'));
  fprintf('LOBO block %d: chi2/n fitted blocks %.3f, predicted block %.3f\n',cfg.lobo_block,lobo.chi2n_fitted,lobo.chi2n_predicted);
  mainY=Y; mainS2=S2; mainEo=Eo; mainMu=F.mu;
 end
 clear cache
end

% ---- tables
frow={}; nrow={};
for v=1:numel(res)
 F=res(v).fit;
 frow(end+1,:)={res(v).name,F.chi2,F.dof,F.chi2dof,F.ho_chi2n,F.ho_mean,F.n_tail_neg,F.exitflag,F.start,res(v).perf.fit_s}; %#ok<AGROW>
 for j=1:2
  for n=1:5, nrow(end+1,:)={res(v).name,j,'E0',intr0.Qn_E0(n),F.E0(j,n),F.sE0(j,n),F.E0(j,n)<=cfg.bounds_E0(j,1)+1|F.E0(j,n)>=cfg.bounds_E0(j,2)-1}; end %#ok<AGROW>
  for n=1:3, nrow(end+1,:)={res(v).name,j,'G',intr0.Qn_G(n),F.G(j,n),F.sG(j,n),F.G(j,n)<=cfg.bounds_G(1)+1|F.G(j,n)>=cfg.bounds_G(2)-1}; end %#ok<AGROW>
 end
end
writetable(cell2table(frow,'VariableNames',{'variant','chi2','dof','chi2dof','holdout_chi2_per_point','holdout_mean_residual_sigma', ...
 'negative_tails','exitflag','best_start','fit_s'}),fullfile(out,'p2b_fits.csv'));
writetable(cell2table(nrow,'VariableNames',{'variant','branch','parameter','Q_node','value_meV','stat_err_meV','at_bound'}),fullfile(out,'p2b_nodes.csv'));
C=res(1).curves; CT=table(C.Q(:),'VariableNames',{'Q'});
for j=1:2
 for f={'E0','G','apex','pole'}
  CT.(sprintf('%s_b%d',f{1},j))=C.(f{1})(j,:).'; CT.(sprintf('%s_b%d_err',f{1},j))=C.(['s' f{1}])(j,:).';
  sys=zeros(numel(C.Q),1);
  for v=2:numel(res), sys=max(sys,abs(res(v).curves.(f{1})(j,:).'-C.(f{1})(j,:).')); end
  CT.(sprintf('%s_b%d_sys',f{1},j))=sys;
 end
end
writetable(CT,fullfile(out,'p2b_curves_main.csv'));
writetable(struct2table(res(1).chan),fullfile(out,'p2b_channel_scores_main.csv'));

% ---- dispersion: P2b nodes (E0 nodes, both measures) + binned fits at 0.0055-0.0155
H=readtable(fullfile(run_root,cfg.highq.file)); H=H(abs(H.q)>=cfg.highq.q_range(1)-1e-9 & abs(H.q)<=cfg.highq.q_range(2)+1e-9,:);
qa=unique(abs(H.q)); hq=struct('q',qa,'E0',zeros(2,numel(qa)),'sE0',zeros(2,numel(qa)),'apex',zeros(2,numel(qa)),'sapex',zeros(2,numel(qa)));
for k=1:numel(qa)
 x=H(abs(abs(H.q)-qa(k))<1e-9,:);
 hq.E0(1,k)=mean(x.E0_low); hq.E0(2,k)=mean(x.E0_high); hq.apex(1,k)=mean(x.apex_low); hq.apex(2,k)=mean(x.apex_high);
 hd=@(y) max(abs(diff(y))/2,cfg.highq.err_floor); if height(x)<2, hd=@(y) 3*cfg.highq.err_floor; end
 hq.sE0(1,k)=hypot(hd(x.E0_low),mean(x.E0_low_sd)); hq.sE0(2,k)=hypot(hd(x.E0_high),mean(x.E0_high_sd));
 hq.sapex(1,k)=hd(x.apex_low); hq.sapex(2,k)=hd(x.apex_high);
end
Qn=intr0.Qn_E0; drow={}; dfit=struct();
for j=1:2
 for meas={'E0','apex'}
  [yn,sn]=local_at_nodes(res,j,meas{1},Qn,intr0);
  x=[Qn hq.q(:).']; y=[yn hq.(meas{1})(j,:)]; s=[sn hq.(['s' meas{1}])(j,:)];
  for model={'gapless','gapped'}
   d=local_dispfit(x,y,s,model{1}); dfit.(sprintf('b%d_%s_%s',j,meas{1},model{1}))=d;
   drow(end+1,:)={j,meas{1},model{1},d.A,d.sA,d.rho,d.srho,d.Delta,d.sDelta,d.Einf,d.chi2,d.dof,d.chi2/max(d.dof,1)}; %#ok<AGROW>
  end
 end
end
writetable(cell2table(drow,'VariableNames',{'branch','measure','model','A_meV2A','A_err','rho_A','rho_err','Delta_meV','Delta_err', ...
 'E_inf_meV','chi2','dof','chi2dof'}),fullfile(out,'p2b_dispersion_fits.csv'));

% ---- figures
local_fig_dispersion(out,res,hq,dfit,cfg,intr0);
local_fig_spectra(out,mainEo,mainY,mainS2,mainMu,cfg);
local_fig_variants(out,res,intr0);
save(fullfile(out,'p2b_results.mat'),'cfg','res','hq','dfit','mainEo','mainY','mainS2','mainMu');
meta=struct('config',cfg,'config_sha256',tag,'code_sha256',b1_sha256({'file',[mfilename('fullpath') '.m']}), ...
 'operator_config',ocfg,'inputs',ind,'response',rd,'runtime_s',toc(t0), ...
 'created',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
fprintf('P2b written to %s (%.0f s)\n',out,toc(t0));
end


% =============================== fitting ===============================
function fo = local_fopt(st, cfg, sw, ex, kernel)
fo = struct('omega', cfg.omega(1):st.step:cfg.omega(3), 'eps', -72:st.step:72, 'variant', kernel, ...
    'slice_weight', {sw}, 'exit_slice', ex);
end

function [Y, S2, Eo] = local_data(cube, Eout, ie, ns)
Eo = Eout(:); E = cube.E(:); J = cube.b6.J; nb = size(J, 3);
Y = nan(numel(Eo), numel(ie), nb); [tf, loc] = ismember(round(Eo), round(E));
for b = 1:nb, Y(tf, :, b) = J(loc(tf), ie, b); end
S2 = max(movmean(Y, ns, 1, 'omitnan'), 1);
end

function F = local_multistart(cache, Y, S2, mask, Eo, intr, cfg, tailOn, opts)
best = [];
for s = 1:numel(cfg.starts)
    th0 = struct('E0sq', cfg.starts(s).E0 .^ 2, 'G', cfg.starts(s).G);
    G = local_fit(cache, Y, S2, mask, Eo, intr, th0, cfg, tailOn, opts); G.start = cfg.starts(s).name;
    if isempty(best) || G.chi2 < best.chi2, best = G; end
end
F = best;
end

function F = local_fit(cache, Y, S2, mask, Eo, intr, th0, cfg, tailOn, opts)
lb = [log(cfg.bounds_E0(1, 1) ^ 2) * ones(1, 5), log(cfg.bounds_E0(2, 1) ^ 2) * ones(1, 5), log(cfg.bounds_G(1)) * ones(1, 6)];
ub = [log(cfg.bounds_E0(1, 2) ^ 2) * ones(1, 5), log(cfg.bounds_E0(2, 2) ^ 2) * ones(1, 5), log(cfg.bounds_G(2)) * ones(1, 6)];
u0 = [reshape(log(th0.E0sq).', 1, []), reshape(log(th0.G).', 1, [])]; u0 = min(max(u0, lb + 1e-6), ub - 1e-6);
D = local_design(Y, S2, mask, Eo, cfg, tailOn);
f = @(u) local_res(u, cache, D, intr, cfg.lambda_curv);
[u, ~, ~, flag, outp, ~, Jac] = lsqnonlin(f, u0, lb, ub, opts);
[r, a, b] = f(u); nd = numel(D.y);
chi2 = sum(r(1:nd) .^ 2); dof = nd - 16 - numel(a) - tailOn * D.ns;
Cov = inv(full(Jac.' * Jac)) * chi2 / dof;
F = struct('u', u, 'a', a, 'b', b, 'chi2', chi2, 'dof', dof, 'chi2dof', chi2 / dof, 'cov', Cov, 'exitflag', flag, ...
    'iterations', outp.iterations, 'mask', mask, 'tailOn', tailOn, 'start', '');
end

function D = local_design(Y, S2, mask, Eo, cfg, tailOn)
[nE, ni, nb] = size(Y); sidx = repmat(reshape(1:ni * nb, 1, ni, nb), nE, 1, 1); Ef = repmat(Eo, 1, ni, nb);
D = struct('y', Y(mask), 'w', 1 ./ sqrt(S2(mask)), 's', sidx(mask), 'T', (Ef(mask) / cfg.tail.E_ref) .^ (-cfg.tail.p), ...
    'mask', mask, 'ns', ni * nb, 'tailOn', tailOn);
end

function [r, a, b] = local_res(u, cache, D, intr, lam)
th = local_decode(u); cols = qe_lowq_joint_forward('eval', cache, th, intr);
A = reshape(permute(cols, [1 3 4 2]), [], size(cols, 2)); A = A(D.mask(:), :);
yw = D.y .* D.w; Aw = A .* D.w; b = zeros(D.ns, 1);
if D.tailOn
    tw = D.T .* D.w; tt = accumarray(D.s, tw .^ 2, [D.ns 1]); ty = accumarray(D.s, tw .* yw, [D.ns 1]);
    tA = zeros(D.ns, size(A, 2)); for k = 1:size(A, 2), tA(:, k) = accumarray(D.s, tw .* Aw(:, k), [D.ns 1]); end
    ywp = yw - tw .* ty(D.s) ./ tt(D.s); Awp = Aw - tw .* tA(D.s, :) ./ tt(D.s);
else
    ywp = yw; Awp = Aw;
end
sc = max(abs(Awp), [], 1); sc(sc == 0) = 1; a = lsqnonneg(Awp ./ sc, ywp) ./ sc(:);
if D.tailOn, b = (ty - tA * a) ./ tt; end
r = [Awp * a - ywp; lam * reshape(diff(log(th.E0sq), 2, 2), [], 1)];
end

function th = local_decode(u)
th = struct('E0sq', exp(reshape(u(1:10), 5, 2).'), 'G', exp(reshape(u(11:16), 3, 2).'));
end

function F = local_post(F, cache, Y, S2, Eo, intr, cfg, V, mwin)
% model on all energies (B1 + extrapolated tails), holdout and per-channel scores
th = local_decode(F.u); F.E0 = sqrt(th.E0sq); F.G = th.G;
J = local_jac_nodes(F.u); s = sqrt(max(diag(J * F.cov * J.'), 0)); F.sE0 = reshape(s(1:10), 5, 2).'; F.sG = reshape(s(11:16), 3, 2).';
cols = qe_lowq_joint_forward('eval', cache, th, intr); mu = squeeze(sum(cols .* reshape(F.a, 1, [], 1, 1), 2));
[nE, ni, nb] = size(Y); T = (Eo / cfg.tail.E_ref) .^ (-cfg.tail.p); T(Eo < cfg.holdout(1)) = 0;
if F.tailOn, mu = mu + T .* reshape(F.b, 1, ni, nb); end
F.mu = mu; ho = Eo >= cfg.holdout(1) & Eo < cfg.holdout(2); z = (Y - mu) ./ sqrt(S2);
F.ho_chi2n = mean(z(ho, :, :) .^ 2, 'all', 'omitnan'); F.ho_mean = mean(z(ho, :, :), 'all', 'omitnan');
F.n_tail_neg = nnz(F.b < 0) * F.tailOn;
chan = struct('q', num2cell(cfg.exit_q(:)), 'chi2dof_window', [], 'holdout_mean_sigma', [], 'holdout_chi2n', []);
for i = 1:ni
    zi = z(mwin, i, :); zh = z(ho, i, :);
    chan(i).chi2dof_window = mean(zi(:) .^ 2, 'omitnan'); chan(i).holdout_mean_sigma = mean(zh(:), 'omitnan'); chan(i).holdout_chi2n = mean(zh(:) .^ 2, 'omitnan');
end
F.chan = chan;
end

function J = local_jac_nodes(u)
% d(E0 nodes, G nodes)/du for u = [log E0^2 (10), log G (6)]
E0 = sqrt(exp(u(1:10))); G = exp(u(11:16)); J = diag([E0 / 2, G]);
end

function C = local_curves(F, intr, cfg)
Q = cfg.Q_curve; f = @(u) local_curve_vals(u, Q, intr);
v0 = f(F.u); g = zeros(numel(v0), 16); h = 1e-5;
for k = 1:16, up = F.u; up(k) = up(k) + h; g(:, k) = (f(up) - v0) / h; end
s = sqrt(max(sum((g * F.cov) .* g, 2), 0)); n = numel(Q);
C = struct('Q', Q); names = {'E0', 'G', 'apex', 'pole'};
for k = 1:4
    C.(names{k}) = reshape(v0((k - 1) * 2 * n + (1:2 * n)), 2, n);
    C.(['s' names{k}]) = reshape(s((k - 1) * 2 * n + (1:2 * n)), 2, n);
end
end

function v = local_curve_vals(u, Q, intr)
th = local_decode(u); lQ = log(Q); n = numel(Q); E0 = zeros(2, n); G = E0;
for j = 1:2
    xn = log(intr.Qn_E0); E0(j, :) = sqrt(exp(interp1(xn, log(th.E0sq(j, :)), min(max(lQ, xn(1)), xn(end)))));
    xg = log(intr.Qn_G); G(j, :) = exp(interp1(xg, log(th.G(j, :)), min(max(lQ, xg(1)), xg(end))));
end
x = 2 * E0 .^ 2 - G .^ 2; apex = sqrt((x + sqrt(x .^ 2 + 12 * E0 .^ 4)) / 6); pole = sqrt(max(E0 .^ 2 - G .^ 2 / 4, 0));
v = [E0(:); G(:); apex(:); pole(:)];      % each 2 x n, column-major (branches interleaved), as read by local_curves
end

function lobo = local_lobo(G, cache, Y, S2, Eo, intr, cfg, V, mwin, mask)
th = local_decode(G.u); cols = qe_lowq_joint_forward('eval', cache, th, intr);
mu = squeeze(sum(cols .* reshape(G.a, 1, [], 1, 1), 2)); [~, ni, nb] = size(Y);
T = (Eo / cfg.tail.E_ref) .^ (-cfg.tail.p); z2 = zeros(nb, 1); n = z2;
for b = 1:nb
    for i = 1:ni
        y = Y(mwin, i, b); m = mu(mwin, i, b); w = 1 ./ S2(mwin, i, b);
        if V.tail, t = T(mwin); c = sum(w .* t .* (y - m)) / sum(w .* t .^ 2); m = m + c * t; end
        z2(b) = z2(b) + sum(w .* (y - m) .^ 2); n(b) = n(b) + numel(y);
    end
end
fit_b = setdiff(1:nb, cfg.lobo_block);
lobo = struct('left_out_block', cfg.lobo_block, 'chi2n_fitted', sum(z2(fit_b)) / sum(n(fit_b)), ...
    'chi2n_predicted', z2(cfg.lobo_block) / n(cfg.lobo_block), 'E0_low_nodes', mat2str(round(sqrt(th.E0sq(1, :)))), ...
    'E0_high_nodes', mat2str(round(sqrt(th.E0sq(2, :)))));
end

% =============================== dispersion ===============================
function [y, s] = local_at_nodes(res, j, meas, Qn, intr)
C = res(1).curves; y = interp1(C.Q, C.(meas)(j, :), Qn); st = interp1(C.Q, C.(['s' meas])(j, :), Qn); sy = zeros(size(Qn));
for v = 2:numel(res), sy = max(sy, abs(interp1(res(v).curves.Q, res(v).curves.(meas)(j, :), Qn) - y)); end
s = hypot(st, sy);                         % statistical and the largest variant deviation
end

function d = local_dispfit(x, y, s, model)
m = isfinite(y) & isfinite(s) & s > 0; x = x(m); y = y(m); s = s(m);
switch model
    case 'gapless', f = @(p, q) sqrt(exp(p(1)) * q ./ (1 + exp(p(2)) * q)); p0 = [log(max(y) ^ 2 / 0.005 * 2) log(200)];
    case 'gapped', f = @(p, q) sqrt(exp(p(3)) ^ 2 + exp(p(1)) * q ./ (1 + exp(p(2)) * q)); p0 = [log(max(y) ^ 2 / 0.005) log(200) log(min(y) * 0.8)];
end
best = []; o = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 2000, 'MaxFunctionEvaluations', 1e4);
for k = 1:9
    pk = p0 + [0.5 * randn(1, 2), zeros(1, numel(p0) - 2)] * (k > 1);
    [p, rn, ~, ~, ~, ~, Jc] = lsqnonlin(@(p) (f(p, x) - y) ./ s, pk, [], [], o);
    if isempty(best) || rn < best.rn, best = struct('p', p, 'rn', rn, 'J', Jc); end
end
p = best.p; dof = numel(y) - numel(p); C = inv(full(best.J.' * best.J)) * max(best.rn / max(dof, 1), 1);
e = sqrt(max(diag(C), 0)).' .* exp(p);       % log parameters -> absolute errors
d = struct('model', model, 'A', exp(p(1)), 'sA', e(1), 'rho', exp(p(2)), 'srho', e(2), 'Delta', NaN, 'sDelta', NaN, ...
    'chi2', best.rn, 'dof', dof, 'f', f, 'p', p, 'Einf', NaN);
if strcmp(model, 'gapped'), d.Delta = exp(p(3)); d.sDelta = e(3); d.Einf = sqrt(exp(p(3)) ^ 2 + exp(p(1)) / exp(p(2)));
else, d.Einf = sqrt(exp(p(1)) / exp(p(2))); end
end

% =============================== figures ===============================
function local_fig_dispersion(out, res, hq, dfit, cfg, intr)
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1700 1050]); tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
C = res(1).curves; cl = [0 0.3 0.8; 0.8 0.2 0.1];
for j = [2 1]
    ax = nexttile(tl); hold(ax, 'on');
    for v = 2:numel(res), plot(ax, res(v).curves.Q * 1e3, res(v).curves.E0(j, :), '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off'); end
    fill(ax, [C.Q fliplr(C.Q)] * 1e3, [C.E0(j, :) + C.sE0(j, :), fliplr(C.E0(j, :) - C.sE0(j, :))], cl(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax, [C.Q fliplr(C.Q)] * 1e3, [C.apex(j, :) + C.sapex(j, :), fliplr(C.apex(j, :) - C.sapex(j, :))], cl(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, C.Q * 1e3, C.E0(j, :), '-', 'Color', cl(1, :), 'LineWidth', 2, 'DisplayName', 'E0, P2b (main; grey: variants)');
    plot(ax, C.Q * 1e3, C.apex(j, :), '--', 'Color', cl(2, :), 'LineWidth', 2, 'DisplayName', 'DL apex, P2b');
    errorbar(ax, hq.q * 1e3, hq.E0(j, :), hq.sE0(j, :), 'o', 'Color', cl(1, :), 'DisplayName', 'E0, binned fits (+/-q mean)');
    errorbar(ax, hq.q * 1e3, hq.apex(j, :), hq.sapex(j, :), '^', 'Color', cl(2, :), 'DisplayName', 'DL apex, binned fits');
    qq = linspace(0, 16e-3, 300);
    for meas = {'E0', 'apex'}
        for model = {'gapless', 'gapped'}
            d = dfit.(sprintf('b%d_%s_%s', j, meas{1}, model{1})); ls = ':'; if strcmp(model{1}, 'gapped'), ls = '-.'; end
            c = cl(1, :); if strcmp(meas{1}, 'apex'), c = cl(2, :); end
            plot(ax, qq * 1e3, d.f(d.p, qq), ls, 'Color', c, 'LineWidth', 1, 'DisplayName', sprintf('%s fit, %s (chi2/dof %.1f)', meas{1}, model{1}, d.chi2 / max(d.dof, 1)));
        end
    end
    grid(ax, 'on'); xlim(ax, [0 16]); ylim(ax, [0 1800]); xlabel(ax, 'Q (10^{-3} A^{-1})'); ylabel(ax, 'energy (meV)');
    nm = {'low branch', 'high branch'}; title(ax, sprintf('%s: E0 and DL apex', nm{j})); legend(ax, 'FontSize', 7, 'Location', 'southeast');
end
for j = [2 1]
    ax = nexttile(tl); hold(ax, 'on');
    for v = 2:numel(res), plot(ax, res(v).curves.Q * 1e3, res(v).curves.G(j, :), '-', 'Color', [0.75 0.75 0.75]); end
    fill(ax, [C.Q fliplr(C.Q)] * 1e3, [C.G(j, :) + C.sG(j, :), fliplr(C.G(j, :) - C.sG(j, :))], [0.2 0.6 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax, C.Q * 1e3, C.G(j, :), '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 2); plot(ax, C.Q * 1e3, C.E0(j, :), ':', 'Color', cl(1, :));
    grid(ax, 'on'); xlim(ax, [0 5.5]); xlabel(ax, 'Q (10^{-3} A^{-1})'); ylabel(ax, 'meV');
    nm = {'low branch', 'high branch'}; title(ax, sprintf('%s: width (green) against E0 (dotted)', nm{j}));
end
exportgraphics(f, fullfile(out, 'fig_p2b_dispersion.png'), 'Resolution', 110); close(f);
end

function local_fig_spectra(out, Eo, Y, S2, mu, cfg)
m = Eo >= 150 & Eo <= 1800; ni = size(Y, 2);
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 2000 1100]); tl = tiledlayout(f, 4, 5, 'TileSpacing', 'compact');
for i = 1:ni
    ax = nexttile(tl); yy = sum(Y(m, i, :), 3); mm = sum(mu(m, i, :), 3); top = max(yy(Eo(m) >= 250)) * 1.3;
    yyaxis(ax, 'left'); plot(ax, Eo(m), yy, 'k.', Eo(m), mm, 'r-', 'MarkerSize', 3); ylim(ax, [0 top]);
    yyaxis(ax, 'right'); plot(ax, Eo(m), (yy - mm) ./ sqrt(sum(S2(m, i, :), 3)), 'b-'); ylim(ax, [-10 10]);
    xline(ax, [150 250], 'k:'); grid(ax, 'on'); xlim(ax, [150 1800]); title(ax, sprintf('q = %+.4f', cfg.exit_q(i)), 'FontSize', 9);
end
exportgraphics(f, fullfile(out, 'fig_p2b_spectra.png'), 'Resolution', 100); close(f);
end

function local_fig_variants(out, res, intr)
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1600 700]); tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
cl = lines(numel(res)); nm = {'low branch', 'high branch'};
for j = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for v = 1:numel(res)
        F = res(v).fit; errorbar(ax, intr.Qn_E0 * 1e3 + (v - 4.5) * 0.02, F.E0(j, :), F.sE0(j, :), 'o-', 'Color', cl(v, :), 'MarkerSize', 4, 'DisplayName', res(v).name);
    end
    grid(ax, 'on'); set(ax, 'XScale', 'log'); xlabel(ax, 'Q node (10^{-3} A^{-1})'); ylabel(ax, 'E0 (meV)'); title(ax, sprintf('%s: E0 nodes, all variants', nm{j}));
    if j == 1, legend(ax, 'FontSize', 7, 'Location', 'best'); end
end
exportgraphics(f, fullfile(out, 'fig_p2b_variants.png'), 'Resolution', 110); close(f);
end
