function out = b1_p2a_synth(input_dir, response_dir, operator_dir)
% P2a step 2.4: injection and recovery with the 590 block operators.
% Generator: slice kernels, measured drift (six blocks), two-Pearson angular
% weights, Gaussian y; finer discretisation (cfg.gen). Three truths: fixed E0
% with Q-dependent widths and weights; a known dispersion; the reference L of the
% operator run (P1b-like). Each is fitted with the slice, integrated and exit
% kernels (main discretisation): 16 nonlinear parameters (log E0^2 at five nodes,
% log width at three nodes, two branches) by lsqnonlin, six residue coefficients
% by NNLS, weak curvature penalty on log E0^2, fit window 250-1800 meV, holdout
% 150-250 meV. Noise-free first (sigma^2 = expected counts), then three Poisson
% (Gaussian-approximated) draws of the reference truth fitted with slice and exit
% kernels. Reported: node values and DL apex energies at |q| = report_q against the
% truth, chi2 in the fit window and the holdout, runtime.
% Output: <operator_dir>/synth_<config hash>/.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_p2a_synth_config(); ocfg=b1_p2a_operator_config(); tag=b1_sha256(cfg); t0=tic;
out=fullfile(operator_dir,sprintf('synth_%s',tag(1:8))); if ~isfolder(out), mkdir(out); end
L=load(fullfile(input_dir,'p2a_cube.mat'),'cube'); cube=L.cube; q=cube.q(:).'; dq=cube.dq;
RR=load(fullfile(response_dir,'response.mat'),'R','J6'); R=RR.R; J6=RR.J6;
O=load(fullfile(operator_dir,'operator_results.mat'),'a_ref','th_ref'); a_ref=O.a_ref; th_ref=O.th_ref;
ie=arrayfun(@(v) find(abs(q-v)<1e-9),ocfg.exit_q); intr=qe_lowq_joint_forward('intrinsic');
base=struct('blocking','b6','variant','slice','ytype','gauss','wide','on','ext','const','n_a',ocfg.n_a, ...
 'n_sx',ocfg.n_sx,'n_sy',ocfg.n_sy,'n_gh',ocfg.n_gh,'Qf',1,'step',ocfg.omega(2),'support',ocfg.support, ...
 'tail_max',ocfg.tail_edges(end),'exits',1:numel(ie),'blocks',[]);
cls=cfg.classes;
for c=1:numel(cls)
 if isempty(cls(c).E0), cls(c).E0=sqrt(th_ref.E0sq); cls(c).G=th_ref.G; end
 if isempty(cls(c).a), cls(c).a=a_ref(:).'; end
end

% ---- generator (finer discretisation)
g=base; g.n_a=cfg.gen.n_a; g.n_sx=cfg.gen.n_s; g.n_sy=cfg.gen.n_s; g.n_gh=cfg.gen.n_gh; g.Qf=cfg.gen.Qf;
[ops,K,sw,ex]=b1_p2a_build_ops(g,ocfg,cube,R,J6,q,dq,ie);
gc=qe_lowq_joint_forward('build',ops,K,local_fopt(g,ocfg,sw,ex,'slice')); Eo=gc.E_out(:);
mu=cell(1,numel(cls));
for c=1:numel(cls)
 th=struct('E0sq',cls(c).E0.^2,'G',cls(c).G); mu{c}=local_combine(qe_lowq_joint_forward('eval',gc,th,intr),cls(c).a(:));
end
clear gc ops K
mfit=Eo>=cfg.fit_window(1)&Eo<=cfg.fit_window(2); mho=Eo>=cfg.holdout(1)&Eo<cfg.holdout(2);

% ---- fit caches (main discretisation)
fc=struct();
[ops,K,sw,ex]=b1_p2a_build_ops(base,ocfg,cube,R,J6,q,dq,ie);
for v=cfg.variants, fc.(v{1})=qe_lowq_joint_forward('build',ops,K,local_fopt(base,ocfg,sw,ex,v{1})); end
clear ops K

% ---- fits
rows={}; nodes={}; apex={};
opts=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10,'StepTolerance',1e-8,'OptimalityTolerance',1e-8, ...
 'MaxIterations',200,'MaxFunctionEvaluations',6000);
jobs={};
for c=1:numel(cls), for v=cfg.variants, jobs(end+1,:)={c,v{1},0,'expected'}; end, end %#ok<AGROW>
for wt=cfg.noise_weights
 for s=cfg.noise_seeds, for v=cfg.noise_variants, jobs(end+1,:)={cfg.noise_class,v{1},s,wt{1}}; end, end %#ok<AGROW>
end
jobs(cellfun(@isempty,jobs(:,end)),4)={'expected'};
for s=cfg.noise_seeds(1:3), jobs(end+1,:)={cfg.noise_class,'slice',s,[cfg.check.weights '-check']}; end %#ok<AGROW>
for jb=1:size(jobs,1)
 c=jobs{jb,1}; v=jobs{jb,2}; seed=jobs{jb,3}; wt=jobs{jb,4}; chk=endsWith(wt,'-check'); wt0=erase(wt,'-check');
 if seed==0
  Y=mu{c}; S2=max(mu{c},1);
 else
  rng(seed); Y=max(round(mu{c}+sqrt(max(mu{c},0)).*randn(size(mu{c}))),0);
  switch wt0
   case 'data', S2=max(Y,1);                                       % Neyman: biased towards low counts
   case 'smoothed', S2=max(movmean(Y,cfg.smooth_points,1),1);       % as the real noise model (smoothed variance)
   case 'oracle', S2=max(mu{c},1);                                  % expected counts
  end
 end
 th0=struct('E0sq',(cls(c).E0*cfg.start_perturb.E0).^2,'G',cls(c).G*cfg.start_perturb.G); o2=opts;
 if chk
  th0=struct('E0sq',(cls(c).E0*cfg.check.start_perturb(1)).^2,'G',cls(c).G*cfg.check.start_perturb(2));
  o2=optimoptions(opts,'FunctionTolerance',cfg.check.FunctionTolerance,'StepTolerance',cfg.check.StepTolerance,'MaxIterations',cfg.check.MaxIterations);
 end
 tt=tic; [th,a,info]=local_fitL(fc.(v),Y,S2,mfit,intr,th0,cfg,o2); rt=toc(tt);
 muf=local_combine(qe_lowq_joint_forward('eval',fc.(v),th,intr),a);
 chi_fit=sum(((Y(mfit,:,:)-muf(mfit,:,:)).^2)./S2(mfit,:,:),'all'); n_fit=nnz(Y(mfit,:,:)>=0);
 chi_ho=sum(((Y(mho,:,:)-muf(mho,:,:)).^2)./S2(mho,:,:),'all'); n_ho=nnz(Y(mho,:,:)>=0);
 dmax=max(abs(muf(mfit,:,:)-mu{c}(mfit,:,:))./sqrt(max(mu{c}(mfit,:,:),1)),[],'all');
 rows(end+1,:)={c,cls(c).name,v,seed,wt,chi_fit,chi_fit/(n_fit-22),chi_ho/n_ho,dmax,info.exitflag,info.iterations,rt}; %#ok<AGROW>
 E0f=sqrt(th.E0sq); E0t=cls(c).E0;
 for j=1:2
  for n=1:5, nodes(end+1,:)={c,v,seed,wt,j,'E0',intr.Qn_E0(n),E0t(j,n),E0f(j,n)}; end %#ok<AGROW>
  for n=1:3, nodes(end+1,:)={c,v,seed,wt,j,'G',intr.Qn_G(n),cls(c).G(j,n),th.G(j,n)}; end %#ok<AGROW>
  for qq=cfg.report_q
   [et,gt]=local_at(qq,struct('E0sq',E0t.^2,'G',cls(c).G),intr,j); [ef,gf]=local_at(qq,th,intr,j);
   apex(end+1,:)={c,v,seed,wt,j,qq,sqrt(et),local_apex(sqrt(et),gt),sqrt(ef),local_apex(sqrt(ef),gf),gt,gf}; %#ok<AGROW>
  end
 end
 fprintf('job %d/%d: class %d %s seed %d (%s weights): chi2/dof %.3f, holdout chi2/n %.3f, %.0f s\n',jb,size(jobs,1),c,v,seed,wt,chi_fit/(n_fit-22),chi_ho/n_ho,rt);
end
FT=cell2table(rows,'VariableNames',{'class','truth','fit_kernel','noise_seed','weights','chi2_fit','chi2dof_fit','chi2_per_point_holdout', ...
 'max_abs_forward_error_sigma','exitflag','iterations','runtime_s'});
writetable(FT,fullfile(out,'synth_fits.csv'));
NT=cell2table(nodes,'VariableNames',{'class','fit_kernel','noise_seed','weights','branch','parameter','Q_node','true_meV','fitted_meV'});
writetable(NT,fullfile(out,'synth_nodes.csv'));
AT=cell2table(apex,'VariableNames',{'class','fit_kernel','noise_seed','weights','branch','Q','E0_true','apex_true','E0_fit','apex_fit','G_true','G_fit'});
writetable(AT,fullfile(out,'synth_apex.csv'));
local_fig(out,AT,cfg);
meta=struct('config',cfg,'config_sha256',tag,'code_sha256',b1_sha256({'file',[mfilename('fullpath') '.m']}), ...
 'operator_dir',operator_dir,'runtime_s',toc(t0),'created',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
fprintf('synthetic tests written to %s (%.0f s)\n',out,toc(t0));
end


function fo = local_fopt(st, cfg, sw, ex, variant)
fo = struct('omega', cfg.omega(1):st.step:cfg.omega(3), 'eps', -72:st.step:72, 'variant', variant, ...
    'slice_weight', {sw}, 'exit_slice', ex);
end

function mu = local_combine(cols, a)
mu = squeeze(sum(cols .* reshape(a, 1, [], 1, 1), 2));
if ndims(mu) == 2, mu = reshape(mu, size(cols, 1), size(cols, 3), size(cols, 4)); end
end

function [th, a, info] = local_fitL(cache, Y, S2, m, intr, th0, cfg, opts)
lb = [log(cfg.bounds_E0(1, 1) ^ 2) * ones(1, 5), log(cfg.bounds_E0(2, 1) ^ 2) * ones(1, 5), log(cfg.bounds_G(1)) * ones(1, 6)];
ub = [log(cfg.bounds_E0(1, 2) ^ 2) * ones(1, 5), log(cfg.bounds_E0(2, 2) ^ 2) * ones(1, 5), log(cfg.bounds_G(2)) * ones(1, 6)];
u0 = [reshape(log(th0.E0sq).', 1, []), reshape(log(th0.G).', 1, [])]; u0 = min(max(u0, lb + 1e-6), ub - 1e-6);
y = Y(m, :, :); w = 1 ./ sqrt(S2(m, :, :)); y = y(:); w = w(:);
f = @(u) local_res(u, cache, y, w, m, intr, cfg.lambda_curv);
[u, ~, ~, flag, outp] = lsqnonlin(f, u0, lb, ub, opts);
[~, a] = f(u); th = local_decode(u);
info = struct('exitflag', flag, 'iterations', outp.iterations);
end

function [r, a] = local_res(u, cache, y, w, m, intr, lam)
th = local_decode(u); cols = qe_lowq_joint_forward('eval', cache, th, intr);
A = reshape(permute(cols(m, :, :, :), [1 3 4 2]), [], size(cols, 2));
sc = max(abs(A), [], 1); sc(sc == 0) = 1; a = lsqnonneg((A ./ sc) .* w, y .* w) ./ sc(:);
r = [(A * a - y) .* w; lam * reshape(diff(log(th.E0sq), 2, 2), [], 1)];
end

function th = local_decode(u)
th = struct('E0sq', exp(reshape(u(1:10), 5, 2).'), 'G', exp(reshape(u(11:16), 3, 2).'));
end

function [E0sq, G] = local_at(Q, th, intr, j)
lQ = log(Q); xn = log(intr.Qn_E0); E0sq = exp(interp1(xn, log(th.E0sq(j, :)), min(max(lQ, xn(1)), xn(end))));
xg = log(intr.Qn_G); G = exp(interp1(xg, log(th.G(j, :)), min(max(lQ, xg(1)), xg(end))));
end

function Em = local_apex(E0, G)
x = 2 * E0 .^ 2 - G .^ 2; Em = sqrt((x + sqrt(x .^ 2 + 12 * E0 .^ 4)) / 6);
end

function local_fig(out, AT, cfg)
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1800 900]); tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact');
mk = struct('slice', 'o', 'integrated', 's', 'exit', '^'); cl = lines(3);
for c = 1:3
    for j = 1:2
        ax = nexttile(tl, (j - 1) * 3 + c); hold(ax, 'on');
        x = AT(AT.class == c & AT.branch == j & AT.noise_seed == 0, :);
        if isempty(x), continue; end
        t = x(strcmp(x.fit_kernel, 'slice'), :); plot(ax, t.Q * 1e3, t.E0_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'true E0');
        plot(ax, t.Q * 1e3, t.apex_true, 'k--', 'LineWidth', 1.2, 'DisplayName', 'true apex');
        for v = 1:numel(cfg.variants)
            y = x(strcmp(x.fit_kernel, cfg.variants{v}), :);
            plot(ax, y.Q * 1e3, y.E0_fit, mk.(cfg.variants{v}), 'Color', cl(v, :), 'DisplayName', ['E0, ' cfg.variants{v} ' kernel']);
            plot(ax, y.Q * 1e3, y.apex_fit, [mk.(cfg.variants{v}) ':'], 'Color', cl(v, :), 'DisplayName', ['apex, ' cfg.variants{v}]);
        end
        n = x(x.noise_seed > 0, :);
        grid(ax, 'on'); xlabel(ax, 'Q (10^{-3} A^{-1})'); ylabel(ax, 'meV'); title(ax, sprintf('class %d, branch %d', c, j));
        if c == 1 && j == 1, legend(ax, 'FontSize', 7, 'Location', 'best'); end
    end
end
exportgraphics(f, fullfile(out, 'fig_synth_recovery.png'), 'Resolution', 110); close(f);
end
