function out = b1_p2a_operator(input_dir, response_dir)
% P2a step 2.2-2.5 on 590: build the block operators from the slice response,
% fit the reference amplitudes, and run the fixed-L forward diagnostics.
%  operator: incident slices = the 27 measured channels (weights = elastic-line
%   areas pi_br, sub-slice shape from the joint two-Pearson density) + tail slices
%   from the density beyond |q| = 0.00675; y: Gaussian sigma_y (main), symmetric
%   narrow Pearson, two-component 2D mixture (slot transmission corrected); exit
%   boxes [q_i +/- dq/2] x [-h, h]; radial grid Qm (qe_acceptance_operator);
%   energy kernels per slice from the line models (qe_elastic_line_model).
%  reference L: two DL branches with the E0 / width nodes of the config (from P1b
%   KB), residue coefficients fitted by NNLS to the 14 x 6 block spectra (300-1800
%   meV, Poisson weights) with the main operator (slice kernels).
%  diagnostics at fixed L: slice / integrated / exit kernels; y variants; broad
%   component on/off for B1; 6 versus 12 blocks; const versus gapless low-Q
%   extension; the old model (centred Gaussian beam, K averaged, L at |q_i|, exit
%   kernel, free scale); numerical convergence (sub-slice nodes, exit samples, y
%   nodes, radial grid, energy step, kernel support, tail range); Q-contribution
%   maps and the share inside the light cone Q < E/(hbar c); timings.
% Output: <response_dir>/operator_<config hash>/.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_p2a_operator_config(); tag=b1_sha256(cfg); t0=tic; perf=struct();
out=fullfile(response_dir,sprintf('operator_%s',tag(1:8))); if ~isfolder(out), mkdir(out); end
L=load(fullfile(input_dir,'p2a_cube.mat'),'cube'); cube=L.cube; E=cube.E(:); q=cube.q(:).'; dq=cube.dq;
RR=load(fullfile(response_dir,'response.mat'),'R','J6'); R=RR.R; J6=RR.J6;
ie=arrayfun(@(v) find(abs(q-v)<1e-9),cfg.exit_q); ni=numel(ie);
intr=qe_lowq_joint_forward('intrinsic'); th_ref=struct('E0sq',cfg.ref_E0.^2,'G',cfg.ref_G);
base=struct('blocking','b6','variant','slice','ytype','gauss','wide','on','ext','const','n_a',cfg.n_a, ...
 'n_sx',cfg.n_sx,'n_sy',cfg.n_sy,'n_gh',cfg.n_gh,'Qf',1,'step',cfg.omega(2),'support',cfg.support,'tail_max',cfg.tail_edges(end), ...
 'exits',1:ni,'blocks',[]);

% ---- main operator and reference amplitudes
tic; [ops,K,sw,ex]=b1_p2a_build_ops(base,cfg,cube,R,J6,q,dq,ie); perf.build_geometry_s=toc;
tic; cache=qe_lowq_joint_forward('build',ops,K,local_fopt(base,cfg,sw,ex)); perf.build_cache_s=toc;
w=whos('cache'); perf.cache_MB=w.bytes/2^20;
tic; cols=qe_lowq_joint_forward('eval',cache,th_ref,intr); perf.eval_s=toc;
[Y,S2,Eo]=local_data(cube,cache.E_out,ie);                         % data on the forward grid
m=Eo>=cfg.ref_fit(1)&Eo<=cfg.ref_fit(2); nb=size(Y,3); ncol=size(cols,2);
A=zeros(nnz(m)*ni*nb,ncol); y=zeros(nnz(m)*ni*nb,1); wv=y; o=0;
for b=1:nb, for i=1:ni
 rr=o+(1:nnz(m)); A(rr,:)=cols(m,:,i,b); y(rr)=Y(m,i,b); wv(rr)=1./sqrt(S2(m,i,b)); o=o+nnz(m);
end, end
sc=max(abs(A),[],1); sc(sc==0)=1; a_ref=lsqnonneg((A./sc).*wv,y.*wv)./sc(:);
mu_ref=local_combine(cols,a_ref); clear cache
fprintf('reference amplitudes: %s\n',mat2str(a_ref.',4));

% ---- fixed-L comparisons
cmp=struct('name',{},'mu',{},'cols',{},'note',{});
for v={'integrated','exit'}
 st=base; st.variant=v{1}; [mu,~,cc]=local_mu(st,cfg,cube,R,J6,q,dq,ie,th_ref,a_ref,intr);
 cmp(end+1)=struct('name',['kernel: ' v{1}],'mu',mu,'cols',cc,'note','same L, geometry, exposure'); %#ok<AGROW>
end
for v={'narrow_pearson','mixture'}
 st=base; st.ytype=v{1}; [mu,~,cc]=local_mu(st,cfg,cube,R,J6,q,dq,ie,th_ref,a_ref,intr);
 cmp(end+1)=struct('name',['y: ' v{1}],'mu',mu,'cols',cc,'note','same L; y distribution changed'); %#ok<AGROW>
end
st=base; st.wide='off'; [mu,~,cc]=local_mu(st,cfg,cube,R,J6,q,dq,ie,th_ref,a_ref,intr);
cmp(end+1)=struct('name','broad component off for B1','mu',mu,'cols',cc,'note','B1 sees the narrow component only');
st=base; st.ext='gapless'; [mu,~,cc]=local_mu(st,cfg,cube,R,J6,q,dq,ie,th_ref,a_ref,intr);
cmp(end+1)=struct('name','low-Q extension: gapless','mu',mu,'cols',cc,'note','E0^2 proportional to Q below 0.00025');
st=base; st.blocking='b12'; [mu12,~,c12]=local_mu(st,cfg,cube,R,J6,q,dq,ie,th_ref,a_ref,intr);
cmp(end+1)=struct('name','12 blocks (summed in pairs)','mu',mu12(:,:,1:2:end)+mu12(:,:,2:2:end),'cols',c12(:,:,:,1:2:end)+c12(:,:,:,2:2:end),'note','same L; time blocking');
mo=local_old(base,cfg,R,dq,ie,th_ref,a_ref,intr,Eo); mm=Eo>=300&Eo<=1800;
mo=mo.*(sum(mu_ref(mm,:,:),1)./sum(mo(mm,:,:),1));                 % free scale per spectrum (as the old single-spectrum fits)
cmp(end+1)=struct('name','old model (centred Gaussian, K averaged, L at |q|)','mu',mo,'cols',[],'note','scale per spectrum matched to the main model, 300-1800 meV');
% shared residues refitted for every variant (to the main-model prediction, data sigma):
% the part of a difference that the six linear amplitudes cannot absorb
for c=1:numel(cmp)
 if isempty(cmp(c).cols), cmp(c).mu_refit=cmp(c).mu; continue; end
 A2=zeros(nnz(m)*ni*nb,ncol); y2=zeros(nnz(m)*ni*nb,1); o=0;
 for b=1:nb, for i=1:ni
  rr=o+(1:nnz(m)); A2(rr,:)=cmp(c).cols(m,:,i,b); y2(rr)=mu_ref(m,i,b); o=o+nnz(m);
 end, end
 sc2=max(abs(A2),[],1); sc2(sc2==0)=1; a2=lsqnonneg((A2./sc2).*wv,y2.*wv)./sc2(:);
 cmp(c).mu_refit=local_combine(cmp(c).cols,a2); cmp(c).a_refit=a2;
end
reg=cfg.regions; crow={};
for c=1:numel(cmp)
 d=(cmp(c).mu-mu_ref)./sqrt(S2); d2=(cmp(c).mu_refit-mu_ref)./sqrt(S2);
 for i=1:ni
  for g=1:size(reg,1)
   mm=Eo>=reg(g,1)&Eo<=reg(g,2); x=d(mm,i,:); x2=d2(mm,i,:);
   crow(end+1,:)={cmp(c).name,cfg.exit_q(i),sprintf('%d-%d',reg(g,1),reg(g,2)),max(abs(x(:))),sqrt(mean(x(:).^2)),mean(x(:)),sum(x(:).^2), ...
    max(abs(x2(:))),sum(x2(:).^2),numel(x)}; %#ok<AGROW>
  end
 end
end
CT=cell2table(crow,'VariableNames',{'comparison','q','region_meV','max_abs_delta_sigma','rms_delta_sigma','mean_delta_sigma','delta_chi2', ...
 'max_abs_delta_sigma_after_residue_refit','delta_chi2_after_residue_refit','n_points'});
writetable(CT,fullfile(out,'fixedL_comparisons.csv'));
% apparent energy shift of each variant against the main model (300-1800 meV, all blocks summed)
srow={};
for c=1:numel(cmp)
 for i=1:ni
  a=squeeze(sum(mu_ref(:,i,:),3)); bq=squeeze(sum(cmp(c).mu(:,i,:),3)); mm=Eo>=300&Eo<=1800;
  srow(end+1,:)={cmp(c).name,cfg.exit_q(i),local_shift(Eo,a,bq,mm),sum(bq(mm))/sum(a(mm))}; %#ok<AGROW>
 end
end
writetable(cell2table(srow,'VariableNames',{'comparison','q','apparent_shift_meV','intensity_ratio_300_1800'}),fullfile(out,'fixedL_shifts.csv'));

% ---- numerical convergence (main model at fixed L)
conv={'sub-slice nodes 8','n_a',cfg.conv.n_a; 'exit samples 17x17','n_s',cfg.conv.n_s; 'y nodes 61','n_gh',cfg.conv.n_gh; ...
 'radial grid x2','Qf',cfg.conv.Qm_factor; 'kernel support 40','support',cfg.conv.support(1); 'kernel support 100','support',cfg.conv.support(2); ...
 'tail slices to 0.08','tail_max',cfg.conv.tail_max};
vrow={};
for c=1:size(conv,1)
 st=base;
 switch conv{c,2}
  case 'n_s', st.n_sx=conv{c,3}; st.n_sy=conv{c,3};
  otherwise, st.(conv{c,2})=conv{c,3};
 end
 tic; muc=local_mu(st,cfg,cube,R,J6,q,dq,ie,th_ref,a_ref,intr); tt=toc; d=(muc-mu_ref)./sqrt(S2);
 for g=1:size(reg,1)
  mm=Eo>=reg(g,1)&Eo<=reg(g,2); x=d(mm,:,:); vrow(end+1,:)={conv{c,1},sprintf('%d-%d',reg(g,1),reg(g,2)),max(abs(x(:))),sqrt(mean(x(:).^2)),tt}; %#ok<AGROW>
 end
end
% energy step 1 meV (subset: four exit channels, blocks 1 and 6), integrated over the 4 meV channels
st=base; st.step=cfg.conv.fine_step; st.exits=find(ismember(cfg.exit_q,[-0.002 -0.0005 0.0005 0.002])); st.blocks=[1 6];
tic; [muf,Ef]=local_mu(st,cfg,cube,R,J6,q,dq,ie,th_ref,a_ref,intr); tt=toc;
mub=zeros(numel(Eo),numel(st.exits),2);
for k=1:numel(Eo)
 j=find(abs(Ef-Eo(k))<=2+1e-9); if numel(j)~=5, mub(k,:,:)=NaN; continue; end
 mub(k,:,:)=sum(muf(j,:,:).*reshape([0.5 1 1 1 0.5]/4,[],1,1),1);
end
d=(mub-mu_ref(:,st.exits,[1 6]))./sqrt(S2(:,st.exits,[1 6]));
for g=1:size(reg,1)
 mm=Eo>=reg(g,1)&Eo<=reg(g,2); x=d(mm,:,:); vrow(end+1,:)={'energy step 1 meV, channel integrated (4 exits, 2 blocks)',sprintf('%d-%d',reg(g,1),reg(g,2)),max(abs(x(:)),[],'omitnan'),sqrt(mean(x(:).^2,'omitnan')),tt}; %#ok<AGROW>
end
writetable(cell2table(vrow,'VariableNames',{'refinement','region_meV','max_abs_delta_sigma','rms_delta_sigma','runtime_s'}),fullfile(out,'convergence.csv'));

% ---- Q-contribution maps and light-cone share (block 1, main operator)
st=base; [ops,~,~,~]=b1_p2a_build_ops(st,cfg,cube,R,J6,q,dq,ie); Qm=ops{1}.Qm;
wgrid=cfg.omega(1):cfg.omega(2):cfg.omega(3); C=qe_lowq_joint_forward('C',Qm,wgrid,th_ref,intr);
Ct=reshape(reshape(C,[],size(C,3))*a_ref,size(C,1),size(C,2));    % M x n_omega
hbarc=1973.27e3;                                                    % meV A
lrow={}; maps=struct('q',{},'E',{},'Q',{},'rho',{});
for i=1:ni
 Wi=reshape(sum(ops{1}.W(i,:,:),2),[],1);                           % M x 1 (all slices)
 rho=Wi.*Ct; tot=sum(rho,1);
 for e=[200 300 400 600 800 1000 1200 1500 1800]
  k=find(abs(wgrid-e)<1e-9); inside=Qm<e/hbarc; lrow(end+1,:)={cfg.exit_q(i),e,sum(rho(inside,k))/tot(k),local_quant(Qm,rho(:,k),0.5),local_quant(Qm,rho(:,k),0.1),local_quant(Qm,rho(:,k),0.9)}; %#ok<AGROW>
 end
 if ismember(cfg.exit_q(i),cfg.maps_q)
  for e=cfg.maps_E, k=find(abs(wgrid-e)<1e-9); maps(end+1)=struct('q',cfg.exit_q(i),'E',e,'Q',Qm,'rho',rho(:,k)/tot(k)); end %#ok<AGROW>
 end
end
LT=cell2table(lrow,'VariableNames',{'q','E_meV','share_inside_light_cone','Q_median','Q_p10','Q_p90'});
writetable(LT,fullfile(out,'q_contribution_and_light_cone.csv'));

% ---- figures
local_figs(out,cfg,Eo,Y,S2,mu_ref,cmp,maps,LT,ni);
perf.total_s=toc(t0);
save(fullfile(out,'operator_results.mat'),'cfg','a_ref','th_ref','mu_ref','cmp','Eo','Y','S2','perf','maps');
meta=struct('config',cfg,'config_sha256',tag,'code_sha256',b1_sha256({'file',[mfilename('fullpath') '.m']}), ...
 'operator_sha256',b1_sha256({'file',which('qe_acceptance_operator')}),'forward_sha256',b1_sha256({'file',which('qe_lowq_joint_forward')}), ...
 'response_dir',response_dir,'perf',perf,'reference_amplitudes',a_ref.', ...
 'created',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
fprintf('operator diagnostics written to %s (%.0f s)\n',out,perf.total_s);
end


% =============================== helpers ===============================
function fo = local_fopt(st, cfg, sw, ex)
fo = struct('omega', cfg.omega(1):st.step:cfg.omega(3), 'eps', -72:st.step:72, 'variant', st.variant, ...
    'slice_weight', {sw}, 'exit_slice', ex);
end

function [mu, Eout, cols] = local_mu(st, cfg, cube, R, J6, q, dq, ie, th, a, intr)
intr.ext = st.ext;
[ops, K, sw, ex] = b1_p2a_build_ops(st, cfg, cube, R, J6, q, dq, ie);
% exit_slice indexes rows of K: the measured slice equal to the exit channel
cache = qe_lowq_joint_forward('build', ops, K, local_fopt(st, cfg, sw, ex));
cols = qe_lowq_joint_forward('eval', cache, th, intr); Eout = cache.E_out(:); mu = local_combine(cols, a);
end

function mu = local_combine(cols, a)
mu = squeeze(sum(cols .* reshape(a, 1, [], 1, 1), 2));
if ndims(mu) == 2, mu = reshape(mu, size(cols, 1), size(cols, 3), size(cols, 4)); end
end

function [Y, S2, Eo] = local_data(cube, Eout, ie)
Eo = Eout(:); E = cube.E(:); J = cube.b6.J; nb = size(J, 3);
Y = nan(numel(Eo), numel(ie), nb);
[tf, loc] = ismember(round(Eo), round(E));
for b = 1:nb, Y(tf, :, b) = J(loc(tf), ie, b); end
S2 = max(Y, 1); S2(isnan(S2)) = Inf;
end

function mu = local_old(st, cfg, R, dq, ie, th, a, intr, Eo)
% centred Gaussian beam (sigma = old sigma_probe in x and y), one slice, exit kernel;
% K averaged over the acceptance, L at the nominal |q_i| (the P1/P1b forward model)
Rb = R.b6; nb = size(Rb.area, 1); Qm = b1_p2a_build_ops('Qm', cfg.Qm_spec, 1);
[gx, gw] = b1_p2a_build_ops('gh', 25); s = cfg.old_sigma; ax = sqrt(2) * s * gx; wx = gw / sqrt(pi);
op = qe_acceptance_operator(cfg.exit_q.', ax, wx, ones(numel(ax), 1), ax, wx, Qm, dq=dq, h=cfg.h, n_sx=st.n_sx, n_sy=st.n_sy);
w = cfg.omega(1):cfg.omega(2):cfg.omega(3); eps = -72:4:72; Kq = local_K(Qm, w, intr);
mu = zeros(numel(Eo), numel(ie), nb);
for i = 1:numel(ie)
    Wi = squeeze(op.W(i, 1, :)); Kbar = (Wi.' * Kq) / sum(Wi);
    Qi = abs(cfg.exit_q(i)); Ci = qe_lowq_joint_forward('C', Qi, w, th, intr); Ki = local_K(Qi, w, intr);
    G = Kbar .* (reshape(Ci, numel(w), []) * a).' ./ Ki;            % K-bar x S(|q_i|, w)
    for b = 1:nb
        k = qe_elastic_line_model('kernel', eps, Rb.th{b, ie(i)}, st.support) * 4;
        y = conv(G(:), k(:)); Ey = (w(1) + eps(1)) + 4 * (0:numel(y) - 1).';
        mu(:, i, b) = interp1(Ey, y, Eo, 'linear', 0);
    end
end
end

function K = local_K(Q, w, intr)
aw = max(abs(w(:).'), abs(w(2) - w(1)) / 4); qE = aw * intr.qE_per_meV;
K = Q(:) ./ (Q(:) .^ 2 + qE .^ 2) .^ 2; K(Q == 0, :) = 0;
end

function d = local_shift(E, a, b, m)
% energy shift of b against a (meV) by least squares on the window m
f = @(u) (interp1(E, a, E(m) - u(1), 'pchip') * u(2) - b(m)) ./ sqrt(max(b(m), 1));
u = lsqnonlin(f, [0 sum(b(m)) / sum(a(m))], [-50 0], [50 Inf], optimoptions('lsqnonlin', 'Display', 'off')); d = u(1);
end

function x = local_quant(Qm, rho, p)
c = cumsum(rho(:)); c = c / c(end); i = find(c >= p, 1, 'first'); x = Qm(i);
end

function local_figs(out, cfg, Eo, Y, S2, mu_ref, cmp, maps, LT, ni)
m = Eo >= 40 & Eo <= 1800; cl = lines(numel(cmp));
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1900 1050]); tl = tiledlayout(f, 3, 5, 'TileSpacing', 'compact');
for i = 1:ni
    ax = nexttile(tl); hold(ax, 'on');
    for c = 1:numel(cmp)
        d = sum(cmp(c).mu(m, i, :) - mu_ref(m, i, :), 3) ./ sqrt(sum(S2(m, i, :), 3));
        plot(ax, Eo(m), d, '-', 'Color', cl(c, :), 'DisplayName', cmp(c).name);
    end
    xline(ax, [150 250], 'k:'); grid(ax, 'on'); xlim(ax, [40 1800]); ylim(ax, [-30 30]);
    title(ax, sprintf('q = %+.4f', cfg.exit_q(i))); if i > 10, xlabel(ax, 'E (meV)'); end
    if mod(i, 5) == 1, ylabel(ax, '(variant - main) / \sigma, blocks summed'); end
end
ax = nexttile(tl); axis(ax, 'off'); hold(ax, 'on');
for c = 1:numel(cmp), plot(ax, NaN, NaN, '-', 'Color', cl(c, :), 'DisplayName', cmp(c).name); end
legend(ax, 'Location', 'west', 'FontSize', 8);
exportgraphics(f, fullfile(out, 'fig_fixedL_differences.png'), 'Resolution', 110); close(f);

f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1900 1050]); tl = tiledlayout(f, 3, 5, 'TileSpacing', 'compact');
for i = 1:ni
    ax = nexttile(tl); yy = sum(Y(m, i, :), 3); mm = sum(mu_ref(m, i, :), 3); top = max(yy(Eo(m) >= 250)) * 1.3;
    yyaxis(ax, 'left'); plot(ax, Eo(m), yy, 'k.', Eo(m), mm, 'r-', 'MarkerSize', 4); ylabel(ax, 'counts per 4 meV'); ylim(ax, [0 top]);
    yyaxis(ax, 'right'); plot(ax, Eo(m), (yy - mm) ./ sqrt(sum(S2(m, i, :), 3)), 'b-'); ylim(ax, [-40 40]);
    xline(ax, [150 250 300], 'k:'); grid(ax, 'on'); xlim(ax, [150 1800]); title(ax, sprintf('q = %+.4f: data (dots), reference model; residual / sigma', cfg.exit_q(i)));
end
exportgraphics(f, fullfile(out, 'fig_reference_model.png'), 'Resolution', 110); close(f);

f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1700 950]); tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on'); qs = unique(LT.q); cq = turbo(numel(qs));
for k = 1:numel(qs), x = LT(LT.q == qs(k), :); plot(ax, x.E_meV, x.share_inside_light_cone, 'o-', 'Color', cq(k, :), 'MarkerSize', 3, 'DisplayName', sprintf('%+.4f', qs(k))); end
grid(ax, 'on'); xlabel(ax, 'E (meV)'); ylabel(ax, 'share of W K S inside Q < E/(hbar c)', 'Interpreter', 'none'); legend(ax, 'FontSize', 7, 'NumColumns', 2, 'Location', 'northwest');
title(ax, '(a) retardation risk: weight inside the light cone');
for k = 1:min(3, numel(cfg.maps_q))
    ax = nexttile(tl); hold(ax, 'on'); mk = maps([maps.q] == cfg.maps_q(k)); ce = lines(numel(mk));
    for j = 1:numel(mk), plot(ax, mk(j).Q * 1e3, cumsum(mk(j).rho), '-', 'Color', ce(j, :), 'DisplayName', sprintf('E = %d meV', mk(j).E)); end
    xline(ax, abs(cfg.maps_q(k)) * 1e3, 'k:'); grid(ax, 'on'); xlim(ax, [0 6]); xlabel(ax, 'Q (10^{-3} A^{-1})'); ylabel(ax, 'cumulative share of counts');
    title(ax, sprintf('(%s) true-Q contributions, exit q = %+.4f', char('a' + k), cfg.maps_q(k))); legend(ax, 'FontSize', 7, 'Location', 'southeast');
end
exportgraphics(f, fullfile(out, 'fig_q_contributions.png'), 'Resolution', 110); close(f);
end
