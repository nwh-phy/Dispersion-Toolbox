function out = b1_p2a_escale(input_dir)
% P2a step 2.0c: where does the energy-angle coupling act? Energy shift of the
% +q spectrum against the -q spectrum, at equal |q|, from the A1-aligned block
% sums of p2a_cube.mat.
%   elastic window |E| <= 40 meV : should reproduce c(+q) - c(-q) (estimator check)
%   B1 windows 300-1800, 300-900, 900-1800 meV : inelastic energy reference
% Model: S+(E) = P(E) * S-~(E - delta), P a quadratic polynomial (multiplicative,
% absorbs the smooth +/-q intensity differences), S-~ the -q spectrum smoothed
% with a Gaussian (sigma 8 meV in the B1 windows, none in the elastic window).
% Weights: Poisson variance of both spectra. Uncertainty: jackknife over the six
% blocks. A coupling formed downstream of the sample (exit-trajectory energy
% error, image drift afterwards) shifts the inelastic spectra with the elastic
% line: delta_B1 = c(+q) - c(-q), about +12 meV at |q| >= 0.0015. A coupling in
% the incident beam leaves the inelastic energy reference common to all exit
% channels: delta_B1 near 0 (small, opposite-sign terms only where K*L varies
% strongly across the beam width).
% Output: <input_dir>/escale/.
L=load(fullfile(input_dir,'p2a_cube.mat'),'cube'); c=L.cube; E=c.E(:); q=c.q(:).'; J=c.b6.J; nb=size(J,3);
out=fullfile(input_dir,'escale'); if ~isfolder(out), mkdir(out); end
opts=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-12,'StepTolerance',1e-10,'MaxIterations',400);
win=struct('name',{'elastic |E|<=40','B1 300-1800','B1 300-900','B1 900-1800'},'r',{[-40 40],[300 1800],[300 900],[900 1800]}, ...
 'smooth',{0,8,8,8});
qa=0.0005:0.0005:max(q)+1e-9; rows={};
for a=qa
 kp=find(abs(q-a)<1e-9); km=find(abs(q+a)<1e-9); if isempty(kp)||isempty(km), continue; end
 cdiff=c.c_q(kp)-c.c_q(km);
 for w=1:numel(win)
  d_all=local_shift(sum(J(:,kp,:),3),sum(J(:,km,:),3),E,win(w),opts);
  d_jk=zeros(nb,1);
  for b=1:nb
   o=setdiff(1:nb,b); d_jk(b)=local_shift(sum(J(:,kp,o),3),sum(J(:,km,o),3),E,win(w),opts);
  end
  se=sqrt((nb-1)/nb*sum((d_jk-mean(d_jk)).^2));
  rows(end+1,:)={a,win(w).name,d_all,se,cdiff}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'abs_q','window','delta_plus_minus_meV','jackknife_se_meV','c_plus_minus_c_minus_meV'});
writetable(T,fullfile(out,'escale_shifts.csv'));

f=figure('Visible','off','Position',[0 0 1100 650]); ax=axes(f); hold(ax,'on'); mk={'ko','rs','b^','mv'};
for w=1:numel(win)
 m=strcmp(T.window,win(w).name);
 errorbar(ax,T.abs_q(m)+(w-2.5)*2e-5,T.delta_plus_minus_meV(m),T.jackknife_se_meV(m),mk{w},'DisplayName',win(w).name,'MarkerSize',5);
end
m=strcmp(T.window,win(1).name); plot(ax,T.abs_q(m),T.c_plus_minus_c_minus_meV(m),'k--','DisplayName','c(+q) - c(-q) from the elastic centres');
yline(ax,0,'k:'); grid(ax,'on'); xlabel(ax,'|q| (A^{-1})'); ylabel(ax,'energy shift of +q against -q (meV)');
legend(ax,'Location','best','FontSize',8);
title(ax,{'Energy reference of the inelastic signal: +q against -q at equal |q|', ...
 'downstream coupling: B1 follows the dashed line; coupling in the incident beam: B1 near 0'});
exportgraphics(f,fullfile(out,'escale_shifts.png'),'Resolution',120); close(f);
fprintf('energy-scale test written to %s\n',out);
disp(T);
end


function d = local_shift(Sp, Sm, E, win, opts)
% Shift delta (meV) and a quadratic multiplicative polynomial between the +q and
% -q spectra inside the window.
m = E >= win.r(1) & E <= win.r(2); x = (E(m) - mean(win.r)) / (diff(win.r) / 2);
if win.smooth > 0
    k = -4 * win.smooth:4:4 * win.smooth; g = exp(-k .^ 2 / (2 * win.smooth ^ 2)); g = g(:) / sum(g);
    St = conv(Sm, g, 'same');
else
    St = Sm;
end
y = Sp(m); vp = max(Sp(m), 1);
f = @(u) local_res(u, E, St, Sm, m, x, y, vp);
u0 = [0 log(sum(Sp(m)) / sum(Sm(m))) 0 0];
u = lsqnonlin(f, u0, [-30 -Inf -Inf -Inf], [30 Inf Inf Inf], opts); d = u(1);
end

function r = local_res(u, E, St, Sm, m, x, y, vp)
mod0 = interp1(E, St, E(m) - u(1), 'pchip');
P = exp(u(2)) * (1 + u(3) * x + u(4) * x .^ 2);
vm = interp1(E, max(Sm, 1), E(m) - u(1), 'linear');
r = (y - P .* mod0) ./ sqrt(vp + P .^ 2 .* vm);
end
