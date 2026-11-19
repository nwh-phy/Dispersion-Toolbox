function b1_forward_bi_p2a_stripe(E,Q,J,dose,C,cfg,out)
%B1_FORWARD_BI_P2A_STRIPE  Probe fixed-detector diagonal structures and weak columns.
idx = round(Q/C.dq);
[tf,mirror] = ismember(-idx,idx);
j = find(idx>=cfg.stripe_qmin_channels & tf);
k = idx(j);
plus = J(:,j,:); minus = J(:,mirror(j),:);
window = (-cfg.stripe_half_width:4:cfg.stripe_half_width).';
train = [1 3 5]; test = [2 4 6];
scan = zeros(numel(cfg.stripe_slope_search)*numel(cfg.stripe_anchor_search),5);
n = 0;
for slope = cfg.stripe_slope_search
    for anchor = cfg.stripe_anchor_search
        n = n+1;
        e = slope*k+anchor;
        eligible = e>=85 & e<=350;
        [~,~,z] = local_trace(E,plus,minus,e,window,cfg.stripe_sideband);
        discover = local_score(z(eligible,train));
        validate = local_score(z(eligible,test));
        scan(n,:) = [slope anchor nnz(eligible) discover validate];
    end
end
S = array2table(scan,'VariableNames',{'slope_meV_ch','intercept_meV','n_channels','odd_train_score','odd_test_score'});
writetable(S,fullfile(out,'b1_forward_bi_stripe_scan.csv'));
% Discovery is restricted to slopes near the independent slot observation.
sel = find(S.n_channels>=30 & S.slope_meV_ch>=1.25 & S.slope_meV_ch<=3.25);
[~,pick] = max(S.odd_train_score(sel)); best = sel(pick);
trials = [cfg.stripe_slope cfg.stripe_anchor; S{best,1:2}];
data = {};
cross = {};
for t = 1:size(trials,1)
    e = trials(t,1)*k+trials(t,2);
    [r,v,z,bg,rp,rm] = local_trace(E,plus,minus,e,window,cfg.stripe_sideband);
    [rg,vg,zg,bgg] = local_trace(E,minus,plus,-e,window,cfg.stripe_sideband);
    for n = 1:numel(j)
        v0 = sum(v(n,:)); vg0 = sum(vg(n,:));
        for b = 0:6
            if b==0
                bidx = 1:6;
                z0 = sum(r(n,:))/sqrt(v0); zg0 = sum(rg(n,:))/sqrt(vg0);
            else
                bidx = b; z0 = z(n,b); zg0 = zg(n,b);
            end
            r0 = sum(r(n,bidx)); g0 = sum(rg(n,bidx));
            data(end+1,:) = {t,b,k(n),Q(j(n)),e(n),r0,sum(rp(n,bidx)),sum(rm(n,bidx)), ...
                sum(bg(n,bidx)),r0/sum(bg(n,bidx)),z0,g0,g0/sum(bgg(n,bidx)),zg0, ...
                r0/sum(dose(bidx))}; %#ok<AGROW>
        end
    end
    side = cfg.stripe_sideband;
    p1 = sum(local_sample(E,plus,e-side),2);
    p2 = sum(local_sample(E,plus,e+side),2);
    m1 = sum(local_sample(E,minus,e-side),2);
    m2 = sum(local_sample(E,minus,e+side),2);
    for off = -24:4:24
        a = (1+off/side)/2;
        py = sum(local_sample(E,plus,e+off),2);
        my = sum(local_sample(E,minus,e+off),2);
        bgp = (1-a)*p1+a*p2;
        bgm = (1-a)*m1+a*m2;
        frac = (py-my-(bgp-bgm))./max(bgp,1);
        for band = [100 150;150 250;250 350].'
            i = e(:)>=band(1) & e(:)<band(2);
            cross(end+1,:) = {t,band(1),band(2),off,median(frac(i),'omitnan'), ...
                prctile(frac(i),[10 90])}; %#ok<AGROW>
        end
    end
end
trace = cell2table(data,'VariableNames',{'template','block','native_k','q','line_E_meV','odd_counts', ...
    'plus_ridge_counts','minus_ridge_counts','plus_local_bg_counts','odd_fraction_plus_bg','approx_z', ...
    'negative_E_mirror_odd_counts','negative_E_odd_fraction','negative_E_approx_z','odd_counts_per_dose'});
writetable(trace,fullfile(out,'b1_forward_bi_stripe_profiles.csv'));
X = cell2table(cross(:,1:5),'VariableNames',{'template','Emin','Emax','offset_meV','median_odd_fraction'});
quantiles = vertcat(cross{:,6});
X.p10 = quantiles(:,1); X.p90 = quantiles(:,2);
writetable(X,fullfile(out,'b1_forward_bi_stripe_cross_section.csv'));

% Energy-separated persistent column deficits or excesses. Skip central core
% curvature and the detector boundary when evaluating neighbor interpolation.
bands = [-300 -75; 85 220; 250 600; 650 1600];
rows = {};
for n = 3:numel(Q)-2
    if abs(idx(n))<7 || any(diff(idx(n-2:n+2))~=1), continue; end
    nb = [n-2 n-1 n+1 n+2];
    ratios = nan(4,6);
    for r = 1:4
        v = squeeze(sum(J(E>=bands(r,1)&E<=bands(r,2),[n nb],:),1));
        ratios(r,:) = v(1,:)./median(v(2:end,:),1);
    end
    u = median(ratios,2);
    f = (all(ratios([1 3 4],:)<cfg.bad_col_ratio(1),1) | ...
         all(ratios([1 3 4],:)>cfg.bad_col_ratio(2),1));
    rows(end+1,:) = {C.zero+idx(n),idx(n),Q(n),u(1),u(2),u(3),u(4),nnz(f),nnz(f)>=5}; %#ok<AGROW>
end
cols = cell2table(rows,'VariableNames',{'detector_col','native_k','q','gain_ratio','low_loss_ratio', ...
    'low_B1_ratio','high_B1_ratio','flag_blocks','persistent_candidate'});
writetable(cols,fullfile(out,'b1_forward_bi_columns.csv'));

a = sum(J,3);
p = a(:,j); m = a(:,mirror(j));
mask = abs(E)<=360;
odd = (p-m)./max(p+m,1);
oddHP = odd-movmean(odd,15,1);
f = figure('Visible','off','Theme','light','Position',[0 0 1400 800]);
tl = tiledlayout(f,2,2,'TileSpacing','compact');
ax = nexttile(tl);
imagesc(ax,k,E(mask),oddHP(mask,:)); axis(ax,'xy'); hold(ax,'on');
for t=1:2, plot(ax,k,trials(t,1)*k+trials(t,2),'k--'); end
ylim(ax,[-360 360]); clim(ax,[-.12 .12]); colorbar(ax);
xlabel(ax,'Detector channel from q=0'); ylabel(ax,'Energy (meV)'); title(ax,'High-pass (+q - mirror) / sum');
ax = nexttile(tl);
imagesc(ax,idx,E(mask),log10(max(a(mask,:),1))); axis(ax,'xy'); colorbar(ax);
xlabel(ax,'Detector channel from q=0'); ylabel(ax,'Energy (meV)'); title(ax,'Log summed native counts');
ax = nexttile(tl);
hold(ax,'on');
for t=1:2
    M=trace.template==t & trace.block==0 & trace.line_E_meV>=85 & trace.line_E_meV<=350;
    plot(ax,trace.line_E_meV(M),100*trace.odd_fraction_plus_bg(M),'DisplayName',sprintf('template %d',t));
end
yline(ax,0); legend(ax,'Location','best'); xlabel(ax,'Line energy (meV)'); ylabel(ax,'Odd ridge contrast (%)');
ax = nexttile(tl);
hold(ax,'on');
for b=1:6
    M=trace.template==1 & trace.block==b & trace.line_E_meV>=85 & trace.line_E_meV<=350;
    plot(ax,trace.line_E_meV(M),trace.odd_counts_per_dose(M));
end
yline(ax,0); xlabel(ax,'Line energy (meV)'); ylabel(ax,'Odd ridge counts / common dose');
exportgraphics(f,fullfile(out,'b1_forward_bi_stripe.png'),'Resolution',120);
close(f);

% Inspect the actual fit channels independently of the far-off-axis scan.
core = find(idx>0 & Q<=.005+1e-9 & tf);
pos = a(:,core);
neg = a(:,mirror(core));
odd = (pos-neg)./max(pos+neg,1);
hp = odd-movmean(odd,15,1);
g = E>=80 & E<=360;
f = figure('Visible','off','Theme','light','Position',[0 0 1050 500]);
ax = axes(f);
imagesc(ax,idx(core),E(g),hp(g,:)); axis(ax,'xy');
clim(ax,[-.15 .15]); colorbar(ax);
xlabel(ax,'Native |q| channel'); ylabel(ax,'Energy (meV)');
title(ax,'Low-q high-pass (+q - mirror) / sum');
exportgraphics(f,fullfile(out,'b1_forward_bi_stripe_core.png'),'Resolution',120);
close(f);
end

function [r,v,z,bg,rp,rm] = local_trace(E,P,M,e,win,side)
x = win+e;
p0 = local_sample(E,P,x); p1 = local_sample(E,P,x-side); p2 = local_sample(E,P,x+side);
m0 = local_sample(E,M,x); m1 = local_sample(E,M,x-side); m2 = local_sample(E,M,x+side);
bg = (p1+p2)/2;
rp = p0-bg; rm = m0-(m1+m2)/2;
r = rp-rm;
v = max(p0+m0+.25*(p1+p2+m1+m2),1);
z = r./sqrt(v);
end

function y = local_sample(E,A,x)
[nE,nQ,nB] = size(A);
j = (x-E(1))/(E(2)-E(1))+1;
bad = any(j<1 | j>=nE,1);
k = min(max(floor(j),1),nE-1);
f = j-k;
z = k+reshape((0:nQ*nB-1)*nE,1,nQ,nB);
y = reshape(sum(A(z).*(1-f)+A(z+1).*f,1),nQ,nB);
y(bad,:) = NaN;
end

function s = local_score(z)
% Trim out isolated large detector residuals. The score is descriptive
% because neighboring points and discovery grid cells are correlated.
if isempty(z), s = -Inf; return; end
s = median(z(:),'omitnan');
end
