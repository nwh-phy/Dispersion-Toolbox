function b1_forward_20w_compare(cfg,out)
%B1_FORWARD_20W_COMPARE  Reuse dispersion fits with a common light-cone cut.
S = load(fullfile(out,'b1_forward_20w_fits.mat'));
assert(isempty(S.issues) && all(~cellfun(@isempty,S.fits)),'Incomplete h scan.');
P = load(fullfile(out,'b1_forward_20w_inputs.mat'),'p','o');
metrics = fullfile(P.p.input_dir,'b1_forward_20w_vs_590_inputs.csv');
if ~isfile(metrics)
    M = [b1_forward_20w_input_metrics(b1_forward_config());b1_forward_20w_input_metrics(P.p)];
    writetable(M,metrics);
end
fitDirs = [S.fits fullfile(cfg.paths.p2b_dir,cfg.baseline_fits)];
labels = [compose("20w h=%.4f",S.h), "590 h66", "590 historical"];
modes = ["combined","nodes"];
report = table();
scan = table();
scores = table();
paths = strings(numel(fitDirs),numel(modes));
for k = 1:numel(fitDirs)
    fit = load(fullfile(fitDirs{k},'p2b_results.mat'),'res','hq','cfg');
    if k<=numel(S.h)
        T = readtable(fullfile(fitDirs{k},'p2b_fits.csv'));
        V = readtable(fullfile(fitDirs{k},'b1_forward_20w_block6.csv'));
        assert(isequal(T.variant,V.variant),'Prediction and training variant order differ.');
        T.h = repmat(S.h(k),height(T),1);
        T.dchi2_within_h = T.chi2-T.chi2(1);
        T.block6_chi2_per_point = V.chi2_per_point;
        T.block6_low_energy_chi2_per_point = V.low_energy_chi2_per_point;
        scores = [scores;T]; %#ok<AGROW>
    end
    q = [0.001 0.0015 0.002 0.0025 0.0035 0.005];
    for v = 1:numel(fit.res)
        C = fit.res(v).curves;
        c = @(f,j) interp1(C.Q,C.(f)(j,:),min(q,C.Q(end))).';
        F = fit.res(v).fit;
        step = table(repmat(labels(k),numel(q),1),repmat(string(fit.res(v).name),numel(q),1),q.', ...
            c('E0',1),c('apex',1),c('E0',2),c('apex',2),c('G',1),c('G',2), ...
            repmat(F.chi2,numel(q),1),repmat(F.exitflag,numel(q),1), ...
            'VariableNames',{'fit','variant','Q','E0_low','apex_low','E0_high','apex_high', ...
            'G_low','G_high','chi2_train','exitflag'});
        scan = [scan;step]; %#ok<AGROW>
    end
    for mode = 1:numel(modes)
        folder = fullfile(out,sprintf('dispersion_%d_%s',k,modes(mode)));
        assert(~isfolder(folder),'Dispersion output already exists: %s',folder);
        mkdir(folder);
        res = struct('curves',{fit.res.curves});
        hq = fit.hq;
        fcfg = fit.cfg;
        fcfg.nodes.E0 = fcfg.nodes.E0(fcfg.nodes.E0>=cfg.q_min);
        if mode==2
            hq.q = hq.q([]);
            for name = {'E0','sE0','apex','sapex'}
                f = name{1};
                hq.(f) = hq.(f)(:,[]);
            end
        end
        % The existing dispersion postprocessor reads these three fields.
        % It uses cfg.nodes.E0 to choose evaluation points of the stored curves.
        local_stage(folder,res,hq,fcfg);
        paths(k,mode) = b1_p2b_dispersion(folder,struct('output_dir',folder));
        A = load(fullfile(folder,'dispersion.mat'),'D');
        T = readtable(fullfile(folder,'dispersion_summary.csv'));
        for r = 1:height(T)
            x = A.D.(sprintf('%s_%s',T.branch{r},T.measure{r}));
            T.chi2_const(r) = x.c.chi2;
            T.chi2_gapless(r) = x.g.chi2;
            T.chi2_gapped(r) = x.d.chi2;
        end
        T.fit = repmat(labels(k),height(T),1);
        T.points = repmat(modes(mode),height(T),1);
        report = [report;T]; %#ok<AGROW>
    end
end
writetable(report,fullfile(out,'b1_forward_20w_dispersion_summary.csv'));
writetable(scan,fullfile(out,'b1_forward_20w_hscan.csv'));
writetable(scores,fullfile(out,'b1_forward_20w_fit_scores.csv'));
B = load(fullfile(paths(1,1),'dispersion.mat'),'D');
A = load(fullfile(paths(end-1,1),'dispersion.mat'),'D');
b1_forward_20w_figure(A.D,B.D,out);
local_layers(cfg,A.D,B.D,out);
local_low_energy(P.p,out);
% Fixed best-fit L check of exit-y midpoint quadrature. Use the same h in
% all rows. The shared routine's legacy column names identify rows 2 and 3.
for k = 1:numel(S.h)
    p = P.p;
    p.fit_dir = S.fits{k};
    c = struct('paths',p,'output_dir',fullfile(S.fits{k},'b1_forward_20w_nsy'));
    assert(~isfolder(c.output_dir),'Convergence output exists.');
    h = S.h(k);
    if k==1
        n = P.o.n_sy;
    else
        n = max(9,2*ceil(h/2.5e-4)+1);
    end
    T = b1_p2b_h66_conv([],c,P.o,[h 9;h 9;h n;h 2*n-1]);
    fprintf('h %.4f: n_sy %d vs %d max %.6f sigma (sum of six blocks)\n', ...
        h,n,2*n-1,max(T.conv_nsy53));
end
end

function local_stage(folder,res,hq,cfg)
save(fullfile(folder,'p2b_results.mat'),'res','hq','cfg');
end

function local_layers(cfg,A,B,out)
rows = {};
for j = 1:2
    branch = {'low','high'};
    for measure = {'E0','apex'}
        key = [branch{j} '_' measure{1}];
        a = A.(key);
        b = B.(key);
        for source = {'P2b','binned'}
            sa = strcmp(a.src,source{1});
            sb = strcmp(b.src,source{1});
            if ~any(sa) || ~any(sb), continue; end
            lo = max([min(a.x(sa)) min(b.x(sb)) cfg.layer_qmin(j)]);
            hi = min([max(a.x(sa)) max(b.x(sb))]);
            q = b.x(sb & b.x>=lo-1e-10 & b.x<=hi+1e-10);
            if isempty(q), continue; end
            x1 = local_values(a.x(sa),a.y(sa),q);
            x2 = local_values(b.x(sb),b.y(sb),q);
            e1 = local_values(a.x(sa),a.s(sa),q);
            e2 = local_values(b.x(sb),b.s(sb),q);
            for i = 1:numel(q)
                rows(end+1,:) = {branch{j},measure{1},source{1},q(i), ...
                    x1(i),x2(i),x2(i)-x1(i),hypot(e1(i),e2(i))}; %#ok<AGROW>
            end
        end
    end
end
names = {'branch','measure','source','Q','energy_590','energy_20w','delta_20w_minus_590','err_scale_meV'};
writetable(cell2table(rows,'VariableNames',names),fullfile(out,'b1_forward_20w_layers.csv'));
base = b1_layer_compare(cfg.paths.run_root);
writetable(base,fullfile(out,'b1_forward_20w_layers_2d_h5.csv'));
end

function y = local_values(x,v,q)
if numel(x)==1
    y = repmat(v,size(q));
else
    y = interp1(x,v,min(max(q,min(x)),max(x)));
end
end

function local_low_energy(p,out)
L = load(fullfile(p.input_dir,'p2a_cube.mat'),'cube');
C = L.cube;
e = C.E(C.E>=40 & C.E<=250);
[~,ip] = min(abs(C.E-e.'),[],1);
[~,im] = min(abs(C.E+e.'),[],1);
% Gain-side mirror is a descriptive reference. It retains any elastic asymmetry.
loss = C.b6.J(ip,:,:);
gain = C.b6.J(im,:,:);
contrast = movmean((loss-gain)./max(gain,1),3,1);
avg = median(contrast,3);
pos = C.q>0 & C.q<=0.005;
neg = C.q<0 & C.q>=-0.005;
sg = avg(:,pos)-fliplr(avg(:,neg));
f = figure('Visible','off','Theme','light','Position',[0 0 1350 560]);
tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl);
imagesc(ax,C.q*1e3,e,avg); axis(ax,'xy'); colorbar(ax);
clim(ax,local_quantiles(avg,[0.02 0.98]));
xlabel(ax,'q (10^{-3} Å^{-1})'); ylabel(ax,'Energy loss (meV)');
title(ax,'Median block (loss − mirrored gain) / gain');
ax = nexttile(tl);
imagesc(ax,C.q(pos)*1e3,e,sg); axis(ax,'xy'); colorbar(ax);
lim = local_quantiles(abs(sg),0.98);
clim(ax,[-lim lim]);
xlabel(ax,'|q| (10^{-3} Å^{-1})'); ylabel(ax,'Energy loss (meV)');
title(ax,'Positive q minus negative q contrast');
exportgraphics(f,fullfile(out,'b1_forward_20w_40_250meV.png'),'Resolution',135);
close(f);
save(fullfile(out,'b1_forward_20w_40_250meV.mat'),'e','contrast','avg','sg');
end

function v = local_quantiles(x,p)
x = sort(x(:));
v = x(round(1+p*(numel(x)-1)));
end
