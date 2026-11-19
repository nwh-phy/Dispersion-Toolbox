function b1_forward_bi_p2b_figures(T,chosen,S,out)
%B1_FORWARD_BI_P2B_FIGURES  Data/model residuals, h dependence and BiSb comparison.
for j = 1:numel(S)
    c = S{j}.cfg;
    ids = char(T.id(chosen(j)));
    data = S{j};
    train = setdiff(1:size(data.mainY,3),c.lobo_block);
    E = data.mainEo;
    m = E>=150 & E<=1800;
    Y = sum(data.mainY(m,:,train),3);
    mu = sum(data.mainMu(m,:,train),3);
    sigma = sqrt(sum(data.mainS2(m,:,train),3));
    z = (Y-mu)./sigma;
    figureFile = figure('Visible','off','Theme','light','Position',[0 0 1450 540]);
    lay = tiledlayout(figureFile,1,3,'TileSpacing','compact');
    x = {log10(max(Y,1)),log10(max(mu,1)),z};
    names = {'Data, log_{10} counts','Model, log_{10} counts','Data − model, σ'};
    lim = [min(x{1},[],'all','omitnan'),max(x{1},[],'all','omitnan')];
    for k = 1:3
        ax = nexttile(lay);
        level = lim;
        if k==3, level=[-8 8]; end
        b1_forward_bi_p2b_heatmap(ax,c.exit_q,E(m),x{k},level);
        title(ax,names{k});
        yline(ax,250,'k--');
    end
    title(lay,strrep(ids,'_','\_'));
    exportgraphics(figureFile,fullfile(out,['fig_' ids '_data_model_residual.png']),'Resolution',145);
    close(figureFile);
end
local_hscan(T,chosen,S,out);
local_overlay(T(chosen,:),S,out);
end

function local_hscan(T,chosen,S,out)
f = figure('Visible','off','Theme','light','Position',[0 0 1400 900]);
tl = tiledlayout(f,2,2,'TileSpacing','compact');
for j = 1:numel(S)
    ax = nexttile(tl);
    hold(ax,'on');
    k0 = chosen(j);
    isScan = T.id==T.id(k0) & startsWith(T.label,"hscan") & T.variant=="main";
    isScan(k0) = true;
    p = find(isScan);
    [~,order] = sort(T.h(p));
    p = p(order);
    cl = parula(numel(p)+1);
    for k = 1:numel(p)
        A = load(fullfile(T.folder(p(k)),'p2b_results.mat'),'res','cfg');
        C = A.res(1).curves;
        q = A.cfg.nodes.E0;
        y = interp1(C.Q,C.E0(1,:),q);
        width = 1.4+1.2*(p(k)==k0);
        plot(ax,1e3*q,y,'o-','Color',cl(k,:),'LineWidth',width,'MarkerSize',4, ...
            'DisplayName',sprintf('h = %.5g',T.h(p(k))));
    end
    set(ax,'XScale','log'); grid(ax,'on');
    xticks(ax,[.25 .5 1 1.5 2.5 3.5 5]);
    xtickformat(ax,'%.2g');
    title(ax,strrep(char(T.id(k0)),'_','\_'));
    xlabel(ax,'Q (10^{-3} Å^{-1})'); ylabel(ax,'Low-branch E_0 (meV)');
    legend(ax,'Location','best','FontSize',8);
end
exportgraphics(f,fullfile(out,'fig_hscan_low_E0.png'),'Resolution',155);
close(f);
end

function local_overlay(T,S,out)
base = b1_forward_config();
archive = {'fit_590_3fb982f2','fit_no_PL2_20w_2film_6be42333'};
names = cellstr(T.id);
names = [names(:).' {'BiSb 590 (1L)','BiSb 20w (2L)'}];
for k = 1:numel(archive)
    baseDir = fullfile(base.p2b_dir,archive{k});
    S{end+1} = load(fullfile(baseDir,'p2b_results.mat'),'cfg','res'); %#ok<AGROW>
end
h = [T.h(:).' .0066 .0033];
measure = {'E0','apex'};
titles = {'E_0','loss-function apex'};
branch = {'Low branch','High branch'};
f = figure('Visible','off','Theme','light','Position',[0 0 1450 980]);
tl = tiledlayout(f,2,2,'TileSpacing','compact');
cl = lines(numel(S));
for b = 1:2
    for k = 1:numel(measure)
        ax = nexttile(tl);
        hold(ax,'on');
        for j = 1:numel(S)
            C = S{j}.res(1).curves;
            q = S{j}.cfg.nodes.E0;
            y = interp1(C.Q,C.(measure{k})(b,:),min(max(q,min(C.Q)),max(C.Q)));
            use = q>=min(C.Q)-1e-12 & q<=max(C.Q)+1e-12;
            q=q(use); y=y(use);
            plot(ax,1e3*q,y,'o-','Color',cl(j,:),'LineWidth',1.5,'MarkerSize',4, ...
                'DisplayName',sprintf('%s, h=%.4g',strrep(names{j},'_','\_'),h(j)));
        end
        q = logspace(log10(.00019),log10(.006),200);
        plot(ax,1e3*q,1973e3*q,'k--','DisplayName','Light line');
        set(ax,'XScale','log'); grid(ax,'on'); xlim(ax,[.19 6]);
        xticks(ax,[.25 .5 1 1.5 2.5 3.5 5]);
        xtickformat(ax,'%.2g');
        if b==1, ylim(ax,[100 1350]); else, ylim(ax,[150 2575]); end
        xlabel(ax,'Q (10^{-3} Å^{-1})'); ylabel(ax,'Energy (meV)');
        title(ax,[branch{b} ' ' titles{k}]);
        if b==2 && k==2, legend(ax,'Location','best','FontSize',8); end
    end
end
exportgraphics(f,fullfile(out,'fig_Bi_BiSb_dispersion_overlay.png'),'Resolution',160);
close(f);
end
