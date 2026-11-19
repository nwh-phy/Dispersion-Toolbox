function b1_forward_20w_figure(A,B,out)
%B1_FORWARD_20W_FIGURE  Compare the 1L and 2L dispersion on log Q.
f = figure('Visible','off','Theme','light','Position',[0 0 1240 775]);
tl = tiledlayout(f,2,2,'TileSpacing','compact');
keys = {'low_E0','low_apex','high_E0','high_apex'};
titles = {'Low branch E_0','Low branch apex','High branch E_0','High branch apex'};
cl = [0.10 0.30 0.63;0.77 0.29 0.11];
for k = 1:numel(keys)
    ax = nexttile(tl);
    hold(ax,'on');
    for j = 1:2
        D = {A,B};
        C = D{j}.(keys{k});
        P = strcmp(C.src,'P2b');
        nm = {'590, 1L','20w, 2L'};
        hs = {'0.0066','0.0033'};
        errorbar(ax,1e3*C.x(P),C.y(P),C.s(P),'s-','Color',cl(j,:), ...
            'MarkerFaceColor',cl(j,:),'DisplayName',[nm{j} ', P2b h=' hs{j}]);
        errorbar(ax,1e3*C.x(~P),C.y(~P),C.s(~P),'o:','Color',cl(j,:), ...
            'DisplayName',[nm{j} ', bins h=0.0005']);
    end
    lim = [350 1100;350 1100;650 2300;650 2000];
    q = logspace(log10(0.00025),log10(0.016),120);
    plot(ax,q*1e3,1973e3*q,'k--','DisplayName','Light line');
    title(ax,titles{k}); ylabel(ax,'Energy (meV)');
    xlabel(ax,'Q (10^{-3} Å^{-1})');
    set(ax,'XScale','log'); grid(ax,'on');
    xlim(ax,[0.25 16]); ylim(ax,lim(k,:));
    if k==1, legend(ax,'Location','southeast','FontSize',8); end
end
name = fullfile(out,'b1_forward_20w_dispersion_590_overlay_labeled.png');
assert(~isfile(name),'Dispersion overlay already exists.');
exportgraphics(f,name,'Resolution',160);
close(f);
end
