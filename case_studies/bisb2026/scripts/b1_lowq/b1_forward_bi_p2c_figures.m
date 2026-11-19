function b1_forward_bi_p2c_figures(P,B,M,out)
%B1_FORWARD_BI_P2C_FIGURES  Native held-block spectra with model components.
e=B.E>=100 & B.E<=400;
m=M{1};
b=m.held;
keep=find(any(m.mask(:,:,b),1));
q=B.q(keep);
inn=min(abs(q));
targets=[inn -.0025 min(.004,P.qmax)];
[~,i]=min(abs(q(:)-targets),[],1);
cols=keep(i);
f=figure('Visible','off','Theme','light','Position',[0 0 1450 1000]);
tiledlayout(f,numel(cols),3,'TileSpacing','compact');
cl=lines(5);
nam={'Two DL branches','Interband','Drude'};
for j=cols
    for k=1:3
        ax=nexttile; hold(ax,'on');
        F=M{k};
        ee=B.E(e); y=B.Y(e,j,b);
        plot(ax,ee,y,'k.','MarkerSize',9,'DisplayName','Measured');
        plot(ax,ee,F.mu(e,j,b),'k-','LineWidth',1.5,'DisplayName','Total');
        plot(ax,ee,F.parts(e,j,b,1),'Color',cl(1,:),'DisplayName','Low DL');
        if size(F.parts,4)>1
            plot(ax,ee,F.parts(e,j,b,2),'Color',cl(2,:),'DisplayName','High DL');
        end
        plot(ax,ee,F.tail(e,j,b),'Color',cl(3,:),'DisplayName','Long tail');
        plot(ax,ee,F.phonon(e,j,b),'Color',cl(4,:),'DisplayName','130 meV');
        if k>1, plot(ax,ee,F.mat(e,j,b),'Color',cl(5,:),'DisplayName','Material'); end
        xline(ax,[P.window(1) 300],':','HandleVisibility','off');
        title(ax,sprintf('%s, q = %+.4f',nam{k},B.q(j)));
        xlabel(ax,'E (meV)'); ylabel(ax,'Counts / 4 meV');
        grid(ax,'on'); xlim(ax,[100 400]); set(ax,'YScale','log');
        ylim(ax,[max(.3,min(y(y>0))*.08) max(y)*2]);
        if j==cols(1) && k==2, legend(ax,'Location','northeast','FontSize',8); end
    end
end
sgtitle(sprintf('%s · held block %d',P.id,b),'Interpreter','none');
exportgraphics(f,fullfile(out,'fig_low_components.png'),'Resolution',125);
close(f);
end
