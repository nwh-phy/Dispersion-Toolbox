function b1_forward_bi_p2c_asymfit(P,out)
%B1_FORWARD_BI_P2C_ASYMFIT  Compare joint odd response and single-side DL fits.
P.seed=P.main;
P.window=[250 1800];
P.qmax=max(abs(P.c.exit_q));
B=b1_forward_bi_p2c_prepare(P,'odd',P.o.h);
modes={'none','gainE','shift','plus','minus'};
fits=cell(numel(modes),numel(P.held));
rows={}; fails=0; t0=tic;
for k=1:numel(modes)
    V=P;
    if k<=3, V.gain=modes{k}; else, V.gain='none'; V.side=9-2*k; end
    for j=1:numel(P.held)
        held=P.held(j);
        try
            t=tic;
            F=b1_forward_bi_p2c_fit(B,V,held,B.seed);
            fits{k,j}=F;
            save(fullfile(out,sprintf('fit_%s_block%d.mat',modes{k},held)),'F','V','held');
            th=NaN; d0=NaN; d1=NaN;
            if strcmp(V.gain,'shift')
                d0=F.gain(1)*1e-3; d1=F.gain(2); th=asind(d1);
            end
            rows(end+1,:)={modes{k},held,F.chi2,F.dof,F.chi2dof,F.held_chi2n,F.held_valley, ...
                d0,d1,th,F.flag,F.iterations,toc(t),mat2str(F.gain,7)}; %#ok<AGROW>
            fprintf('%s %s held %d: train %.5f held %.5f, flag %d, %.2fs\n',P.id,modes{k},held, ...
                F.chi2dof,F.held_chi2n,F.flag,toc(t));
            if k==1 && j==1
                fprintf('Pilot %.2f s optimizer + %.2f s geometry. 40 fits at this rate: %.2f h serial; low-window pilot will time its own profiles.\n', ...
                    toc(t),B.build_s,(40*toc(t)+4*B.build_s)/3600);
            end
        catch err
            fails=fails+1;
            fprintf(2,'%s held %d: %s\n',modes{k},held,getReport(err,'extended','hyperlinks','off'));
        end
    end
end
T=cell2table(rows,'VariableNames',{'model','held_block','chi2','dof','train_chi2dof','held_chi2n','held_250_300_chi2n', ...
    'delta0_Ainv','delta1','effective_angle_deg','exitflag','iterations','elapsed_s','odd_coef'});
writetable(T,fullfile(out,'asymmetry_fits.csv'));
assert(fails==0,'%d odd fit(s) failed.',fails);
local_dispersion(P,out,fits,modes);
local_residuals(B,P,out,fits,modes);
fprintf('%s odd fit with figures %.1f s\n',P.id,toc(t0));
end

function local_dispersion(P,out,fits,modes)
archive=load(fullfile(P.archive,sprintf('fit_%s_%s',P.id,P.main),'p2b_results.mat'),'res','cfg','hq');
dirs=[{'P2b'} modes];
f=figure('Visible','off','Theme','light','Position',[0 0 900 590]);
ax=axes(f); hold(ax,'on'); cl=lines(numel(dirs));
rows={};
for k=1:numel(dirs)
    cfg=archive.cfg; hq=archive.hq; res=archive.res(1);
    if k>1
        F=fits{k-1,1};
        res.curves=local_curve_grid(F,cfg);
        res.fit.E0=F.E0;
    end
    path=fullfile(out,['disp_' dirs{k}]);
    mkdir(path);
    save(fullfile(path,'p2b_results.mat'),'cfg','hq','res');
    Q=res.curves.Q; EE=res.curves.E0(1,:);
    ok=Q>=cfg.nodes.E0(1)-1e-10 & Q<=max(cfg.nodes.E0)+1e-10;
    plot(ax,Q(ok)*1e3,EE(ok),'LineWidth',1.6,'Color',cl(k,:),'DisplayName',dirs{k});
    for j=1:2
        dq=fullfile(path,sprintf('qcut_%d',j));
        v=struct('node_analysis',true,'make_figure',false,'output_dir',dq);
        if j==2, v.min_Q=.001; end
        b1_p2b_dispersion(path,v);
        D=readtable(fullfile(dq,'dispersion_details.csv'));
        d=D(strcmp(D.branch,'low')&strcmp(D.measure,'E0'),:);
        e=interp1(Q,EE,[.001 .0025]);
        rows(end+1,:)={dirs{k},j,e(1),e(2),d.Delta_meV,d.Delta_profile_lo_meV,d.Delta_profile_hi_meV, ...
            d.rho_lower_A,d.chi2_constant,d.chi2_gapless,d.chi2_gapped,d.n}; %#ok<AGROW>
    end
end
set(ax,'XScale','log'); xlabel(ax,'Q (10^{-3} Å^{-1})'); ylabel(ax,'Low-branch E_0 (meV)');
ylim(ax,[0 1850]); xx=logspace(log10(.12),log10(5),200);
plot(ax,xx,1973.27*xx,'k:','DisplayName','light line');
grid(ax,'on'); legend(ax,'Location','northwest'); title(ax,P.id,'Interpreter','none');
exportgraphics(f,fullfile(out,'fig_dispersion.png'),'Resolution',125); close(f);
T=cell2table(rows,'VariableNames',{'model','qcut','E0_at_001','E0_at_0025','Delta','Delta_lo','Delta_hi', ...
    'rho_lower','chi2_const','chi2_gapless','chi2_gapped','nodes'});
writetable(T,fullfile(out,'dispersion.csv'));
end

function C=local_curve_grid(F,cfg)
Q=unique([cfg.Q_curve cfg.nodes.E0]); q=log(Q); C=struct('Q',Q);
nn={'E0','G'};
for k=1:numel(nn)
    v=nn{k}; x=log(cfg.nodes.(v));
    xx=min(max(q,x(1)),x(end));
    C.(v)=exp(interp1(x,log(F.(v)).',xx).');
    C.(['s' v])=interp1(x,F.(['s' v]).',xx).';
end
E=C.E0; G=C.G;
a=@(e,g) sqrt((2*e.^2-g.^2+sqrt((2*e.^2-g.^2).^2+12*e.^4))/6);
C.apex=a(E,G);
h=.001;
C.sapex=hypot((a(E+h,G)-a(E-h,G))./(2*h).*C.sE0,(a(E,G+h)-a(E,G-h))./(2*h).*C.sG);
C.pole=sqrt(max(E.^2-G.^2/4,0));
C.spole=hypot(E.*C.sE0,G.*C.sG/4)./max(C.pole,1);
end

function local_residuals(B,P,out,fits,modes)
f=figure('Visible','off','Theme','light','Position',[0 0 1450 510]);
tiledlayout(f,1,3,'TileSpacing','compact'); e=B.E>=250 & B.E<=1800;
for k=1:3
    ax=nexttile;
    F=fits{k,1}; b=F.held;
    z=(B.Y(e,:,b)-F.mu(e,:,b))./sqrt(B.S2(e,:,b));
    b1_forward_bi_p2b_heatmap(ax,B.q,B.E(e),z,[-5 5]);
    title(ax,[modes{k} ' · held block 6']);
end
sgtitle(P.id,'Interpreter','none');
exportgraphics(f,fullfile(out,'fig_residuals.png'),'Resolution',125); close(f);
end
