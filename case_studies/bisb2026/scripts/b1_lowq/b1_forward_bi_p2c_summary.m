function b1_forward_bi_p2c_summary(run,label)
%B1_FORWARD_BI_P2C_SUMMARY  Assemble comparative scores and profile envelopes.
root=bisb_find_project_root(fileparts(mfilename('fullpath')));
run=fullfile(root,'paper_results','b1_forward_bi_20261119','p2c',run);
out=fullfile(run,label);
assert(~isfolder(out),'Summary already exists.');
mkdir(out);
ids={'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'};
P=dir(fullfile(run,'low_*'));
P=P([P.isdir] & ~startsWith({P.name},'low_pilot'));
if isfile(fullfile(run,'first_low_profile_superseded_utc.txt'))
    P=P(~strcmp({P.name},'low_primary'));
end
[~,order]=sort(string({P.name}),'descend');
P=P(order);
S=table(); U={}; Profiles=table();
for j=1:numel(ids)
    id=ids{j};
    seen={};
    for k=1:numel(P)
        folders=dir(fullfile(run,P(k).name,id));
        folders=folders([folders.isdir] & ~ismember({folders.name},{'.','..'}));
        for l=1:numel(folders)
            name=folders(l).name;
            folder=fullfile(folders(l).folder,name);
            f=fullfile(folder,'scores_6.csv');
            if ismember(name,seen) || ~isfile(f), continue; end
            seen{end+1}=name; %#ok<AGROW>
            M=local_tag(readtable(f,'TextType','string'),id,name,folder);
            s1=fullfile(folder,'scores_1.csv');
            S=[S;M]; %#ok<AGROW>
            if isfile(s1), S=[S;local_tag(readtable(s1,'TextType','string'),id,name,folder)]; end %#ok<AGROW>
            pr=local_tag(readtable(fullfile(folder,'profile_6.csv'),'TextType','string'),id,name,folder);
            Profiles=[Profiles;pr]; %#ok<AGROW>
            q=readtable(fullfile(folder,'limits_6.csv'),'TextType','string');
            for fam=["interband","drude"]
                qf=q(q.family==fam,:);
                L=pr(pr.family==fam & pr.refit,:);
                vf=M(M.family==fam,:);
                w=qf(qf.metric=="W",:);
                r=qf(qf.metric=="ratio",:);
                a=L(L.delta_scaled<=3.841459,:);
                null=vf.delta_objective/max(vf.train_chi2dof,1)<=3.841459;
                U(end+1,:)={id,name,char(fam),vf.h,vf.Elo,vf.n_q, ...
                    max(w.upper,[],'omitnan'),max(r.upper,[],'omitnan'),max(w.upper_common,[],'omitnan'), ...
                    max(r.upper_common,[],'omitnan'),all(w.closed)&&all(r.closed),vf.parameter,vf.width, ...
                    min(a.parameter,[],'omitnan'),max(a.parameter,[],'omitnan'),null,sum(a.flag==0),folder}; %#ok<AGROW>
            end
        end
    end
end
T=cell2table(U,'VariableNames',{'id','variant','family','h','Elo','n_q','W_upper','ratio_upper', ...
    'W_common_upper','ratio_common_upper','closed','best_parameter','best_width','prof_min','prof_max', ...
    'null_in_profile','flag0_within_profile','folder'});
writetable(S,fullfile(out,'scores_all.csv'));
writetable(T,fullfile(out,'limits_all.csv'));
writetable(Profiles,fullfile(out,'shape_profiles_all.csv'));
F=local_envelope(T,ids);
writetable(F,fullfile(out,'envelopes.csv'));
local_figures(out,F,Profiles,ids);
end

function T=local_tag(T,id,name,folder)
n=height(T);
T.id=repmat(string(id),n,1);
T.variant=repmat(string(name),n,1);
T.folder=repmat(string(folder),n,1);
end

function F=local_envelope(T,ids)
A=string(T.variant);
r={};
for j=1:numel(ids)
    id=ids{j};
    for fam={'interband','drude'}
        D=T(strcmp(T.id,id) & strcmp(T.family,fam{1}),:);
        n=D(strcmp(D.variant,'main'),:);
        if isempty(n), continue; end
        a=A(strcmp(T.id,id) & strcmp(T.family,fam{1}));
        bg=a=="main" | startsWith(a,"power") | a=="anchor" | startsWith(a,"no_phonon");
        vh=a=="main" | startsWith(a,"h_");
        eq=D.Elo==n.Elo; qeq=D.n_q==n.n_q;
        for kind={'nominal','background','h','same_window_and_mask','all_common'}
            switch kind{1}
                case 'nominal', m=a=="main";
                case 'background', m=bg;
                case 'h', m=vh;
                case 'same_window_and_mask', m=eq & qeq;
                case 'all_common', m=D.Elo<300;
            end
            X=D(m,:);
            [w,iw]=max(X.W_upper); [v,iv]=max(X.ratio_upper);
            [z,iz]=max(X.ratio_common_upper);
            if strcmp(kind{1},'all_common')
                w=NaN; v=NaN;
            end
            r(end+1,:)={id,fam{1},kind{1},height(X),n.Elo,w,v,max(X.W_common_upper),z, ...
                X.variant{iw},X.variant{iv},X.variant{iz},min(X.prof_min),max(X.prof_max), ...
                any(X.null_in_profile),all(X.closed)}; %#ok<AGROW>
        end
    end
end
F=cell2table(r,'VariableNames',{'id','family','envelope','n_variants','Elo','W_upper','ratio_upper','W_common_upper', ...
    'ratio_common_upper','W_max_variant','ratio_max_variant','common_ratio_max_variant','prof_min','prof_max','any_null','all_closed'});
end

function local_figures(out,F,P,ids)
f=figure('Visible','off','Theme','light','Position',[0 0 1350 675]);
tiledlayout(f,1,2,'TileSpacing','compact');
legendNames={'Nominal','Tail and 130 meV','Same window / channels'};
groups={'nominal','background','same_window_and_mask'};
col=lines(numel(groups));
for k=1:2
    ax=nexttile; hold(ax,'on'); fam={'interband','drude'};
    a=strcmp(F.family,fam{k});
    for j=1:numel(groups)
        v=F(a & strcmp(F.envelope,groups{j}),:);
        [~,idx]=ismember(v.id,ids);
        plot(ax,idx+(j-2)*.13,100*v.ratio_common_upper,'o','MarkerSize',8, ...
            'LineWidth',1.7,'Color',col(j,:),'DisplayName',legendNames{j});
    end
    set(ax,'YScale','log','XTick',1:4,'XTickLabel',ids,'TickLabelInterpreter','none');
    xtickangle(ax,20); grid(ax,'on'); xlim(ax,[.55 4.45]);
    ylabel(ax,'Conditional upper weight / low DL (%)'); title(ax,[fam{k} '   250-300 meV']);
    legend(ax,'Location','best');
end
exportgraphics(f,fullfile(out,'fig_W_upper.png'),'Resolution',130); close(f);

f=figure('Visible','off','Theme','light','Position',[0 0 1250 860]);
tiledlayout(f,2,2,'TileSpacing','compact');
for j=1:numel(ids)
    ax=nexttile; hold(ax,'on');
    u=P(P.id==string(ids{j}) & P.family=="interband" & P.refit,:);
    vs=unique(u.variant);
    nominal="main";
    if any(vs=="threshold_refined"), nominal="threshold_refined"; end
    for k=1:numel(vs)
        v=u(u.variant==vs(k),:);
        g=groupsummary(v,'parameter','min','delta_scaled');
        if vs(k)==nominal
            plot(ax,g.parameter,g.min_delta_scaled,'o-','LineWidth',2,'Color',[.08 .32 .65], ...
                'DisplayName','Nominal');
        else
            plot(ax,g.parameter,g.min_delta_scaled,'.-','Color',[.65 .69 .74],'HandleVisibility','off');
        end
    end
    yline(ax,3.841459,'--','DisplayName','95% grid criterion');
    grid(ax,'on'); xlim(ax,[0 800]); ylim(ax,[-.2 30]);
    title(ax,ids{j},'Interpreter','none'); xlabel(ax,'Interband threshold (meV)'); ylabel(ax,'Scaled profile difference');
    if j==1, legend(ax,'Location','northwest'); end
end
exportgraphics(f,fullfile(out,'fig_threshold_profiles.png'),'Resolution',130);
close(f);
end
