function out = b1_forward_bi_p2b_compare(phase,reference)
%B1_FORWARD_BI_P2B_COMPARE  Assemble dispersion, held predictions and geometry sensitivity.
% A final comparison may reuse validated cross-prediction and quadrature
% files from reference when the selected main fits are identical.
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
c = b1_forward_bi_p2b_config('Bi1107_10w');
p = c.paths.p2b_dir;
ids = ["Bi1107_10w","Bi1107_20w","Bi1117_10w","Bi1117_20w"];
T = b1_forward_bi_p2b_catalog(p);
assert(height(T)>=49,'Expected all 49 first-tier fits.');
tag = b1_sha256(struct('phase',phase,'fits',{cellstr(T.sha256)}));
out = fullfile(p,['compare_' tag(1:8)]);
assert(~isfolder(out),'Comparison directory already exists: %s',out);
mkdir(out);
writetable(T,fullfile(out,'fit_inventory.csv'));
S = cell(size(ids));
chosen = zeros(size(ids));
t = tic;
for j = 1:numel(ids)
    candidates = find(T.id==ids(j) & (T.label=="main" | startsWith(T.label,"start-")));
    [~,k] = min(T.train_chi2(candidates));
    chosen(j) = candidates(k);
    S{j} = load(fullfile(T.folder(chosen(j)),'p2b_results.mat'));
end
writetable(T(chosen,:),fullfile(out,'selected_main_fits.csv'));
tasks = {@() local_dispersion(T,chosen,out), @() local_transfer(T(chosen,:),S,out), ...
    @() local_nodes(T,chosen,S,out), @() b1_forward_bi_p2b_figures(T,chosen,S,out), ...
    @() local_convergence(T(chosen,:),S,out)};
names = {'dispersion','cross prediction','nodes and differences','summary figures','quadrature'};
if nargin>1
    P=readtable(fullfile(reference,'selected_main_fits.csv'),'TextType','string');
    assert(isequal(P.sha256,T.sha256(chosen(:))), ...
        'Reusing diagnostics requires identical selected main fits.');
    tasks{2}=@() local_copy(reference,out,{'cross_*.mat','cross_scores.csv','fig_1117_cross_prediction.png'});
    tasks{5}=@() local_copy(reference,out,{'nsy_convergence.csv'});
end
if strcmp(phase,'final')
    tasks{end+1} = @() local_twelve(T.folder(chosen(1)),S{1},out);
    names{end+1} = '12-block drift';
end
failed = {};
for k = 1:numel(tasks)
    try
        tk = tic;
        fprintf('START %s\n',names{k});
        tasks{k}();
        fprintf('DONE %s %.2f s\n',names{k},toc(tk));
    catch err
        failed{end+1} = names{k}; %#ok<AGROW>
        fprintf(2,'FAILED %s: %s\n',names{k},getReport(err,'extended','hyperlinks','off'));
    end
end
assert(isempty(failed),'Failed comparison tasks: %s.',strjoin(failed,', '));
fprintf('Compare %s complete %.2f s %s\n',phase,toc(t),out);
end

function local_copy(reference,out,patterns)
for k=1:numel(patterns)
    files=dir(fullfile(reference,patterns{k}));
    assert(~isempty(files),'Missing validated file %s.',patterns{k});
    for j=1:numel(files)
        f=files(j);
        target=fullfile(out,f.name);
        assert(~isfile(target),'Comparison file already exists: %s',target);
        copyfile(fullfile(f.folder,f.name),target);
    end
end
end

function local_dispersion(T,chosen,out)
parts = cell(height(T)+numel(chosen),1);
count = 0;
for k = 1:height(T)
    if T.branches(k)~=2, continue; end
    try
        id = char(T.id(k));
        path = fullfile(out,['dispersion_' id '_' char(extractBefore(T.sha256(k),9))]);
        cfg = struct('output_dir',path,'node_analysis',true,'make_figure',false);
        b1_p2b_dispersion(char(T.folder(k)),cfg);
        count = count+1;
        parts{count} = local_details(path,T(k,:),"individual");
    catch err
        fprintf(2,'Dispersion %s %s failed: %s\n',T.id(k),T.label(k),getReport(err,'extended','hyperlinks','off'));
    end
end
expected = nnz(T.branches==2);
assert(count==expected && all(~cellfun(@isempty,parts(1:count))),'Some per-run dispersions failed.');
for k = chosen
    variant = T.id==T.id(k) & ismember(T.label, ...
        ["y: Gaussian","broad component off for B1","kernel: exit","low-Q extension: gapless","empirical"]);
    tag = char(extractBefore(T.sha256(k),9));
    path = fullfile(out,['dispersion_' char(T.id(k)) '_' tag '_systematics']);
    cfg = struct('output_dir',path,'node_analysis',true,'variant_dirs',{cellstr(T.folder(variant))});
    b1_p2b_dispersion(char(T.folder(k)),cfg);
    count = count+1;
    parts{count} = local_details(path,T(k,:),"instrument spreads");
end
D = vertcat(parts{1:count});
writetable(D,fullfile(out,'dispersion_all.csv'));
local_envelopes(D,T(chosen,:),out);
local_reference_dispersion(T(chosen,:),out);
end

function D = local_details(path,row,kind)
D = readtable(fullfile(path,'dispersion_details.csv'),'TextType','string');
N = height(D);
meta = repmat(row(:,{'id','label','variant','start','sha256','h'}),N,1);
meta.kind = repmat(kind,N,1);
meta.dispersion_folder = repmat(string(path),N,1);
D = [meta D];
end

function local_reference_dispersion(T,out)
base = b1_forward_config();
locations = [T.folder;string(fullfile(base.p2b_dir,'fit_590_3fb982f2')); ...
    string(fullfile(base.p2b_dir,'fit_no_PL2_20w_2film_6be42333'))];
ids = [T.id;"BiSb_590_h66";"BiSb_20w_2L"];
parts = cell(numel(ids),1);
for k = 1:numel(ids)
    path = fullfile(out,['dispersion_commonQ_' char(ids(k))]);
    cfg = struct('output_dir',path,'node_analysis',true,'min_Q',.001,'make_figure',false);
    b1_p2b_dispersion(char(locations(k)),cfg);
    D = readtable(fullfile(path,'dispersion_details.csv'),'TextType','string');
    D = addvars(D,repmat(ids(k),height(D),1),'Before',1,'NewVariableNames','id');
    parts{k} = D;
end
writetable(vertcat(parts{:}),fullfile(out,'dispersion_commonQ.csv'));
end

function local_envelopes(D,T,out)
rows = cell(0,10);
for k = 1:height(T)
    d = D(D.id==T.id(k)&D.kind=="individual"&D.branch=="low",:);
    primary = d.sha256==T.sha256(k);
    hscan = startsWith(d.label,"hscan") & d.variant=="main";
    gaussian = d.variant=="y: Gaussian";
    nuisance = ismember(d.label,["y: Gaussian","broad component off for B1", ...
        "kernel: exit","low-Q extension: gapless","empirical"]);
    sets = {primary|hscan,primary|nuisance,primary|hscan|nuisance,gaussian};
    strata = ["main h scan","nominal instrument variants","main h and nominal variants","Gaussian h scan"];
    for m = ["E0","apex"]
        for j = 1:numel(sets)
            use = sets{j} & d.measure==m & isfinite(d.Delta_profile_lo_meV) & isfinite(d.Delta_profile_hi_meV);
            B = d(use,:);
            if isempty(B), continue; end
            rows(end+1,:) = {T.id(k),m,strata(j),height(B), ...
                min(B.Delta_meV),max(B.Delta_meV),min(B.Delta_profile_lo_meV),max(B.Delta_profile_hi_meV), ...
                min(B.rho_lower_A,[],'omitnan'),max(B.rho_lower_A,[],'omitnan')}; %#ok<AGROW>
        end
    end
end
E = cell2table(rows,'VariableNames',{'id','measure','variants','n_fits','Delta_estimate_min', ...
    'Delta_estimate_max','Delta_profile_envelope_lo','Delta_profile_envelope_hi','rho_lower_min','rho_lower_max'});
writetable(E,fullfile(out,'dispersion_sensitivity.csv'));
end

function local_transfer(T,S,out)
rows = cell(2,1);
Q = max(S{3}.cfg.nodes.E0(1),S{4}.cfg.nodes.E0(1));
f = figure('Visible','off','Theme','light','Position',[0 0 1500 870]);
lay = tiledlayout(f,2,3,'TileSpacing','compact');
for j = 3:4
    donor = 7-j;
    c = S{j}.cfg;
    o = local_operator(T.folder(j));
    M = b1_forward_bi_p2b_material(S{donor},Q);
    P = b1_p2b_fit([],c,o,struct('material',M));
    F = S{j}.res(1).fit;
    X = b1_forward_bi_p2b_transfer(P,F,c);
    assert(abs(X.table.own_chi2(1)-F.chi2)<1e-6, ...
        'Recipient training score does not match archive.');
    assert(abs(X.table.own_chi2dof(2)-F.prediction.chi2dof)<1e-10, ...
        'Recipient held score does not match archive.');
    n = height(X.table);
    header = table(repmat(T.id(donor),n,1),repmat(T.id(j),n,1),repmat(X.scale,n,1), ...
        repmat(Q,n,1),repmat(T.sha256(donor),n,1),repmat(T.sha256(j),n,1), ...
        'VariableNames',{'donor','recipient','scale','common_Qmin','donor_sha256','recipient_sha256'});
    rows{j-2} = [header X.table];
    name = ['cross_' char(T.id(donor)) '_to_' char(T.id(j))];
    save(fullfile(out,[name '.mat']),'P','X','M','c','o');
    ax = nexttile(lay,3*(j-3)+1);
    hold(ax,'on');
    curve = S{j}.res(1).curves;
    from = S{donor}.res(1).curves;
    branch = {'low','high'};
    for b = 1:2
        y = curve.Q>=Q & curve.Q<=.005;
        z = from.Q>=Q & from.Q<=.005;
        plot(ax,1e3*curve.Q(y),curve.E0(b,y),'b-','LineWidth',1.6, ...
            'DisplayName',['Recipient ' branch{b}]);
        plot(ax,1e3*from.Q(z),from.E0(b,z),'r--','LineWidth',1.6, ...
            'DisplayName',['Donor ' branch{b}]);
    end
    set(ax,'XScale','log');
    xticks(ax,[.5 1 1.5 2.5 3.5 5]);
    xtickformat(ax,'%.2g');
    title(ax,strrep([char(T.id(donor)) ' → ' char(T.id(j))],'_','\_'));
    xlabel(ax,'Q (10^{-3} Å^{-1})'); ylabel(ax,'E_0 (meV)'); grid(ax,'on');
    legend(ax,'Location','best');
    g = P.E>=c.variants(1).window(1) & P.E<=c.variants(1).window(2);
    ax = nexttile(lay,3*(j-3)+2);
    b1_forward_bi_p2b_heatmap(ax,c.exit_q,P.E(g),X.own_z(g,:,c.lobo_block),[-6 6]);
    title(ax,sprintf('Own held block, χ²/n %.3f',X.table.own_chi2dof(2)));
    ax = nexttile(lay,3*(j-3)+3);
    b1_forward_bi_p2b_heatmap(ax,c.exit_q,P.E(g),X.z(g,:,c.lobo_block),[-6 6]);
    title(ax,sprintf('Cross held block, χ²/n %.3f',X.table.cross_chi2dof(2)));
end
exportgraphics(f,fullfile(out,'fig_1117_cross_prediction.png'),'Resolution',150);
close(f);
writetable(vertcat(rows{:}),fullfile(out,'cross_scores.csv'));
end

function local_nodes(T,chosen,S,out)
rows = cell(0,14);
for k = 1:height(T)
    A = load(fullfile(T.folder(k),'p2b_results.mat'),'cfg','res');
    C = A.res(1).curves;
    F = A.res(1).fit;
    q = A.cfg.nodes.E0;
    a = reshape(F.a,numel(A.cfg.nodes.A),[]).';
    for b = 1:size(F.E0,1)
        gamma = local_at(C.Q,C.G(b,:),q);
        apex = local_at(C.Q,C.apex(b,:),q);
        sapex = local_at(C.Q,C.sapex(b,:),q);
        aq = local_at(A.cfg.nodes.A,a(b,:),q).*(q/.001);
        for j = 1:numel(q)
            rows(end+1,:) = {T.id(k),T.label(k),T.sha256(k),T.h(k),b,q(j),F.E0(b,j),F.sE0(b,j), ...
                apex(j),sapex(j),gamma(j),aq(j),F.E0(b,j)>=A.cfg.bounds_E0(b,2)-1, ...
                F.E0(b,j)<=A.cfg.bounds_E0(b,1)+1}; %#ok<AGROW>
        end
    end
end
N = cell2table(rows,'VariableNames',{'id','label','sha256','h','branch','Q','E0','sE0','apex','sapex', ...
    'Gamma','A_at_Q','E0_at_upper','E0_at_lower'});
writetable(N,fullfile(out,'all_nodes.csv'));
q = S{1}.cfg.nodes.E0;
E = nan(3,numel(q));
Ap = E;
St = E;
j = 0;
for i = [1 3 4]
    j = j+1;
    C = S{i}.res(1).curves;
    E(j,:) = local_at(C.Q,C.E0(1,:),q);
    Ap(j,:) = local_at(C.Q,C.apex(1,:),q);
    St(j,:) = local_at(C.Q,C.sE0(1,:),q);
end
R = table(q.',E(1,:).',E(2,:).',E(3,:).',(E(1,:)-E(2,:)).',(E(1,:)-E(3,:)).', ...
    hypot(St(1,:),St(2,:)).',hypot(St(1,:),St(3,:)).',(Ap(1,:)-Ap(2,:)).',(Ap(1,:)-Ap(3,:)).', ...
    'VariableNames',{'Q','A10_E0','same10_E0','same20_E0','A_minus_same10_E0','A_minus_same20_E0', ...
    'A_minus_same10_conditional_SE','A_minus_same20_conditional_SE','A_minus_same10_apex','A_minus_same20_apex'});
writetable(R,fullfile(out,'area_differences.csv'));
R = local_paired_h(T,N,chosen,q);
writetable(R,fullfile(out,'area_differences_paired_h.csv'));
end

function R = local_paired_h(T,N,chosen,q)
nominal = false(height(T),1);
nominal(chosen) = true;
mainScan = (startsWith(T.label,"hscan") & T.variant=="main") | nominal;
rows = cell(0,8);
for j = find(mainScan & T.id=="Bi1107_10w").'
    for k = 3:4
        ref = chosen(k);
        hRatio = T.h(j)/T.h(chosen(1));
        match = find(mainScan & T.id==T.id(ref) & abs(T.h/T.h(ref)-hRatio)<1e-8,1);
        if isempty(match), continue; end
        A=N(N.sha256==T.sha256(j) & N.branch==1,:);
        B=load(fullfile(T.folder(match),'p2b_results.mat'),'res');
        C=B.res(1).curves;
        e=local_at(C.Q,C.E0(1,:),q);
        ae=local_at(C.Q,C.apex(1,:),q);
        se=local_at(C.Q,C.sE0(1,:),q);
        for u = 1:numel(q)
            rows(end+1,:) = {T.id(ref),hRatio,T.h(j),T.h(match),q(u), ...
                A.E0(u)-e(u),A.apex(u)-ae(u),hypot(A.sE0(u),se(u))}; %#ok<AGROW>
        end
    end
end
R=cell2table(rows,'VariableNames',{'reference','h_ratio','h_A','h_reference','Q', ...
    'A_minus_reference_E0','A_minus_reference_apex','conditional_E0_SE'});
end

function local_convergence(T,S,out)
rows = cell(0,8);
for k = 1:height(T)
    c = S{k}.cfg;
    o = local_operator(T.folder(k));
    o.n_sy = 2*o.n_sy-1;
    M = b1_forward_bi_p2b_material(S{k});
    P = b1_p2b_fit([],c,o,struct('material',M));
    F = S{k}.res(1).fit;
    tail = zeros(size(P.E));
    m = P.E>=c.holdout(1);
    tail(m) = (P.E(m)/c.tail.E_ref).^(-c.tail.p);
    nuisance = tail.*reshape(F.b,1,numel(c.exit_q),[]);
    dz = (P.mu+nuisance-F.mu)./sqrt(P.S2);
    windows = [100 150;150 250;250 1800];
    for j = 1:size(windows,1)
        g = P.E>=windows(j,1) & P.E<windows(j,2);
        if j==3, g = g | P.E==windows(j,2); end
        v = dz(g,:,:);
        rows(end+1,:) = {T.id(k),T.h(k),T.n_sy(k),o.n_sy,windows(j,1),windows(j,2), ...
            max(abs(v),[],'all','omitnan'),sqrt(mean(v.^2,'all','omitnan'))}; %#ok<AGROW>
    end
end
Q = cell2table(rows,'VariableNames',{'id','h','n_sy','n_sy_dense','Emin','Emax','max_sigma','rms_sigma'});
writetable(Q,fullfile(out,'nsy_convergence.csv'));
end

function o = local_operator(path)
M = jsondecode(fileread(fullfile(path,'config.json')));
o = M.operator_config;
o.tail_edges = o.tail_edges(:).';
end

function local_twelve(path,S,out)
c = S.cfg;
c.blocking = 'b12';
o = local_operator(path);
M = b1_forward_bi_p2b_material(S);
P = b1_p2b_fit([],c,o,struct('material',M));
Y = P.Y(:,:,1:2:end)+P.Y(:,:,2:2:end);
assert(isequaln(Y,S.mainY),'Pair sums of archived 12 blocks differ from six blocks.');
F = S.res(1).fit;
tail = zeros(size(P.E));
m = P.E>=c.holdout(1);
tail(m)=(P.E(m)/c.tail.E_ref).^(-c.tail.p);
background=tail.*reshape(F.b,1,numel(c.exit_q),[]);
mu = P.mu(:,:,1:2:end)+P.mu(:,:,2:2:end)+background;
sigma = sqrt(S.mainS2);
z = (Y-mu)./sigma;
z6 = (Y-F.mu)./sigma;
dz = (mu-F.mu)./sigma;
hold = P.E>=c.holdout(1)&P.E<c.holdout(2);
window = P.E>=c.variants(1).window(1)&P.E<=c.variants(1).window(2);
B = [true(1,6);1:6==c.lobo_block];
rows = cell(4,5);
k = 0;
for m = {window,hold}
    for j = 1:size(B,1)
        k=k+1;
        rows(k,:) = {P.E(find(m{1},1)),P.E(find(m{1},1,'last')),j==2, ...
            mean(z(m{1},:,B(j,:)).^2,'all','omitnan')-mean(z6(m{1},:,B(j,:)).^2,'all','omitnan'), ...
            max(abs(dz(m{1},:,B(j,:))),[],'all','omitnan')};
    end
end
D = cell2table(rows,'VariableNames',{'Emin','Emax','held_only','delta_chi2_per_point','max_delta_sigma'});
writetable(D,fullfile(out,'Bi1107_10w_12block_drift.csv'));
end

function y = local_at(q,v,x)
[u,i] = unique(log(q),'stable');
y = interp1(u,v(i),log(min(max(x,q(1)),q(end))));
end
