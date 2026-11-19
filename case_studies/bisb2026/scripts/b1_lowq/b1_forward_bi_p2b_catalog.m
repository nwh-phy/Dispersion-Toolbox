function T = b1_forward_bi_p2b_catalog(out)
%B1_FORWARD_BI_P2B_CATALOG  Index completed fits and their independent block scores.
files = dir(fullfile(out,'fit_Bi*','config.json'));
rows = cell(numel(files),22);
for k = 1:numel(files)
    p = files(k).folder;
    M = jsondecode(fileread(fullfile(p,files(k).name)));
    S = load(fullfile(p,'p2b_results.mat'),'cfg','res');
    c = S.cfg;
    R = S.res(1);
    F = R.fit;
    label = "main";
    version = regexprep(c.version,'\+retry\d+$','');
    first = find(version=='+',1);
    if ~isempty(first), label = string(version(first+1:end)); end
    q = c.nodes.E0;
    a = reshape(F.a,numel(c.nodes.A),[]).';
    highFraction = NaN;
    nHighBounds = 0;
    if size(a,1)==2
        highFraction = sum(a(2,:))/sum(a,'all');
        nHighBounds = nnz(F.E0(2,:)>=c.bounds_E0(2,2)-1);
    end
    rows(k,:) = {string(c.paths.session),label,string(R.name),string(F.start), ...
        string(p),string(M.config_sha256),M.operator_config.h,M.operator_config.n_sy, ...
        numel(c.exit_q),numel(q),size(a,1),c.lobo_block,F.exitflag,F.chi2,F.dof,F.chi2dof, ...
        F.prediction.chi2dof,F.prediction.low_chi2_per_point,M.runtime_s, ...
        highFraction,nHighBounds,F.ho_chi2n};
end
T = cell2table(rows,'VariableNames',{'id','label','variant','start','folder','sha256','h','n_sy', ...
    'channels','energy_nodes','branches','held_block','exitflag','train_chi2','train_dof', ...
    'train_chi2dof','held_chi2dof','held_150_250_chi2_per_point','runtime_s', ...
    'high_residue_fraction','high_E0_at_upper','all_150_250_chi2_per_point'});
T = sortrows(T,{'id','label'});
end
