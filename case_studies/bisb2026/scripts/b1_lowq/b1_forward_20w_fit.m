function b1_forward_20w_fit(cfg, out)
%B1_FORWARD_20W_FIT  Fit blocks 1–5 and predict block 6 at each h.
S = load(fullfile(out,'b1_forward_20w_inputs.mat'),'p');
p = S.p;
o = b1_p2a_operator_config(p);
inputTrain = fullfile(p.input_dir,'b1_forward_20w_training');
responseTrain = fullfile(inputTrain,'response');
assert(~isfolder(inputTrain),'Training data already exist.');
mkdir(responseTrain);
L = load(fullfile(p.input_dir,'p2a_cube.mat'),'cube');
cube = L.cube;
exposure = cube.b6.d;
cube.b6.J = cube.b6.J(:,:,1:5);
save(fullfile(inputTrain,'p2a_cube.mat'),'cube','-v7.3');
L = load(fullfile(p.response_dir,'response.mat'),'R','J6');
R = L.R;
R.b6.area = R.b6.area(1:5,:);
J6 = L.J6;
save(fullfile(responseTrain,'response.mat'),'R','J6');
train = p;
train.input_dir = inputTrain;
train.response_dir = responseTrain;
c = b1_p2b_config(train,cfg.nodes);
c.version = cfg.version;
c.Q_curve = unique([c.Q_curve c.nodes.E0]);
c.lobo_block = 5; % Existing fitter's internal check within the training set.
assert(numel(c.exit_q)==40 && numel(exposure)==6);
c.highq.q_range(2) = 0.01175; % maximum paired 20w bin centre in the archive
H = readtable(c.highq.file);
m = abs(H.q)>=c.highq.q_range(1)-1e-9 & abs(H.q)<=c.highq.q_range(2)+1e-9;
assert(all(strcmp(H.status(m),'converged')) && all(ismember(-H.q(m),H.q(m))), ...
    'High-Q source has unpaired or unconverged rows.');
h = [o.h cfg.h_scan];
fits = cell(size(h));
issues = [];
for k = 1:numel(h)
    try
        tm = tic;
        if k==1
            fitCfg = c;
            fitCfg.variants = c.variants(ismember({c.variants.name},cfg.variants));
            fitOp = o;
        else
            [fitCfg,fitOp] = b1_p2b_slot_config('hscan',h(k),c,o);
        end
        fits{k} = b1_p2b_fit([],fitCfg,fitOp);
        F = load(fullfile(fits{k},'p2b_results.mat'),'res','cfg');
        local_predict(p,F,fitOp,exposure,fits{k});
        fprintf('h %.4f complete %.2f s: %s\n',h(k),toc(tm),fits{k});
    catch err
        issues(end+1) = k; %#ok<AGROW>
        fprintf(2,'h %.4f failed: %s\n',h(k),getReport(err,'extended'));
    end
end
save(fullfile(out,'b1_forward_20w_fits.mat'),'h','fits','issues');
assert(isempty(issues),'A fit or prediction failed.');
end

function local_predict(p,S,o,exposure,fitDir)
cubeFile = load(fullfile(p.input_dir,'p2a_cube.mat'),'cube');
cube = cubeFile.cube;
q = cube.q(:).';
response = load(fullfile(p.response_dir,'response.mat'),'R','J6');
c = S.cfg;
o.exit_q = c.exit_q;
ie = arrayfun(@(v)find(abs(q-v)<1e-9),c.exit_q);
intr0 = qe_lowq_joint_forward('intrinsic','Qn_E0',c.nodes.E0, ...
    'Qn_G',c.nodes.G,'Qn_A',c.nodes.A);
rows = cell(numel(S.res),6);
for k = 1:numel(S.res)
    V = c.variants(k);
    intr = intr0;
    intr.ext = V.ext;
    st = struct('blocking','b6','variant',V.kernel,'ytype',V.ytype,'wide',V.wide, ...
        'n_a',o.n_a,'n_sx',o.n_sx,'n_sy',o.n_sy,'n_gh',o.n_gh,'Qf',1, ...
        'step',o.omega(2),'support',o.support,'tail_max',o.tail_edges(end), ...
        'exits',1:numel(ie),'blocks',6);
    [ops,K,sw,ex] = b1_p2a_build_ops(st,o,cube,response.R,response.J6,q,cube.dq,ie);
    fopt = struct('omega',o.omega(1):st.step:o.omega(3),'eps',-72:st.step:72, ...
        'variant',V.kernel,'slice_weight',{sw},'exit_slice',ex);
    cache = qe_lowq_joint_forward('build',ops,K,fopt);
    F = S.res(k).fit;
    cols = qe_lowq_joint_forward('eval',cache,b1_p2b_decode(F.u,c.nodes),intr);
    mu = squeeze(sum(cols.*reshape(F.a,1,[],1,1),2));
    E = cache.E_out(:);
    [tf,loc] = ismember(round(E),round(cube.E));
    Y = nan(size(mu));
    Y(tf,:) = cube.b6.J(loc(tf),ie,6);
    S2 = max(movmean(Y,c.noise_smooth,1,'omitnan'),1);
    tails = reshape(F.b,numel(ie),5);
    tail6 = sum(tails,2).' * exposure(6)/sum(exposure(1:5));
    eligible = E>=c.holdout(1);
    mu(eligible,:) = mu(eligible,:) + (E(eligible)/c.tail.E_ref).^(-c.tail.p).*tail6;
    z = (Y-mu)./sqrt(S2);
    win = E>=V.window(1) & E<=V.window(2);
    hold = E>=c.holdout(1) & E<c.holdout(2);
    rows(k,:) = {V.name,nnz(isfinite(z(win,:))),mean(z(win,:).^2,'all','omitnan'), ...
        mean(z(win,:),'all','omitnan'),mean(z(hold,:).^2,'all','omitnan'), ...
        mean(z(hold,:),'all','omitnan')};
    if k==1
        save(fullfile(fitDir,'b1_forward_20w_block6_prediction.mat'),'E','Y','mu','S2','tail6');
    end
end
T = cell2table(rows,'VariableNames',{'variant','n_window','chi2_per_point','mean_residual_sigma', ...
    'low_energy_chi2_per_point','low_energy_mean_residual_sigma'});
writetable(T,fullfile(fitDir,'b1_forward_20w_block6.csv'));
disp(T);
end
