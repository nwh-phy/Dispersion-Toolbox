function b1_forward_bi_p2a_reference(p,cube,R,J6,N,out)
%B1_FORWARD_BI_P2A_REFERENCE  Fixed-L forward wiring and numerical sensitivity.
cfg = b1_p2a_operator_config(p);
fit = b1_p2b_config(p);
cfg.exit_q = fit.exit_q;
ref = b1_p2a_operator([],[],cfg,true);
E = ref.Eo(:);
q = cube.q;
dq = cube.dq;
[tf,ie] = ismembertol(cfg.exit_q,q,1e-10,'DataScale',1);
assert(all(tf) && numel(ref.a_ref)==6 && all(isfinite(ref.a_ref)));
assert(any(ref.a_ref>0) && all(isfinite(ref.mu_ref(:))));
[eOk,ix] = ismember(round(E),round(cube.E));
assert(isequal(ref.Y(eOk,:,:),cube.b6.J(ix(eOk),ie,:)));

% Check the same residues at successively finer midpoint quadrature, including
% the two lowest nonzero native channels and the first and last six-block sums.
qa = unique([dq .001 .002 .005]);
key = find(ismembertol(abs(cfg.exit_q),qa,1e-10,'DataScale',1));
blk = [1 6];
st = struct('blocking','b6','variant','slice','ytype','gauss','wide','on', ...
    'n_a',cfg.n_a,'n_sx',cfg.n_sx,'n_sy',cfg.n_sy,'n_gh',cfg.n_gh,'Qf',1, ...
    'step',cfg.omega(2),'support',cfg.support,'tail_max',cfg.tail_edges(end), ...
    'exits',key,'blocks',blk);
base = ref.mu_ref(:,key,blk);
Y = ref.Y(:,key,blk);
G = interp1(N.E,N.g,E,'linear',1);
[~,gn] = min(abs(cfg.exit_q(key).'-N.q_noise),[],2);
sigma = sqrt(max(Y,1).*G(:,gn));
intr = qe_lowq_joint_forward('intrinsic');
theta = struct('E0sq',cfg.ref_E0.^2,'G',cfg.ref_G);
pairs = {'rebuild',st; '2x n_sy',st; '4x n_sy',st; 'exit kernel',st; ...
    'integrated kernel',st; 'wide off',st; '12 blocks pair sum',st; ...
    'narrow Pearson y',st; 'mixture y',st; 'remove q0 block zero from all kernels',st; ...
    'empirical mode for each slice',st; 'mixture y 2x n_sy',st};
pairs{2,2}.n_sy = 2*st.n_sy-1;
pairs{3,2}.n_sy = 4*st.n_sy-3;
pairs{4,2}.variant = 'exit';
pairs{5,2}.variant = 'integrated';
pairs{6,2}.wide = 'off';
pairs{7,2}.blocking = 'b12';
pairs{7,2}.blocks = [1 2 11 12];
pairs{8,2}.ytype = 'narrow_pearson';
pairs{9,2}.ytype = 'mixture';
% The 2D y alternatives use the same two blocks and four exits.
pairs{8,2}.exits = key(ismembertol(abs(cfg.exit_q(key)),[dq .002],1e-10,'DataScale',1));
pairs{9,2}.exits = pairs{8,2}.exits;
pairs{12,2} = pairs{9,2};
pairs{12,2}.n_sy = pairs{2,2}.n_sy;
rows = {};
allRows = {};
last = base;
for a = 1:size(pairs,1)
    s = pairs{a,2};
    t = tic;
    response = R;
    if strcmp(pairs{a,1},'remove q0 block zero from all kernels')
        for b = blk
            z0 = R.b6.centre(b,q==0);
            for r = 1:numel(q)
                response.b6.th{b,r}.c = R.b6.th{b,r}.c-z0;
            end
        end
    end
    if strcmp(pairs{a,1},'empirical mode for each slice')
        for b = blk
            for r = 1:numel(q)
                c = qe_zlp_center(cube.b6.J(:,r,b),cube.E);
                response.b6.th{b,r}.c = c;
            end
        end
    end
    [ops,K,sw,ex] = b1_p2a_build_ops(s,cfg,cube,response,J6,q,dq,ie);
    radial = max(cellfun(@(o)o.out_of_range,ops));
    o = struct('omega',cfg.omega(1):s.step:cfg.omega(3),'eps',-72:s.step:72, ...
        'variant',s.variant,'exit_slice',ex,'slice_weight',{sw});
    H = qe_lowq_joint_forward('build',ops,K,o);
    Z = qe_lowq_joint_forward('eval',H,theta,intr);
    z = reshape(sum(Z.*reshape(ref.a_ref,1,[],1,1),2),numel(E),numel(s.exits),[]);
    clear H Z
    if strcmp(s.blocking,'b12'), z = z(:,:,1:2:end)+z(:,:,2:2:end); end
    [~,idx] = ismember(s.exits,key);
    v = sigma(:,idx,:);
    u = (z-base(:,idx,:))./v;
    if a<=3
        w = (z-last)./v;
        last = z;
    elseif strcmp(pairs{a,1},'mixture y')
        previousMixture = z;
        w = u;
    elseif strcmp(pairs{a,1},'mixture y 2x n_sy')
        w = (z-previousMixture)./v;
    else
        w = u;
    end
    if a==1, assert(max(abs(u(:)),[],'omitnan')<1e-10,'Forward values differ on exact subset rebuild.'); end
    assert(all(isfinite(z(:)))&&radial==0,'Invalid reference forward output.');
    elapsed = toc(t);
    windows = [100 150;150 250;250 1800];
    for k = 1:size(windows,1)
        m = E>=windows(k,1)&E<=windows(k,2)&eOk;
        dd = u(m,:,:); dc = w(m,:,:);
        rows(end+1,:) = {pairs{a,1},s.n_sy,windows(k,1),windows(k,2), ...
            max(abs(dd(:))),sqrt(mean(dd(:).^2)),max(abs(dc(:))),sqrt(mean(dc(:).^2)),elapsed}; %#ok<AGROW>
        for j = 1:numel(s.exits)
            rr = dd(:,j,:);
            allRows(end+1,:) = {pairs{a,1},cfg.exit_q(s.exits(j)),windows(k,1),windows(k,2), ...
                max(abs(rr(:))),sqrt(mean(rr(:).^2)),sum(rr(:).^2)}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows,'VariableNames',{'test','n_sy','Emin','Emax','max_sigma_vs_base','rms_sigma_vs_base', ...
    'max_sigma_vs_previous','rms_sigma_vs_previous','elapsed_s'});
U = cell2table(allRows,'VariableNames',{'test','q','Emin','Emax','max_sigma','rms_sigma','sum_delta_chi2'});
writetable(T,fullfile(out,'b1_forward_bi_reference_comparisons.csv'));
writetable(U,fullfile(out,'b1_forward_bi_reference_by_q.csv'));
save(fullfile(out,'b1_forward_bi_reference.mat'),'ref','cfg','fit','theta','intr','key','blk');
end
