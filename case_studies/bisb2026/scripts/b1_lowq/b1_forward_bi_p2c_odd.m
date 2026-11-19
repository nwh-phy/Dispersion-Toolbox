function b1_forward_bi_p2c_odd(P,out)
%B1_FORWARD_BI_P2C_ODD  Native paired spectra and held-block source diagnostics.
t0 = tic;
S = load(fullfile(P.archive,sprintf('fit_%s_%s',P.id,P.main),'p2b_results.mat'));
p = S.cfg.paths;
D = load(fullfile(p.input_dir,'p2a_cube.mat'),'cube');
R = load(fullfile(p.response_dir,'response.mat'),'R','J6');
cube = D.cube;
E = S.mainEo(:);
q = S.cfg.exit_q(:).';
jp = find(q>0);
[tf,jm] = ismembertol(-q(jp),q,1e-10,'DataScale',1);
assert(all(tf),'Signed channels need one-to-one partners.');
qp = q(jp);
nb = size(S.mainY,3);
[~,ip] = ismembertol(qp,cube.q,1e-10,'DataScale',1);
[~,im] = ismembertol(-qp,cube.q,1e-10,'DataScale',1);
Yp = S.mainY(:,jp,:); Ym = S.mainY(:,jm,:);
Mp = S.mainMu(:,jp,:); Mm = S.mainMu(:,jm,:);
V = S.mainS2;
sd = sqrt(V(:,jp,:)+V(:,jm,:));
logR = log(max(Yp,1)./max(Ym,1));
logM = log(max(Mp,1)./max(Mm,1));
zr = (Yp-Ym-Mp+Mm)./sd;
ratio = Yp./Ym; ratio(Ym<=0) = NaN;
[Egrid,Qgrid,Bgrid] = ndgrid(E,qp,1:nb);
w = Egrid>=250 & Egrid<=1800 & isfinite(zr);
T = table(Egrid(w),Qgrid(w),Bgrid(w),Yp(w),Ym(w),ratio(w), ...
    Yp(w)-Ym(w)-Mp(w)+Mm(w),zr(w), (Yp(w)./max(Mp(w),1))./(Ym(w)./max(Mm(w),1)), ...
    'VariableNames',{'E_meV','abs_q','block','Iplus','Iminus','raw_ratio','residual_difference','odd_sigma','ratio_over_model'});
writetable(T,fullfile(out,'odd_native.csv'));

[Qall,J] = local_full(P,cube,p);
[tf,jFar] = ismembertol(-Qall,Qall,1e-10,'DataScale',1);
far = find(Qall>=.055 & Qall<=.06 & tf);
far2 = find(Qall>=.065 & Qall<=.085 & tf);
Band = [P.bands; -500 -250; -240 -100; -84 -68; -52 -24; 24 52; 68 84];
rows = {};
for b = 1:nb
    for k = 1:size(Band,1)
        m = cube.E>=Band(k,1) & cube.E<Band(k,2)+(k==3);
        pp = sum(cube.b6.J(m,ip,b),1); pm = sum(cube.b6.J(m,im,b),1);
        fd = [sum(J(m,far,b),'all')/sum(J(m,jFar(far),b),'all'), ...
              sum(J(m,far2,b),'all')/sum(J(m,jFar(far2),b),'all')];
        x = E>=Band(k,1) & E<Band(k,2)+(k==3);
        t = sum(Yp(x,:,b)-Ym(x,:,b)-Mp(x,:,b)+Mm(x,:,b),1);
        if k<=3
            v = sum(sd(x,:,b).^2,1);
            rp = sum(Yp(x,:,b)-Mp(x,:,b),1); rm = sum(Ym(x,:,b)-Mm(x,:,b),1);
        else
            t(:)=NaN; v=nan(size(t)); rp=v; rm=v;
        end
        for i = 1:numel(qp)
            rows(end+1,:) = {b,k,Band(k,1),Band(k,2),qp(i),pp(i),pm(i),pp(i)/pm(i), ...
                (pp(i)/pm(i))/(R.R.b6.area(b,ip(i))/R.R.b6.area(b,im(i))), ...
                t(i),t(i)/sqrt(v(i)),rp(i),rm(i),rp(i)/rm(i),fd(1),fd(2),R.J6.pos(b), ...
                R.R.b6.centre(b,ip(i))-R.R.b6.centre(b,im(i))}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows,'VariableNames',{'block','band','lo_meV','hi_meV','abs_q', ...
    'Iplus','Iminus','raw_ratio','ratio_over_elastic_area_ratio','residual_difference','odd_sigma', ...
    'resid_plus','resid_minus','resid_ratio','far_R_1','far_R_2','beam_pos','elastic_center_diff'});
writetable(T,fullfile(out,'odd_bands.csv'));

% Small-displacement tests on the symmetric model. delta denotes the full
% +/- pair displacement, so log R changes by -delta d(log I)/d|q|.
mid = log(max((Mp+Mm)/2,1));
dg = zeros(size(mid));
for i = 1:numel(qp)
    l = max(1,i-1); r = min(numel(qp),i+1);
    dg(:,i,:) = (mid(:,r,:)-mid(:,l,:))/(qp(r)-qp(l));
end
dg = movmean(dg,9,1);
ded = movmean(mid,9,1);
ded = (ded([2:end end],:,:)-ded([1 1:end-1],:,:))./(E([2:end end])-E([1 1:end-1]));
cDiff = reshape((R.R.b6.centre(:,ip)-R.R.b6.centre(:,im)).',1,numel(qp),nb);
intr = qe_lowq_joint_forward('intrinsic');
xc = -dg*1e-3; xt = -dg.*(E*intr.qE_per_meV);
ze = -ded.*cDiff;
u = qp/max(qp);
g = cat(4,repmat(reshape(u,1,[]),size(E,1),1,nb), ...
    repmat(reshape(u.*(1-u),1,[]),size(E,1),1,nb), ...
    repmat(reshape(u.*(1-u).^2,1,[]),size(E,1),1,nb));
gainE = g.*reshape(log(max(E,250)/700),[],1,1,1);
xx = {zeros(numel(w),0),xc(:),[xc(:) xt(:)], ...
    reshape(g,[],3),[reshape(g,[],3) reshape(gainE,[],3)],ze(:),[xc(:) xt(:) ze(:)]};
names = {'archive','constant_q_shift','shift_plus_tilt','energy_independent_gain','energy_dependent_gain','energy_coupling','shift_tilt_coupling'};
pr = {};
rt = logR-logM;
iv = 1./(V(:,jp,:)./max(Yp,1).^2+V(:,jm,:)./max(Ym,1).^2);
for b = [6 1]
    tr = w & Bgrid~=b; ho = w & Bgrid==b;
    for k = 1:numel(xx)
        X = xx{k};
        coef = (X(tr(:),:).*sqrt(iv(tr)))\(rt(tr).*sqrt(iv(tr)));
        diff = rt(:)-X*coef;
        pr(end+1,:) = {names{k},b,mean(diff(tr).^2.*iv(tr)), ...
            mean(diff(ho).^2.*iv(ho)),mat2str(coef',6)}; %#ok<AGROW>
    end
end
A = cell2table(pr,'VariableNames',{'model','held_block','train_odd_chi2_per_pair','held_odd_chi2_per_pair','coef'});
writetable(A,fullfile(out,'source_predictors.csv'));

f = figure('Visible','off','Theme','light','Position',[0 0 1220 740]);
tiledlayout(f,2,3,'TileSpacing','compact');
for b = 1:nb
    ax = nexttile;
    use = E>=250 & E<=1800;
    z = log2(movsum(Yp(:, :, b),9,1)./movsum(Ym(:, :, b),9,1));
    imagesc(ax,qp*1e3,E(use),z(use,:),[-1 1]); axis(ax,'xy'); colorbar(ax);
    title(ax,sprintf('Block %d',b)); xlabel(ax,'|q| (10^{-3} Å^{-1})'); ylabel(ax,'E (meV)');
end
sgtitle(sprintf('%s   log_2 [I(+q)/I(-q)]',P.id),'Interpreter','none');
exportgraphics(f,fullfile(out,'fig_raw_ratio.png'),'Resolution',120); close(f);
fprintf('%s source diagnostics %.1f s, hbar v %.2f eV A\n',P.id,toc(t0),1e-3/intr.qE_per_meV);
disp(A);
end

function [Q,J] = local_full(P,cube,p)
S = load(p.stage1_file,'C','cfg');
d = bisb_datasets(P.id);
A = qe_align_sequence_frames(d.path,d.zero,d.dq,q_select=Inf, ...
    reference_q=S.cfg.reference_abs_q,alignment_window=S.cfg.zlp_window, ...
    energy_zero="detector_peak",raw_profile_window=S.cfg.raw_profile_window);
assert(isequal(A.E,cube.E) && isequal(A.alignment.measured_offset_pixels(:),S.C.offsets(:)));
Q = A.Q;
J = zeros(numel(A.E),numel(Q),6);
for b = 1:6
    J(:,:,b) = squeeze(sum(A.F(:,cube.b6.blk==b,:),2));
end
[~,i] = ismembertol(cube.q,Q,1e-10,'DataScale',1);
assert(isequal(J(:,i,:),cube.b6.J),'Re-aligned native counts differ from archive.');
end
