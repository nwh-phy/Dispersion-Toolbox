function b1_forward_bi_p2a_check(cfg)
%B1_FORWARD_BI_P2A_CHECK  Verify the saved response against acquisition diagnostics.
t0 = tic;
p = cfg.paths;
input = cfg.input.output_dir;
hash = b1_sha256(cfg.response);
p.input_dir = input;
p.response_dir = fullfile(input,['response_' hash(1:8)]);
out = fullfile(input,'b1_forward_bi_checks');
if isfield(cfg,'diagnostic_output_dir'), out = cfg.diagnostic_output_dir; end
assert(isfile(fullfile(p.response_dir,'response.mat')),'Missing P2a response.');
assert(~isfolder(out),'Diagnostic output already exists: %s',out);
mkdir(out);

S = load(p.stage1_file,'C','cfg');
C = S.C;
d = bisb_datasets(p.session);
assert(strcmp(b1_sha256({'file',d.path}),C.raw_sha256),'Raw acquisition differs from stage1.');
D = load(fullfile(input,'p2a_cube.mat'),'cube'); cube = D.cube;
L = load(fullfile(p.response_dir,'response.mat'),'R','J6');
I = load(fullfile(input,'inputs_p2a.mat'),'ct','fwt','pf_mu','cq','tie','call');
N = load(fullfile(p.diag_dir,'noise_model.mat'),'E','g','q_noise');
F = load(fullfile(p.diag_dir,'common_ref.mat'),'ref');
[A, raw] = qe_align_sequence_frames(d.path,d.zero,d.dq,q_select=Inf, ...
    reference_q=S.cfg.reference_abs_q, alignment_window=S.cfg.zlp_window, ...
    energy_zero="detector_peak",raw_profile_window=S.cfg.raw_profile_window);
assert(isequal(A.E,cube.E) && isequal(A.alignment.measured_offset_pixels(:),C.offsets(:)));
assert(isequal(raw.valid_channels,C.valid_channels) && isequal(A.selected,raw.valid_channels));
fprintf('Full detector loaded after %.1f s, valid columns %d\n',toc(t0),numel(A.selected));
[ok, j] = ismember(cube.q,A.Q);
assert(all(ok));
assert(isequal(sum(A.F(:,cube.b6.blk==1,j),2),permute(cube.b6.J(:,:,1),[1 3 2])));
J = zeros(numel(A.E),numel(A.Q),6);
for b = 1:6
    J(:,:,b) = squeeze(sum(A.F(:,cube.b6.blk==b,:),2));
end
local_alignment(A.E,A.F(:,:,j),cube,I,C,S.cfg,cfg.diagnostic,out);
local_noise(cube,A.F(:,:,j),raw.counts_by_energy_frame,C,F.ref,out);
local_band_sensitivity(cube,A.F(:,:,j),I,C,out);
local_response_windows(cube,L.R,out);
Q = A.Q;
clear A
b1_forward_bi_p2a_stripe(cube.E,Q,J,cube.b6.d,C,cfg.diagnostic,out);
b1_forward_bi_p2a_low_energy(cube,L.R,d,cfg.diagnostic,out);
fprintf('Acquisition and low-energy checks completed after %.1f s; starting reference\n',toc(t0));
b1_forward_bi_p2a_reference(p,cube,L.R,L.J6,N,out);
b1_forward_bi_p2a_stats(p,out,cube,L.R);

m = abs(cube.q)<1e-9;
c = I.ct(:,m);
g = N.g(N.E>=250 & N.E<=1800,:);
a0 = sum(J(:,abs(Q)<1e-9,:),3);
[centre,fwhm] = qe_zlp_center(a0,cube.E);
centroid = I.pf_mu(:,1)/C.dq;
serial = corrcoef(centroid(1:end-1),centroid(2:end));
summary = array2table([C.T,C.dq,C.h,C.zero,numel(Q),numel(C.q_axis)-numel(Q), ...
    C.dE,C.raw_peak_pixel,C.reference_pixel, std(F.ref)/mean(F.ref), ...
    min(c+C.offsets(:)*C.dE), max(c+C.offsets(:)*C.dE), ...
    min(c),max(c),mean(c),std(c),centre,fwhm, ...
    std(centroid),max(centroid)-min(centroid),std(diff(centroid)),serial(1,2), ...
    (max(cube.b6.pos_narrow_mode)-min(cube.b6.pos_narrow_mode))/C.dq, ...
    (max(cube.b12.pos_narrow_mode)-min(cube.b12.pos_narrow_mode))/C.dq, ...
    (max(L.J6.pos)-min(L.J6.pos))/C.dq,min(L.J6.logrms),max(L.J6.logrms), ...
    min(I.tie(:,5)),max(I.tie(:,5)),median(g(:),'omitnan')], ...
    'VariableNames',{'frames','dq','h','zero_col','valid_columns','invalid_columns', ...
    'dE','raw_peak_pixel','reference_pixel','beam_cv','raw_c0_min','raw_c0_max', ...
    'c0_min','c0_max','c0_mean','c0_sd','sum_c0','sum_fwhm0','position_sd_ch', ...
    'position_range_ch','position_step_sd_ch','position_lag1','core_range6_ch', ...
    'core_range12_ch','area_range6_ch','area_logrms_min','area_logrms_max', ...
    'coupling_slope_min','coupling_slope_max','fast_g_median'});
writetable(summary,fullfile(out,'b1_forward_bi_summary.csv'));
scripts = {'config','run','check','stripe','low_energy','reference','stats'};
code = struct();
for k = 1:numel(scripts)
    code.(scripts{k}) = b1_sha256({'file',which(['b1_forward_bi_p2a_' scripts{k}])});
end
record = struct('config',cfg,'config_sha256',b1_sha256(cfg),'code_sha256',code,'stage1_sha256',C.config_sha256, ...
    'input_sha256',b1_sha256(cfg.input),'response_sha256',hash, ...
    'raw_sha256',C.raw_sha256,'elapsed_s',toc(t0));
fid = fopen(fullfile(out,'b1_forward_bi_config.json'),'w');
assert(fid>=0); fwrite(fid,jsonencode(record,'PrettyPrint',true)); fclose(fid);
end

function local_response_windows(cube,R,out)
opts = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-12, ...
    'StepTolerance',1e-10,'MaxIterations',800,'MaxFunctionEvaluations',30000);
test = [-60 20;-40 32;-28 28];
rows = {};
for b = 1:6
    for j = find(abs(cube.q)<=cube.dq+1e-12)
        y = cube.b6.J(:,j,b);
        [c,w] = qe_zlp_center(y,cube.E);
        rows(end+1,:) = {b,cube.q(j),test(1,1),test(1,2),c,w,R.b6.centre(b,j), ...
            R.b6.chi2dof(b,j),R.b6.at_bound(b,j)}; %#ok<AGROW>
        for k = 2:size(test,1)
            [t,fit] = qe_elastic_line_model('fit',cube.E,y,test(k,:),opts);
            rows(end+1,:) = {b,cube.q(j),test(k,1),test(k,2),c,w,t.c,fit.chi2dof,fit.at_bound}; %#ok<AGROW>
        end
    end
end
writetable(cell2table(rows,'VariableNames',{'block','q','window_lo','window_hi', ...
    'native_mode_meV','native_fwhm_meV','fit_mode_meV','chi2dof','at_bound'}), ...
    fullfile(out,'b1_forward_bi_elastic_response_windows.csv'));
end

function local_alignment(E,F,cube,I,C,prep,cfg,out)
q = cube.q;
i0 = abs(q)<1e-9;
Y = F(:,:,i0);
offset = C.offsets(:);
step = C.dE;
u = E>=-140 & E<=140;
X = E(u);
beam = I.pf_mu(:,1)/C.dq;
aligned = I.ct(:,i0);
raw0 = nan(nnz(u),C.T);
sets = {i0,abs(q)<=prep.reference_abs_q+1e-12,true(size(q))};
names = {'q0','P0 core','full P2a q'};
rr = {};
for k = 1:numel(sets)
    ref = sum(F(:,:,sets{k}),3);
    for frame = 1:C.T
        y = interp1(E+offset(frame)*step,ref(:,frame),X,'linear');
        if k==1, raw0(:,frame) = y; end
        refRaw(:,frame) = y; %#ok<AGROW>
    end
    for w = cfg.elastic_windows
        try
            a = qe_zlp_integer_align(X,refRaw,refRaw,[-w w]);
        catch err
            if ~strcmp(err.identifier,'qe_zlp_integer_align:NoZLP'), rethrow(err); end
            rr(end+1,:) = {names{k},w,0,NaN,NaN,NaN,NaN,NaN,NaN}; %#ok<AGROW>
            continue
        end
        peak = X(a.peak_pixel(a.valid));
        alt = a.measured_offset_pixels(:);
        m = a.valid(:);
        c = aligned(m)+step*(offset(m)-alt(m));
        piv = X(a.reference_pixel);
        rr(end+1,:) = {names{k},w,nnz(m),piv,sum(offset(m)~=alt(m)), ...
            mean(c),std(c),mean(c)-piv,max(peak)-min(peak)}; %#ok<AGROW>
    end
end
writetable(cell2table(rr,'VariableNames',{'reference','half_window_meV','valid_frames', ...
    'reference_pixel_energy','different_offsets','c0_mean','c0_sd','c0_mean_vs_own_ref','raw_peak_range'}), ...
    fullfile(out,'b1_forward_bi_alignment_windows.csv'));

% Look for secondary, separated maxima in the frame core; half-width and
% interpolation are evaluated on the same native 4 meV grid for all controls.
Ys = movmean(Y(u,:),3,1);
peak = islocalmax(Ys,1,'MinSeparation',2);
peak = peak & Ys>=0.2*max(Ys,[],1);
nr = sum(peak,1);
W = zeros(size(Y));
for frame = 1:C.T
    W(:,frame) = interp1(E-aligned(frame),Y(:,frame),E,'linear',0);
end
rows = {};
shape = readtable(fullfile(fileparts(out),'block_shape_consistency.csv'));
for nb = [6 12]
    b = cube.(sprintf('b%d',nb));
    pos = b.pos_narrow_mode./C.dq;
    for j = 1:nb
        f = b.blk==j;
        y = sum(Y(:,f),2);
        [c,w] = qe_zlp_center(y,E);
        [~,wc] = qe_zlp_center(sum(W(:,f),2),E);
        yr = sum(raw0(:,f),2);
        [~,rawPeak] = max(yr);
        [~,wr] = qe_zlp_center(yr,X-X(rawPeak));
        meanFrame = median(I.fwt(f,i0));
        region = strcmp(shape.region_meV,'400-1800') & shape.block==ceil(j*6/nb) & abs(shape.q)>=C.dq/2;
        % Shape checks are on the original six blocks, even for a child block.
        sh = median(shape.chi2dof_shape_only_after_shift(region),'omitnan');
        rows(end+1,:) = {nb,j,nnz(f),pos(j),std(beam(f)),c,w,wc,wr,meanFrame,std(aligned(f)), ...
            b.c_allq(j),b.fwhm_allq(j),nnz(nr(f)>1),sh}; %#ok<AGROW>
    end
end
writetable(cell2table(rows,'VariableNames',{'n_blocks','block','frames','narrow_mode_ch','position_sd_ch', ...
    'c0_meV','fwhm_q0','fwhm_q0_subpixel','fwhm_q0_raw','median_frame_fwhm','c0_frame_sd', ...
    'c_allq','fwhm_allq','multi_max_frames','B1_shape_chi2_median'}), ...
    fullfile(out,'b1_forward_bi_elastic_blocks.csv'));
frame = (1:C.T).';
writetable(table(frame,beam,aligned,nr(:),I.fwt(:,i0),I.ct(:,i0)+step*offset, ...
    'VariableNames',{'frame','position_ch','aligned_c0_meV','core_maxima','fwhm_meV','raw_c0_meV'}), ...
    fullfile(out,'b1_forward_bi_elastic_frames.csv'));
end

function local_noise(cube,F,detector,C,ref,out)
E = cube.E(:);
ref = ref(:);
ref2 = sum(detector,1).';
ho = C.support(E>=100 & E<=300);
ho0 = C.support(E>=150 & E<=250);
check = ref2;
for t = 1:C.T
    ref2(t) = ref2(t)-sum(detector(ho+C.offsets(t),t));
    check(t) = check(t)-sum(detector(ho0+C.offsets(t),t));
end
assert(max(abs(check-ref))<=1e-12*max(ref) && all(ref2>0),'Common reference differs from P0.');
r = ref2./mean(ref2)-ref./mean(ref);
fprintf('P2c exclusion 100-300: P0 minus candidate normalized max %.5g, RMS %.5g\n',max(abs(r)),sqrt(mean(r.^2)));
rows = {};
bands = [100 150;150 250;250 400;400 1800];
for j = 1:numel(cube.q)
    M = F(:,:,j);
    avg = max(mean(M,2),1);
    raw = var(M,0,2)./avg;
    [~,common] = b1_p1b_noise(M,E,ref,[150 250],'common',21);
    [~,fast] = b1_p1b_noise(M,E,ref,[150 250],'first_diff',21);
    [~,broad] = b1_p1b_noise(M,E,ref2,[100 300],'first_diff',21);
    for k = 1:size(bands,1)
        e = E>=bands(k,1)&E<bands(k,2);
        rows(end+1,:) = {cube.q(j),bands(k,1),bands(k,2),median(raw(e)), ...
            median(common(e)),median(fast(e)),median(broad(e)), ...
            max(abs(fast(e)-broad(e)))}; %#ok<AGROW>
    end
end
writetable(cell2table(rows,'VariableNames',{'q','Emin','Emax','g_raw','g_common', ...
    'g_first_diff_P0','g_first_diff_100_300_holdout','max_g_difference'}), ...
    fullfile(out,'b1_forward_bi_noise.csv'));
end

function local_band_sensitivity(cube,F,I,C,out)
E = cube.E(:);
dE = C.dE;
offsets = C.offsets(:);
J = cube.b6.J;
J2 = cube.b12.J(:,:,1:2:end)+cube.b12.J(:,:,2:2:end);
assert(max(abs(J(:)-J2(:)))<=1e-12*max(J(:)),'Six-block and twelve-block sums differ.');
a0 = I.ct(:,cube.q==0);
bands = [120 150;150 250;250 400;400 1800];
rows = {};
for b = 1:6
    frames = find(cube.b6.blk==b);
    a0b = mean(a0(frames));
    ab = mean(I.call(frames));
    for s = 1:size(bands,1)
        W = @(x) max(0,min(x+dE/2,bands(s,2))-max(x-dE/2,bands(s,1)))/dE;
        wt = W(E);
        ea = E+offsets(frames).'*dE;
        ew = E-(a0(frames)-a0b).';
        ei = E-(I.call(frames)-ab).';
        base = reshape(sum(F(:,frames,:).*wt,[1 2]),[],1);
        raw = reshape(sum(F(:,frames,:).*W(ea),[1 2]),[],1);
        q0 = reshape(sum(F(:,frames,:).*W(ew),[1 2]),[],1);
        int = reshape(sum(F(:,frames,:).*W(ei),[1 2]),[],1);
        for j = 1:numel(cube.q)
            rows(end+1,:) = {b,cube.q(j),bands(s,1),bands(s,2),base(j), ...
                raw(j)/base(j)-1,q0(j)/base(j)-1,int(j)/base(j)-1}; %#ok<AGROW>
        end
    end
end
writetable(cell2table(rows,'VariableNames',{'block','q','Emin','Emax','aligned_counts', ...
    'raw_over_aligned_minus1','q0_recentered_over_aligned_minus1', ...
    'integrated_recentered_over_aligned_minus1'}),fullfile(out,'b1_forward_bi_elastic_band_sensitivity.csv'));
end
