function b1_forward_bi_p2a_stats(p,out,cube,R)
%B1_FORWARD_BI_P2A_STATS  Compact, reproducible tables for the P2a report.
T = readtable(fullfile(out,'b1_forward_bi_stripe_profiles.csv'));
S = readtable(fullfile(out,'b1_forward_bi_stripe_scan.csv'));
E = readtable(fullfile(out,'b1_forward_bi_tail.csv'));
P = readtable(fullfile(out,'b1_forward_bi_phonon_blocks.csv'));
R0 = readtable(fullfile(out,'b1_forward_bi_columns.csv'));
A = readtable(fullfile(out,'b1_forward_bi_elastic_blocks.csv'));
B = readtable(fullfile(out,'b1_forward_bi_elastic_band_sensitivity.csv'));
F = readtable(fullfile(p.input_dir,'block_shape_consistency.csv'));
J = readtable(fullfile(p.response_dir,'response_joint_fit.csv'));
for e = unique(E.E_meV).'
    m = E.q==0 & E.E_meV==e;
    M = E.E_meV==e & abs(E.q)<=.002;
    v = E.envelope_over_loss(m);
    fprintf('tail E=%g q0 proxy/loss=[%.3f %.3f] width=%.3f loss counts %.0f lowq width median/max %.3f %.3f\n', ...
        e,median(E.proxy_min(m)./E.measured_loss(m)),median(E.proxy_max(m)./E.measured_loss(m)), ...
        median(v),median(E.measured_loss(m)),median(E.envelope_over_loss(M)),max(E.envelope_over_loss(M)));
end
fprintf('phonon c130 q0 by block %s\n',mat2str(P.legacy_c130(P.q==0).',4));
fprintf('phonon area q0 by block %s\n',mat2str(P.gauss_area_counts(P.q==0).',5));
meanF = mean(P.gauss_over_elastic_area(P.q==0));
fprintf('phonon q0 mean area/elastic=%.5g\n',meanF);
for nb = [6 12]
    m = A.n_blocks==nb;
    fprintf('elastic n=%d q0 cent %s FWHM %s raw %s subpix %s frame %s integrated %s %s raw multi %d\n', ...
        nb,local_span(A.c0_meV(m)),local_span(A.fwhm_q0(m)),local_span(A.fwhm_q0_raw(m)), ...
        local_span(A.fwhm_q0_subpixel(m)),local_span(A.median_frame_fwhm(m)), ...
        local_span(A.c_allq(m)),local_span(A.fwhm_allq(m)),sum(A.multi_max_frames(m)));
end
for e = unique(B.Emin).'
    n = B.Emin==e & abs(B.q)<=.005 & abs(B.q)>0;
    V = abs(B{n,6:8});
    fprintf('alignment band %g raw, q0 subpix, integrated subpix max %s\n',e,mat2str(max(V),4));
end
z = abs(cube.q)<1e-9;
SL = readtable(fullfile(p.response_dir,'response_slices.csv'));
SW = readtable(fullfile(out,'b1_forward_bi_elastic_response_windows.csv'));
for wh = unique(SW.window_hi).'
    i = SW.q==0 & SW.window_hi==wh;
    fprintf('response window upper %g q0 centre %s fit minus native %s\n',wh, ...
        local_span(SW.fit_mode_meV(i)),local_span(SW.fit_mode_meV(i)-SW.native_mode_meV(i)));
end
extra = median(SL.extra_mass_60_to_100meV(SL.q==0 & SL.n_blocks==6));
fprintf('slice q0 6-block centre %s FWHM %s chi %s bound %d; all bound %d/%d +60..100 extra median %.5g\n', ...
    local_span(R.b6.centre(:,z)),local_span(R.b6.fwhm(:,z)),local_span(R.b6.chi2dof(:,z)), ...
    nnz(R.b6.at_bound(:,z)),nnz(R.b6.at_bound),numel(R.b6.at_bound), ...
    extra);
% Response capture and line-fit statistics by block and core-window fit.
for k = [1 7]
    fprintf('angular fit row %d: wide %.4f narrow FWHM %.6g 2m %.3f wide 2m %.3f wide dQ %.6g narrow scales %.6g %.6g\n', ...
        k,J.broad_fraction(k),J.narrow_fwhm_Ainv(k),J.narrow_2m(k),J.broad_2m(k), ...
        J.broad_offset_Ainv(k),J.narrow_scale_minus(k),J.narrow_scale_plus(k));
end
wings = sum(R.b6.area(:,abs(cube.q)>=.0015-1e-9),2)./sum(R.b6.area,2);
fprintf('empirical slice area fraction |q|>=.0015 in measured support %s\n',local_span(wings));
for group = {'400-1800','250-400','150-250','40-150'}
    m = strcmp(F.region_meV,group{1}) & abs(F.q)<=.005 & abs(F.q)>0;
    fprintf('shape %s q nonzero median %.3f -> %.3f, worst %.3f, shifted RMS %.4f\n',group{1}, ...
        median(F.chi2dof_shape_only(m)),median(F.chi2dof_shape_only_after_shift(m)), ...
        max(F.chi2dof_shape_only_after_shift(m)),median(F.rms_shape_dev_after_shift(m)));
end

bad = R0.flag_blocks>=5;
fprintf('persistent columns %s, q %s\n',mat2str(R0.detector_col(bad).'), ...
    mat2str(R0.q(bad).',5));
i = S.n_channels>=30 & S.slope_meV_ch>=1.25 & S.slope_meV_ch<=3.25;
K = S(i,:); [~,j] = max(K.odd_train_score);
fprintf('stripe nearest selected %.2f %.2f train %.3f test %.3f\n',K{j,1:2},K.odd_train_score(j),K.odd_test_score(j));
i = S.n_channels>=30; K = S(i,:); [~,j] = max(K.odd_train_score);
fprintf('stripe best exploratory %.2f %.2f train %.3f test %.3f\n',K{j,1:2},K.odd_train_score(j),K.odd_test_score(j));
for t = 1:2
    for w = [100 150;150 250;250 350].'
        m = T.template==t&T.block==0&T.line_E_meV>=w(1)&T.line_E_meV<w(2);
        H = T(m,:);
        fprintf('stripe %d %g-%g q %s frac_med %.4f ridgecounts_med %.1f ratio90abs %.4f z_med %.3f gainfrac %.4f same-sign z>3 %d/%d\n', ...
            t,w(1),w(2),local_span(H.q),median(H.odd_fraction_plus_bg),median(H.odd_counts), ...
            prctile(abs(H.odd_fraction_plus_bg),90),median(H.approx_z),median(H.negative_E_odd_fraction,'omitnan'), ...
            nnz(H.approx_z>3),height(H));
    end
end
end

function s = local_span(x)
s = sprintf('%.5g..%.5g',min(x),max(x));
end
