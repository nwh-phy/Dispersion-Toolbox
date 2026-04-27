function tests = test_qe_background_subtraction
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testPowerBackgroundSubtractionRecoversSyntheticPeakHeight(testCase)
[qe, peak_only] = makeSyntheticQeDataset('power', [1.0]);
opts = makeBgOpts('Power');

[qe_out, bg_diag] = qe_preprocess(qe, opts);

[peak_true, idx_true] = max(peak_only(:, 1));
peak_measured = qe_out.intensity(idx_true, 1);

verifyGreaterThan(testCase, peak_measured, 0);
verifyLessThan(testCase, abs(peak_measured - peak_true) / peak_true, 0.20);
verifyLessThan(testCase, bg_diag(1).neg_fraction, 0.10);
verifyLessThan(testCase, bg_diag(1).neg_area_fraction, 0.05);
end


function testAutoBackgroundSelectionPrefersPowerForPowerLawTail(testCase)
[qe, ~] = makeSyntheticQeDataset('power', [1.0]);
opts = makeBgOpts('Auto');
opts.bg_candidate_methods = {'Power', 'Exp2', 'ExpPoly3'};

[~, bg_diag] = qe_preprocess(qe, opts);

verifyEqual(testCase, bg_diag(1).selected_method, 'Power');
verifyGreaterThanOrEqual(testCase, numel(bg_diag(1).candidate_methods), 3);
verifyGreaterThanOrEqual(testCase, numel(bg_diag(1).candidate_scores), 3);
end


function testAutoBackgroundSelectionSupportsLogDomainCandidates(testCase)
[qe, ~] = makeSyntheticQeDataset('power', [1.0]);
opts = makeBgOpts('Auto');
opts.bg_candidate_methods = {'Power', 'Pearson', 'ExpPoly3'};

[~, bg_diag] = qe_preprocess(qe, opts);

[best_score, best_idx] = min(bg_diag(1).candidate_scores);
verifyEqual(testCase, bg_diag(1).selected_method, bg_diag(1).candidate_methods{best_idx});
verifyEqual(testCase, {bg_diag(1).candidate_details.method}, ...
    {'Power', 'Pearson', 'ExpPoly3'});
verifyTrue(testCase, all(isfinite([bg_diag(1).candidate_details.linear_rmse])));
verifyTrue(testCase, isfinite(best_score));
verifyEqual(testCase, bg_diag(1).candidate_scores, ...
    [bg_diag(1).candidate_details.score]);
end


function testBackgroundDiagnosticsExposeFractionsAndSelectedMethod(testCase)
[qe, ~] = makeSyntheticQeDataset('power', [1.0, 0.85]);
opts = makeBgOpts('Auto');
opts.bg_candidate_methods = {'Power', 'Exp2'};

[~, bg_diag] = qe_preprocess(qe, opts);

verifyTrue(testCase, isfield(bg_diag, 'selected_method'));
verifyTrue(testCase, isfield(bg_diag, 'candidate_methods'));
verifyTrue(testCase, isfield(bg_diag, 'candidate_scores'));
verifyTrue(testCase, isfield(bg_diag, 'candidate_details'));
verifyTrue(testCase, isfield(bg_diag, 'linear_rmse'));
verifyTrue(testCase, isfield(bg_diag, 'neg_fraction'));
verifyTrue(testCase, isfield(bg_diag, 'neg_area_fraction'));
verifyTrue(testCase, isfield(bg_diag, 'bg_fraction'));

for qi = 1:numel(bg_diag)
    verifyGreaterThanOrEqual(testCase, bg_diag(qi).neg_fraction, 0);
    verifyLessThanOrEqual(testCase, bg_diag(qi).neg_fraction, 1);
    verifyGreaterThanOrEqual(testCase, bg_diag(qi).neg_area_fraction, 0);
    verifyLessThanOrEqual(testCase, bg_diag(qi).neg_area_fraction, 1);
    verifyGreaterThanOrEqual(testCase, bg_diag(qi).bg_fraction, 0);
    verifyLessThanOrEqual(testCase, bg_diag(qi).bg_fraction, 1);
    verifyGreaterThan(testCase, bg_diag(qi).linear_rmse, 0);
    verifyEqual(testCase, numel(bg_diag(qi).candidate_details), numel(bg_diag(qi).candidate_methods));
end
end


function testPearsonCandidateDetailsUseLinearDomainScore(testCase)
[qe, ~] = makeSyntheticQeDataset('pearson_like', [1.0]);
opts = makeBgOpts('Auto');
opts.bg_candidate_methods = {'Pearson', 'Power'};

[~, bg_diag] = qe_preprocess(qe, opts);

pearson_idx = find(strcmp(bg_diag(1).candidate_methods, 'Pearson'), 1);
power_idx = find(strcmp(bg_diag(1).candidate_methods, 'Power'), 1);
verifyNotEmpty(testCase, pearson_idx);
verifyNotEmpty(testCase, power_idx);
verifyGreaterThan(testCase, bg_diag(1).candidate_details(pearson_idx).linear_rmse, 0);
verifyGreaterThan(testCase, bg_diag(1).candidate_details(power_idx).linear_rmse, 0);
verifyTrue(testCase, isfinite(bg_diag(1).candidate_scores(pearson_idx)));
verifyTrue(testCase, isfinite(bg_diag(1).candidate_scores(power_idx)));
end


function testDualWindowBackgroundSubtractionRecoversLogLinearPeakHeight(testCase)
[qe, peak_only] = makeDualWindowSyntheticQeDataset([1.0]);
opts = makeBgOpts('DualWindow');
opts.bg_win_lo = [120, 360];
opts.bg_win_hi = [2600, 3200];

[qe_out, bg_diag] = qe_preprocess(qe, opts);

[peak_true, idx_true] = max(peak_only(:, 1));
peak_measured = qe_out.intensity(idx_true, 1);

verifyEqual(testCase, bg_diag(1).selected_method, 'DualWindow');
verifyLessThan(testCase, abs(peak_measured - peak_true) / peak_true, 0.12);
verifyLessThan(testCase, bg_diag(1).neg_area_fraction, 0.08);
end


function testAutoBackgroundSelectionCanSelectExplicitDualWindowCandidate(testCase)
[qe, ~] = makeDualWindowSyntheticQeDataset([1.0]);
opts = makeBgOpts('Auto');
opts.bg_win_lo = [120, 360];
opts.bg_win_hi = [2600, 3200];
opts.bg_candidate_methods = {'Power', 'DualWindow'};

[~, bg_diag] = qe_preprocess(qe, opts);

verifyTrue(testCase, any(strcmp(bg_diag(1).candidate_methods, 'DualWindow')));
verifyEqual(testCase, bg_diag(1).selected_method, 'DualWindow');
[~, best_idx] = min(bg_diag(1).candidate_scores);
verifyEqual(testCase, bg_diag(1).candidate_methods{best_idx}, 'DualWindow');
end


function testDefaultAutoCandidatesExcludeDualWindow(testCase)
[qe, ~] = makeDualWindowSyntheticQeDataset([1.0]);
opts = makeBgOpts('Auto');
opts.bg_win_lo = [120, 360];
opts.bg_win_hi = [2600, 3200];

[~, bg_diag] = qe_preprocess(qe, opts);

verifyFalse(testCase, any(strcmp(bg_diag(1).candidate_methods, 'DualWindow')));
end


function testDualWindowMissingSecondWindowFallsBackConservatively(testCase)
[qe, ~] = makeDualWindowSyntheticQeDataset([1.0]);
opts = makeBgOpts('DualWindow');
opts.bg_win_lo = [120, 360];
opts.bg_win_hi = [];

[qe_out, bg_diag] = qe_preprocess(qe, opts);

verifyEqual(testCase, bg_diag(1).selected_method, 'DualWindow');
verifyTrue(testCase, all(isfinite(qe_out.intensity(:))));
verifyTrue(testCase, all(isfinite(bg_diag(1).bg_curve(:))));
verifyGreaterThanOrEqual(testCase, bg_diag(1).bg_fraction, 0);
verifyLessThanOrEqual(testCase, bg_diag(1).bg_fraction, 1);
end


function testAutoDualWindowCandidateAcceptsStringArray(testCase)
[qe, ~] = makeDualWindowSyntheticQeDataset([1.0]);
opts = makeBgOpts('Auto');
opts.bg_win_lo = [120, 360];
opts.bg_win_hi = [2600, 3200];
opts.bg_candidate_methods = ["Power", "DualWindow"];

[~, bg_diag] = qe_preprocess(qe, opts);

verifyEqual(testCase, bg_diag(1).candidate_methods, {'Power', 'DualWindow'});
verifyEqual(testCase, bg_diag(1).selected_method, 'DualWindow');
end


function opts = makeBgOpts(method)
opts = struct();
opts.do_despike = false;
opts.do_normalize = false;
opts.do_denoise = false;
opts.do_bg_sub = true;
opts.bg_method = method;
opts.bg_win_lo = [50, 300];
opts.bg_win_hi = [];
opts.bg_iterative = false;
opts.do_deconv = false;
end


function [qe, peak_only] = makeSyntheticQeDataset(background_kind, scales)
energy_meV = (-100:10:4000)';
q_channel = 1:numel(scales);
intensity = zeros(numel(energy_meV), numel(scales));
peak_only = zeros(numel(energy_meV), numel(scales));

for qi = 1:numel(scales)
    scale = scales(qi);
    switch char(background_kind)
        case 'power'
            background = syntheticPowerBackground(energy_meV, scale);
        case 'pearson_like'
            background = syntheticPearsonLikeBackground(energy_meV, scale);
        otherwise
            error('Unknown background kind: %s', background_kind);
    end
    peak = syntheticLorentzPeak(energy_meV, 1200, 180, 0.24 * scale);
    shoulder = syntheticLorentzPeak(energy_meV, 1850, 260, 0.05 * scale);
    intensity(:, qi) = background + peak + shoulder;
    peak_only(:, qi) = peak + shoulder;
end

qe = make_qe_struct(intensity, energy_meV, 0.005, q_channel, ...
    source_kind="synthetic", label="synthetic-bg-test");
end


function [qe, peak_only] = makeDualWindowSyntheticQeDataset(scales)
energy_meV = (-100:10:4000)';
q_channel = 1:numel(scales);
intensity = zeros(numel(energy_meV), numel(scales));
peak_only = zeros(numel(energy_meV), numel(scales));

for qi = 1:numel(scales)
    scale = scales(qi);
    background = syntheticLogLinearBackground(energy_meV, scale);
    peak = syntheticLorentzPeak(energy_meV, 1350, 170, 0.30 * scale);
    shoulder = syntheticLorentzPeak(energy_meV, 1900, 220, 0.035 * scale);
    intensity(:, qi) = background + peak + shoulder;
    peak_only(:, qi) = peak + shoulder;
end

qe = make_qe_struct(intensity, energy_meV, 0.005, q_channel, ...
    source_kind="synthetic", label="synthetic-dualwindow-bg-test");
end


function bg = syntheticPowerBackground(energy_meV, scale)
energy_pos = max(energy_meV, 1) + 25;
bg = 18 * scale * energy_pos .^ (-1.32);
bg(energy_meV < 0) = max(bg(energy_meV < 0), 0.02 * scale);
end


function bg = syntheticLogLinearBackground(energy_meV, scale)
anchor1_e = 240;
anchor2_e = 2900;
anchor1_log_i = log(0.16 * scale);
anchor2_log_i = log(0.030 * scale);
slope = (anchor2_log_i - anchor1_log_i) / (anchor2_e - anchor1_e);
bg = exp(anchor1_log_i + slope .* (energy_meV - anchor1_e));
bg(energy_meV < 0) = max(bg(energy_meV < 0), 0.18 * scale);
end


function bg = syntheticPearsonLikeBackground(energy_meV, scale)
energy_pos = max(energy_meV, 1) + 30;
log_e = log(energy_pos);
log_bg = 2.7 - 0.95 * log_e - 0.12 * (log_e - mean(log_e)).^2;
bg = scale * exp(log_bg);
bg(energy_meV < 0) = max(bg(energy_meV < 0), 0.015 * scale);
end


function peak = syntheticLorentzPeak(energy_meV, center, gamma, amplitude)
x = energy_meV - center;
peak = amplitude ./ (1 + (x ./ gamma) .^ 2);
end
