function tests = test_qe_summarize_lowq_background_eval
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
addpath(fullfile(project_root, 'scripts'));
testCase.TestData.project_root = project_root;
end


function testSeparatesCenterAndSignedShoulderRegions(testCase)
[qe_focus, qe_bg, bg_diag] = makeSyntheticBackgroundEvaluationFixture();

summary = qe_summarize_lowq_background_eval(qe_focus, qe_bg, bg_diag);

verifyEqual(testCase, summary.regions.full_low_q.channel_count, 9);
verifyEqual(testCase, summary.regions.center_core.channel_count, 3);
verifyEqual(testCase, summary.regions.positive_shoulder.channel_count, 2);
verifyEqual(testCase, summary.regions.negative_shoulder.channel_count, 2);
verifyEqual(testCase, summary.regions.center_core.q_min, -0.005, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.center_core.q_max, 0.005, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.positive_shoulder.q_min, 0.03, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.positive_shoulder.q_max, 0.04, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.negative_shoulder.q_min, -0.04, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.negative_shoulder.q_max, -0.02, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.center_core.selected_methods, {'PearsonVII'});
verifyEqual(testCase, summary.regions.center_core.selected_method_counts, 3);
verifyEqual(testCase, lookupMethodCount(summary.regions.positive_shoulder, 'Pearson'), 1);
verifyEqual(testCase, lookupMethodCount(summary.regions.positive_shoulder, 'PearsonVII'), 1);
verifyEqual(testCase, lookupMethodCount(summary.regions.negative_shoulder, 'Exp2'), 1);
verifyEqual(testCase, lookupMethodCount(summary.regions.negative_shoulder, 'Pearson'), 1);
end


function testComputesBranchPreservationMetricsPerRegion(testCase)
[qe_focus, qe_bg, bg_diag] = makeSyntheticBackgroundEvaluationFixture();

summary = qe_summarize_lowq_background_eval(qe_focus, qe_bg, bg_diag);

verifyEqual(testCase, summary.regions.full_low_q.branch1_peak_mean, 30, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.full_low_q.branch2_peak_mean, 25, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.full_low_q.branch3_peak_mean, 22, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.center_core.branch1_peak_mean, 45, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.center_core.branch2_peak_mean, 40, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.center_core.branch3_peak_mean, 37, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.center_core.branch1_peak_min, 40, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.positive_shoulder.branch1_peak_mean, 30, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.positive_shoulder.branch2_peak_min, 20, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.negative_shoulder.branch3_peak_mean, 17, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.center_core.neg_fraction_mean, 0.02, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.positive_shoulder.neg_peak_fraction_mean, 0.125, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.regions.negative_shoulder.neg_area_fraction_mean, 0.065, 'AbsTol', 1e-12);
end


function testCapturesQZeroSummaryAndAutoCandidateSnapshot(testCase)
[qe_focus, qe_bg, bg_diag] = makeSyntheticBackgroundEvaluationFixture();

summary = qe_summarize_lowq_background_eval(qe_focus, qe_bg, bg_diag);

verifyEqual(testCase, summary.q0.q_index, qe_focus.q_zero_index);
verifyEqual(testCase, summary.q0.q_Ainv, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.q0.selected_method, 'PearsonVII');
verifyEqual(testCase, summary.q0.neg_fraction, 0.02, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.q0.processed_min_after50, -3, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.q0.negative_points_after50, 1);
verifyEqual(testCase, summary.q0.branch1_peak, 50, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.q0.branch2_peak, 45, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.q0.branch3_peak, 42, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.q0.candidate_methods, {'Power', 'PearsonVII'});
verifyEqual(testCase, summary.q0.candidate_scores, [1.2, 0.8], 'AbsTol', 1e-12);
verifyEqual(testCase, summary.branch_windows.branch1_meV, [600 1700]);
verifyEqual(testCase, summary.evaluation_ranges.center_qmax, 0.01, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.evaluation_ranges.positive_shoulder_qrange, [0.02 0.05], 'AbsTol', 1e-12);
end


function testCompareBgDualWindow590ReturnsStructuredLowQSummary(testCase)
dataset_path = fullfile(testCase.TestData.project_root, '20260120 Bi', ...
    '590 PL2 10w 0.004 10sx300', 'eq3D.mat');
assumeTrue(testCase, exist(dataset_path, 'file') == 2, ...
    sprintf('Representative dataset missing: %s', dataset_path));

results = compare_bg_dual_window_590();
first = results{1};

verifyTrue(testCase, isfield(first, 'evaluation_ranges'));
verifyTrue(testCase, isfield(first, 'branch_windows'));
verifyTrue(testCase, isfield(first, 'regions'));
verifyTrue(testCase, isfield(first, 'q0_summary'));
verifyTrue(testCase, isfield(first.regions, 'center_core'));
verifyTrue(testCase, isfield(first.regions, 'positive_shoulder'));
verifyTrue(testCase, isfield(first.regions.center_core, 'branch1_peak_mean'));
verifyTrue(testCase, isfield(first.q0_summary, 'processed_min_after50'));
end


function testFallsBackToNearestZeroWhenQZeroIndexIsInvalid(testCase)
[qe_focus, qe_bg, bg_diag] = makeSyntheticBackgroundEvaluationFixture();
qe_focus.q_zero_index = 999;

summary = qe_summarize_lowq_background_eval(qe_focus, qe_bg, bg_diag);

verifyEqual(testCase, summary.q0.q_index, 5);
verifyEqual(testCase, summary.q0.q_Ainv, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, summary.q0.selected_method, 'PearsonVII');
end


function [qe_focus, qe_bg, bg_diag] = makeSyntheticBackgroundEvaluationFixture()
qe_focus = struct();
qe_focus.energy_meV = [0; 100; 1000; 2000; 3200];
qe_focus.q_Ainv = [-0.12 -0.04 -0.02 -0.005 0 0.005 0.03 0.04 0.12];
qe_focus.q_abs_Ainv = abs(qe_focus.q_Ainv);
qe_focus.q_zero_index = 5;

qe_bg = struct();
qe_bg.intensity = [
    0   0   0   0   0   0   0   0   0;
    1  -1  -2  -1  -3  -1  -2  -4  -1;
    10 20 30 40 50 45 35 25 15;
    5  15 25 35 45 40 30 20 10;
    2  12 22 32 42 37 27 17 7];

empty_diag = struct( ...
    'selected_method', '', ...
    'neg_fraction', NaN, ...
    'neg_area_fraction', NaN, ...
    'neg_peak_fraction', NaN, ...
    'bg_fraction', NaN, ...
    'linear_rmse', NaN, ...
    'candidate_methods', {{}}, ...
    'candidate_scores', [], ...
    'candidate_details', struct('linear_rmse', {}, 'score', {}));
bg_diag = repmat(empty_diag, 1, numel(qe_focus.q_Ainv));

methods = {'Power', 'Exp2', 'Pearson', 'PearsonVII', 'PearsonVII', 'PearsonVII', 'Pearson', 'PearsonVII', 'Power'};
neg_fraction = [0.11 0.07 0.08 0.01 0.02 0.03 0.05 0.06 0.12];
neg_area_fraction = [0.12 0.06 0.07 0.02 0.03 0.04 0.08 0.09 0.13];
neg_peak_fraction = [0.20 0.10 0.11 0.04 0.05 0.06 0.12 0.13 0.21];
bg_fraction = [0.90 0.88 0.89 0.86 0.85 0.87 0.91 0.92 0.93];
linear_rmse = [4.0 3.0 3.2 1.0 1.1 1.2 2.4 2.5 4.1];

for i = 1:numel(methods)
    bg_diag(i).selected_method = methods{i};
    bg_diag(i).neg_fraction = neg_fraction(i);
    bg_diag(i).neg_area_fraction = neg_area_fraction(i);
    bg_diag(i).neg_peak_fraction = neg_peak_fraction(i);
    bg_diag(i).bg_fraction = bg_fraction(i);
    bg_diag(i).linear_rmse = linear_rmse(i);
end

bg_diag(5).candidate_methods = {'Power', 'PearsonVII'};
bg_diag(5).candidate_scores = [1.2 0.8];
bg_diag(5).candidate_details = struct( ...
    'linear_rmse', {1.4, 1.0}, ...
    'score', {1.2, 0.8});
end


function count = lookupMethodCount(region, method)
idx = find(strcmp(region.selected_methods, method), 1);
if isempty(idx)
    count = 0;
else
    count = region.selected_method_counts(idx);
end
end
