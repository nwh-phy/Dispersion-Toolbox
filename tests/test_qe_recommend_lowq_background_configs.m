function tests = test_qe_recommend_lowq_background_configs
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
addpath(fullfile(project_root, 'scripts'));
testCase.TestData.project_root = project_root;
end


function testRanksCenterConservativeAndChoosesSafeDefault(testCase)
comparisons = makeSyntheticComparisons(false);

screening = qe_recommend_lowq_background_configs(comparisons);

verifyEqual(testCase, screening.rankings.center_conservative_names, ...
    {'dual_auto_plus_vii_group_0p01', 'dual_auto_plus_vii_group_0p01_pos_shoulder', 'dual_auto_plus_vii', 'dual_exppoly3'});
verifyEqual(testCase, screening.safe_center_names, ...
    {'dual_auto_plus_vii_group_0p01', 'dual_auto_plus_vii_group_0p01_pos_shoulder'});
verifyEqual(testCase, screening.recommended_center_name, 'dual_auto_plus_vii_group_0p01');
verifyEqual(testCase, screening.rankings.center_branch1_names, ...
    {'dual_exppoly3', 'dual_auto_plus_vii_group_0p01', 'dual_auto_plus_vii_group_0p01_pos_shoulder', 'dual_auto_plus_vii'});
end


function testEnablesCenterAndPositiveShoulderGroupingWhenMetricsImprove(testCase)
comparisons = makeSyntheticComparisons(false);

screening = qe_recommend_lowq_background_configs(comparisons);

verifyTrue(testCase, screening.decisions.center_grouping.available);
verifyTrue(testCase, screening.decisions.center_grouping.enabled);
verifyEqual(testCase, screening.decisions.center_grouping.baseline_name, 'dual_auto_plus_vii');
verifyEqual(testCase, screening.decisions.center_grouping.candidate_name, 'dual_auto_plus_vii_group_0p01');
verifyTrue(testCase, screening.decisions.center_grouping.center_neg_area_improvement_frac > 0.5);
verifyTrue(testCase, screening.decisions.positive_shoulder_grouping.available);
verifyTrue(testCase, screening.decisions.positive_shoulder_grouping.enabled);
verifyEqual(testCase, screening.recommended_positive_shoulder_name, 'dual_auto_plus_vii_group_0p01_pos_shoulder');
end


function testDisablesPositiveShoulderGroupingWhenShoulderGetsWorse(testCase)
comparisons = makeSyntheticComparisons(true);

screening = qe_recommend_lowq_background_configs(comparisons);

verifyTrue(testCase, screening.decisions.center_grouping.enabled);
verifyFalse(testCase, screening.decisions.positive_shoulder_grouping.enabled);
verifyEqual(testCase, screening.recommended_positive_shoulder_name, 'dual_auto_plus_vii_group_0p01');
end


function testHandlesZeroBaselineMetricsWithoutExplodingFractions(testCase)
comparisons = makeSyntheticComparisons(false);
comparisons{1}.regions.center_core.neg_area_fraction_mean = 0;
comparisons{2}.regions.center_core.neg_area_fraction_mean = 0.001;

screening = qe_recommend_lowq_background_configs(comparisons);

verifyFalse(testCase, screening.decisions.center_grouping.enabled);
verifyTrue(testCase, isnan(screening.decisions.center_grouping.center_neg_area_improvement_frac));
verifyTrue(testCase, isnan(screening.decisions.center_grouping.primary_neg_area_improvement_frac));
end


function testShoulderRecommendationStaysCoherentWhenCenterDefaultIsDifferent(testCase)
comparisons = makeSyntheticComparisons(false);
comparisons{4}.name = 'dual_pearsonvii';
comparisons{4}.regions.center_core.neg_area_fraction_mean = 0.0035;
comparisons{4}.regions.center_core.neg_peak_fraction_mean = 0.09;
comparisons{4}.regions.center_core.branch1_peak_mean = 130;
comparisons{4}.q0_summary.neg_area_fraction = 0.0035;

screening = qe_recommend_lowq_background_configs(comparisons);

verifyEqual(testCase, screening.recommended_center_name, 'dual_pearsonvii');
verifyEqual(testCase, screening.recommended_positive_shoulder_name, 'dual_pearsonvii');
end


function testCenterGroupingIsNotEnabledWhenAnotherConfigWinsOverall(testCase)
comparisons = makeSyntheticComparisons(false);
comparisons{4}.name = 'dual_pearsonvii';
comparisons{4}.regions.center_core.neg_area_fraction_mean = 0.0035;
comparisons{4}.regions.center_core.neg_peak_fraction_mean = 0.09;
comparisons{4}.regions.center_core.branch1_peak_mean = 130;
comparisons{4}.q0_summary.neg_area_fraction = 0.0035;

screening = qe_recommend_lowq_background_configs(comparisons);

verifyEqual(testCase, screening.recommended_center_name, 'dual_pearsonvii');
verifyFalse(testCase, screening.decisions.center_grouping.enabled);
end


function testShoulderGroupingRequiresShoulderBranchPreservation(testCase)
comparisons = makeSyntheticComparisons(false);
comparisons{3}.regions.positive_shoulder.branch1_peak_mean = 20;

screening = qe_recommend_lowq_background_configs(comparisons);

verifyFalse(testCase, screening.decisions.positive_shoulder_grouping.enabled);
verifyEqual(testCase, screening.recommended_positive_shoulder_name, 'dual_auto_plus_vii_group_0p01');
end


function testCompareBgDualWindowAllReturnsDatasetScreeningSummary(testCase)
project_root = testCase.TestData.project_root;
required_paths = { ...
    fullfile(project_root, '20260120 Bi', '590 PL2 10w 0.004 10sx300', 'eq3D.mat'), ...
    fullfile(project_root, '20260120 Bi', 'n0 pl2 10w 0.004 10s x300', 'eq3D.mat'), ...
    fullfile(project_root, '20260120 Bi', 'no pl2 20w 0.004 10sx300 2film', 'eq3D.mat') ...
    };
assumeTrue(testCase, all(cellfun(@(p) exist(p, 'file') == 2, required_paths)), ...
    'Representative low-q datasets missing for integration smoke test.');

results = compare_bg_dual_window_all();
first = results{1};

verifyTrue(testCase, isfield(first, 'screening'));
verifyTrue(testCase, isfield(first.screening, 'recommended_center_name'));
verifyTrue(testCase, isfield(first.screening, 'decisions'));
verifyTrue(testCase, isfield(first.screening.decisions, 'center_grouping'));
verifyTrue(testCase, isfield(first.screening.decisions, 'positive_shoulder_grouping'));
end


function comparisons = makeSyntheticComparisons(shoulder_worse)
comparisons = { ...
    makeComparison('dual_auto_plus_vii', 0.010, 0.20, 90, 0.020, 0.25, 0.018, 0.14, 0.010), ...
    makeComparison('dual_auto_plus_vii_group_0p01', 0.004, 0.10, 120, 0.021, 0.23, 0.017, 0.13, 0.004), ...
    makeComparison('dual_auto_plus_vii_group_0p01_pos_shoulder', 0.0042, 0.11, 119, 0.015, 0.20, 0.0175, 0.13, 0.0042), ...
    makeComparison('dual_exppoly3', 0.020, 0.60, 150, 0.030, 0.55, 0.010, 0.08, 0.020) ...
    };
if shoulder_worse
    comparisons{3}.regions.positive_shoulder.neg_area_fraction_mean = 0.024;
    comparisons{3}.regions.positive_shoulder.neg_peak_fraction_mean = 0.27;
end
end


function item = makeComparison(name, centerNegArea, centerNegPeak, centerBranch1, posNegArea, posNegPeak, negNegArea, negNegPeak, q0NegArea)
item = struct();
item.name = name;
item.regions = struct();
item.regions.center_core = struct( ...
    'neg_area_fraction_mean', centerNegArea, ...
    'neg_peak_fraction_mean', centerNegPeak, ...
    'branch1_peak_mean', centerBranch1);
item.regions.positive_shoulder = struct( ...
    'neg_area_fraction_mean', posNegArea, ...
    'neg_peak_fraction_mean', posNegPeak, ...
    'branch1_peak_mean', 80);
item.regions.negative_shoulder = struct( ...
    'neg_area_fraction_mean', negNegArea, ...
    'neg_peak_fraction_mean', negNegPeak, ...
    'branch1_peak_mean', 70);
item.q0_summary = struct('neg_area_fraction', q0NegArea);
end
