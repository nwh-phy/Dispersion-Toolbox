function tests = test_b1_double_peak_binning_workflow
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testMandatoryDoublePeakExtractionCoversLowQ(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(false);
opts = local_default_opts(old_points);
opts.noise_threshold = Inf;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

verifyEqual(testCase, height(result.lower_points), numel(qe.q_Ainv));
verifyEqual(testCase, height(result.upper_points), numel(qe.q_Ainv));
verifyEqual(testCase, height(result.fit_failures), 0);
verifyTrue(testCase, all(strcmp(result.lower_points.source_mode, ...
    'single_q_direct')));
verifyTrue(testCase, all(result.lower_points.q_abs_Ainv <= 0.05));
verifyTrue(testCase, all(result.lower_points.energy_meV < ...
    result.upper_points.energy_meV));
end


function testLowSnrRegionBinsBeforeDoublePeakFit(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(true);
opts = local_default_opts(old_points);
opts.noise_threshold = 0.08;
opts.bin_size = 3;
opts.low_q_no_bin_abs_Ainv = 0;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

combined = result.binning_map(strcmp(result.binning_map.source_mode, ...
    'combined_q_binning_3'), :);
verifyEqual(testCase, height(combined), 1);
verifyEqual(testCase, combined.source_q_count(1), 3);
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.03'));
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.05'));
verifyEqual(testCase, height(result.fit_failures), 0);

combined_points = result.combined_points(strcmp( ...
    result.combined_points.source_mode, 'combined_q_binning_3'), :);
verifyEqual(testCase, height(combined_points), 2);
verifyTrue(testCase, all(combined_points.source_q_count == 3));
end


function testLowQNoiseIsProtectedFromBinning(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(true);
opts = local_default_opts(old_points);
opts.noise_threshold = 0.08;
opts.bin_size = 3;
opts.low_q_no_bin_abs_Ainv = 0.03;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

protected = result.noise_profile(result.noise_profile.low_q_no_bin_protected, :);
verifyNotEmpty(testCase, protected);
verifyTrue(testCase, all(~protected.use_binning));

combined = result.binning_map(strcmp(result.binning_map.source_mode, ...
    'combined_q_binning_3'), :);
verifyEqual(testCase, height(combined), 1);
verifyEqual(testCase, combined.source_q_count(1), 2);
verifyEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.03'));
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.04'));
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.05'));
end


function testNoiseThresholdIsAutoPerDatasetUnlessOverridden(testCase)
[qe_clean, qe_clean_raw, old_clean] = local_synthetic_qe(false);
clean_opts = local_default_opts(old_clean);
clean_opts.noise_threshold = NaN;
clean_opts.low_q_no_bin_abs_Ainv = 0;
clean = b1_double_peak_binning_extract(qe_clean, qe_clean_raw, clean_opts);

[qe_noisy, qe_noisy_raw, old_noisy] = local_synthetic_qe(true);
noisy_opts = local_default_opts(old_noisy);
noisy_opts.noise_threshold = NaN;
noisy_opts.low_q_no_bin_abs_Ainv = 0;
noisy = b1_double_peak_binning_extract(qe_noisy, qe_noisy_raw, noisy_opts);

verifyTrue(testCase, ismember('noise_threshold_source', ...
    clean.noise_profile.Properties.VariableNames));
verifyTrue(testCase, ismember('noise_threshold_source', ...
    noisy.noise_profile.Properties.VariableNames));
verifyTrue(testCase, all(strcmp(clean.noise_profile.noise_threshold_source, ...
    'auto_per_dataset')));
verifyTrue(testCase, all(strcmp(noisy.noise_profile.noise_threshold_source, ...
    'auto_per_dataset')));
verifyEqual(testCase, numel(unique(clean.noise_profile.noise_threshold)), 1);
verifyEqual(testCase, numel(unique(noisy.noise_profile.noise_threshold)), 1);
verifyNotEqual(testCase, clean.noise_profile.noise_threshold(1), ...
    noisy.noise_profile.noise_threshold(1));

manual_opts = local_default_opts(old_noisy);
manual_opts.noise_threshold = 0.123;
manual = b1_double_peak_binning_extract(qe_noisy, qe_noisy_raw, manual_opts);
verifyTrue(testCase, all(strcmp(manual.noise_profile.noise_threshold_source, ...
    'manual_override')));
verifyEqual(testCase, unique(manual.noise_profile.noise_threshold), 0.123);
end


function testDenoisedFitPathRecordsFitSpectrumSource(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(true);
opts = local_default_opts(old_points);
opts.noise_threshold = Inf;
opts.fit_denoise_method = 'sgolay';
opts.fit_denoise_window = 11;
opts.fit_denoise_order = 3;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

verifyEqual(testCase, height(result.fit_failures), 0);
verifyTrue(testCase, ismember('fit_spectrum_source', ...
    result.lower_points.Properties.VariableNames));
verifyTrue(testCase, all(strcmp(result.lower_points.fit_spectrum_source, ...
    'denoised')));
verifyTrue(testCase, all(strcmp(result.upper_points.fit_denoise_method, ...
    'sgolay')));
verifyTrue(testCase, all(result.binning_map.fit_denoise_window == 11));

first_detail = result.fit_details{1};
verifyEqual(testCase, first_detail.fit_spectrum_source, 'denoised');
verifyEqual(testCase, first_detail.fit_denoise_method, 'sgolay');
verifyGreaterThan(testCase, first_detail.fit_spectrum_delta_rms, 0);
end


function testAdaptiveAbsQDenoiseUsesStrongerHighQWindow(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(true);
opts = local_default_opts(old_points);
opts.noise_threshold = Inf;
opts.fit_denoise_method = 'sgolay';
opts.fit_denoise_profile = 'adaptive_absq';
opts.fit_denoise_low_window = 11;
opts.fit_denoise_high_window = 31;
opts.fit_denoise_q_start_Ainv = 0.02;
opts.fit_denoise_q_end_Ainv = 0.05;
opts.fit_denoise_order = 3;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

low_q = result.binning_map.q_abs_Ainv <= 0.02;
high_q = result.binning_map.q_abs_Ainv >= 0.05;
verifyTrue(testCase, any(low_q));
verifyTrue(testCase, any(high_q));
verifyTrue(testCase, all(result.binning_map.fit_denoise_window(low_q) == 11));
verifyTrue(testCase, all(result.binning_map.fit_denoise_window(high_q) == 31));
verifyTrue(testCase, all(strcmp(result.binning_map.fit_denoise_profile, ...
    'adaptive_absq')));
end


function testHighQForceBinningOverridesConservativeNoiseThreshold(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(true);
opts = local_default_opts(old_points);
opts.noise_threshold = Inf;
opts.low_q_no_bin_abs_Ainv = 0.02;
opts.high_q_force_bin_abs_Ainv = 0.03;
opts.bin_size = 3;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

verifyTrue(testCase, ismember('high_q_force_binning', ...
    result.noise_profile.Properties.VariableNames));
protected = result.noise_profile(result.noise_profile.q_abs_Ainv <= 0.02, :);
verifyTrue(testCase, all(~protected.use_binning));

forced = result.noise_profile(result.noise_profile.q_abs_Ainv >= 0.03, :);
verifyTrue(testCase, all(forced.high_q_force_binning));
verifyTrue(testCase, all(forced.use_binning));

combined = result.binning_map(strcmp(result.binning_map.source_mode, ...
    'combined_q_binning_3'), :);
verifyEqual(testCase, height(combined), 1);
verifyEqual(testCase, combined.source_q_count(1), 3);
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.03'));
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.05'));
end


function testBalancedHighQBin7KeepsLowQProtectedAndRecordsProvenance(testCase)
[qe, qe_raw, old_points] = local_synthetic_highq_binning_qe();
opts = local_default_opts(old_points);
opts.energy_window_meV = [500 1700];
opts.q_range_Ainv = [0 0.15];
opts.noise_threshold = Inf;
opts.low_q_no_bin_abs_Ainv = 0.05;
opts.high_q_force_bin_abs_Ainv = 0.08;
opts.bin_size = 7;
opts.peak_model = 'lorentz';
opts.tracking_mode = 'propagated_double_peak';
opts.fit_denoise_method = 'sgolay';
opts.fit_denoise_profile = 'adaptive_absq';
opts.fit_denoise_low_window = 11;
opts.fit_denoise_high_window = 71;
opts.fit_denoise_q_start_Ainv = 0.07;
opts.fit_denoise_q_end_Ainv = 0.15;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

protected = result.noise_profile(result.noise_profile.q_abs_Ainv <= 0.05, :);
verifyNotEmpty(testCase, protected);
verifyTrue(testCase, all(protected.low_q_no_bin_protected));
verifyTrue(testCase, all(~protected.use_binning));

forced = result.noise_profile(result.noise_profile.q_abs_Ainv >= 0.08, :);
verifyNotEmpty(testCase, forced);
verifyTrue(testCase, all(forced.high_q_force_binning));
verifyTrue(testCase, all(forced.use_binning));

combined = result.binning_map(strcmp(result.binning_map.source_mode, ...
    'combined_q_binning_7'), :);
verifyEqual(testCase, height(combined), 1);
verifyEqual(testCase, combined.source_q_count(1), 7);
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.08'));
verifyNotEmpty(testCase, strfind(combined.source_q_Ainv{1}, '0.14'));
verifyEqual(testCase, combined.bin_size_requested(1), 7);
verifyGreaterThan(testCase, combined.fit_denoise_window(1), 11);
verifyLessThanOrEqual(testCase, combined.fit_denoise_window(1), 71);

verifyTrue(testCase, all(strcmp(result.combined_points.peak_model, ...
    'lorentz')));
verifyTrue(testCase, all(strcmp(result.binning_map.tracking_mode, ...
    'propagated_double_peak')));
end


function testLorentzPropagatedTrackingRecordsModeAndModel(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(true);
opts = local_default_opts(old_points);
opts.noise_threshold = Inf;
opts.peak_model = 'lorentz';
opts.tracking_mode = 'propagated_double_peak';
opts.max_tracking_shift_meV = 180;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

verifyEqual(testCase, height(result.fit_failures), 0);
verifyTrue(testCase, ismember('peak_model', ...
    result.combined_points.Properties.VariableNames));
verifyTrue(testCase, ismember('tracking_mode', ...
    result.binning_map.Properties.VariableNames));
verifyTrue(testCase, all(strcmp(result.combined_points.peak_model, ...
    'lorentz')));
verifyTrue(testCase, all(strcmp(result.binning_map.tracking_mode, ...
    'propagated_double_peak')));
end


function testPropagatedTrackingReducesSpuriousJump(testCase)
[qe, qe_raw, old_points] = local_synthetic_tracking_challenge();
independent_opts = local_default_opts(old_points);
independent_opts.noise_threshold = Inf;
independent_opts.peak_model = 'lorentz';
independent_opts.tracking_mode = 'independent_double_peak';
independent_opts.min_prominence = 0.01;
independent_opts.energy_window_meV = [500 2000];

tracked_opts = independent_opts;
tracked_opts.tracking_mode = 'propagated_double_peak';
tracked_opts.max_tracking_shift_meV = 180;

independent = b1_double_peak_binning_extract(qe, qe_raw, independent_opts);
tracked = b1_double_peak_binning_extract(qe, qe_raw, tracked_opts);

independent_jumps = local_large_jump_count(independent.upper_points);
tracked_jumps = local_large_jump_count(tracked.upper_points);

verifyGreaterThan(testCase, independent_jumps, 0);
verifyLessThanOrEqual(testCase, tracked_jumps, independent_jumps);
end


function testWorkflowScriptDeclaresRequiredOutputsAndNoFallback(testCase)
script_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'scripts', ...
    'run_b1_double_peak_binning_analysis.m');

verifyTrue(testCase, isfile(script_path));
src = fileread(script_path);

required = { ...
    'b1_double_peak_combined_q_points.csv', ...
    'b1_double_peak_lower_points.csv', ...
    'b1_double_peak_upper_points.csv', ...
    'b1_double_peak_binning_map.csv', ...
    'b1_double_peak_fit_failures.csv', ...
    'b1_double_peak_binning_spectrum_comparison.png', ...
    'b1_double_peak_heatmap_overlay.png', ...
    'b1_double_peak_stacked_spectra_signed_q.png', ...
    'b1_double_peak_waterfall_signed_q', ...
    'b1_double_peak_fit_spectrum_denoise_comparison.png', ...
    'b1_double_peak_manual_window_seed.csv', ...
    'b1_double_peak_plot_q_binning_map.csv', ...
    'b1_single_peak_vs_double_peak_comparison.png'};

for i = 1:numel(required)
    verifyNotEmpty(testCase, strfind(src, required{i}), ...
        sprintf('Missing required workflow reference: %s', required{i}));
end
verifyEmpty(testCase, strfind(src, 'fallback_to_single'));
verifyNotEmpty(testCase, strfind(src, ...
    'options.qRangeOverride_Ainv (1,2) double = [-0.15 0.15]'));
verifyNotEmpty(testCase, strfind(src, ...
    'options.b1EnergyWindowOverrideMeV (1,2) double = [300 2100]'));
verifyNotEmpty(testCase, strfind(src, ...
    'opts.energy_window_meV = sort(run_options.b1EnergyWindowOverrideMeV)'));
verifyNotEmpty(testCase, strfind(src, 'local_binning_map_q_groups'));
verifyNotEmpty(testCase, strfind(src, 'single_q_direct'));
verifyNotEmpty(testCase, strfind(src, 'combined_q_binning_3'));
verifyNotEmpty(testCase, strfind(src, 'low_q_no_bin_abs_Ainv'));
verifyNotEmpty(testCase, strfind(src, 'source_q_count'));
verifyNotEmpty(testCase, strfind(src, 'source_mode'));
verifyNotEmpty(testCase, strfind(src, 'options.fitDenoiseMethod'));
verifyNotEmpty(testCase, strfind(src, 'options.fitDenoiseProfile'));
verifyNotEmpty(testCase, strfind(src, 'options.fitDenoiseHighWindow'));
verifyNotEmpty(testCase, strfind(src, 'options.highQForceBinAbsAinv'));
verifyNotEmpty(testCase, strfind(src, 'options.binSize (1,1) double = 3'));
verifyNotEmpty(testCase, strfind(src, 'opts.bin_size = run_options.binSize'));
verifyNotEmpty(testCase, strfind(src, 'options.peakModelOverride'));
verifyNotEmpty(testCase, strfind(src, 'options.trackingMode'));
verifyNotEmpty(testCase, strfind(src, 'options.fallbackSplitCandidatesMeV'));
verifyNotEmpty(testCase, strfind(src, 'options.maxTrackingShiftMeV'));
verifyNotEmpty(testCase, strfind(src, 'options.trackingWindowHalfWidthMeV'));
verifyNotEmpty(testCase, strfind(src, 'fit_denoise_method'));
end


function testResultIndexDocumentsCurrent300WindowAndDeprecated500Outputs(testCase)
index_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'RESULTS_INDEX.md');

verifyTrue(testCase, isfile(index_path));
src = fileread(index_path);

verifyNotEmpty(testCase, strfind(src, ...
    '00_DEPRECATED_B1_500MEV_WINDOW_260510'));
verifyNotEmpty(testCase, strfind(src, ...
    'old 500 meV lower-bound window'));
verifyNotEmpty(testCase, strfind(src, ...
    'Current clean 300 meV-window upper-stability diagnostic'));
verifyNotEmpty(testCase, strfind(src, ...
    'U1 builds a fresh propagated Lorentz seed'));
verifyNotEmpty(testCase, strfind(src, ...
    '300-1800 meV'));
verifyNotEmpty(testCase, strfind(src, ...
    'No new physical-fit directory is generated'));
end


function testB1FitScriptsExposeCustomBranchInputsWithoutChangingDefaults(testCase)
analysis_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'scripts', ...
    'run_b1_physical_fit_analysis.m');
enhancement_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'scripts', ...
    'run_b1_physical_fit_enhancements.m');

analysis_src = fileread(analysis_path);
enhancement_src = fileread(enhancement_path);

verifyNotEmpty(testCase, strfind(analysis_src, 'options.branchFileName'));
verifyNotEmpty(testCase, strfind(analysis_src, 'options.outputTag'));
verifyNotEmpty(testCase, strfind(analysis_src, 'options.filePrefix'));
verifyNotEmpty(testCase, strfind(analysis_src, 'branch1_points.csv'));
verifyNotEmpty(testCase, strfind(analysis_src, ...
    'b1_physical_fit_points_qabs.csv'));

verifyNotEmpty(testCase, strfind(enhancement_src, 'options.branchFileName'));
verifyNotEmpty(testCase, strfind(enhancement_src, 'options.outputTag'));
verifyNotEmpty(testCase, strfind(enhancement_src, 'options.filePrefix'));
verifyNotEmpty(testCase, strfind(enhancement_src, 'branch1_points.csv'));
verifyNotEmpty(testCase, strfind(enhancement_src, ...
    'b1_enhancement_points_qabs.csv'));
end


function testBalancedHighQBin7LorentzTrackingScriptDeclaresPlan(testCase)
script_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'scripts', ...
    'run_b1_lorentz_tracking_balanced_highqbin7_sg71.m');

verifyTrue(testCase, isfile(script_path));
src = fileread(script_path);

required = { ...
    '260510_lorentz_tracking_balanced_highqbin7_sg71', ...
    'runFits=false', ...
    'peakModelOverride=''lorentz''', ...
    'trackingMode=''propagated_double_peak''', ...
    'trackingMode=''windowed_branch_tracking''', ...
    'binSize=7', ...
    'highQForceBinAbsAinv=0.08', ...
    'fitDenoiseHighWindow=71', ...
    'fitDenoiseLowWindow=11', ...
    'waterfallAreaNormWindowMeV=[50 3800]', ...
    'waterfallStartMeV=250', ...
    'waterfallEndMeV=1600'};
for i = 1:numel(required)
    verifyNotEmpty(testCase, strfind(src, required{i}), ...
        sprintf('Missing V5 plan token: %s', required{i}));
end
verifyEmpty(testCase, strfind(src, 'run_b1_physical_fit'));
end


function testWindowInvalidFallbackUsesDoublePeakAndLogsRepair(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(false);
opts = local_default_opts(old_points);
opts.noise_threshold = Inf;
opts.peak_model = 'lorentz';
opts.tracking_mode = 'windowed_branch_tracking';
opts.tracking_window_half_width_meV = 10;
opts.tracking_window_highq_half_width_meV = 10;
opts.reference_lower_points = table(qe.q_Ainv, repmat(100, numel(qe.q_Ainv), 1), ...
    'VariableNames', {'q_Ainv', 'energy_meV'});
opts.reference_upper_points = table(qe.q_Ainv, repmat(180, numel(qe.q_Ainv), 1), ...
    'VariableNames', {'q_Ainv', 'energy_meV'});
opts.tracking_window_invalid_fallback = 'independent_double_peak';

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

verifyEqual(testCase, height(result.fit_failures), 0);
verifyEqual(testCase, height(result.lower_points), numel(qe.q_Ainv));
verifyTrue(testCase, ismember('repair_source', ...
    result.lower_points.Properties.VariableNames));
verifyTrue(testCase, any(strcmp(result.lower_points.repair_source, ...
    'tracking_window_invalid_independent_double_peak')));
verifyTrue(testCase, any(strcmp(result.binning_map.fit_status, ...
    'ok_repaired_tracking_window_invalid')));
verifyNotEmpty(testCase, result.repair_log);
verifyTrue(testCase, all(strcmp(result.repair_log.repair_source, ...
    'tracking_window_invalid_independent_double_peak')));
verifyTrue(testCase, all(strcmp(result.combined_points.peak_model, ...
    'lorentz')));
end


function testJumpRepairKeepsPointAndLogsOriginalEnergy(testCase)
points = local_repair_points_table([0.01; 0.02; 0.03; 0.04], ...
    [820; 840; 1450; 880], 1, 'b1_double_peak_lower');
failures = table();

repair = b1_double_peak_repair_tracking_points(points, table(), failures, ...
    largeJumpThresholdMeV=250);

verifyEqual(testCase, height(repair.lower_points), 4);
verifyEqual(testCase, repair.lower_points.energy_meV(3), 860, 'AbsTol', 1e-9);
verifyEqual(testCase, repair.lower_points.original_energy_meV(3), 1450, ...
    'AbsTol', 1e-9);
verifyEqual(testCase, repair.lower_points.repair_source{3}, ...
    'jump_repair_interpolated_energy');
verifyEqual(testCase, height(repair.repair_log), 1);
verifyEqual(testCase, repair.repair_log.original_energy_meV(1), 1450, ...
    'AbsTol', 1e-9);
verifyEqual(testCase, repair.repair_log.repaired_energy_meV(1), 860, ...
    'AbsTol', 1e-9);
verifyEqual(testCase, height(repair.exclusion_points), 0);
end


function testUntil300ScriptUsesHardGateAndDoesNotFitAboveGate(testCase)
script_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'scripts', ...
    'run_b1_lorentz_tracking_until300.m');

verifyTrue(testCase, isfile(script_path));
src = fileread(script_path);

required = { ...
    'targetScore (1,1) double = 300', ...
    'b1_lorentz_tracking_optimization_260510_until300', ...
    '260510_lorentz_tracking_until300', ...
    'runFits=false', ...
    'peakModelOverride=''lorentz''', ...
    'trackingWindowInvalidFallback=''independent_double_peak''', ...
    'enableJumpRepair=true', ...
    'largeJumpThresholdMeV=250', ...
    'final_candidate_fit', ...
    'if best_score <= options.targetScore', ...
    'run_b1_physical_fit_analysis', ...
    'run_b1_physical_fit_enhancements'};
for i = 1:numel(required)
    verifyNotEmpty(testCase, strfind(src, required{i}), ...
        sprintf('Missing until300 token: %s', required{i}));
end
verifyEmpty(testCase, strfind(src, '432'));
end


function testUpperStabilityMetricsDetectMediumJumpsAndOverbroadUpper(testCase)
points = local_repair_points_table([0.01; 0.02; 0.03; 0.04], ...
    [1000; 1130; 1255; 1440], 2, 'b1_double_peak_upper');
points.gamma_meV = [900; 1700; 900; 2100];

metrics_20w = b1_double_peak_upper_stability_metrics(points, ...
    'no_PL2_20w_2film');
verifyEqual(testCase, metrics_20w.upper_medium_jump_count, 3);
verifyEqual(testCase, metrics_20w.upper_large_jump_count, 1);
verifyEqual(testCase, metrics_20w.upper_overbroad_count, 2);
verifyEqual(testCase, metrics_20w.upper_medium_jump_threshold_meV, 120);

metrics_590 = b1_double_peak_upper_stability_metrics(points, ...
    '590_PL2_10w');
verifyEqual(testCase, metrics_590.upper_medium_jump_count, 1);
verifyEqual(testCase, metrics_590.upper_large_jump_count, 1);
verifyEqual(testCase, metrics_590.upper_medium_jump_threshold_meV, 150);
end


function testUpperQualityRetryRejectsOverbroadLorentzWithoutSingleFallback(testCase)
[qe, qe_raw, old_points] = local_synthetic_qe(false);
opts = local_default_opts(old_points);
opts.noise_threshold = Inf;
opts.peak_model = 'lorentz';
opts.tracking_mode = 'windowed_branch_tracking';
opts.reference_lower_points = table(old_points.q_Ainv, ...
    old_points.energy_meV - 95, ...
    'VariableNames', {'q_Ainv', 'energy_meV'});
opts.reference_upper_points = table(old_points.q_Ainv, ...
    old_points.energy_meV + 95, ...
    'VariableNames', {'q_Ainv', 'energy_meV'});
opts.tracking_window_half_width_meV = 180;
opts.tracking_window_highq_half_width_meV = 180;
opts.upper_tracking_window_half_width_meV = 120;
opts.upper_tracking_window_highq_half_width_meV = 120;
opts.upper_quality_retry = true;
opts.upper_max_gamma_meV = 5;
opts.upper_max_gamma_over_E = 0.001;
opts.upper_retry_window_half_width_meV = 80;
opts.upper_retry_highq_half_width_meV = 80;

result = b1_double_peak_binning_extract(qe, qe_raw, opts);

verifyGreaterThan(testCase, height(result.fit_failures), 0);
verifyTrue(testCase, all(strcmp(result.fit_failures.status, ...
    'overbroad_upper_peak')));
verifyTrue(testCase, all(strcmp(result.fit_failures.peak_model, ...
    'lorentz')));
verifyTrue(testCase, all(strcmp(result.fit_failures.tracking_mode, ...
    'windowed_branch_tracking')));
verifyEqual(testCase, height(result.lower_points), 0);
verifyEqual(testCase, height(result.upper_points), 0);
verifyEqual(testCase, height(result.combined_points), 0);
end


function testUpperStabilityScriptDeclaresGateAndNoFitAboveGate(testCase)
script_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'scripts', ...
    'run_b1_lorentz_tracking_upper_stability.m');

verifyTrue(testCase, isfile(script_path));
src = fileread(script_path);

required = { ...
    'b1_lorentz_tracking_optimization_260510_upper_stability', ...
    'u1_window300_2100_propagated_seed', ...
    '260510_lorentz_tracking_upper_stability_u1_window300_2100_propagated_seed', ...
    'propagated_double_peak', ...
    'tracking_mode', ...
    'trackingMode=variant.tracking_mode', ...
    'Old 500 meV lower-bound results are deprecated', ...
    'candidate_for_recommendation', ...
    'b1_energy_window', ...
    'b1EnergyWindowOverrideMeV=variant.b1_energy_window', ...
    '[300 2100]', ...
    '[300 2000]', ...
    '[300 1800]', ...
    'u6_window300_2100_repair180', ...
    'largeJumpThresholdMeV=variant.large_jump_threshold_meV', ...
    'B1 double-peak fit candidate energy windows: ', ...
    'local_copy_variant_overlays', ...
    'upper_medium_jump_count', ...
    'upperMaxGammaOverE=1.4', ...
    'upperMaxGammaMeV=1600', ...
    'upperTrackingWindowHalfWidthMeV=180', ...
    'upperTrackingWindowHighQHalfWidthMeV=240', ...
    'visualGateApproved', ...
    'if gate_passed && options.visualGateApproved && options.runSandboxFit', ...
    'upper_stability_final_candidate_fit', ...
    'b1_double_peak_upper_manual_anchor_template.csv'};
for i = 1:numel(required)
    verifyNotEmpty(testCase, strfind(src, required{i}), ...
        sprintf('Missing upper-stability token: %s', required{i}));
end
verifyEmpty(testCase, strfind(src, 'peakModelOverride=''fano'''));
verifyEmpty(testCase, strfind(src, 'u1_baseline_diagnostics'));
end


function [qe, qe_raw, old_points] = local_synthetic_qe(with_highq_noise)
energy = (400:10:2200).';
q = (0.01:0.01:0.05).';
intensity = zeros(numel(energy), numel(q));

for i = 1:numel(q)
    center = 850 + 4200 * q(i);
    lower = center - 95;
    upper = center + 95;
    y = 0.010 + 0.95 * exp(-0.5 * ((energy - lower) ./ 45) .^ 2) + ...
        0.78 * exp(-0.5 * ((energy - upper) ./ 55) .^ 2);
    if with_highq_noise && q(i) >= 0.03
        edge_mask = energy < 650 | energy > 1350;
        ripple = 0.16 * sign(sin((1:nnz(edge_mask)).'));
        y(edge_mask) = y(edge_mask) + ripple;
    end
    intensity(:, i) = y;
end

qe = struct('energy_meV', energy, 'q_Ainv', q, 'intensity', intensity);
qe_raw = qe;
old_points = table(q, 850 + 4200 * q, ...
    'VariableNames', {'q_Ainv', 'energy_meV'});
end


function [qe, qe_raw, old_points] = local_synthetic_tracking_challenge()
energy = (400:10:2200).';
q = (0.01:0.01:0.08).';
intensity = zeros(numel(energy), numel(q));

for i = 1:numel(q)
    lower = 760 + 1700 * q(i);
    upper = 1210 + 1200 * q(i);
    y = 0.010 + 0.80 * local_gaussian(energy, lower, 45) + ...
        0.40 * local_gaussian(energy, upper, 55);
    if i == 5
        y = y + 2.00 * local_gaussian(energy, 1850, 35);
    end
    intensity(:, i) = y;
end

qe = struct('energy_meV', energy, 'q_Ainv', q, 'intensity', intensity);
qe_raw = qe;
old_points = table(q, 980 + 1300 * q, ...
    'VariableNames', {'q_Ainv', 'energy_meV'});
end


function [qe, qe_raw, old_points] = local_synthetic_highq_binning_qe()
energy = (400:10:2200).';
q = (0.01:0.01:0.14).';
intensity = zeros(numel(energy), numel(q));

for i = 1:numel(q)
    center = 820 + 3600 * q(i);
    lower = center - 90;
    upper = center + 105;
    y = 0.010 + 0.95 * local_gaussian(energy, lower, 45) + ...
        0.82 * local_gaussian(energy, upper, 55);
    if q(i) >= 0.08
        edge_mask = energy < 650 | energy > 1600;
        y(edge_mask) = y(edge_mask) + 0.08 * sin((1:nnz(edge_mask)).' .* 0.9);
    end
    intensity(:, i) = y;
end

qe = struct('energy_meV', energy, 'q_Ainv', q, 'intensity', intensity);
qe_raw = qe;
old_points = table(q, 820 + 3600 * q, ...
    'VariableNames', {'q_Ainv', 'energy_meV'});
end


function y = local_gaussian(x, center, sigma)
y = exp(-0.5 * ((x - center) ./ sigma) .^ 2);
end


function n = local_large_jump_count(points)
if isempty(points) || height(points) < 2
    n = 0;
    return
end
[~, order] = sort(points.q_Ainv);
energy = points.energy_meV(order);
n = sum(abs(diff(energy)) > 250);
end


function opts = local_default_opts(old_points)
opts = struct();
opts.energy_window_meV = [500 1500];
opts.q_range_Ainv = [0 0.06];
opts.q_skip_Ainv = 0;
opts.bin_size = 3;
opts.peak_model = 'gaussian';
opts.pre_subtracted = true;
opts.min_prominence = 0.01;
opts.smooth_width = 1;
opts.old_branch_points = old_points;
opts.fallback_split_meV = 190;
opts.min_peak_separation_meV = 20;
end


function tbl = local_repair_points_table(q, energy, branch_id, branch_label)
n = numel(q);
tbl = table(q(:), abs(q(:)), energy(:), repmat(30, n, 1), ...
    repmat(0.95, n, 1), repmat(1, n, 1), energy(:)-1, energy(:)+1, ...
    repmat(29, n, 1), repmat(31, n, 1), repmat(0.9, n, 1), ...
    repmat(1.1, n, 1), repmat(1, n, 1), repmat(branch_id, n, 1), ...
    repmat({branch_label}, n, 1), repmat(1, n, 1), ...
    repmat({'single_q_direct'}, n, 1), ones(n, 1), ...
    arrayfun(@(x) sprintf('%.12g', x), q(:), 'UniformOutput', false), ...
    arrayfun(@(x) sprintf('%.12g', x), (1:n).', 'UniformOutput', false), ...
    repmat(3, n, 1), repmat({'lorentz'}, n, 1), ...
    repmat({'windowed_branch_tracking'}, n, 1), ...
    repmat({'denoised'}, n, 1), repmat({'sgolay'}, n, 1), ...
    repmat(71, n, 1), repmat(3, n, 1), ...
    repmat({'adaptive_absq'}, n, 1), repmat(0.07, n, 1), ...
    repmat(0.15, n, 1), repmat(0.01, n, 1), ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'energy_meV', ...
    'gamma_meV', 'R2', 'amplitude_fit', 'E_ci_lo', 'E_ci_hi', ...
    'gamma_ci_lo', 'gamma_ci_hi', 'A_ci_lo', 'A_ci_hi', ...
    'raw_height', 'branch', 'branch_label', 'E_ci_half_meV', ...
    'source_mode', 'source_q_count', 'source_q_Ainv', ...
    'source_q_index', 'bin_size_requested', 'peak_model', ...
    'tracking_mode', 'fit_spectrum_source', 'fit_denoise_method', ...
    'fit_denoise_window', 'fit_denoise_order', 'fit_denoise_profile', ...
    'fit_denoise_q_start_Ainv', 'fit_denoise_q_end_Ainv', ...
    'fit_spectrum_delta_rms'});
end
