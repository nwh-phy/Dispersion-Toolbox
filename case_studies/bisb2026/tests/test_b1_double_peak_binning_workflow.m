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
    'b1_double_peak_stacked_spectra_positive_q.png', ...
    'b1_double_peak_stacked_spectra_negative_q.png', ...
    'b1_double_peak_stacked_spectra_absq_combined.png', ...
    'b1_double_peak_combined_q_stacked_spectra.png', ...
    'b1_double_peak_waterfall_signed_q.png', ...
    'b1_double_peak_waterfall_absq_combined.png', ...
    'b1_double_peak_waterfall_combined_q.png', ...
    'b1_double_peak_manual_window_seed.csv', ...
    'b1_double_peak_plot_q_binning_map.csv', ...
    'b1_single_peak_vs_double_peak_comparison.png', ...
    'b1_double_peak_lower_physical_fit_binning_260508', ...
    'b1_double_peak_upper_physical_fit_binning_260508'};

for i = 1:numel(required)
    verifyNotEmpty(testCase, strfind(src, required{i}), ...
        sprintf('Missing required workflow reference: %s', required{i}));
end
verifyEmpty(testCase, strfind(src, 'fallback_to_single'));
verifyNotEmpty(testCase, strfind(src, ...
    'options.qRangeOverride_Ainv (1,2) double = [-0.15 0.15]'));
verifyNotEmpty(testCase, strfind(src, 'local_binning_map_q_groups'));
verifyNotEmpty(testCase, strfind(src, 'single_q_direct'));
verifyNotEmpty(testCase, strfind(src, 'combined_q_binning_3'));
verifyNotEmpty(testCase, strfind(src, 'low_q_no_bin_abs_Ainv'));
verifyNotEmpty(testCase, strfind(src, 'source_q_count'));
verifyNotEmpty(testCase, strfind(src, 'source_mode'));
end


function testResultIndexDocumentsStackedManualWindowOutputs(testCase)
index_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'RESULTS_INDEX.md');

verifyTrue(testCase, isfile(index_path));
src = fileread(index_path);

verifyNotEmpty(testCase, strfind(src, ...
    'b1_double_peak_stacked_spectra_positive_q.png'));
verifyNotEmpty(testCase, strfind(src, ...
    'b1_double_peak_combined_q_stacked_spectra.png'));
verifyNotEmpty(testCase, strfind(src, ...
    'b1_double_peak_waterfall_signed_q.png'));
verifyNotEmpty(testCase, strfind(src, ...
    'b1_double_peak_manual_window_seed.csv'));
verifyNotEmpty(testCase, strfind(src, ...
    'b1_double_peak_plot_q_binning_map.csv'));
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
