function tests = test_thesis_pipeline_infra
%TEST_THESIS_PIPELINE_INFRA Tests for the thesis-specific reproducible pipeline helpers.
tests = functiontests(localfunctions);
end


function setupOnce(~)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testThesisConfigDefinesConservativeDefaults(testCase)
cfg = thesis_config();

verifyEqual(testCase, cfg.pipeline.name, 'thesis_qeels_pipeline');
verifyEqual(testCase, cfg.preprocess.norm_method, 'ZLP Peak');
verifyTrue(testCase, cfg.preprocess.do_normalize);
verifyTrue(testCase, cfg.preprocess.do_denoise);
verifyTrue(testCase, cfg.preprocess.do_bg_sub);
verifyFalse(testCase, cfg.preprocess.do_deconv);
verifyEqual(testCase, cfg.fit.pre_subtracted, true);
verifyEqual(testCase, cfg.fit.max_peaks, 2);
verifyEqual(testCase, cfg.branch.low.energy_window_meV, [500 2100]);
verifyEqual(testCase, cfg.branch.high.energy_window_meV, [2800 3800]);
verifyEqual(testCase, cfg.models.low, 'quasi2d_plasmon');
verifyEqual(testCase, cfg.models.high, 'optical_constant');
end


function testAssignQeBranchesUsesEnergyWindowsAndContinuity(testCase)
peaks = makePeakMatrix();
opts = thesis_config().branch;
opts.low.max_delta_meV = 250;
opts.high.max_delta_meV = 250;

branches = assign_qe_branches(peaks, opts);

verifyEqual(testCase, size(branches.low.accepted, 1), 4);
verifyEqual(testCase, size(branches.high.accepted, 1), 4);
verifyTrue(testCase, all(branches.low.accepted(:,2) >= 500 & branches.low.accepted(:,2) <= 2100));
verifyTrue(testCase, all(branches.high.accepted(:,2) >= 2800 & branches.high.accepted(:,2) <= 3800));
verifyTrue(testCase, any(strcmp(branches.rejected.reason, 'outside_energy_window')));
verifyTrue(testCase, any(strcmp(branches.rejected.reason, 'continuity_jump')));
verifyTrue(testCase, all(diff(branches.low.accepted(:,1)) >= 0));
end


function testAssignQeBranchesKeepsBestCandidatePerSignedQ(testCase)
peaks = makePeakMatrix();
% Add a weaker low-energy candidate at the same q. It is inside the low window
% but should lose against the higher-confidence candidate at q = 0.02.
weak = peaks(4,:);
weak(2) = 1180;
weak(3) = 900;
weak(4) = 0.35;
weak(5) = 5;
weak(12) = 4;
peaks = [peaks; weak];

branches = assign_qe_branches(peaks, thesis_config().branch);
q002 = abs(branches.low.accepted(:,1) - 0.02) < 1e-12;

verifyEqual(testCase, sum(q002), 1);
verifyEqual(testCase, branches.low.accepted(q002, 2), 1010, 'AbsTol', 1e-12);
verifyTrue(testCase, any(strcmp(branches.rejected.reason, 'duplicate_signed_q')));
end


function testSymmetrizeQeBranchesWeightedAveragesSignedPairs(testCase)
branches = assign_qe_branches(makePeakMatrix(), thesis_config().branch);
sym = symmetrize_qe_branches(branches, thesis_config().symmetry);

verifyEqual(testCase, size(sym.low.symmetrized, 1), 2);
verifyEqual(testCase, size(sym.high.symmetrized, 1), 2);
verifyEqual(testCase, sym.low.symmetrized(:,1), [0.01; 0.02], 'AbsTol', 1e-12);
verifyGreaterThan(testCase, sym.low.symmetrized(2,2), 990);
verifyLessThan(testCase, sym.low.symmetrized(2,2), 1010);
verifyTrue(testCase, all(sym.low.symmetrized(:,4) > 0));
end


function testThesisSensitivityConfigsExposeSmallAuditMatrix(testCase)
cfg = thesis_config();
variants = thesis_sensitivity_configs(cfg);
variant_names = {variants.name};

verifyGreaterThanOrEqual(testCase, numel(variants), 4);
verifyTrue(testCase, any(strcmp(variant_names, 'baseline')));
verifyTrue(testCase, any(strcmp(variant_names, 'deconv_iter5')));
verifyTrue(testCase, any(strcmp(variant_names, 'bg_auto')));
verifyTrue(testCase, any(strcmp(variant_names, 'prominence_strict')));

baseline = variants(strcmp(variant_names, 'baseline')).cfg;
deconv = variants(strcmp(variant_names, 'deconv_iter5')).cfg;
verifyFalse(testCase, baseline.preprocess.do_deconv);
verifyTrue(testCase, deconv.preprocess.do_deconv);
verifyEqual(testCase, deconv.preprocess.deconv_iter, 5);
end


function testBuildThesisManifestRecordsConfigSessionsAndResults(testCase)
cfg = thesis_config();
sessions = thesis_sessions(cfg);
result = struct();
result.session_name = 'synthetic';
result.success = true;
result.branch_summary = struct('low_n', 2, 'high_n', 2);
manifest = build_thesis_manifest(cfg, sessions(1), result, struct('git_commit', 'abc123'));

verifyEqual(testCase, manifest.pipeline.name, cfg.pipeline.name);
verifyEqual(testCase, manifest.git.commit, 'abc123');
verifyEqual(testCase, manifest.sessions(1).name, sessions(1).name);
verifyEqual(testCase, manifest.results(1).session_name, 'synthetic');
verifyTrue(testCase, isfield(manifest, 'generated_at'));
verifyTrue(testCase, isfield(manifest, 'caution'));
end


function testSummarizeThesisResultsIncludesModelParameters(testCase)
result = struct();
result.session_name = 'synthetic';
result.success = true;
result.error = '';
result.n_raw_peaks = 5;
result.n_fit_success = 4;
result.branch_summary = struct('low_n', 2, 'high_n', 3, 'rejected_n', 1);
result.models.low.fit = struct('rho0', 21, 'E_flat_meV', 1800, ...
    'R_squared', 0.93, 'RMSE_meV', 40);
result.models.high.fit = struct('params', 3400, 'R_squared', 0.01, ...
    'RMSE_meV', 100);

summary = summarize_thesis_results(result, 'unit_variant');

verifyEqual(testCase, height(summary), 1);
verifyEqual(testCase, summary.variant(1), "unit_variant");
verifyEqual(testCase, summary.session(1), "synthetic");
verifyEqual(testCase, summary.low_n(1), 2);
verifyEqual(testCase, summary.high_n(1), 3);
verifyEqual(testCase, summary.low_rho0_A(1), 21);
verifyEqual(testCase, summary.low_Eflat_meV(1), 1800);
verifyEqual(testCase, summary.high_constant_meV(1), 3400);
end


function peaks = makePeakMatrix()
% Columns follow qe_auto_fit blind output convention:
% [q, E, Gamma, R2, A, E_ci_lo, E_ci_hi, G_ci_lo, G_ci_hi, A_ci_lo, A_ci_hi, raw_h]
peaks = nan(11, 12);
peaks(:,1) = [-0.02; -0.01; 0.01; 0.02; -0.02; -0.01; 0.01; 0.02; 0.03; 0.04; -0.04];
peaks(:,2) = [980; 820; 830; 1010; 3200; 3100; 3120; 3210; 2500; 1800; 900];
peaks(:,3) = [120; 110; 115; 125; 350; 340; 360; 355; 300; 120; 100];
peaks(:,4) = [0.90; 0.88; 0.87; 0.91; 0.85; 0.86; 0.84; 0.87; 0.80; 0.95; 0.90];
peaks(:,5) = [100; 90; 92; 110; 80; 75; 78; 82; 10; 95; 90];
peaks(:,12) = peaks(:,5);
end
