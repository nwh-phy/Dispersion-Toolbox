function tests = test_b1_summary_tables
%TEST_B1_SUMMARY_TABLES  b1_layer_compare and b1_p2b_hscan_summary on small made-up runs.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(case_root, 'scripts'));
b1_simple_root();
testCase.TestData.run_root = tempname;
mkdir(fullfile(testCase.TestData.run_root, 'stage5'));
end


function teardownOnce(testCase)
rmdir(testCase.TestData.run_root, 's');
end


function testLayerCompareAveragesConvergedBinsPerBandAndSide(testCase)
run_root = testCase.TestData.run_root;
ids = {'590_PL2_10w', 'n0_PL2_10w_repeat', 'no_PL2_20w_2film'};
q = [-0.0025; 0.0025; 0.003; -0.004; 0.004; 0.02];        % band 1: +-0.0025, 0.003; band 2: +-0.004
status = {'converged'; 'converged'; 'boundary'; 'converged'; 'converged'; 'converged'};
for s = 1:3
    v = 100 * s + [10; 20; 999; 40; 60; 999];             % the boundary bin and q = 0.02 must not count
    T = table(q, status, v, v + 1, v + 2, v + 3, v + 4, 'VariableNames', ...
        {'q', 'status', 'apex_low', 'apex_high', 'apex_measured_LF', 'E0_low', 'E0_high'});
    writetable(T, fullfile(run_root, 'stage5', sprintf('e0_stage2_2d_h5_%s.csv', ids{s})));
end

[~, T] = evalc('b1_layer_compare(run_root)');

verifyEqual(testCase, height(T), 5 * 5 * 3);
row = @(m, lo, side) T(strcmp(T.measure, m) & abs(T.q_lo - lo) < 1e-12 & strcmp(T.side, side), :);
r = row('apex_low', 0.0021, 'both');
verifyEqual(testCase, [r.mean_590 r.sem_590 r.n_590], [115 5 2], 'AbsTol', 1e-9);
verifyEqual(testCase, [r.n0_minus_590 r.w20_minus_590 r.w20_minus_n0], [100 200 100], 'AbsTol', 1e-9);
r = row('E0_high', 0.0031, 'both');
verifyEqual(testCase, [r.mean_20w r.sem_20w r.n_20w], [354 10 2], 'AbsTol', 1e-9);
r = row('apex_low', 0.0021, '+q');
verifyEqual(testCase, [r.mean_590 r.n_590], [120 1]);
verifyTrue(testCase, isnan(r.sem_590));
r = row('apex_low', 0.0021, '-q');
verifyEqual(testCase, [r.mean_590 r.n_590], [110 1]);
r = row('apex_low', 0.0051, 'both');
verifyEqual(testCase, r.n_590, 0);
verifyTrue(testCase, isnan(r.mean_590));
end


function testHscanSummaryFindsRunsByConfigHash(testCase)
run_root = testCase.TestData.run_root; h = [0.001 0.004]; chi2 = [50 40; 70 90];
for k = 1:2
    tag = b1_sha256(b1_p2b_slot_config('hscan', h(k)));
    d = fullfile(run_root, 'p2b', ['fit_590_' tag(1:8)]); mkdir(d);
    writetable(table({'main'; 'y: Gaussian'}, chi2(k, :).', [1.1; 1.2], [9; 8], [0.5; -0.5], [3; 4], 'VariableNames', ...
        {'variant', 'chi2', 'chi2dof', 'holdout_chi2_per_point', 'holdout_mean_residual_sigma', 'negative_tails'}), ...
        fullfile(d, 'p2b_fits.csv'));
    writetable(table(6, 1.3 + k, 'VariableNames', {'left_out_block', 'chi2n_predicted'}), fullfile(d, 'p2b_lobo.csv'));
    Q = [0.001; 0.002; 0.004]; z = 1000 * k + (1:3).';
    writetable(table(Q, z, z / 100, z + 1, z / 100, z + 2, z / 100, z + 3, z / 100, z + 4, z / 100, 'VariableNames', ...
        {'Q', 'apex_b1', 'apex_b1_err', 'E0_b1', 'E0_b1_err', 'apex_b2', 'apex_b2_err', 'E0_b2', 'E0_b2_err', 'G_b2', 'G_b2_err'}), ...
        fullfile(d, 'p2b_curves_main.csv'));
end

[~, F, C] = evalc('b1_p2b_hscan_summary(run_root, h, [0.001 0.0038])');

verifyEqual(testCase, F.h.', [0.001 0.001 0.004 0.004]);
verifyEqual(testCase, F.dchi2.', [10 0 30 50]);
verifyEqual(testCase, F.lobo_chi2n_predicted([1 3]).', [2.3 3.3], 'AbsTol', 1e-12);
verifyTrue(testCase, all(isnan(F.lobo_chi2n_predicted([2 4]))));
c = C(abs(C.h - 0.004) < 1e-12 & strcmp(C.measure, 'E0_b1'), :);
verifyEqual(testCase, c.Q_curve.', [0.001 0.004]);
verifyEqual(testCase, c.value_meV.', [2002 2004]);
verifyEqual(testCase, c.stat_err_meV.', [20.01 20.03], 'AbsTol', 1e-12);
verifyError(testCase, @() b1_p2b_hscan_summary(run_root, 0.002), 'b1_p2b_hscan_summary:NoRun');
end
