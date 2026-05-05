function tests = test_qe_preprocess_preset
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testBiSbStablePresetEncodesLowQBackgroundContract(testCase)
opts = qe_preprocess_preset('bisb_lowq_stable_v1');

verifyTrue(testCase, opts.do_normalize);
verifyEqual(testCase, opts.norm_method, 'ZLP Peak');
verifyTrue(testCase, opts.do_denoise);
verifyEqual(testCase, opts.denoise_method, 'Wiener2D');
verifyTrue(testCase, opts.do_bg_sub);
verifyEqual(testCase, opts.bg_method, 'Auto');
verifyEqual(testCase, opts.bg_win_lo, [50 300]);
verifyEqual(testCase, opts.bg_win_hi, [2800 3400]);
verifyEqual(testCase, opts.bg_auto_group_qmax, 0.01, 'AbsTol', 1e-12);
verifyTrue(testCase, any(strcmp(opts.bg_candidate_methods, 'PearsonVII')));
verifyFalse(testCase, opts.do_deconv);
end


function testPresetAllowsExplicitOverrides(testCase)
opts = qe_preprocess_preset('bisb_lowq_stable_v1', ...
    struct('do_denoise', false, 'bg_auto_group_qmax', 0.02));

verifyFalse(testCase, opts.do_denoise);
verifyEqual(testCase, opts.bg_auto_group_qmax, 0.02, 'AbsTol', 1e-12);
verifyEqual(testCase, opts.bg_method, 'Auto');
end
