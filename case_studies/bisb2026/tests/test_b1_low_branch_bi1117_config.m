function tests = test_b1_low_branch_bi1117_config
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(case_root, 'scripts'), fullfile(case_root, 'scripts', 'b1_lowq'), ...
    fullfile(case_root, 'scripts', 'b1_binned'));
root = b1_simple_root();
run(fullfile(root, 'startup.m'));
testCase.TestData.root = root;
end


function testOnlyTrainFramesAndVersionDifferFromTheBaseConfig(testCase)
root = testCase.TestData.root;
base = b1_low_branch_matlab_config(root, {'Bi1117_20w'});
cfg = b1_low_branch_bi1117_config(root);
d = bisb_datasets('Bi1117_20w');
verifyEqual(testCase, cfg.preparation.train_frames, 300);
verifyLessThan(testCase, cfg.preparation.train_frames, d.n_frames);
verifyEqual(testCase, cfg.session, base.session);
verifyEqual(testCase, cfg.session.h, d.h);
verifyEqual(testCase, cfg.fit, base.fit);
verifyEqual(testCase, cfg.bootstrap, base.bootstrap);
verifyEqual(testCase, rmfield(cfg.preparation, 'train_frames'), rmfield(base.preparation, 'train_frames'));
end


function testHOverride(testCase)
root = testCase.TestData.root;
cfg = b1_low_branch_bi1117_config(root, 0);
verifyEqual(testCase, cfg.session.h, 0);
verifyEqual(testCase, cfg.session.id, bisb_datasets('Bi1117_20w').short);
verifyError(testCase, @() b1_low_branch_bi1117_config(root, -0.001), 'b1_low_branch_bi1117_config:H');
end
