function tests = test_b1_forward_paths
tests = functiontests(localfunctions);
end

function testConfiguredPathsTakePrecedence(t)
p = b1_forward_config('Bi1117_20w',tempname);
p.input_dir = fullfile(p.run_root,'elsewhere','custom_input');
p.response_dir = fullfile(p.run_root,'alternate_response');
cfg = b1_p2a_operator_config(p);
verifyEqual(t,b1_forward_paths(cfg),p);
verifyEqual(t,b1_forward_paths(cfg,p.run_root),p);
verifyError(t,@() b1_forward_paths(cfg,tempname),'b1_forward_paths:Root');
end

function testLegacyRootOverride(t)
p = b1_forward_paths(struct(),tempname);
verifyEqual(t,p.stage1_file,fullfile(p.run_root,'stage1',['channels_' p.session '.mat']));
verifyEqual(t,p.legacy_fit_prefix,'fit_590_');
verifyEqual(t,p.fit_prefix,'fit_590_PL2_10w_');
end
