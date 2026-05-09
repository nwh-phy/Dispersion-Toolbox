function tests = test_gui_history_highq_refinement_profile
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.script_src = fileread(fullfile(case_root, ...
    'scripts', 'run_590_gui_history_area_analysis.m'));
end

function testTwentyWHighQRefinementIsSessionScoped(testCase)
src = testCase.TestData.script_src;

verifyHasText(testCase, src, ...
    'sessions(1).refinement_profile = local_empty_refinement_profile();');
verifyHasText(testCase, src, ...
    'sessions(2).refinement_profile = local_empty_refinement_profile();');
verifyHasText(testCase, src, ...
    'sessions(3).refinement_profile = local_20w_highq_refinement_profile();');
verifyHasText(testCase, src, ...
    'function profile = local_20w_highq_refinement_profile()');
end

function testTwentyWProfileTargetsOnlyB1LargeQ(testCase)
src = testCase.TestData.script_src;

verifyHasText(testCase, src, 'profile.branch_index = 1;');
verifyHasText(testCase, src, 'profile.q_min_Ainv = 0.10;');
verifyHasText(testCase, src, 'profile.refit_window_meV = [1000 1700];');
verifyHasText(testCase, src, 'profile.max_energy_shift_meV = 250;');
verifyHasText(testCase, src, 'profile.window_edge_margin_meV = 50;');
verifyHasText(testCase, src, ...
    'function [branches, refinement_log] = local_apply_session_refinement');
end

function verifyHasText(testCase, text, pattern)
verifyTrue(testCase, contains(string(text), string(pattern)), ...
    sprintf('Expected source to contain: %s', pattern));
end
