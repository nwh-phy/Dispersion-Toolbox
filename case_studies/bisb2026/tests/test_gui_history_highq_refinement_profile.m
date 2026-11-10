function tests = test_gui_history_highq_refinement_profile
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.script_src = fileread(fullfile(case_root, ...
    'scripts', 'legacy_202605', 'run_590_gui_history_area_analysis.m'));
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
verifyHasText(testCase, src, 'profile.q_min_Ainv = 0.010;');
verifyHasText(testCase, src, 'profile.refit_window_meV = [1000 1700];');
verifyHasText(testCase, src, 'profile.max_energy_shift_meV = 250;');
verifyHasText(testCase, src, 'profile.window_edge_margin_meV = 50;');
verifyHasText(testCase, src, ...
    'function [branches, refinement_log] = local_apply_session_refinement');
end

function testDispersionPlotsFitLineOnlyForB1(testCase)
src = testCase.TestData.script_src;

verifyHasText(testCase, src, 'if b == 1');
verifyHasText(testCase, src, 'local_plot_supported_fit_curve(ax, best.fit, col, b);');
end

function testSupportsLegacyB2PointOverride(testCase)
src = testCase.TestData.script_src;

verifyHasText(testCase, src, 'local_apply_branch_point_overrides');
verifyHasText(testCase, src, 'legacy_branch_points_260521');
verifyHasText(testCase, src, 'branch_point_override_log.csv');
end

function verifyHasText(testCase, text, pattern)
verifyTrue(testCase, contains(string(text), string(pattern)), ...
    sprintf('Expected source to contain: %s', pattern));
end
