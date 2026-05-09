function tests = test_20w_candidate_overlay_export_workflow
tests = functiontests(localfunctions);
end


function testOverlayScriptDeclaresEvidenceInputsAndOutputs(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
script_path = fullfile(case_root, 'scripts', ...
    'run_20w_candidate_overlay_export.m');

assertTrue(testCase, isfile(script_path), ...
    'case_studies/bisb2026/scripts/run_20w_candidate_overlay_export.m should exist.');

txt = fileread(script_path);

verifyNotEmpty(testCase, strfind(txt, ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined'));
verifyNotEmpty(testCase, strfind(txt, '20w_B1_highq_audit_260506'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_highq_rescue_candidates.csv'));
verifyNotEmpty(testCase, strfind(txt, 'analysis_results.mat'));
verifyNotEmpty(testCase, strfind(txt, '20w_B1_candidate_overlay_260506'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_candidate_points.csv'));
verifyNotEmpty(testCase, strfind(txt, ...
    '20w_B1_conservative_qe_map.png'));
verifyNotEmpty(testCase, strfind(txt, ...
    '20w_B1_exploratory_qe_map_candidate_overlay.png'));
verifyNotEmpty(testCase, strfind(txt, ...
    '20w_B1_conservative_dispersion.png'));
verifyNotEmpty(testCase, strfind(txt, ...
    '20w_B1_exploratory_dispersion_candidate_overlay.png'));
verifyNotEmpty(testCase, strfind(txt, ...
    '20w_B1_candidate_overlay_report.md'));
end


function testOverlayScriptSeparatesMainFitFromCandidateOnlyLayer(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
script_path = fullfile(case_root, 'scripts', ...
    'run_20w_candidate_overlay_export.m');

assertTrue(testCase, isfile(script_path), ...
    'case_studies/bisb2026/scripts/run_20w_candidate_overlay_export.m should exist.');

txt = fileread(script_path);

verifyNotEmpty(testCase, strfind(txt, 'included_in_main_fit'));
verifyNotEmpty(testCase, strfind(txt, 'candidate_only_not_for_fit'));
verifyNotEmpty(testCase, strfind(txt, 'local_plot_qe_map_overlay'));
verifyNotEmpty(testCase, strfind(txt, 'local_plot_dispersion_overlay'));
verifyNotEmpty(testCase, strfind(txt, 'show_candidates'));
verifyNotEmpty(testCase, strfind(txt, ...
    'candidate-only points are not included in quasi-2D fitting'));
end
