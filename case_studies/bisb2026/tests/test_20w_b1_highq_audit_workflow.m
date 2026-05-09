function tests = test_20w_b1_highq_audit_workflow
tests = functiontests(localfunctions);
end


function testAuditScriptDeclaresEvidenceInputsAndOutputs(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
script_path = fullfile(case_root, 'scripts', 'run_20w_b1_highq_audit.m');

assertTrue(testCase, isfile(script_path), ...
    'case_studies/bisb2026/scripts/run_20w_b1_highq_audit.m should exist.');

txt = fileread(script_path);

verifyNotEmpty(testCase, strfind(txt, ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined'));
verifyNotEmpty(testCase, strfind(txt, 'branch_refinement_log.csv'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_points.csv'));
verifyNotEmpty(testCase, strfind(txt, 'dispersion_model_summary.csv'));
verifyNotEmpty(testCase, strfind(txt, 'analysis_results.mat'));
verifyNotEmpty(testCase, strfind(txt, '20w_B1_highq_audit_260506'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_highq_audit.csv'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_highq_symmetry_audit.csv'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_highq_rescue_candidates.csv'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_model_comparison.csv'));
verifyNotEmpty(testCase, strfind(txt, '20w_B1_highq_audit_report.md'));
end


function testAuditScriptIncludesRescueAndSymmetryLogic(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
script_path = fullfile(case_root, 'scripts', 'run_20w_b1_highq_audit.m');

assertTrue(testCase, isfile(script_path), ...
    'case_studies/bisb2026/scripts/run_20w_b1_highq_audit.m should exist.');

txt = fileread(script_path);

verifyNotEmpty(testCase, strfind(txt, 'local_build_symmetry_audit'));
verifyNotEmpty(testCase, strfind(txt, 'local_classify_rescue_candidates'));
verifyNotEmpty(testCase, strfind(txt, 'model_residual_abs_meV'));
verifyNotEmpty(testCase, strfind(txt, 'symmetry_abs_delta_meV'));
verifyNotEmpty(testCase, strfind(txt, 'candidate_only_not_for_fit'));
verifyNotEmpty(testCase, strfind(txt, 'local_plot_highq_diagnostic_grid'));
verifyNotEmpty(testCase, strfind(txt, 'local_fit_branch1_model'));
end
