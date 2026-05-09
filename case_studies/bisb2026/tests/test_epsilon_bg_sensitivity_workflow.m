function tests = test_epsilon_bg_sensitivity_workflow
tests = functiontests(localfunctions);
end


function testScriptDeclaresPhysicalEpsilonBackgroundSweep(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
script_path = fullfile(case_root, 'scripts', 'run_epsilon_bg_sensitivity.m');

assertTrue(testCase, isfile(script_path), ...
    'case_studies/bisb2026/scripts/run_epsilon_bg_sensitivity.m should exist.');

txt = fileread(script_path);

verifyNotEmpty(testCase, regexp(txt, ...
    'epsilon_bg_values\s*=\s*\[1,\s*4\.5,\s*10,\s*15\]', 'once'), ...
    'The sensitivity sweep should include the baseline and MoS2-relevant bracket.');
verifyNotEmpty(testCase, strfind(txt, 'local_fit_quasi2d_epsilon_bg'), ...
    'The workflow should use a direct epsilon_bg-aware quasi-2D fit.');
verifyNotEmpty(testCase, strfind(txt, 'q_c_Ainv = epsilon_bg / rho0_fit'), ...
    'q_c must be reported as epsilon_bg/rho0 for each background value.');
verifyNotEmpty(testCase, strfind(txt, 'rho0_max_A = 5000'), ...
    'rho0 upper bound should be loose enough for high epsilon_bg B3 fits.');
verifyNotEmpty(testCase, strfind(txt, 'epsilon_bg_sensitivity_summary.csv'));
verifyNotEmpty(testCase, strfind(txt, 'epsilon_bg_sensitivity_report.md'));
end


function testScriptUsesCurrentThreeAreaReports(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
script_path = fullfile(case_root, 'scripts', 'run_epsilon_bg_sensitivity.m');

assertTrue(testCase, isfile(script_path), ...
    'case_studies/bisb2026/scripts/run_epsilon_bg_sensitivity.m should exist.');

txt = fileread(script_path);

verifyNotEmpty(testCase, strfind(txt, '590_gui_history_area_260506'));
verifyNotEmpty(testCase, strfind(txt, 'n0_PL2_10w_gui_history_area_260506'));
verifyNotEmpty(testCase, strfind(txt, 'no_PL2_20w_2film_gui_history_area_260506_highq_refined'));
verifyNotEmpty(testCase, strfind(txt, 'branch1_points.csv'));
verifyNotEmpty(testCase, strfind(txt, 'branch3_points.csv'));
end
