function tests = test_b1_peak_evidence_audit_workflow
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testEvidenceAuditScriptDeclaresV15CleanNoFitWorkflow(testCase)
script_path = fullfile(testCase.TestData.project_root, ...
    'case_studies', 'bisb2026', 'scripts', ...
    'run_b1_lorentz_peak_evidence_audit.m');

verifyTrue(testCase, isfile(script_path));
src = fileread(script_path);

required = { ...
    'b1_lorentz_tracking_peak_evidence_audit_260511', ...
    '260510_lorentz_tracking_v15_upper_rapidrise_plateau', ...
    'b1_peak_evidence_audit_points.csv', ...
    'b1_peak_evidence_constraint_release.csv', ...
    'b1_peak_evidence_robustness.csv', ...
    'b1_peak_evidence_candidate_competition.csv', ...
    'b1_peak_evidence_suspicious_points.csv', ...
    'b1_peak_evidence_three_dataset_summary.png', ...
    'energyWindowMeV=[300 1800]', ...
    'peakModel=''lorentz''', ...
    'runFits=false', ...
    'no physical fit was run'};
for i = 1:numel(required)
    verifyNotEmpty(testCase, strfind(src, required{i}), ...
        sprintf('Missing evidence-audit token: %s', required{i}));
end

forbidden = {'DEPRECATED_B1_500MEV_WINDOW', ...
    'branch1Min'', 500', 'peakModelOverride=''fano''', ...
    'run_b1_physical_fit_analysis', '+q/-q'};
for i = 1:numel(forbidden)
    verifyEmpty(testCase, strfind(src, forbidden{i}), ...
        sprintf('Forbidden evidence-audit token found: %s', forbidden{i}));
end
end


function testEvidenceClassifierSeparatesSupportedTrackingAndSuspicious(testCase)
tbl = table();
tbl.q_abs_Ainv = [0.01; 0.04; 0.10];
tbl.branch_label = {'b1_double_peak_upper'; 'b1_double_peak_upper'; ...
    'b1_double_peak_upper'};
tbl.local_support_score = [0.92; 0.45; 0.15];
tbl.component_sse_increase_fraction = [0.18; 0.06; 0.01];
tbl.constraint_release_delta_meV = [35; 75; 230];
tbl.robustness_max_delta_meV = [40; 90; 260];
tbl.gamma_over_E = [0.50; 0.80; 1.70];
tbl.gamma_meV = [420; 760; 1800];

out = b1_peak_evidence_audit_classify(tbl);

verifyEqual(testCase, out.evidence_class{1}, 'data_supported');
verifyEqual(testCase, out.evidence_class{2}, 'tracking_assisted');
verifyEqual(testCase, out.evidence_class{3}, 'suspicious');
verifyTrue(testCase, out.is_data_supported(1));
verifyTrue(testCase, out.is_tracking_assisted(2));
verifyTrue(testCase, out.is_suspicious(3));
verifyNotEmpty(testCase, out.evidence_reason{3});
end


function testEvidenceClassifierFlagsConstraintOnlyPointAsSuspicious(testCase)
tbl = table();
tbl.q_abs_Ainv = 0.08;
tbl.branch_label = {'b1_double_peak_lower'};
tbl.local_support_score = 0.05;
tbl.component_sse_increase_fraction = 0.00;
tbl.constraint_release_delta_meV = 190;
tbl.robustness_max_delta_meV = 210;
tbl.gamma_over_E = 0.4;
tbl.gamma_meV = 300;

out = b1_peak_evidence_audit_classify(tbl);

verifyEqual(testCase, out.evidence_class{1}, 'suspicious');
verifyNotEmpty(testCase, strfind(out.evidence_reason{1}, ...
    'release_unstable'));
end
