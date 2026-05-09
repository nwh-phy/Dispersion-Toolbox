classdef WideqLinewidthAnalysisWorkflowTest < matlab.unittest.TestCase
    methods (Test)
        function scriptDeclaresLinewidthOnlyOutputs(testCase)
            caseRoot = fileparts(fileparts(mfilename('fullpath')));
            scriptPath = fullfile(caseRoot, 'scripts', ...
                'run_wideq_linewidth_analysis.m');

            testCase.verifyTrue(isfile(scriptPath));
            src = fileread(scriptPath);

            testCase.verifyNotEmpty(strfind(src, ...
                'wideq_linewidth_260506'));
            testCase.verifyNotEmpty(strfind(src, ...
                'qRangeOverride_Ainv=[-0.30 0.30]'));
            testCase.verifyNotEmpty(strfind(src, ...
                'gamma_summary.csv'));
            testCase.verifyNotEmpty(strfind(src, ...
                'do_style_linewidth_summary.csv'));
            testCase.verifyNotEmpty(strfind(src, ...
                'linewidth_uncertainty_budget.csv'));
            testCase.verifyNotEmpty(strfind(src, ...
                'b1_linewidth_q.png'));
            testCase.verifyNotEmpty(strfind(src, ...
                'b3_linewidth_q.png'));
            testCase.verifyNotEmpty(strfind(src, ...
                'b1_b3_fwhm_lorentz_comparison.png'));
            testCase.verifyNotEmpty(strfind(src, ...
                'b1_b3_quality_factor.png'));
            testCase.verifyNotEmpty(strfind(src, ...
                'wideq_boundary_map.png'));
            testCase.verifyNotEmpty(strfind(src, ...
                'wideq_linewidth_addendum.md'));
        end

        function scriptDoesNotIntroduceStrengthAnalysis(testCase)
            caseRoot = fileparts(fileparts(mfilename('fullpath')));
            scriptPath = fullfile(caseRoot, 'scripts', ...
                'run_wideq_linewidth_analysis.m');

            testCase.verifyTrue(isfile(scriptPath));
            src = fileread(scriptPath);

            testCase.verifyEmpty(strfind(src, 'I_kin'));
            testCase.verifyEmpty(strfind(src, 'local_area'));
            testCase.verifyEmpty(strfind(src, 'raw_height'));
            testCase.verifyNotEmpty(strfind(src, 'gamma_over_E'));
            testCase.verifyNotEmpty(strfind(src, 'quality_factor'));
            testCase.verifyNotEmpty(strfind(src, 'fano_fwhm_meV'));
            testCase.verifyNotEmpty(strfind(src, 'lorentz_gamma_meV'));
        end
    end
end
