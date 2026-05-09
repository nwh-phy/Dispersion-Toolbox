classdef LorentzFanoComparisonWorkflowTest < matlab.unittest.TestCase
    methods (Test)
        function scriptDeclaresRequestedComparisonOutputs(testCase)
            caseRoot = fileparts(fileparts(mfilename('fullpath')));
            scriptPath = fullfile(caseRoot, 'scripts', ...
                'run_lorentz_fano_comparison.m');

            testCase.verifyTrue(isfile(scriptPath));
            src = fileread(scriptPath);

            testCase.verifyNotEmpty(strfind(src, ...
                'lorentz_fano_compare_260506'));
            testCase.verifyNotEmpty(strfind(src, ...
                'qRangeOverride_Ainv=[-0.15 0.15]'));
            testCase.verifyNotEmpty(strfind(src, ...
                'peakModelOverride="lorentz"'));
            testCase.verifyNotEmpty(strfind(src, ...
                'peak_extraction_fano_vs_lorentz.png'));
            testCase.verifyNotEmpty(strfind(src, ...
                'gamma_fano_vs_lorentz.png'));
            testCase.verifyNotEmpty(strfind(src, ...
                'lorentz_fano_branch_comparison.csv'));
        end

        function scriptKeepsComparisonFocusedOnPeakAndGamma(testCase)
            caseRoot = fileparts(fileparts(mfilename('fullpath')));
            scriptPath = fullfile(caseRoot, 'scripts', ...
                'run_lorentz_fano_comparison.m');

            testCase.verifyTrue(isfile(scriptPath));
            src = fileread(scriptPath);

            testCase.verifyEmpty(strfind(src, 'I_kin'));
            testCase.verifyNotEmpty(strfind(src, 'local_plot_peak_extraction'));
            testCase.verifyNotEmpty(strfind(src, 'local_plot_gamma_comparison'));
        end
    end
end
