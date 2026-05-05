classdef QeBrowserCuratedWorkflowTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPath(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            run(fullfile(projectRoot, 'startup.m'));
        end
    end

    methods (Test)
        function testAutoFitBuildsEditableBranchesWithoutDispersionFit(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            body = extractLocalFunctionBody(src, ...
                'local_on_auto_fit_dispersion', 'local_on_reassign_points');

            testCase.verifyTrue(contains(body, 'qe_assign_peak_branches_by_windows'));
            testCase.verifyFalse(contains(body, 'fit_dispersion_generic'));
            testCase.verifyFalse(contains(body, 'plot_fit_curve'));
            testCase.verifyTrue(contains(body, 'state.fitResults = {}'));
        end

        function testAutoFitPassesQualityFiltersToBranchAssignment(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            body = extractLocalFunctionBody(src, ...
                'local_on_auto_fit_dispersion', 'local_on_reassign_points');

            testCase.verifyTrue(contains(body, 'branch_filter_opts.q_skip_Ainv'));
            testCase.verifyTrue(contains(body, 'branch_filter_opts.min_R2'));
            testCase.verifyTrue(contains(body, 'branch_filter_opts.max_gamma_ratio'));
            testCase.verifyTrue(contains(body, ...
                'qe_assign_peak_branches_by_windows(peaks, branch_specs, branch_filter_opts)'));
        end

        function testSingleSpectrumFitDisplaysPeakQualityMetadata(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            body = extractLocalFunctionBody(src, ...
                'local_on_fit_spectrum', 'local_on_accept_fit');

            testCase.verifyTrue(contains(body, 'apex_energy_meV'));
            testCase.verifyTrue(contains(body, 'gamma_ratio'));
            testCase.verifyTrue(contains(body, 'peak_valid'));
            testCase.verifyTrue(contains(body, 'qF=%.2f'));
        end

        function testFitCuratedRefreshesSingleSpectrumBranchOverlay(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            singleBody = extractLocalFunctionBody(src, ...
                'local_plot_single_spectrum', 'local_plot_waterfall');
            fitBody = extractLocalFunctionBody(src, ...
                'local_on_fit_dispersion', 'local_on_export_dispersion');

            testCase.verifyTrue(contains(singleBody, 'local_plot_single_branch_fit_overlay'));
            testCase.verifyTrue(contains(src, 'function local_plot_single_branch_fit_overlay'));
            testCase.verifyTrue(contains(src, 'local_eval_branch_fit_at_q'));
            testCase.verifyTrue(contains(fitBody, 'local_update_selected_q_views'));
        end

        function testMainUiExposesCuratedWorkflowAndBranchWindows(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'ui', 'qe_browser_ui.m'));

            testCase.verifyTrue(contains(src, 'Auto Fit Branches'));
            testCase.verifyTrue(contains(src, 'Fit Curated'));
            testCase.verifyTrue(contains(src, 'Branch1MinField'));
            testCase.verifyTrue(contains(src, 'Branch2MinField'));
            testCase.verifyTrue(contains(src, 'Branch3MinField'));
            testCase.verifyTrue(contains(src, '''Fano'''));
            testCase.verifyTrue(contains(src, '''fano'''));
            testCase.verifyFalse(contains(src, 'Auto Fit 蠅(q)'));
            testCase.verifyFalse(contains(src, 'Fit Dispersion'));
        end
    end
end


function body = extractLocalFunctionBody(src, startName, nextName)
startPattern = sprintf('function %s', startName);
nextPattern = sprintf('function %s', nextName);
startIndex = strfind(src, startPattern);
nextIndex = strfind(src, nextPattern);

assert(~isempty(startIndex), 'Missing local function %s', startName);
assert(~isempty(nextIndex), 'Missing local function %s', nextName);

startIndex = startIndex(1);
nextIndex = nextIndex(find(nextIndex > startIndex, 1));
assert(~isempty(nextIndex), 'Missing following local function %s', nextName);

body = src(startIndex:nextIndex-1);
end
