classdef QeBrowserCuratedWorkflowTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPath(~)
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

            testCase.verifyTrue(contains(body, ...
                'auto_opts.bootstrap_ci_samples = ui.BootstrapCiSamplesField.Value'));
            testCase.verifyTrue(contains(body, 'auto_opts.branch_specs = branch_specs'));
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
            testCase.verifyTrue(contains(body, 'local_fit_peak_energy'));
            testCase.verifyTrue(contains(body, ...
                '''bootstrap_ci_samples'', ui.BootstrapCiSamplesField.Value'));
        end

        function testAcceptedSingleSpectrumFitsUseApexEnergy(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            body = extractLocalFunctionBody(src, ...
                'local_on_accept_fit', 'local_on_auto_fit_dispersion');

            testCase.verifyTrue(contains(src, 'function peak_energy = local_fit_peak_energy'));
            testCase.verifyTrue(contains(body, 'peak_energy = local_fit_peak_energy(result)'));
            testCase.verifyTrue(contains(body, 'peak_energy(p)'));
            testCase.verifyFalse(contains(body, 'state.manual_branches, q_value, result.omega_p(p)'));
            testCase.verifyFalse(contains(body, 'new_pt = [q_value, result.omega_p(p)]'));
        end

        function testBranchPlotFallbackProvidesEveryPointErrorbar(testCase)
            branch_without_ci = [0.01, 700; 0.02, 760; 0.03, 820];
            err = qe_plot_helpers.energy_error_half_width(branch_without_ci);

            testCase.verifySize(err, [3, 1]);
            testCase.verifyTrue(all(isfinite(err)));
            testCase.verifyTrue(all(err > 0));

            branch_with_ci = [0.01, 700, 50, 0.9, 1, 690, 715; ...
                              0.02, 760, 50, 0.9, 1, NaN, NaN; ...
                              0.03, 820, 50, 0.9, 1, 830, 840];
            err = qe_plot_helpers.energy_error_half_width(branch_with_ci);

            testCase.verifyEqual(err(1), 12.5, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isfinite(err)));
            testCase.verifyTrue(all(err > 0));
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
            testCase.verifyTrue(contains(src, 'CI Samples:'));
            testCase.verifyTrue(contains(src, 'BootstrapCiSamplesField'));
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
