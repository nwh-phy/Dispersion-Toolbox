classdef QeSelectedBranchFitWorkflowTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPath(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            run(fullfile(projectRoot, 'startup.m'));
        end
    end

    methods (Test)
        function testBranchSelectionParserSupportsAllAndSingleBranch(testCase)
            allIdx = qe_selected_branch_indices(3, 'All');
            oneIdx = qe_selected_branch_indices(3, 'Branch 2');

            testCase.verifyEqual(allIdx, 1:3);
            testCase.verifyEqual(oneIdx, 2);
        end

        function testBranchSelectionParserRejectsOutOfRangeBranch(testCase)
            idx = qe_selected_branch_indices(2, 'Branch 3');

            testCase.verifyEmpty(idx);
        end

        function testFitUiExposesDedicatedBranchSelector(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'ui', 'qe_browser_ui.m'));

            testCase.verifyTrue(contains(src, 'FitBranchDropdown'));
            testCase.verifyTrue(contains(src, 'fit_branch_dropdown'));
            testCase.verifyTrue(contains(src, 'Fit Branch'));
        end

        function testFitCuratedUsesSelectedBranchIndices(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            body = extractLocalFunctionBody(src, ...
                'local_on_fit_dispersion', 'local_on_export_dispersion');

            testCase.verifyTrue(contains(body, 'qe_selected_branch_indices'));
            testCase.verifyTrue(contains(body, 'ui.FitBranchDropdown.Value'));
            testCase.verifyFalse(contains(body, 'for b = 1:n_branches'));
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
