classdef QeViewModeNamingTest < matlab.unittest.TestCase
    methods (Test)
        function testViewModeDropdownUsesComparisonName(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'ui', 'qe_browser_ui.m'));

            testCase.verifyTrue(contains(src, '"Items", {''Physical'', ''Comparison''}'));
            testCase.verifyFalse(contains(src, '"Items", {''Physical'', ''Normalized''}'));
        end

        function testInteractiveBrowserAcceptsComparisonAndLegacyNormalized(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            testCase.verifyTrue(contains(src, '"comparison"'));
            testCase.verifyTrue(contains(src, '"normalized"'));
            testCase.verifyTrue(contains(src, 'local_normalize_view_mode_value'));
        end

        function testUserFacingComparisonTextAvoidsNormalizedViewName(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));
            comparisonSrc = fileread(fullfile(projectRoot, 'src', 'build_comparison_qe.m'));

            testCase.verifyTrue(contains(src, 'Comparison view'));
            testCase.verifyTrue(contains(comparisonSrc, 'comparison off-axis component'));
            testCase.verifyFalse(contains(src, 'Normalized view data'));
            testCase.verifyFalse(contains(comparisonSrc, 'normalized off-axis component'));
        end
    end
end
