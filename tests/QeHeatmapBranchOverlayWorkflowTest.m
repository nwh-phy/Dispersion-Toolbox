classdef QeHeatmapBranchOverlayWorkflowTest < matlab.unittest.TestCase
    methods (Test)
        function testAutoFitRefreshesHeatmapBranchOverlay(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            body = extractLocalFunctionBody(src, ...
                'local_on_auto_fit_dispersion', 'local_on_reassign_points');

            testCase.verifyTrue(contains(body, 'local_update_heatmap_branch_overlays'));
        end

        function testHeatmapOverlayUsesStableTag(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            src = fileread(fullfile(projectRoot, 'src', 'interactive_qe_browser.m'));

            body = extractLocalFunctionBody(src, ...
                'local_plot_dispersion_overlay', 'local_plot_dispersion_result');

            testCase.verifyTrue(contains(src, 'function local_update_heatmap_branch_overlays'));
            testCase.verifyTrue(contains(body, 'branch_fit_overlay'));
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
