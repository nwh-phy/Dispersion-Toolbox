classdef QeAssignPeakBranchesByWindowsTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPath(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            run(fullfile(projectRoot, 'startup.m'));
        end
    end

    methods (Test)
        function testAssignsThreeEditableEnergyWindows(testCase)
            peaks = [
                makePeakRows([-0.02; 0.02], [900; 930], 0.90, 10)
                makePeakRows([-0.03; 0.03], [1900; 1940], 0.91, 12)
                makePeakRows([-0.04; 0.04], [3200; 3220], 0.92, 14)
                ];
            specs = makeBranchSpecs();

            result = qe_assign_peak_branches_by_windows(peaks, specs);

            testCase.verifyEqual(numel(result.branches), 3);
            testCase.verifyEqual(size(result.branches{1}, 1), 2);
            testCase.verifyEqual(size(result.branches{2}, 1), 2);
            testCase.verifyEqual(size(result.branches{3}, 1), 2);
            testCase.verifyTrue(all(result.branches{1}(:,2) >= 500 & result.branches{1}(:,2) <= 1200));
            testCase.verifyTrue(all(result.branches{2}(:,2) >= 1700 & result.branches{2}(:,2) <= 2400));
            testCase.verifyTrue(all(result.branches{3}(:,2) >= 3000 & result.branches{3}(:,2) <= 3600));
        end

        function testOverlappingWindowsAssignEachPeakOnlyOnce(testCase)
            peaks = [
                makePeakRows(0.02, 1720, 0.95, 10)
                makePeakRows(0.03, 2050, 0.95, 10)
                ];
            specs = [
                struct('name', 'Low', 'energy_window_meV', [500 2100], 'enabled', true)
                struct('name', 'Mid', 'energy_window_meV', [1700 2600], 'enabled', true)
                ];

            result = qe_assign_peak_branches_by_windows(peaks, specs);

            testCase.verifyEqual(size(result.branches{1}, 1), 1);
            testCase.verifyEqual(size(result.branches{2}, 1), 1);
            testCase.verifyEqual(size(qe_flatten_branch_points(result.branches), 1), 2);
            testCase.verifyEqual(result.branches{1}(1,2), 1720);
            testCase.verifyEqual(result.branches{2}(1,2), 2050);
        end

        function testKeepsBestPeakPerSignedQWithinBranch(testCase)
            weaker = makePeakRows(0.02, 900, 0.40, 2);
            stronger = makePeakRows(0.02, 940, 0.95, 10);
            peaks = [weaker; stronger];
            specs = struct('name', 'Low', 'energy_window_meV', [500 1200], 'enabled', true);

            result = qe_assign_peak_branches_by_windows(peaks, specs);

            testCase.verifyEqual(size(result.branches{1}, 1), 1);
            testCase.verifyEqual(result.branches{1}(1,2), 940);
            testCase.verifyEqual(height(result.rejected), 1);
            testCase.verifyEqual(result.rejected.reason{1}, 'duplicate_signed_q');
        end
    end
end


function specs = makeBranchSpecs()
specs = [
    struct('name', 'Low', 'energy_window_meV', [500 1200], 'enabled', true)
    struct('name', 'Mid', 'energy_window_meV', [1700 2400], 'enabled', true)
    struct('name', 'High', 'energy_window_meV', [3000 3600], 'enabled', true)
    ];
end


function rows = makePeakRows(q, energy, r2, rawHeight)
q = q(:);
energy = energy(:);
n = numel(q);
if isscalar(energy)
    energy = repmat(energy, n, 1);
end
if isscalar(r2)
    r2 = repmat(r2, n, 1);
end
if isscalar(rawHeight)
    rawHeight = repmat(rawHeight, n, 1);
end

rows = nan(n, 12);
rows(:,1) = q;
rows(:,2) = energy;
rows(:,3) = 0.15 .* energy;
rows(:,4) = r2(:);
rows(:,5) = rawHeight(:);
rows(:,12) = rawHeight(:);
end
