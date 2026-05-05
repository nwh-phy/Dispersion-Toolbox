classdef QeSelectLowEnergyBranchTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPath(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            run(fullfile(projectRoot, 'startup.m'));
        end
    end

    methods (Test)
        function testRejectsMismatchedPositiveSideAfterOutlier(testCase)
            peaks = makeMismatchedLowBranchPeaks();

            result = qe_select_low_energy_branch(peaks, makeBranchOptions());

            testCase.verifyEqual(result.quality.selected_side, -1);
            testCase.verifyFalse(result.quality.symmetry_trusted);
            testCase.verifyFalse(result.symmetry_average);
            testCase.verifyTrue(all(result.physics_branch(:,1) < 0));
            testCase.verifyGreaterThanOrEqual(size(result.physics_branch, 1), 8);
            testCase.verifyGreaterThan(result.quality.selected_fit_R2, 0.95);
        end

        function testUsesSymmetryAverageWhenSidesMatch(testCase)
            peaks = makeSymmetricLowBranchPeaks();

            result = qe_select_low_energy_branch(peaks, makeBranchOptions());

            testCase.verifyEqual(result.quality.selected_side, 0);
            testCase.verifyTrue(result.quality.symmetry_trusted);
            testCase.verifyTrue(result.symmetry_average);
            testCase.verifyTrue(any(result.physics_branch(:,1) < 0));
            testCase.verifyTrue(any(result.physics_branch(:,1) > 0));
            testCase.verifyGreaterThan(result.quality.combined_fit_R2, 0.98);
        end
    end
end


function opts = makeBranchOptions()
opts = thesis_config().branch;
opts.low.energy_window_meV = [500 2100];
opts.low.max_delta_meV = 350;
opts.low.max_q_gap_Ainv = 0.015;
opts.selection.min_segment_points = 4;
opts.selection.min_fit_points = 3;
opts.selection.q_min_fit_Ainv = 0.01;
opts.selection.min_side_fit_R2 = 0.70;
opts.selection.min_symmetry_overlap = 4;
opts.selection.max_symmetry_abs_delta_meV = 120;
opts.selection.max_symmetry_rel_delta = 0.12;
opts.selection.max_combined_R2_drop = 0.10;
end


function peaks = makeMismatchedLowBranchPeaks()
q_neg = -(0.005:0.005:0.08)';
E_neg = quasi2dEnergy(abs(q_neg), 6.2e7, 24);
E_neg = E_neg .* (1 + 0.01 * sin((1:numel(E_neg))'));

q_pos_outlier = 0.005;
E_pos_outlier = 1380;
q_pos = (0.010:0.005:0.050)';
E_pos = [900; 835; 890; 945; 1085; 1045; 755; 765; 775];

peaks = [
    makeRows(q_neg, E_neg, 0.95)
    makeRows(q_pos_outlier, E_pos_outlier, 0.97)
    makeRows(q_pos, E_pos, 0.93)
    makeRows([0.11; -0.11], [3100; 3150], 0.90)
    ];
end


function peaks = makeSymmetricLowBranchPeaks()
q_abs = (0.01:0.01:0.08)';
E = quasi2dEnergy(q_abs, 5.8e7, 22);
peaks = [
    makeRows(-q_abs, E .* 1.01, 0.96)
    makeRows(q_abs, E .* 0.99, 0.95)
    ];
end


function rows = makeRows(q, energy, r2)
q = q(:);
energy = energy(:);
if isscalar(r2)
    r2 = repmat(r2, size(q));
end
gamma = 0.18 .* energy;
amp = 10 ./ max(abs(q), 0.005).^2;
rawHeight = amp .* 1.2;
rows = nan(numel(q), 12);
rows(:,1) = q;
rows(:,2) = energy;
rows(:,3) = gamma;
rows(:,4) = r2(:);
rows(:,5) = amp;
rows(:,12) = rawHeight;
end


function energy = quasi2dEnergy(q, A, rho0)
energy = sqrt(A .* q ./ (1 + rho0 .* q));
end
