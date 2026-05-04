function tests = test_qe_branch_manual_corrections
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end

function testAutoCorrectionReplacesNearestSameQPeak(testCase)
branches = makeBranches();

[updated, correction] = qe_apply_branch_correction(branches, 0.0202, 735, ...
    'Branch', 'auto', ...
    'QTolerance', 0.001);

verifyEqual(testCase, correction.action, 'replace');
verifyEqual(testCase, correction.branch_index, 1);
verifyEqual(testCase, correction.point_index, 2);
verifyEqual(testCase, correction.old_q_Ainv, 0.0200, 'AbsTol', 1e-12);
verifyEqual(testCase, correction.old_energy_meV, 700, 'AbsTol', 1e-12);
verifyEqual(testCase, updated{1}(2, 1), 0.0202, 'AbsTol', 1e-12);
verifyEqual(testCase, updated{1}(2, 2), 735, 'AbsTol', 1e-12);
verifyEqual(testCase, updated{2}, branches{2});
end

function testSpecificBranchAddsWhenNoSameQPointExists(testCase)
branches = makeBranches();

[updated, correction] = qe_apply_branch_correction(branches, 0.0450, 1620, ...
    'Branch', 2, ...
    'QTolerance', 0.001);

verifyEqual(testCase, correction.action, 'add');
verifyEqual(testCase, correction.branch_index, 2);
verifyEqual(testCase, size(updated{2}, 1), size(branches{2}, 1) + 1);
verifyEqual(testCase, updated{2}(end, 1:2), [0.0450 1620], 'AbsTol', 1e-12);
end

function testExplicitBranchTargetOverridesAutoChoice(testCase)
branches = makeBranches();

[updated, correction] = qe_apply_branch_correction(branches, 0.0200, 735, ...
    'Branch', 'Branch 2', ...
    'QTolerance', 0.001);

verifyEqual(testCase, correction.action, 'replace');
verifyEqual(testCase, correction.branch_index, 2);
verifyEqual(testCase, correction.old_energy_meV, 1500, 'AbsTol', 1e-12);
verifyEqual(testCase, updated{2}(2, 2), 735, 'AbsTol', 1e-12);
verifyEqual(testCase, updated{1}, branches{1});
end

function testNewBranchTargetAddsTrailingBranch(testCase)
branches = makeBranches();

[updated, correction] = qe_apply_branch_correction(branches, 0.0250, 990, ...
    'Branch', 'New Branch 3', ...
    'QTolerance', 0.001);

verifyEqual(testCase, correction.action, 'add');
verifyEqual(testCase, correction.branch_index, 3);
verifyEqual(testCase, numel(updated), 3);
verifyEqual(testCase, updated{3}(:, 1:2), [0.0250 990], 'AbsTol', 1e-12);
end

function testUndoReplaceCorrectionRestoresPreviousPeak(testCase)
branches = makeBranches();

[updated, correction] = qe_apply_branch_correction(branches, 0.0202, 735, ...
    'Branch', 'auto', ...
    'QTolerance', 0.001);
restored = qe_revert_branch_correction(updated, correction);

verifyEqual(testCase, restored, branches);
end

function testUndoAddCorrectionRemovesInsertedPeak(testCase)
branches = makeBranches();

[updated, correction] = qe_apply_branch_correction(branches, 0.0450, 1620, ...
    'Branch', 2, ...
    'QTolerance', 0.001);
restored = qe_revert_branch_correction(updated, correction);

verifyEqual(testCase, restored, branches);
end

function testFlattenBranchesUsesAllBranchPoints(testCase)
branches = makeBranches();

pts = qe_flatten_branch_points(branches);

verifySize(testCase, pts, [6 2]);
verifyEqual(testCase, pts(1, :), [0.0100 650], 'AbsTol', 1e-12);
verifyEqual(testCase, pts(end, :), [0.0300 1550], 'AbsTol', 1e-12);
end

function branches = makeBranches()
branches = cell(2, 1);
branches{1} = [
    0.0100 650 100 0.95 1.0;
    0.0200 700 110 0.94 1.1;
    0.0300 760 115 0.93 1.2];
branches{2} = [
    0.0100 1450 300 0.90 0.8;
    0.0200 1500 310 0.91 0.85;
    0.0300 1550 320 0.92 0.9];
end
