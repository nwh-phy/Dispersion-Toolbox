function tests = test_bisb_datasets
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(case_root, 'scripts'));
root = b1_simple_root();
run(fullfile(root, 'startup.m'));
testCase.TestData.root = root;
end


function testNamesAreUniqueAndResolve(testCase)
D = bisb_datasets();
names = [{D.id}, {D.short}, [D.aliases]];
verifyEqual(testCase, numel(unique(lower(names))), numel(names), 'Ambiguous dataset names');
for k = 1:numel(D)
    for n = [{D(k).id, D(k).short}, D(k).aliases]
        verifyEqual(testCase, bisb_datasets(n{1}).id, D(k).id);
    end
end
verifyError(testCase, @() bisb_datasets('no_such_dataset'), 'bisb_datasets:Unknown');
end


function testEntriesAreComplete(testCase)
D = bisb_datasets();
for k = 1:numel(D)
    d = D(k);
    verifyGreaterThan(testCase, d.dq, 0, d.id);
    verifyEqual(testCase, d.zero, round(d.zero), d.id);
    verifyTrue(testCase, d.zero >= 1 && d.zero <= 512, d.id);
    verifyMatches(testCase, d.sha256, '^[0-9a-f]{64}$');
    verifyEqual(testCase, d.path, fullfile(testCase.TestData.root, d.set, d.folder, d.file));
    verifyTrue(testCase, endsWith(d.file, '.npy'), d.id);
    verifyTrue(testCase, isnan(d.h) || d.h > 0, d.id);
end
end


function testLowBranchConfigTakesItsSessionsFromTheRegistry(testCase)
cfg = b1_low_branch_matlab_config(testCase.TestData.root);
verifyEqual(testCase, {cfg.session.id}, {'590', 'n0', '20w'});
for s = cfg.session
    d = bisb_datasets(s.id);
    verifyEqual(testCase, [s.dq s.zero s.h], [d.dq d.zero d.h]);
    verifyEqual(testCase, {s.folder, s.file, s.sha256}, {d.folder, d.file, d.sha256});
end
verifyError(testCase, @() b1_low_branch_matlab_config(testCase.TestData.root, {'Bi1107_10w'}), ...
    'b1_low_branch_matlab_config:NoSlot');
end


function testFrozenC20SessionListAgreesWithTheRegistry(testCase)
% legacy_c20/thesis_sessions.m is the session list the frozen C20 scripts read.
case_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(case_root, 'scripts', 'legacy_c20'), fullfile(case_root, 'scripts', 'legacy_202605'));
S = thesis_sessions();
verifyEqual(testCase, sort({S.name}), sort({'590_PL2_10w', 'n0_PL2_10w_repeat', 'no_PL2_20w_2film'}));
for s = S
    d = bisb_datasets(s.name);
    verifyEqual(testCase, s.path, fullfile(testCase.TestData.root, d.set, d.folder), s.name);
    verifyEqual(testCase, s.dq_Ainv, d.dq, s.name);
end
end


function testRawFilesMatchTheirHeaderWhenPresent(testCase)
D = bisb_datasets();
D = D(arrayfun(@(d) isfile(d.path), D));
assumeNotEmpty(testCase, D, 'No raw data on this machine');
for d = D
    fid = fopen(d.path, 'r'); header = fread(fid, 160, '*char')'; fclose(fid);
    shape = regexp(header, '''shape'': \((\d+), (\d+), (\d+)\)', 'tokens', 'once');
    verifyEqual(testCase, str2double(shape), [d.n_frames 512 1028], d.id);
    verifyTrue(testCase, isfile(strrep(d.path, '.npy', '.json')), d.id);
end
end
