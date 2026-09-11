function tests = test_qe_dq_calibration
tests = functiontests(localfunctions);
end


function setupOnce(~)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testTenWPathInfersCorrectMomentumStep(testCase)
work_dir = tempname;
mkdir(work_dir);
cleanup = onCleanup(@() local_remove_dir(work_dir));

session_dir = fullfile(work_dir, '590 PL2 10w 0.004 10sx300');
mkdir(session_dir);
local_write_eq3d(fullfile(session_dir, 'eq3D.mat'));

dataset = load_qe_dataset(session_dir);

verifyEqual(testCase, dataset.dq_Ainv, 0.0005, 'AbsTol', 1e-12);
verifyEqual(testCase, median(diff(dataset.qe.q_Ainv)), 0.0005, 'AbsTol', 1e-12);
end


function testTwentyWPathInfersCorrectMomentumStep(testCase)
work_dir = tempname;
mkdir(work_dir);
cleanup = onCleanup(@() local_remove_dir(work_dir));

session_dir = fullfile(work_dir, 'no pl2 20w 0.004 10sx300 2film');
mkdir(session_dir);
local_write_eq3d(fullfile(session_dir, 'eq3D.mat'));

dataset = load_qe_dataset(session_dir);

verifyEqual(testCase, dataset.dq_Ainv, 0.00025, 'AbsTol', 1e-12);
verifyEqual(testCase, median(diff(dataset.qe.q_Ainv)), 0.00025, 'AbsTol', 1e-12);
end


function testExplicitDqOverrideStillWins(testCase)
work_dir = tempname;
mkdir(work_dir);
cleanup = onCleanup(@() local_remove_dir(work_dir));

session_dir = fullfile(work_dir, '590 PL2 10w 0.004 10sx300');
mkdir(session_dir);
local_write_eq3d(fullfile(session_dir, 'eq3D.mat'));

dataset = load_qe_dataset(session_dir, 0.00125);

verifyEqual(testCase, dataset.dq_Ainv, 0.00125, 'AbsTol', 1e-12);
verifyEqual(testCase, median(diff(dataset.qe.q_Ainv)), 0.00125, 'AbsTol', 1e-12);
end


function local_write_eq3d(path_name)
n_E = 64;
n_q = 21;
e = ((1:n_E) - 1) * 4;
a3 = zeros(n_E, n_q);
for q_idx = 1:n_q
    center = 20 + 0.1 * (q_idx - 11);
    a3(:, q_idx) = exp(-(((1:n_E)' - center).^2) ./ (2 * 4^2));
end
save(path_name, 'a3', 'e');
end


function local_remove_dir(path_name)
if isfolder(path_name)
    rmdir(path_name, 's');
end
end
