function tests = test_load_qe_dataset_cache
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testRawCacheInvalidatesWhenQCropChanges(testCase)
work_dir = tempname;
mkdir(work_dir);
cleanup = onCleanup(@() local_remove_dir(work_dir));

raw_path = fullfile(work_dir, 'synthetic_10w_raw.mat');
local_write_synthetic_raw_mat(raw_path);

first = load_qe_dataset(raw_path, 0.005, q_crop=[20 70]);
verifyEqual(testCase, size(first.qe.intensity, 2), 51);

cache_path = fullfile(work_dir, 'eq3D_processed.mat');
verifyEqual(testCase, exist(cache_path, 'file'), 2);

second = load_qe_dataset(raw_path, 0.005, q_crop=[30 50]);
verifyEqual(testCase, size(second.qe.intensity, 2), 21);
verifyEqual(testCase, second.qe.q_channel, 1:21);
end


function local_write_synthetic_raw_mat(raw_path)
n_frames = 3;
n_E = 512;
n_q = 128;
energy = ((1:n_E)' - 1) * 4;
raw = zeros(n_frames, n_E, n_q);

for frame_idx = 1:n_frames
    for q_idx = 1:n_q
        zlp_center = 80 + round(0.03 * (q_idx - 64));
        raw(frame_idx, :, q_idx) = ...
            exp(-((energy - zlp_center).^2) ./ (2 * 7^2)) + 0.01;
    end
end

save(raw_path, 'raw', 'energy', '-v7.3');
end


function local_remove_dir(path_name)
if isfolder(path_name)
    rmdir(path_name, 's');
end
end
