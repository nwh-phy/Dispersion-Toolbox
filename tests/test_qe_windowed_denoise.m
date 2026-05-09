function tests = test_qe_windowed_denoise
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testAbsQWindowsOnlyDenoiseSelectedColumns(testCase)
qe = local_synthetic_qe();
input_intensity = qe.intensity;

opts = local_base_opts();
opts.denoise_method = 'Wiener2D';
opts.denoise_windows = [
    local_window('SavGol', [], [0 0.025]), ...
    local_window('SavGol', [], [0.075 0.085])];

qe_out = qe_preprocess(qe, opts);

unselected_cols = [2 6];
selected_cols = [1 3 4 5 7];
verifyEqual(testCase, qe_out.intensity(:, unselected_cols), ...
    input_intensity(:, unselected_cols), 'AbsTol', 0);
verifyGreaterThan(testCase, max(abs(qe_out.intensity(:, selected_cols) - ...
    input_intensity(:, selected_cols)), [], 'all'), 1e-8);
end


function testSignedQWindowCanTargetOneSide(testCase)
qe = local_synthetic_qe();
input_intensity = qe.intensity;

opts = local_base_opts();
opts.denoise_method = 'Wiener2D';
opts.denoise_windows = local_window('SavGol', [0.01 0.09], []);

qe_out = qe_preprocess(qe, opts);

positive_cols = find(qe.q_Ainv >= 0.01 & qe.q_Ainv <= 0.09);
other_cols = setdiff(1:numel(qe.q_Ainv), positive_cols);
verifyEqual(testCase, qe_out.intensity(:, other_cols), ...
    input_intensity(:, other_cols), 'AbsTol', 0);
verifyGreaterThan(testCase, max(abs(qe_out.intensity(:, positive_cols) - ...
    input_intensity(:, positive_cols)), [], 'all'), 1e-8);
end


function testBm3dZeroUsesWeakDefaultFactor(testCase)
tmp_dir = tempname;
mkdir(tmp_dir);
cleanup = onCleanup(@() local_cleanup_path(tmp_dir));
mock_file = fullfile(tmp_dir, 'BM3D_QRS.m');
fid = fopen(mock_file, 'w');
fprintf(fid, 'function denoised = BM3D_QRS(raw, BMfactor, varargin)\n');
fprintf(fid, 'denoised = raw + BMfactor;\n');
fprintf(fid, 'end\n');
fclose(fid);
addpath(tmp_dir, '-begin');
clear BM3D_QRS

qe = local_synthetic_qe();
qe.intensity = zeros(size(qe.intensity));
opts = local_base_opts();
opts.denoise_method = 'BM3D';
opts.denoise_sigma = 0;

qe_out = qe_preprocess(qe, opts);
verifyEqual(testCase, qe_out.intensity, 0.5 * ones(size(qe.intensity)), 'AbsTol', 1e-12);

opts.denoise_sigma = 0.25;
qe_out = qe_preprocess(qe, opts);
verifyEqual(testCase, qe_out.intensity, 0.25 * ones(size(qe.intensity)), 'AbsTol', 1e-12);
end


function testAdaptiveStrengthProfileUsesSameMethodWithSmoothQWeights(testCase)
tmp_dir = tempname;
mkdir(tmp_dir);
cleanup = onCleanup(@() local_cleanup_path(tmp_dir));
mock_file = fullfile(tmp_dir, 'BM3D_QRS.m');
fid = fopen(mock_file, 'w');
fprintf(fid, 'function denoised = BM3D_QRS(raw, BMfactor, varargin)\n');
fprintf(fid, 'denoised = raw + BMfactor;\n');
fprintf(fid, 'end\n');
fclose(fid);
addpath(tmp_dir, '-begin');
clear BM3D_QRS

qe = struct();
qe.energy_meV = (1:5)';
qe.q_Ainv = [0.01 0.02 0.04 0.06 0.08];
qe.intensity = zeros(numel(qe.energy_meV), numel(qe.q_Ainv));

opts = local_base_opts();
opts.denoise_method = 'BM3D';
opts.denoise_q_ramp = struct( ...
    'q_start_Ainv', 0.02, ...
    'q_end_Ainv', 0.06, ...
    'low_sigma', 1, ...
    'high_sigma', 3);

qe_out = qe_preprocess(qe, opts);

verifyEqual(testCase, qe_out.intensity(:, 1), ones(5, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, qe_out.intensity(:, 3), 2 * ones(5, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, qe_out.intensity(:, 5), 3 * ones(5, 1), 'AbsTol', 1e-12);
verifyLessThan(testCase, abs(qe_out.intensity(1, 3) - qe_out.intensity(1, 2)), 1.01);
verifyLessThan(testCase, abs(qe_out.intensity(1, 4) - qe_out.intensity(1, 3)), 1.01);
end


function qe = local_synthetic_qe()
energy = linspace(0, 200, 31)';
q_axis = [-0.08 -0.05 -0.02 0 0.02 0.05 0.08];
base = sin(energy / 19) + 0.02 * energy;
intensity = repmat(base, 1, numel(q_axis));

for qi = 1:numel(q_axis)
    intensity(8 + mod(qi, 5), qi) = intensity(8 + mod(qi, 5), qi) + 3 + qi / 10;
    intensity(19 - mod(qi, 4), qi) = intensity(19 - mod(qi, 4), qi) - 2 - qi / 20;
end

qe = struct();
qe.energy_meV = energy;
qe.q_Ainv = q_axis;
qe.intensity = intensity;
end


function opts = local_base_opts()
opts = struct();
opts.do_normalize = false;
opts.do_denoise = true;
opts.denoise_method = 'Wiener2D';
opts.denoise_sigma = 0;
opts.sg_order = 2;
opts.sg_framelen = 5;
opts.do_bg_sub = false;
opts.do_deconv = false;
end


function window = local_window(method, q_range, abs_q_range)
window = struct();
window.method = method;
window.q_range_Ainv = q_range;
window.abs_q_range_Ainv = abs_q_range;
window.denoise_sigma = 0;
window.sg_order = 2;
window.sg_framelen = 5;
end


function local_cleanup_path(tmp_dir)
if exist(tmp_dir, 'dir')
    rmpath(tmp_dir);
    clear BM3D_QRS
    try
        rmdir(tmp_dir, 's');
    catch
    end
end
end
