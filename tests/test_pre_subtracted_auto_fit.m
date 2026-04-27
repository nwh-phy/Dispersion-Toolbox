function tests = test_pre_subtracted_auto_fit
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testFitLossFunctionPreSubtractedKeepsBackgroundZero(testCase)
[energy_meV, spectrum, center_meV, baseline] = makeSingleSpectrum();

result = fit_loss_function(energy_meV, spectrum, ...
    'E_min', 80, 'E_max', 2500, ...
    'max_peaks', 1, ...
    'min_prominence', 0.05, ...
    'smooth_width', 9, ...
    'peak_model', 'gaussian', ...
    'pre_subtracted', true);

verifyTrue(testCase, result.pre_subtracted);
verifyEqual(testCase, result.B0, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, result.alpha, 0, 'AbsTol', 1e-12);
verifyLessThan(testCase, abs(result.offset - baseline), 0.02);
verifyEqual(testCase, result.n_peaks, 1);
verifyLessThan(testCase, abs(result.omega_p(1) - center_meV), 40);
end


function testFitLossFunctionRawModeStillUsesPowerLawBackground(testCase)
[energy_meV, spectrum] = makeRawPowerLawSpectrum();

result = fit_loss_function(energy_meV, spectrum, ...
    'E_min', 80, 'E_max', 2500, ...
    'max_peaks', 1, ...
    'min_prominence', 0.05, ...
    'smooth_width', 9, ...
    'peak_model', 'gaussian');

verifyFalse(testCase, result.pre_subtracted);
verifyGreaterThan(testCase, result.B0, 0);
verifyGreaterThan(testCase, result.alpha, 0);
verifyEqual(testCase, result.offset, 0, 'AbsTol', 1e-12);
end


function testQeAutoFitPropagatesPreSubtractedBlindMode(testCase)
[qe, qe_raw, opts] = makeAutoFitCase();
opts.guesses = [];

results = qe_auto_fit(qe, qe_raw, opts);

verifyFalse(testCase, results.used_seed);
verifyPreSubtractedFitDetails(testCase, results.fit_details, 2);
end


function testPropagateSeedPeaksPropagatesPreSubtractedSeedMode(testCase)
[qe, ~, opts, centers_meV] = makeAutoFitCase();
seed_idx = 2;

results = propagate_seed_peaks(qe.intensity, qe.energy_meV, qe.q_Ainv, ...
    'seed_guesses', centers_meV(seed_idx), ...
    'seed_idx', seed_idx, ...
    'direction', 'both', ...
    'max_shift', 120, ...
    'min_R2', 0, ...
    'peak_model', 'gaussian', ...
    'pre_subtracted', true, ...
    'E_min', opts.E_min, 'E_max', opts.E_max, ...
    'smooth_width', opts.smooth_width, ...
    'verbose', false);

verifyPreSubtractedFitDetails(testCase, results.fit_details, 2);
end


function verifyPreSubtractedFitDetails(testCase, fit_details, min_count)
populated = fit_details(~cellfun(@isempty, fit_details));
verifyGreaterThanOrEqual(testCase, numel(populated), min_count);

for i = 1:numel(populated)
    detail = populated{i};
    verifyTrue(testCase, detail.pre_subtracted);
    verifyEqual(testCase, detail.B0, 0, 'AbsTol', 1e-12);
    verifyEqual(testCase, detail.alpha, 0, 'AbsTol', 1e-12);
    verifyTrue(testCase, isfinite(detail.offset));
end
end


function [qe, qe_raw, opts, centers_meV] = makeAutoFitCase()
energy_meV = linspace(80, 2500, 400)';
q_axis = [0.010; 0.020; 0.030; 0.040];
centers_meV = [650; 700; 750; 800];
baseline = 0.020;
sigma = 55;

intensity = zeros(numel(energy_meV), numel(q_axis));
for qi = 1:numel(q_axis)
    amp = 0.23 + 0.015 * qi;
    ripple = 0.0015 * sin(energy_meV / 180 + qi);
    intensity(:, qi) = baseline + ripple + gaussianPeak(energy_meV, centers_meV(qi), sigma, amp);
end

qe = struct();
qe.intensity = intensity;
qe.energy_meV = energy_meV;
qe.q_Ainv = q_axis;
qe_raw = qe;

opts = struct();
opts.E_min = 80;
opts.E_max = 2500;
opts.q_start = min(q_axis);
opts.q_end = max(q_axis);
opts.prominence = 0.05;
opts.smooth_width = 9;
opts.max_peaks = 1;
opts.peak_model = 'gaussian';
opts.guesses = [];
opts.seed_idx = 2;
opts.max_shift = 120;
opts.energy_mask = true(size(energy_meV));
opts.energy_axis = energy_meV;
opts.R2_threshold = 0;
opts.verbose = false;
opts.progress_fn = [];
opts.pre_subtracted = true;
end


function [energy_meV, spectrum, center_meV, baseline] = makeSingleSpectrum()
energy_meV = linspace(80, 2500, 400)';
center_meV = 720;
baseline = 0.020;
spectrum = baseline ...
    + 0.0015 * sin(energy_meV / 170) ...
    + gaussianPeak(energy_meV, center_meV, 55, 0.24);
end


function [energy_meV, spectrum] = makeRawPowerLawSpectrum()
energy_meV = linspace(80, 2500, 400)';
center_meV = 720;
alpha = 0.75;
b0 = 0.020 * 80 ^ alpha;
background = b0 .* energy_meV .^ (-alpha);
spectrum = background + gaussianPeak(energy_meV, center_meV, 55, 0.24);
end


function y = gaussianPeak(x, center, sigma, amplitude)
y = amplitude .* exp(-0.5 .* ((x - center) ./ sigma) .^ 2);
end
