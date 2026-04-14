function tests = test_qe_background_subtraction
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testPowerBackgroundSubtractionRecoversSyntheticPeakHeight(testCase)
[qe, peak_only] = makeSyntheticQeDataset([1.0]);
opts = makeBgOpts('Power');

[qe_out, bg_diag] = qe_preprocess(qe, opts);

[peak_true, idx_true] = max(peak_only(:, 1));
peak_measured = qe_out.intensity(idx_true, 1);

verifyGreaterThan(testCase, peak_measured, 0);
verifyLessThan(testCase, abs(peak_measured - peak_true) / peak_true, 0.20);
verifyLessThan(testCase, bg_diag(1).neg_fraction, 0.10);
verifyLessThan(testCase, bg_diag(1).neg_area_fraction, 0.05);
end


function testAutoBackgroundSelectionPrefersPowerForPowerLawTail(testCase)
[qe, ~] = makeSyntheticQeDataset([1.0]);
opts = makeBgOpts('Auto');
opts.bg_candidate_methods = {'Power', 'Exp2'};

[~, bg_diag] = qe_preprocess(qe, opts);

verifyEqual(testCase, bg_diag(1).selected_method, 'Power');
verifyGreaterThanOrEqual(testCase, numel(bg_diag(1).candidate_methods), 2);
verifyGreaterThanOrEqual(testCase, numel(bg_diag(1).candidate_scores), 2);
end


function testBackgroundDiagnosticsExposeFractionsAndSelectedMethod(testCase)
[qe, ~] = makeSyntheticQeDataset([1.0, 0.85]);
opts = makeBgOpts('Auto');
opts.bg_candidate_methods = {'Power', 'Exp2'};

[~, bg_diag] = qe_preprocess(qe, opts);

verifyTrue(testCase, isfield(bg_diag, 'selected_method'));
verifyTrue(testCase, isfield(bg_diag, 'candidate_methods'));
verifyTrue(testCase, isfield(bg_diag, 'candidate_scores'));
verifyTrue(testCase, isfield(bg_diag, 'neg_fraction'));
verifyTrue(testCase, isfield(bg_diag, 'neg_area_fraction'));
verifyTrue(testCase, isfield(bg_diag, 'bg_fraction'));

for qi = 1:numel(bg_diag)
    verifyGreaterThanOrEqual(testCase, bg_diag(qi).neg_fraction, 0);
    verifyLessThanOrEqual(testCase, bg_diag(qi).neg_fraction, 1);
    verifyGreaterThanOrEqual(testCase, bg_diag(qi).neg_area_fraction, 0);
    verifyLessThanOrEqual(testCase, bg_diag(qi).neg_area_fraction, 1);
    verifyGreaterThanOrEqual(testCase, bg_diag(qi).bg_fraction, 0);
    verifyLessThanOrEqual(testCase, bg_diag(qi).bg_fraction, 1);
end
end


function opts = makeBgOpts(method)
opts = struct();
opts.do_despike = false;
opts.do_normalize = false;
opts.do_denoise = false;
opts.do_bg_sub = true;
opts.bg_method = method;
opts.bg_win_lo = [50, 300];
opts.bg_win_hi = [];
opts.bg_iterative = false;
opts.do_deconv = false;
end


function [qe, peak_only] = makeSyntheticQeDataset(scales)
energy_meV = (-100:10:4000)';
q_channel = 1:numel(scales);
intensity = zeros(numel(energy_meV), numel(scales));
peak_only = zeros(numel(energy_meV), numel(scales));

for qi = 1:numel(scales)
    scale = scales(qi);
    background = syntheticPowerBackground(energy_meV, scale);
    peak = syntheticLorentzPeak(energy_meV, 1200, 180, 0.24 * scale);
    shoulder = syntheticLorentzPeak(energy_meV, 1850, 260, 0.05 * scale);
    intensity(:, qi) = background + peak + shoulder;
    peak_only(:, qi) = peak + shoulder;
end

qe = make_qe_struct(intensity, energy_meV, 0.005, q_channel, ...
    source_kind="synthetic", label="synthetic-bg-test");
end


function bg = syntheticPowerBackground(energy_meV, scale)
energy_pos = max(energy_meV, 1) + 25;
bg = 18 * scale * energy_pos .^ (-1.32);
bg(energy_meV < 0) = max(bg(energy_meV < 0), 0.02 * scale);
end


function peak = syntheticLorentzPeak(energy_meV, center, gamma, amplitude)
x = energy_meV - center;
peak = amplitude ./ (1 + (x ./ gamma) .^ 2);
end
