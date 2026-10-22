function tests = test_qe_fit_spectrum_joint
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.q = 0.005;
testCase.TestData.dq = 0.0005;
end


function testRecoversDrudeLorentzPeakWithKinematicFactor(testCase)
[E, Y] = synthSpectrum(testCase, 1000, 400, true);
r = qe_fit_spectrum_joint(E, Y, testCase.TestData.q, dq=testCase.TestData.dq, ...
    n_peaks=1, signal_window=[300 1800]);
verifyTrue(testCase, r.success);
verifyLessThan(testCase, abs(r.omega_p - 1000), 25);
verifyLessThan(testCase, abs(r.gamma - 400) / 400, 0.15);
verifyTrue(testCase, r.peak_valid);
end


function testWithoutKinematicFactorE0IsBiasedLow(testCase)
% At q = 0.002 the prefactor falls by ~4x across the peak (qE = 1.54e-3 per eV at 30 kV):
% fitting without it puts E0 ~70 meV low (Gamma 400), with it within a few meV.
q = 0.002;
[E, Y] = synthSpectrum(testCase, 1000, 400, true, q);
with_k = qe_fit_spectrum_joint(E, Y, q, dq=testCase.TestData.dq, n_peaks=1);
without_k = qe_fit_spectrum_joint(E, Y, q, dq=testCase.TestData.dq, n_peaks=1, kinematic=false);
verifyFalse(testCase, without_k.kinematic.used);
verifyLessThan(testCase, abs(with_k.omega_p - 1000), 25);
verifyLessThan(testCase, without_k.omega_p, 1000 - 40);
end


function testTwoPeakResultFollowsTheGuiContract(testCase)
[E, Y1] = synthSpectrum(testCase, 750, 350, true);
[~, Y2] = synthSpectrum(testCase, 1450, 500, true);
Y = Y1 + 0.6 * (Y2 - zlpCounts(E));
r = qe_fit_spectrum_joint(E, Y, testCase.TestData.q, dq=testCase.TestData.dq, n_peaks=2);

verifyEqual(testCase, r.method, 'zlp_joint');
verifyEqual(testCase, r.n_peaks, 2);
verifySize(testCase, r.omega_p, [2 1]);
verifyTrue(testCase, issorted(r.omega_p));
verifySize(testCase, r.apex_energy_meV, [2 1]);
verifySize(testCase, r.peak_valid, [2 1]);
verifyClass(testCase, r.peak_valid, 'logical');
verifySize(testCase, r.peak_curves, [1 2]);
n = numel(r.energy_fit);
verifyEqual(testCase, numel(r.curve_fit), n);
verifyEqual(testCase, numel(r.bg_curve), n);
verifyEqual(testCase, numel(r.peak_curves{1}), n);
verifyEqual(testCase, r.curve_fit, r.bg_curve + r.peak_curves{1} + r.peak_curves{2}, ...
    'RelTol', 1e-9, 'AbsTol', 1e-6 * max(abs(r.curve_fit)));
verifyGreaterThanOrEqual(testCase, min(r.energy_fit), 50);
verifyGreaterThan(testCase, r.R_squared, 0.9);
verifyLessThan(testCase, abs(r.omega_p - [750; 1450]), [40; 80]);
end


function testApexIsTheLossFunctionMaximum(testCase)
[E, Y] = synthSpectrum(testCase, 1000, 600, true);
dl = qe_fit_spectrum_joint(E, Y, testCase.TestData.q, dq=testCase.TestData.dq, n_peaks=1);
[~, k] = max(dl.joint.loss_function(:, 1));
verifyLessThanOrEqual(testCase, abs(dl.apex_energy_meV - E(k)), E(2) - E(1));
verifyLessThan(testCase, dl.apex_energy_meV, dl.omega_p);   % damped DL loss function: apex below E0

sym = qe_fit_spectrum_joint(E, Y, testCase.TestData.q, dq=testCase.TestData.dq, n_peaks=1, ...
    peak_model='lorentz_symmetric');
verifyEqual(testCase, sym.apex_energy_meV, sym.omega_p);
end


function [E, Y] = synthSpectrum(testCase, E0, G, use_kinematic, q)
% ZLP (Pearson VII) + one Drude-Lorentz loss peak convolved with the unit-area ZLP,
% optionally times the 2D kinematic prefactor, with Gaussian-approximated Poisson noise.
if nargin < 5, q = testCase.TestData.q; end
E = (-300:4:2000).';
zlp = zlpCounts(E);
peak = zeros(size(E));
pos = E > 0;
peak(pos) = E(pos) * G ./ ((E(pos) .^ 2 - E0 ^ 2) .^ 2 + E(pos) .^ 2 * G ^ 2);
if use_kinematic
    K = qe_kinematic_prefactor(q, dq=testCase.TestData.dq, form='2d', beam_kV=30);
    peak(pos) = peak(pos) .* K(E(pos));
end
kernel = zlp / sum(zlp);
peak = conv(peak, kernel(E >= -120 & E <= 120), 'same');
peak = peak / max(peak) * 2e3;
mu = zlp + peak + 20;
rng(20261022, 'twister');
Y = mu + sqrt(mu) .* randn(size(mu));
end


function zlp = zlpCounts(E)
zlp = 1e6 * (1 + (E / 18) .^ 2 / 1.5) .^ (-1.5);
end
