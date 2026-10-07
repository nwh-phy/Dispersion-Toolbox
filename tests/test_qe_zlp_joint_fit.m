function tests = test_qe_zlp_joint_fit
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testRecoversPlasmonOnZlpTail(testCase)
[E, Y] = makeSpectrum(struct('E0', 700, 'W', 400, 'height', 300));
fit = qe_zlp_joint_fit(E, Y, signal_window=[300 1800], n_peaks=1, n_starts=4);

verifyTrue(testCase, fit.success);
verifyEqual(testCase, fit.numerical_status, 'converged');
verifyLessThan(testCase, abs(fit.parameters(1) - 700), 15);
verifyLessThan(testCase, abs(fit.parameters(2) - 400) / 400, 0.08);
verifyLessThan(testCase, abs(fit.zlp_parameters.center_meV - 1.5), 0.5);
verifyLessThan(testCase, fit.chi2_red_signal, 1.5);
end


function testKernelRemovesInstrumentBroadening(testCase)
zlp = struct('h', [5e5 2e4], 's', [20 40; 20 40], 'm', [3 1; 3 1], 'c', 0);
[E, Y] = makeSpectrum(struct('E0', 400, 'W', 30, 'height', 2000), zlp);
deconv = qe_zlp_joint_fit(E, Y, signal_window=[250 600], n_peaks=1, n_starts=4);
plain = qe_zlp_joint_fit(E, Y, signal_window=[250 600], n_peaks=1, n_starts=4, kernel_mode='none');

verifyLessThan(testCase, abs(deconv.parameters(2) - 30), 6);
verifyGreaterThan(testCase, plain.parameters(2), 45);
end


function testGainSidePrefitExtrapolatesSymmetricLossTail(testCase)
% Stage A sees only E <= 2*core_halfwidth; the loss tail at 300-1000 meV is
% pure extrapolation through the gain-side shape.
zlp = struct('h', [5e5 2e4], 's', [8 25; 8 25], 'm', [3 1; 3 1], 'c', 1.5);
[E, Y, truth] = makeSpectrum(struct('E0', 700, 'W', 400, 'height', 0), zlp);
fit = qe_zlp_joint_fit(E, Y, n_peaks=0);

at = E >= 300 & E <= 1000;
rel = abs(fit.prefit_zlp(at) - truth.zlp(at)) ./ truth.zlp(at);
verifyLessThan(testCase, max(rel), 0.15);
end


function testGainSideOfPeakFollowsDetailedBalance(testCase)
[E, Y] = makeSpectrum(struct('E0', 80, 'W', 20, 'height', 3000));
fit = qe_zlp_joint_fit(E, Y, signal_window=[40 160], n_peaks=1, n_starts=2, temperature_K=300);

[~, ip] = min(abs(E - fit.parameters(1)));
[~, in] = min(abs(E + fit.parameters(1)));
ratio = fit.peaks_intrinsic(in) / fit.peaks_intrinsic(ip);
verifyEqual(testCase, ratio, exp(-E(ip) / (0.08617333262 * 300)), 'RelTol', 1e-6);
end


function testAuxPeakSeparatesLowEnergyLossFromPlasmon(testCase)
peaks = struct('E0', {700, 90}, 'W', {400, 60}, 'height', {300, 400}, 'shape', {'sym', 'dl'});
[E, Y] = makeSpectrum(peaks);
fit = qe_zlp_joint_fit(E, Y, signal_window=[300 1800], n_peaks=1, n_starts=4, aux_windows=[30 300]);

verifyEqual(testCase, fit.numerical_status, 'converged');
verifyLessThan(testCase, abs(fit.parameters(1) - 700), 15);
verifyLessThan(testCase, abs(fit.parameters(2) - 400) / 400, 0.10);
verifyLessThan(testCase, abs(fit.aux_parameters(1) - 90), 10);
verifySize(testCase, fit.aux_peaks, [numel(fit.energy_meV) 1]);
end


function testTwoAuxPeaksResolvePhononAndHigherBand(testCase)
peaks = struct('E0', {700, 50, 120}, 'W', {400, 10, 50}, 'height', {300, 2000, 300}, ...
    'shape', {'sym', 'dl', 'dl'});
[E, Y] = makeSpectrum(peaks);
fit = qe_zlp_joint_fit(E, Y, signal_window=[300 1800], n_peaks=1, n_starts=4, ...
    aux_windows=[30 80; 80 300]);

verifyEqual(testCase, size(fit.aux_parameters), [2 3]);
verifyLessThan(testCase, abs(fit.aux_parameters(1, 1) - 50), 6);
verifyLessThan(testCase, abs(fit.aux_parameters(2, 1) - 120), 15);
verifyLessThan(testCase, abs(fit.parameters(1) - 700), 15);
end


function testKinematicPrefactorMatchesAnalyticLimit(testCase)
[K3, info] = qe_kinematic_prefactor(0.01, dq=1e-9, form='3d');
K2 = qe_kinematic_prefactor(0.01, dq=1e-9, form='2d');
qE = info.qE_per_meV * 1000;

verifyEqual(testCase, info.qE_per_meV, 1.5434e-6, 'RelTol', 1e-3);
verifyEqual(testCase, K3(1000), 1 / (0.01 ^ 2 + qE ^ 2), 'RelTol', 1e-6);
verifyEqual(testCase, K2(1000), 0.01 / (0.01 ^ 2 + qE ^ 2) ^ 2, 'RelTol', 1e-6);
end


function testPrefactorRemovesKinematicTilt(testCase)
% A 2D prefactor at |q| = 0.0025 A^-1 (30 kV) drops by ~55% between 0.5 and
% 1.3 eV; fitting without it pulls a broad peak to lower energy.
Kfun = qe_kinematic_prefactor([0.002 0.0025 0.003], form='2d', sigma_probe=4e-4);
[E, Y] = makeSpectrum(struct('E0', 900, 'W', 800, 'height', 300), [], Kfun);
with = qe_zlp_joint_fit(E, Y, signal_window=[300 1800], n_peaks=1, n_starts=4, prefactor=Kfun);
without = qe_zlp_joint_fit(E, Y, signal_window=[300 1800], n_peaks=1, n_starts=4);

verifyLessThan(testCase, abs(with.parameters(1) - 900), 20);
verifyLessThan(testCase, abs(with.parameters(2) - 800) / 800, 0.08);
verifyLessThan(testCase, without.parameters(1), 870);
verifySize(testCase, with.loss_function, [numel(with.energy_meV) 1]);
end


function testRejectsNonUniformGrid(testCase)
[E, Y] = makeSpectrum(struct('E0', 700, 'W', 400, 'height', 300));
E(300:end) = E(300:end) + 1;
verifyError(testCase, @() qe_zlp_joint_fit(E, Y), 'qe_zlp_joint_fit:NonUniform');
end


function [E, Y, truth] = makeSpectrum(peaks, zlp, prefactor)
% Two-component asymmetric Pearson ZLP plus loss peaks (symmetric Lorentzian
% by default, 'dl' = Drude-Lorentz) with detailed-balance gain side, each
% convolved with the area-normalised ZLP; Gaussian-approximated Poisson noise.
if nargin < 3, prefactor = []; end
if nargin < 2 || isempty(zlp)
    zlp = struct('h', [5e5 2e4], 's', [8 25; 8 27], 'm', [3 1.0; 3 0.95], 'c', 1.5);
end
E = (-180:4:1800).';
pear = @(x, sL, sR, mL, mR) (x < 0) .* (1 + (x / sL) .^ 2 / mL) .^ (-mL) + ...
    (x >= 0) .* (1 + (x / sR) .^ 2 / mR) .^ (-mR);
Z = zeros(size(E));
lag = (-300:4:300).';
K = zeros(size(lag));
for k = 1:numel(zlp.h)
    Z = Z + zlp.h(k) * pear(E - zlp.c, zlp.s(1, k), zlp.s(2, k), zlp.m(1, k), zlp.m(2, k));
    K = K + zlp.h(k) * pear(lag - zlp.c, zlp.s(1, k), zlp.s(2, k), zlp.m(1, k), zlp.m(2, k));
end
K = K / sum(K);
Ex = (E(1) - 300:4:E(end) + 300).';
Ea = abs(Ex);
kT = 0.08617333262 * 300;
S = zeros(size(E));
for j = 1:numel(peaks)
    pk = peaks(j);
    if isfield(pk, 'shape') && strcmp(pk.shape, 'dl')
        f = Ea * pk.W ./ ((Ea .^ 2 - pk.E0 ^ 2) .^ 2 + Ea .^ 2 * pk.W ^ 2);
    else
        f = (pk.W / 2) ./ ((Ea - pk.E0) .^ 2 + (pk.W / 2) ^ 2);
    end
    f(Ex < 0) = f(Ex < 0) .* exp(-Ea(Ex < 0) / kT);
    if ~isempty(prefactor), f = f .* prefactor(Ea); end
    Sj = conv(f, K, 'valid');
    S = S + pk.height * Sj / max(Sj);
end
mu = Z + S;
rng(20261007);
Y = mu + sqrt(mu) .* randn(size(mu));
truth = struct('zlp', Z, 'signal', S, 'mu', mu);
end
