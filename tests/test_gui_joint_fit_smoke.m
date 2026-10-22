function tests = test_gui_joint_fit_smoke
% Runs the real GUI on a synthetic 10w dataset (ZLP + one dispersing Drude-Lorentz
% plasmon with the 2D kinematic prefactor) and drives the single-spectrum fit
% through the test hooks of interactive_qe_browser.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
data_dir = fullfile(folder, 'synth 10w smoke');   % '10w' in the path gives dq = 0.0005 1/A
mkdir(data_dir);
[a3, e] = synthDataset();
save(fullfile(data_dir, 'eq3D.mat'), 'a3', 'e');
testCase.TestData.data_dir = data_dir;
end


function setup(testCase)
[app, hooks] = interactive_qe_browser(testCase.TestData.data_dir);
testCase.addTeardown(@delete, app);
ui = hooks.ui;
ui.EnergyMinField.Value = 300;
ui.EnergyMaxField.Value = 1800;
ui.MaxPeaksField.Value = 1;
ui.PeakModelDropdown.Value = 'lorentz';
ui.GuessField.Value = '';
hooks.select_q(channelAt(hooks.state(), 0.005));
testCase.TestData.hooks = hooks;
end


function testZlpJointFitFindsThePlasmonAndAcceptsIt(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();

st = h.state();
verifyNotEmpty(testCase, st.pendingFit);
r = st.pendingFit.result;
verifyEqual(testCase, r.method, 'zlp_joint');
verifyTrue(testCase, r.kinematic.used);
verifyEqual(testCase, st.pendingFit.q_Ainv, 0.005, 'AbsTol', 1e-9);
verifyLessThan(testCase, abs(r.omega_p - truthE0(0.005)), 30);
verifySubstring(testCase, char(h.ui.FitInfoLabel.Text), 'ZLP joint');
verifySubstring(testCase, char(h.ui.FitInfoLabel.Text), 'chi2_red');
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'on');
names = get(findobj(h.ui.SingleAxes, 'Tag', 'spectrum_fit_overlay', 'Type', 'line'), 'DisplayName');
verifyTrue(testCase, any(strcmp(names, 'Fit input (raw counts)')));
peak_line = findobj(h.ui.SingleAxes, 'Tag', 'spectrum_fit_overlay', 'DisplayName', 'Peak 1');
verifyNumElements(testCase, peak_line, 1);
xl = xlim(h.ui.SingleAxes); yl = ylim(h.ui.SingleAxes);
in_view = peak_line.XData >= xl(1) & peak_line.XData <= xl(2);
verifyGreaterThanOrEqual(testCase, min(peak_line.YData(in_view)), yl(1));   % Auto Y keeps the component visible
verifyLessThanOrEqual(testCase, max(peak_line.YData(in_view)), yl(2));

h.accept_fit();
st = h.state();
verifySize(testCase, st.manual_points, [1 2]);
verifyEqual(testCase, st.manual_points(1, :), [0.005 r.apex_energy_meV], 'AbsTol', 1e-9);
end


function testBackgroundPlusPeaksFitStillRuns(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'bg_peaks';
h.fit_spectrum();

st = h.state();
verifyNotEmpty(testCase, st.pendingFit);
verifyFalse(testCase, isfield(st.pendingFit.result, 'method'));
verifyEmpty(testCase, strfind(char(h.ui.FitInfoLabel.Text), 'ZLP joint'));
end


function testZlpJointFitRefusesOtherPeakModels(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.ui.PeakModelDropdown.Value = 'gaussian';
h.fit_spectrum();

verifyEmpty(testCase, h.state().pendingFit);
verifySubstring(testCase, char(h.ui.FitInfoLabel.Text), 'Drude-Lorentz and Sym. Lorentz only');
end


function k = channelAt(st, q)
[~, k] = min(abs(st.physicalQE.q_Ainv - q));
end


function E0 = truthE0(q)
E0 = 700 + 6e4 * abs(q);
end


function [a3, e] = synthDataset()
% 101 channels around q = 0, |q| <= 0.025 (covers the default comparison reference
% range 0.010-0.020); the ZLP is strongest in the centre channel, where
% make_qe_struct puts q = 0. Gaussian-approximated Poisson noise, fixed seed.
e = (-300:4:2400).';
nq = 101; q0 = 51; dq = 0.0005;
rng(20261022, 'twister');
a3 = zeros(numel(e), nq);
for j = 1:nq
    q = (j - q0) * dq;
    zlp = 1e6 * (1 + 0.2 * exp(-(j - q0) ^ 2 / 8)) * (1 + (e / 18) .^ 2 / 1.5) .^ (-1.5);
    E0 = truthE0(q); G = 400;
    pos = e > 0;
    peak = zeros(size(e));
    peak(pos) = e(pos) * G ./ ((e(pos) .^ 2 - E0 ^ 2) .^ 2 + e(pos) .^ 2 * G ^ 2);
    K = qe_kinematic_prefactor(abs(q), dq=dq, form='2d', beam_kV=30);
    peak(pos) = peak(pos) .* K(e(pos));
    kernel = zlp / sum(zlp);
    peak = conv(peak, kernel(e >= -120 & e <= 120), 'same');
    mu = zlp + peak / max(peak) * 2e3 + 20;
    a3(:, j) = mu + sqrt(mu) .* randn(size(mu));
end
end
