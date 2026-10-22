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


function testFailedRefitCannotAcceptThePreviousResult(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();
verifyNotEmpty(testCase, h.state().pendingFit);
h.ui.PeakModelDropdown.Value = 'gaussian';
h.fit_spectrum();
verifyEmpty(testCase, h.state().pendingFit);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'off');
verifyEmpty(testCase, findobj(h.ui.SingleAxes, 'Tag', 'spectrum_fit_overlay'));
h.accept_fit();
verifyEmpty(testCase, h.state().manual_points);
end


function testChangingChannelOrFitSettingsInvalidatesPendingFit(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();
h.select_q(channelAt(h.state(), -0.005));
verifyEmpty(testCase, h.state().pendingFit);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'off');
h.accept_fit();
verifyEmpty(testCase, h.state().manual_points);
h.fit_spectrum();
verifyNotEmpty(testCase, h.state().pendingFit);
h.ui.BeamKVField.Value = 60;
h.ui.BeamKVField.ValueChangedFcn(h.ui.BeamKVField, []);
verifyEmpty(testCase, h.state().pendingFit);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'off');
end


function testHistoryRestoresNormalizationCalibrationAndSelectedChannel(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
history_path = fullfile(folder, 'history.mat');
h.ui.DqOverrideField.Value = 0.00025;
h.ui.DqOverrideField.ValueChangedFcn(h.ui.DqOverrideField, []);
h.select_q(channelAt(h.state(), -0.0025));
h.ui.NormMethodDropdown.Value = 'ZLP Peak';
h.ui.AreaNormCheckbox.Value = true;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.ui.KinematicCheckbox.Value = false;
h.ui.BeamKVField.Value = 60;
h.ui.PoissonWeightsCheckbox.Value = false;
h.save_history(history_path);
saved = load(history_path);
snapshot = saved.opHistory{end}.snapshot;
delete(h.ui.Figure);
[app, restored] = interactive_qe_browser(testCase.TestData.data_dir, history_path);
testCase.addTeardown(@delete, app);
verifyEqual(testCase, restored.ui.NormMethodDropdown.Value, 'ZLP Peak');
verifyTrue(testCase, restored.ui.AreaNormCheckbox.Value);
verifyEqual(testCase, restored.ui.FitMethodDropdown.Value, 'zlp_joint');
verifyFalse(testCase, restored.ui.KinematicCheckbox.Value);
verifyEqual(testCase, restored.ui.BeamKVField.Value, 60);
verifyFalse(testCase, restored.ui.PoissonWeightsCheckbox.Value);
st = restored.state();
verifyEqual(testCase, st.physicalQE.dq_Ainv, 0.00025);
verifyEqual(testCase, st.dataset.dq_Ainv, st.physicalQE.dq_Ainv);
verifyEqual(testCase, st.physicalQE.q_Ainv(st.selectedQIndex), -0.0025, 'AbsTol', 1e-12);
verifyEqual(testCase, st.selectedQIndex, snapshot.selectedQIndex);
verifyEqual(testCase, st.stackOffsetMode, string(snapshot.stackOffsetMode));
end


function testOldHistoryRestoresLegacyMethodWithoutNewFields(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
history_path = fullfile(folder, 'old_history.mat');
h.save_history(history_path);
saved = load(history_path);
opHistory = saved.opHistory;
new_fields = {'fitMethod', 'kinematicFactor', 'beamKV', 'poissonWeights', ...
    'normMethod', 'stackOffsetMode'};
for k = 1:numel(opHistory)
    fields = intersect(fieldnames(opHistory{k}.snapshot), new_fields);
    opHistory{k}.snapshot = rmfield(opHistory{k}.snapshot, fields);
end
save(history_path, 'opHistory');
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.ui.NormMethodDropdown.Value = 'ZLP Peak';
h.load_history(history_path);
verifyEqual(testCase, h.ui.FitMethodDropdown.Value, 'bg_peaks');
verifyEqual(testCase, h.ui.NormMethodDropdown.Value, 'Area');
verifySubstring(testCase, char(h.ui.InfoLabel.Text), 'History loaded');
end


function testBrokenHistoryDoesNotClaimSuccessfulRestoration(testCase)
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
history_path = fullfile(folder, 'broken_history.mat');
opHistory = {struct('label', 'broken', 'timestamp', '00:00:00', 'snapshot', struct())};
save(history_path, 'opHistory');
[app, restored] = interactive_qe_browser(testCase.TestData.data_dir, history_path);
testCase.addTeardown(@delete, app);
verifySubstring(testCase, char(restored.ui.InfoLabel.Text), 'History restore failed');
verifyEmpty(testCase, strfind(char(restored.ui.InfoLabel.Text), 'History loaded'));
end


function testExportSettingsPreservePendingFitAndOverlay(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();
before = h.state().pendingFit;
assertNotEmpty(testCase, before);
lines = findobj(h.ui.SingleAxes, 'Tag', 'spectrum_fit_overlay');
assertNotEmpty(testCase, lines);
h.ui.ExportModeDropdown.Value = 'Publication';
h.ui.ExportModeDropdown.ValueChangedFcn(h.ui.ExportModeDropdown, []);
h.ui.ExportRatioDropdown.Value = '4:3';
h.ui.ExportRatioDropdown.ValueChangedFcn(h.ui.ExportRatioDropdown, []);
verifyEqual(testCase, h.state().pendingFit, before);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'on');
verifyEqual(testCase, findobj(h.ui.SingleAxes, 'Tag', 'spectrum_fit_overlay'), lines);
snapshot = h.state().opHistory{end}.snapshot;
verifyEqual(testCase, snapshot.exportMode, 'Publication');
verifyEqual(testCase, snapshot.exportRatio, '4:3');
h.accept_fit();
verifyEqual(testCase, h.state().manual_points, ...
    [before.q_Ainv before.result.apex_energy_meV], 'AbsTol', 1e-12);
end


function testIncompleteHistoryIsRejectedBeforeChangingLiveState(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'incomplete_history.mat');
h.save_history(p);
loaded = load(p);
opHistory = loaded.opHistory;
opHistory{end}.snapshot.dqOverride = 0.00025;
opHistory{end}.snapshot = rmfield(opHistory{end}.snapshot, 'energyMin');
save(p, 'opHistory');
before = h.state();
verifyError(testCase, @() h.load_history(p), 'qe_read_browser_history:InvalidSnapshot');
after = h.state();
verifyEqual(testCase, after.physicalQE, before.physicalQE);
verifyEqual(testCase, after.dataset, before.dataset);
verifyEqual(testCase, after.selectedQIndex, before.selectedQIndex);
verifyEqual(testCase, after.pendingFit, before.pendingFit);
verifyEqual(testCase, after.opHistory, before.opHistory);
verifyEqual(testCase, h.ui.DqOverrideField.Value, before.physicalQE.dq_Ainv);
h.accept_fit();
verifyEqual(testCase, h.state().manual_points(1, 1), ...
    before.physicalQE.q_Ainv(before.selectedQIndex), 'AbsTol', 1e-12);
end


function testHistoryValueErrorRollsBackCalibrationAndInvalidatesFit(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'invalid_value_history.mat');
h.save_history(p);
loaded = load(p);
opHistory = loaded.opHistory;
opHistory{end}.snapshot.dqOverride = 0.00025;
opHistory{end}.snapshot.yScale = 'unsupported';
save(p, 'opHistory');
before = h.state();
old_y_scale = h.ui.YScaleDropdown.Value;
threw = false;
try
    h.load_history(p);
catch
    threw = true;
end
verifyTrue(testCase, threw, 'Invalid UI values must fail restoration.');
after = h.state();
verifyEqual(testCase, after.physicalQE, before.physicalQE);
verifyEqual(testCase, after.dataset, before.dataset);
verifyEqual(testCase, after.selectedQIndex, before.selectedQIndex);
verifyEqual(testCase, after.opHistory, before.opHistory);
verifyEqual(testCase, after.manual_points, before.manual_points);
verifyEqual(testCase, h.ui.DqOverrideField.Value, before.physicalQE.dq_Ainv);
verifyEqual(testCase, h.ui.YScaleDropdown.Value, old_y_scale);
verifyEmpty(testCase, after.pendingFit);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'off');
verifyEmpty(testCase, findobj(h.ui.SingleAxes, 'Tag', 'spectrum_fit_overlay'));
h.accept_fit();
verifyEqual(testCase, h.state().manual_points, before.manual_points);
end


function testCuratedPointsCorrectionsAndTargetsSurviveReopening(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
points_path = fullfile(folder, 'points.mat');
history_path = fullfile(folder, 'history.mat');
manual_dispersion = struct('branches', {{[-0.005 600; 0.005 700], ...
    [-0.005 900; 0.005 1000]}});
save(points_path, 'manual_dispersion');
h.load_points(points_path);
h.select_q(channelAt(h.state(), -0.005));
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.ui.CorrectionBranchDropdown.Value = 'Branch 2';
h.fit_spectrum();
result = h.state().pendingFit.result;
h.accept_fit();
before = h.state();
verifyEqual(testCase, before.manual_branches{2}(1, 1:2), ...
    [-0.005 result.apex_energy_meV], 'AbsTol', 1e-12);
verifyNumElements(testCase, before.branchCorrectionLog, 1);
h.save_points(points_path);
h.save_history(history_path);
saved = load(points_path);
verifyEqual(testCase, saved.manual_dispersion.q_Ainv, before.manual_points(:, 1));
verifyEqual(testCase, saved.manual_dispersion.energy_meV, before.manual_points(:, 2));
verifyEqual(testCase, saved.manual_dispersion.abs_q_Ainv, abs(before.manual_points(:, 1)));
delete(h.ui.Figure);
[app, restored] = interactive_qe_browser(testCase.TestData.data_dir);
testCase.addTeardown(@delete, app);
restored.load_points(points_path);
restored.load_history(history_path);
after = restored.state();
verifyEqual(testCase, after.manual_points, before.manual_points);
verifyEqual(testCase, after.manual_branches, before.manual_branches);
verifyEqual(testCase, after.branchCorrectionLog, before.branchCorrectionLog);
verifyEqual(testCase, restored.ui.CorrectionBranchDropdown.Value, 'Branch 2');
verifyEqual(testCase, after.selectedQIndex, before.selectedQIndex);
end


function testDataExportUsesPhysicalAxesAndPreprocessingInComparisonView(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
export_path = fullfile(folder, 'export.mat');
h.ui.ViewModeDropdown.Value = 'Comparison';
h.ui.AreaNormCheckbox.Value = true;
h.ui.NormMethodDropdown.Value = 'ZLP Peak';
h.export_data(export_path);
out = load(export_path);
qe = h.state().physicalQE;
[expected, ~] = qe_preprocess(qe, out.preprocessing);
verifyEqual(testCase, out.data_kind, 'physical');
verifyEqual(testCase, out.view_mode_at_export, 'Comparison');
verifyEqual(testCase, out.intensity_unprocessed, qe.intensity);
verifyEqual(testCase, out.intensity, expected.intensity, 'AbsTol', 1e-12);
verifyEqual(testCase, out.q_Ainv, qe.q_Ainv(:));
verifyEqual(testCase, out.energy_meV, qe.energy_meV(:));
verifyEqual(testCase, out.dq_Ainv, qe.dq_Ainv);
verifyEqual(testCase, out.preprocessing.norm_method, 'ZLP Peak');
verifyEqual(testCase, out.dimensions.n_q, numel(qe.q_Ainv));
verifyEqual(testCase, out.dimensions.n_energy, numel(qe.energy_meV));
end


function testCancelledFileOperationsDoNotChangeCuratedState(testCase)
h = testCase.TestData.hooks;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();
h.accept_fit();
before = h.state();
h.save_history('');
h.load_history('');
h.save_points('');
h.load_points('');
h.export_data('');
after = h.state();
verifyEqual(testCase, after.manual_points, before.manual_points);
verifyEqual(testCase, after.manual_branches, before.manual_branches);
verifyEqual(testCase, after.branchCorrectionLog, before.branchCorrectionLog);
verifyEqual(testCase, after.opHistory, before.opHistory);
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
