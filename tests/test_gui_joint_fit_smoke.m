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
peak_label = findobj(h.ui.SingleAxes, 'Tag', 'spectrum_fit_overlay', 'Type', 'text');
verifyLessThan(testCase, peak_label.Position(2), yl(2));

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


function testCompleteSessionRestoresDataBranchesFitsAndPendingOverlay(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'complete_session.mat');
h.ui.DqOverrideField.Value = 0.00025;
h.ui.DqOverrideField.ValueChangedFcn(h.ui.DqOverrideField, []);
points = fullfile(folder, 'points.mat');
q = [-0.005; -0.0025; 0.0025; 0.005; 0.0075; 0.010];
manual_dispersion = struct('branches', {{[q 600+3e4*abs(q)], [q 900+2e4*abs(q)]}});
save(points, 'manual_dispersion');
h.load_points(points);
h.ui.CorrectionBranchDropdown.Value = 'Branch 2';
h.ui.FitBranchDropdown.Value = 'Branch 2';
h.ui.EnergyMinField.Value = 300;
h.ui.EnergyMaxField.Value = 1800;
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.select_q(channelAt(h.state(), -0.0025));
h.fit_spectrum();
h.accept_fit();
h.ui.DispModelDropdown.Value = 'acoustic_linear';
h.ui.FitDispButton.ButtonPushedFcn(h.ui.FitDispButton, []);
assertNotEmpty(testCase, h.state().fitResults);
h.ui.NormMethodDropdown.Value = 'ZLP Peak';
h.ui.AreaNormCheckbox.Value = true;
h.ui.ViewModeDropdown.Value = 'Comparison';
h.ui.BuildViewsButton.ButtonPushedFcn(h.ui.BuildViewsButton, []);
h.select_q(channelAt(h.state(), 0.005));
h.fit_spectrum();
before = h.state();
assertNotEmpty(testCase, before.pendingFit);
h.save_session(p);
verifyEqual(testCase, h.state().pendingFit, before.pendingFit);
expected_curves = overlayData(h.ui.SingleAxes);
delete(h.ui.Figure);
[app, restored] = interactive_qe_browser(p);
testCase.addTeardown(@delete, app);
after = restored.state();
for field = {'dataset', 'eq3dQE', 'physicalQE', 'comparisonQE', 'manual_points', ...
        'manual_branches', 'branchCorrectionLog', 'fitResults', 'acceptedFits', 'pendingFit', 'sourceIdentity', 'opHistory'}
    verifyEqual(testCase, after.(field{1}), before.(field{1}));
end
verifyEqual(testCase, restored.ui.ViewModeDropdown.Value, 'Comparison');
verifyEqual(testCase, restored.ui.NormMethodDropdown.Value, 'ZLP Peak');
verifyEqual(testCase, restored.ui.CorrectionBranchDropdown.Value, 'Branch 2');
verifyEqual(testCase, restored.ui.FitBranchDropdown.Value, 'Branch 2');
verifyEqual(testCase, restored.ui.DqOverrideField.Value, 0.00025);
verifyEqual(testCase, restored.ui.SourceLabel.Text, before.dataset.label);
verifyEqual(testCase, restored.ui.ResetStageDropdown.Value, 'eq3d');
verifyEqual(testCase, after.selectedQIndex, before.selectedQIndex);
verifyEqual(testCase, overlayData(restored.ui.SingleAxes), expected_curves);
verifyEqual(testCase, char(restored.ui.AcceptFitButton.Enable), 'on');
restored.accept_fit();
verifyNumElements(testCase, restored.state().acceptedFits, 2);
record = restored.state().acceptedFits{end};
verifyEqual(testCase, record.accepted_energy_meV, before.pendingFit.result.apex_energy_meV);
verifyEqual(testCase, record.provenance.energy_kind, {'loss_function_apex'});
verifyEqual(testCase, record.provenance.settings.fitMethod, 'zlp_joint');
verifyEqual(testCase, record.provenance.signal_window_meV, [300 1800]);
verifyEqual(testCase, record.provenance.input_kind, 'physical_counts');
verifyEqual(testCase, record.result.omega_p, before.pendingFit.result.omega_p);
verifyEqual(testCase, record.result.peak_valid, before.pendingFit.result.peak_valid);
end


function testBackgroundFitPendingSessionRestoresWithoutRefitting(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'background.mat');
h.ui.FitMethodDropdown.Value = 'bg_peaks';
h.fit_spectrum();
before = h.state().pendingFit;
assertNotEmpty(testCase, before);
h.save_session(p);
h.select_q(channelAt(h.state(), -0.005));
h.load_session(p);
verifyEqual(testCase, h.state().pendingFit, before);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'on');
h.accept_fit();
verifyEqual(testCase, h.state().acceptedFits{end}.provenance.input_kind, 'preprocessed_spectrum');
end


function testAcceptedOverlayRestoresWithoutEnablingAccept(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'accepted.mat');
h.ui.FitMethodDropdown.Value = 'zlp_joint'; h.fit_spectrum(); h.accept_fit();
verifyEmpty(testCase, h.state().pendingFit);
curves = overlayData(h.ui.SingleAxes);
verifyNotEmpty(testCase, curves);
h.save_session(p);
h.select_q(channelAt(h.state(), -0.005));
h.load_session(p);
verifyEmpty(testCase, h.state().pendingFit);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'off');
verifyEqual(testCase, overlayData(h.ui.SingleAxes), curves);
points = h.state().manual_points;
h.accept_fit();
verifyEqual(testCase, h.state().manual_points, points);
end


function testPanelExportContainsTheRestoredSpectrum(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
h.ui.FitMethodDropdown.Value = 'zlp_joint'; h.fit_spectrum();
p = fullfile(folder, 'session.mat'); h.save_session(p); h.load_session(p);
paths = qe_export_browser_panel(h.ui.SingleAxes, fullfile(folder, 'spectrum'), ...
    aspect_ratio=4/3, resolution=100, write_pdf=true);
pixels = imread(paths.png);
nonwhite = any(pixels(:,:,1:3) < 240, 3);
verifyGreaterThan(testCase, mean(nonwhite, 'all'), 0.01);
verifyEqual(testCase, size(pixels,2)/size(pixels,1), 4/3, 'AbsTol', 0.02);
verifyTrue(testCase, isfile(paths.pdf));
verifyNotEmpty(testCase, h.state().pendingFit);
end


function testAutoFitDetailsSurviveSessionRoundTrip(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'auto.mat');
h.ui.QStartField.Value = 0.005;
h.ui.QEndField.Value = 0.006;
h.ui.GuessField.Value = '1000';
h.ui.AutoFitButton.ButtonPushedFcn(h.ui.AutoFitButton, []);
before = h.state();
assertNotEmpty(testCase, before.autoFitResults);
assertGreaterThan(testCase, before.autoFitResults.n_success, 0);
h.save_session(p);
h.ui.ClearPtsButton.ButtonPushedFcn(h.ui.ClearPtsButton, []);
h.load_session(p);
verifyEqual(testCase, h.state().autoFitResults, before.autoFitResults);
verifyEqual(testCase, h.state().manual_branches, before.manual_branches);
verifyEqual(testCase, h.state().manual_points, before.manual_points);
end


function testUnavailableSourceCanRecoverEmbeddedCounts(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'offline.mat');
h.save_session(p);
loaded = load(p); session = loaded.session;
% The unavailable source is intentional; embedded counts stay fingerprinted.
session.sources(1).path = fullfile(folder, 'unavailable', 'eq3D.mat');
save(p, 'session');
[app, restored] = interactive_qe_browser();
testCase.addTeardown(@delete, app);
restored.load_session(p);
verifyEqual(testCase, restored.state().physicalQE, h.state().physicalQE);
verifySubstring(testCase, char(restored.ui.InfoLabel.Text), 'using embedded q-E counts');
end


function testChangedSourceIsRejectedWithoutLosingPendingFit(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'session.mat');
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum(); h.save_session(p);
loaded = load(p); session = loaded.session;
copy = fullfile(folder, 'changed.mat');
copyfile(session.sources(1).path, copy);
session.sources(1).path = copy;
fid = fopen(copy, 'ab'); fwrite(fid, uint8(42)); fclose(fid);
save(p, 'session');
before = h.state();
verifyError(testCase, @() h.load_session(p), 'qe_check_browser_sources:SourceChanged');
verifyEqual(testCase, h.state().pendingFit, before.pendingFit);
verifyEqual(testCase, h.state().physicalQE, before.physicalQE);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'on');
end


function testSessionUIFailureRollsBackDataControlsAndOverlay(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'bad_control.mat');
h.save_session(p);
loaded = load(p); session = loaded.session;
session.snapshot.yScale = 'unsupported';
save(p, 'session');
h.ui.DqOverrideField.Value = 0.00025;
h.ui.DqOverrideField.ValueChangedFcn(h.ui.DqOverrideField, []);
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum();
before = h.state();
curves = overlayData(h.ui.SingleAxes);
try
    h.load_session(p);
    verifyFail(testCase, 'Invalid control value unexpectedly restored');
catch ME
    verifyNotEmpty(testCase, ME.identifier);
end
verifyEqual(testCase, h.state().pendingFit, before.pendingFit);
verifyEqual(testCase, h.state().dataset, before.dataset);
verifyEqual(testCase, h.state().manual_points, before.manual_points);
verifyEqual(testCase, h.state().opHistory, before.opHistory);
verifyEqual(testCase, char(h.ui.AcceptFitButton.Enable), 'on');
verifyEqual(testCase, overlayData(h.ui.SingleAxes), curves);
verifyEqual(testCase, h.ui.YScaleDropdown.Value, 'linear');
end


function testInvalidSessionCountsVersionAndPendingChannelAreRejected(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'invalid.mat');
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.fit_spectrum(); h.save_session(p);
loaded = load(p); original = loaded.session;
before = h.state();
for scenario = 1:5
    session = original;
    switch scenario
        case 1, session.data.physicalQE.intensity(1) = session.data.physicalQE.intensity(1) + 1;
        case 2, session.version = 999;
        case 3, session.pending_fit.q_index = 1;
        case 4, session.snapshot.normMethod = 'unsupported';
        case 5, session.pending_fit.keep = ~session.pending_fit.keep;
    end
    save(p, 'session');
    verifyError(testCase, @() h.load_session(p), 'qe_validate_browser_session:InvalidSession');
    verifyEqual(testCase, h.state().pendingFit, before.pendingFit);
    verifyEqual(testCase, h.state().physicalQE, before.physicalQE);
end
end


function testSessionCancellationAndFailedSavePreserveExistingFile(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'existing.mat');
h.save_session(p);
original_hash = qe_file_sha256(p);
before = h.state();
h.save_session(''); h.load_session('');
verifyEqual(testCase, h.state(), before);
h.ui.FitMethodDropdown.Value = 'zlp_joint'; h.fit_spectrum();
% Changing a fit control without its callback must not persist a stale fit.
h.ui.BeamKVField.Value = 60;
verifyError(testCase, @() h.save_session(p), 'qe_validate_browser_session:InvalidSession');
verifyEqual(testCase, qe_file_sha256(p), original_hash);
verifyEqual(testCase, numel(dir(fullfile(folder, '*.mat'))), 1);
source = h.state().dataset.source_path;
verifyError(testCase, @() qe_write_browser_session(source, loadSession(p)), ...
    'qe_write_browser_session:SourceOverwrite');
end


function testAcceptedFitProvenanceSurvivesPointsExportAndLegacyImport(testCase)
h = testCase.TestData.hooks;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'points.mat');
h.ui.FitMethodDropdown.Value = 'zlp_joint'; h.fit_spectrum(); h.accept_fit();
records = h.state().acceptedFits;
h.save_points(p); h.load_points(p);
verifyEqual(testCase, h.state().acceptedFits, records);
loaded = load(p); manual_dispersion = rmfield(loaded.manual_dispersion, 'accepted_fits');
save(p, 'manual_dispersion'); h.load_points(p);
verifyEmpty(testCase, h.state().acceptedFits);
verifyNotEmpty(testCase, h.state().manual_points);
end


function session = loadSession(path)
loaded = load(path); session = loaded.session;
end


function values = overlayData(ax)
lines = findobj(ax, 'Type', 'line', 'Tag', 'spectrum_fit_overlay');
values = arrayfun(@(line) struct('X', line.XData, 'Y', line.YData, 'name', line.DisplayName), lines);
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
