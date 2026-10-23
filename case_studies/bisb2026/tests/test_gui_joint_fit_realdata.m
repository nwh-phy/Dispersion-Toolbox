function tests = test_gui_joint_fit_realdata
% Bounded numerical acceptance: the GUI must fit the same raw channel as the API.
% These checks do not establish the physical validity of a single-channel model.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
case_root = fileparts(fileparts(mfilename('fullpath')));
scripts_dir = fullfile(case_root, 'scripts');
addpath(scripts_dir);
project_root = bisb_find_project_root(scripts_dir);
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function test590PositiveAndLowNegativeQ(testCase)
checkDataset(testCase, '590 PL2 10w 0.004 10sx300', [-0.002, 0.005]);
end


function testN0RepeatPositiveAndNegativeQ(testCase)
checkDataset(testCase, 'n0 pl2 10w 0.004 10s x300', [-0.005, 0.002]);
end


function test20wLowQAndWeakerHighQ(testCase)
checkDataset(testCase, 'no pl2 20w 0.004 10sx300 2film', [-0.002, 0.01]);
end


function checkDataset(testCase, session, q_values)
dataset_path = fullfile(testCase.TestData.project_root, '20260120 BiSb', session, 'eq3D.mat');
assumeTrue(testCase, isfile(dataset_path), 'Local BiSb eq3D dataset is unavailable.');
[app, h] = interactive_qe_browser(dataset_path);
testCase.addTeardown(@delete, app);
h.ui.FitMethodDropdown.Value = 'zlp_joint';
h.ui.MaxPeaksField.Value = 2;
h.ui.EnergyMinField.Value = 300;
h.ui.EnergyMaxField.Value = 1800;
h.ui.GuessField.Value = '650,1350';
h.ui.AreaNormCheckbox.Value = true;
h.ui.NormMethodDropdown.Value = 'ZLP Peak';
h.ui.ViewModeDropdown.Value = 'Comparison';
qe = h.state().physicalQE;
for q = q_values
    [~, k] = min(abs(qe.q_Ainv - q));
    h.select_q(k);
    h.fit_spectrum();
    pending = h.state().pendingFit;
    assertNotEmpty(testCase, pending, char(h.ui.FitInfoLabel.Text));
    gui = pending.result;
    direct = qe_fit_spectrum_joint(double(qe.energy_meV(:)), double(qe.intensity(:, k)), ...
        qe.q_Ainv(k), n_peaks=2, peak_model='lorentz', signal_window=[300 1800], ...
        dq=median(abs(diff(qe.q_Ainv))), kinematic=true, beam_kV=30, ...
        weights='poisson', initial_E0=[650 1350]);
    verifyTrue(testCase, direct.success);
    verifyEqual(testCase, pending.q_Ainv, qe.q_Ainv(k));
    verifyEqual(testCase, gui.omega_p, direct.omega_p, 'AbsTol', 1e-8);
    verifyEqual(testCase, gui.gamma, direct.gamma, 'AbsTol', 1e-8);
    verifyEqual(testCase, gui.apex_energy_meV, direct.apex_energy_meV, 'AbsTol', 1e-8);
    verifyEqual(testCase, gui.curve_fit, direct.curve_fit, 'RelTol', 1e-12);
    verifyEqual(testCase, gui.peak_valid, direct.peak_valid);
    verifyEqual(testCase, gui.chi2_red, direct.chi2_red, 'RelTol', 1e-12);
    fprintf('GUI/API %s q=%+.6f: E0=%s apex=%s valid=%s\n', session, ...
        pending.q_Ainv, mat2str(gui.omega_p.', 7), mat2str(gui.apex_energy_meV.', 7), ...
        mat2str(gui.peak_valid.'));
end
end
