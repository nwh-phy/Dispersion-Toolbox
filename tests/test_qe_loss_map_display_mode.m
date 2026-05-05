function tests = test_qe_loss_map_display_mode
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testDisplayModeDoesNotApplyIkinCorrectionOrClipSignedData(testCase)
old_visible = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
cleanup = onCleanup(@() set(groot, 'DefaultFigureVisible', old_visible));

qe = struct();
qe.energy_meV = [100; 200; 300];
qe.q_Ainv = [-0.1 0.1];
qe.intensity = [1 -2; 3 4; 5 6];

phys = struct();
phys.q = [0.1];
phys.I_kin = [1000];
phys.omega_p = [200];
phys.rho0 = 25;

fig = qe_loss_map(qe, phys, struct('mode', 'display'));
close_fig = onCleanup(@() close(fig));

images = findobj(fig, 'Type', 'image');
cdatas = arrayfun(@(h) h.CData, images, 'UniformOutput', false);
has_original_map = any(cellfun(@(c) isequal(c, qe.intensity), cdatas));

verifyTrue(testCase, has_original_map);
end
