function tests = test_qe_fit_point_energy
% Which energy of a fitted peak goes into the dispersion.
tests = functiontests(localfunctions);
end


function setupOnce(~)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testZlpJointFitsGiveTheModelE0(testCase)
r = struct('method', 'zlp_joint', 'omega_p', [700; 1600], 'apex_energy_meV', [640; 1400]);
[energy, kind] = qe_fit_point_energy(r);
verifyEqual(testCase, energy, [700; 1600]);
verifyEqual(testCase, kind, {'model_E0'; 'model_E0'});
end


function testBackgroundPeakFitsGiveTheApexWhereFinite(testCase)
r = struct('omega_p', [700; 1600], 'apex_energy_meV', [640; NaN]);
[energy, kind] = qe_fit_point_energy(r);
verifyEqual(testCase, energy, [640; 1600]);
verifyEqual(testCase, kind, {'loss_function_apex'; 'model_E0'});
end


function testRecordedKindsAreKept(testCase)
% A fit saved while ZLP joint points were apexes keeps its apex.
r = struct('method', 'zlp_joint', 'omega_p', [700; 1600], 'apex_energy_meV', [640; 1400]);
verifyEqual(testCase, qe_fit_point_energy(r, {'loss_function_apex'; 'model_E0'}), [640; 1600]);
end
