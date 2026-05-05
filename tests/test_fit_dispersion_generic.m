function tests = test_fit_dispersion_generic
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testQuasi2DUsesSubstrateDielectricInFit(testCase)
q = linspace(0.01, 0.12, 18).';
epsilon_s = 9;
epsilon_bg = (1 + epsilon_s) / 2;
A_true = 2.4e7;
rho0_true = 30;
energy = sqrt(A_true .* q ./ (epsilon_bg + rho0_true .* q));

fit = fit_dispersion_generic(q, energy, ...
    'model', 'quasi2d_plasmon', ...
    'epsilon_s', epsilon_s);

verifyEqual(testCase, fit.epsilon_s, epsilon_s);
verifyEqual(testCase, fit.rho0, rho0_true, 'RelTol', 0.03);
verifyEqual(testCase, fit.q_c_Ainv, epsilon_bg / rho0_true, 'RelTol', 0.03);
verifyGreaterThan(testCase, fit.R_squared, 0.999);
end
