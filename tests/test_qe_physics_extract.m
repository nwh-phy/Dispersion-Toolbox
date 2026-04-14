function tests = test_qe_physics_extract
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testRecoversScreeningLengthFromQuasi2DDispersion(testCase)
q_abs = [0.01; 0.02; 0.03; 0.04; 0.06; 0.08; 0.10];
q_signed = [-flipud(q_abs); q_abs];

drude_weight_true = 1.8e6;  % meV^2 * Angstrom
rho0_true = 25;
epsilon_inter_true = 1 + rho0_true * abs(q_signed);
omega_p = sqrt(drude_weight_true .* abs(q_signed) ./ epsilon_inter_true);

q_crossover_true = 0.04;
I_kin_true = 8 .* abs(q_signed).^(-3) .* (1 + abs(q_signed) ./ q_crossover_true);
A_measured = I_kin_true ./ epsilon_inter_true;
gamma = 0.12 .* omega_p;

branch = nan(numel(q_signed), 12);
branch(:, 1) = q_signed;
branch(:, 2) = omega_p;
branch(:, 3) = gamma;
branch(:, 4) = 0.99;
branch(:, 5) = A_measured * 0.4;  % intentionally different from measured intensity
branch(:,12) = A_measured;

phys = qe_physics_extract({branch}, struct('q_min_fit', 0.01, 'q_max_fit', 0.10));

verifyEqual(testCase, numel(phys.q), numel(q_abs));
verifyEqual(testCase, phys.rho0, rho0_true, 'AbsTol', 2.5);
verifyEqual(testCase, phys.epsilon_inter, 1 + rho0_true .* q_abs, 'RelTol', 0.12);
end


function testDefaultFitUsesFullAvailableQRange(testCase)
q_abs = [0.01; 0.02; 0.03; 0.04; 0.06; 0.08; 0.10; 0.12];
q_signed = [-flipud(q_abs); q_abs];

drude_weight_true = 1.5e6;
rho0_true = 22;
epsilon_inter_true = 1 + rho0_true * abs(q_signed);
omega_p = sqrt(drude_weight_true .* abs(q_signed) ./ epsilon_inter_true);
gamma = 0.10 .* omega_p;
A_measured = 5 ./ epsilon_inter_true;

branch = nan(numel(q_signed), 12);
branch(:, 1) = q_signed;
branch(:, 2) = omega_p;
branch(:, 3) = gamma;
branch(:, 4) = 0.98;
branch(:, 5) = A_measured;
branch(:,12) = A_measured;

phys_default = qe_physics_extract({branch});
phys_limited = qe_physics_extract({branch}, struct('q_max_fit', 0.06));

verifyEqual(testCase, phys_default.Drude_fit.q_range(2), max(q_abs), 'AbsTol', 1e-12);
verifyLessThan(testCase, abs(phys_default.rho0 - rho0_true), abs(phys_limited.rho0 - rho0_true) + 1e-9);
end


function testKeepsStrongestCandidatePerSignedQBeforeSymmetryAverage(testCase)
q_abs = [0.01; 0.03; 0.05];
q_signed = [-flipud(q_abs); q_abs];

drude_weight_true = 9.0e5;
rho0_true = 18;
epsilon_inter_true = 1 + rho0_true * abs(q_signed);
omega_true = sqrt(drude_weight_true .* abs(q_signed) ./ epsilon_inter_true);
gamma_true = 0.08 .* omega_true;
A_true = 4 ./ epsilon_inter_true;

spurious_energy = omega_true .* 0.65;
spurious_gamma = gamma_true .* 2.5;
spurious_A = A_true .* 0.08;

branch_true = nan(numel(q_signed), 12);
branch_true(:, 1) = q_signed;
branch_true(:, 2) = omega_true;
branch_true(:, 3) = gamma_true;
branch_true(:, 4) = 0.99;
branch_true(:, 5) = A_true;
branch_true(:,12) = A_true;

branch_spurious = branch_true;
branch_spurious(:, 2) = spurious_energy;
branch_spurious(:, 3) = spurious_gamma;
branch_spurious(:, 5) = spurious_A;
branch_spurious(:,12) = spurious_A;

branch = [branch_true; branch_spurious];
phys = qe_physics_extract({branch}, struct('q_min_fit', 0.01, 'q_max_fit', 0.05));

verifyEqual(testCase, phys.q, q_abs, 'AbsTol', 1e-12);
verifyEqual(testCase, phys.omega_p, omega_true(numel(q_abs)+1:end), 'RelTol', 0.08);
end
