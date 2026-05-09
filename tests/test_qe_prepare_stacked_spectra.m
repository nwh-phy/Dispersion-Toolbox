function tests = test_qe_prepare_stacked_spectra
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testSelectsNearestQChannelsInStableOrderInsideRequestedRange(testCase)
energy_meV = (0:50:200)';
spectra = [
    1  10  100  1000  10000;
    2  20  200  2000  20000;
    3  30  300  3000  30000;
    4  40  400  4000  40000;
    5  50  500  5000  50000];
q_axis = [-0.20 -0.05 0.00 0.08 0.21];
opts = struct('q_start', -0.20, 'q_end', 0.20, 'q_step', 0.10, 'offset', 2.0);

stack = qe_prepare_stacked_spectra(energy_meV, spectra, q_axis, opts);

verifyEqual(testCase, stack.q_indices, [1 2 3 4]);
verifyEqual(testCase, stack.q_values, q_axis([1 2 3 4]), 'AbsTol', 1e-12);
verifyEqual(testCase, stack.energy_meV, energy_meV);
verifyEqual(testCase, stack.traces, spectra(:, [1 2 3 4]));
end


function testDoesNotSelectNearestChannelsOutsideQRange(testCase)
energy_meV = [100; 200; 300];
spectra = [
    1   10   100;
    2   20   200;
    3   30   300];
q_axis = [0.00 0.04 0.08];
opts = struct('q_start', 0.02, 'q_end', 0.06, 'q_step', 0.02, 'offset', 1.0);

stack = qe_prepare_stacked_spectra(energy_meV, spectra, q_axis, opts);

verifyEqual(testCase, stack.q_indices, 2);
verifyEqual(testCase, stack.q_values, 0.04, 'AbsTol', 1e-12);
verifyEqual(testCase, stack.traces, spectra(:, 2));
end


function testKeepsNumericallyEquivalentBoundaryQChannel(testCase)
energy_meV = [100; 200];
spectra = [
    1   10   100   1000;
    2   20   200   2000];
q_axis = [0.0 0.1 0.2 0.1 + 0.2];
opts = struct('q_start', 0.0, 'q_end', 0.3, 'q_step', 0.1, 'offset', 1.0);

stack = qe_prepare_stacked_spectra(energy_meV, spectra, q_axis, opts);

verifyEqual(testCase, stack.q_indices, [1 2 3 4]);
verifyEqual(testCase, stack.q_values, q_axis, 'AbsTol', 1e-12);
end


function testAppliesOffsetsAndCollapsesDuplicateTargets(testCase)
energy_meV = [0; 10; 20];
spectra = [
    1   10   100;
    2   20   200;
    3   30   300];
q_axis = [-0.10 0.00 0.10];
opts = struct('q_start', 0.10, 'q_end', -0.10, 'q_step', 0.05, 'offset', 5.0);

stack = qe_prepare_stacked_spectra(energy_meV, spectra, q_axis, opts);

verifyEqual(testCase, stack.q_indices, [1 2 3]);
verifyEqual(testCase, stack.offsets, [0 5 10], 'AbsTol', 1e-12);
verifyEqual(testCase, stack.shifted_traces, [
    1  15  110;
    2  25  210;
    3  35  310], 'AbsTol', 1e-12);
end
