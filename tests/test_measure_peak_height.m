function tests = test_measure_peak_height
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testReturnsMeasuredLocalMaximum(testCase)
energy_meV = (0:10:1000)';
spectrum = zeros(size(energy_meV));
spectrum(48:53) = [0.15; 0.55; 1.10; 1.35; 0.90; 0.30];

peak_height = measure_peak_height(energy_meV, spectrum, 500, 20);

verifyEqual(testCase, peak_height, 1.35, 'AbsTol', 1e-12);
end


function testUsesMinimumWindowAroundPeak(testCase)
energy_meV = (0:5:200)';
spectrum = zeros(size(energy_meV));
spectrum(19:23) = [0.2; 0.6; 1.4; 0.8; 0.3];

peak_height = measure_peak_height(energy_meV, spectrum, 100, 1);

verifyEqual(testCase, peak_height, 1.4, 'AbsTol', 1e-12);
end
