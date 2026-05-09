classdef PeakFwhmMeasurementTest < matlab.unittest.TestCase
    methods (Test)
        function lorentzCurveReturnsExpectedFwhm(testCase)
            energy = linspace(500, 1500, 5001)';
            center = 1000;
            fwhm = 240;
            half_width = fwhm / 2;
            curve = half_width^2 ./ ((energy - center).^2 + half_width^2);

            measured = measure_peak_fwhm(energy, curve);

            testCase.verifyEqual(measured.fwhm_meV, fwhm, 'AbsTol', 1.0);
            testCase.verifyEqual(measured.left_meV, center - half_width, 'AbsTol', 1.0);
            testCase.verifyEqual(measured.right_meV, center + half_width, 'AbsTol', 1.0);
        end

        function asymmetricCurveUsesLocalContrastHalfMaximum(testCase)
            energy = linspace(0, 1000, 5001)';
            baseline = 0.25;
            curve = baseline + exp(-((energy - 420) / 75).^2) .* ...
                (1 + 0.35 * (energy - 420) / 200);

            measured = measure_peak_fwhm(energy, curve);

            testCase.verifyTrue(isfinite(measured.fwhm_meV));
            testCase.verifyGreaterThan(measured.fwhm_meV, 50);
            testCase.verifyLessThan(measured.fwhm_meV, 250);
            testCase.verifyGreaterThan(measured.apex_meV, measured.left_meV);
            testCase.verifyLessThan(measured.apex_meV, measured.right_meV);
        end
    end
end
