function tests = test_qe_spectrum_helpers
% qe_zlp_center, qe_pm_spectra, qe_powerlaw_tail, qe_peak_apex, qe_plot_heatmap
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testZlpCenterMatchesTheP2aRule(testCase)
E = (-200:4:600)';
for c0 = [-3.1 0 1.7 2.4]
    y = 1e6 * exp(-(E - c0).^2 / (2 * 9^2));
    [c, fw] = qe_zlp_center(y, E);
    % the old b1_p2a_inputs local_centre, written out
    dE = E(2) - E(1); [~, i0] = max(y .* (abs(E) < 40));
    ly = log(max(y(i0-1:i0+1), 1)); cref = E(i0) + (ly(3) - ly(1)) / (2 * (2 * ly(2) - ly(1) - ly(3))) * dE;
    verifyEqual(testCase, c, cref, 'AbsTol', 1e-12);
    verifyEqual(testCase, c, c0, 'AbsTol', 1e-9);          % exact for a Gaussian
    verifyEqual(testCase, fw, 2 * sqrt(2 * log(2)) * 9, 'RelTol', 0.05);
end
end


function testPmSpectraAveragesBothSides(testCase)
q = (-3:3) * 0.0005; M = repmat(1:7, 4, 1);
[S, found] = qe_pm_spectra(M, q, [0 0.0005 0.0015 0.002]);
verifyEqual(testCase, found, [true true true false]);
verifyEqual(testCase, S(1, 1:3), [4 4 4]);
verifyTrue(testCase, all(isnan(S(:, 4))));
verifyEqual(testCase, qe_pm_spectra(M, q, 0.001, 'plus'), 6 * ones(4, 1));
verifyEqual(testCase, qe_pm_spectra(M, q, 0.001, 'minus'), 2 * ones(4, 1));
end


function testPowerlawTailRecoversTheExponent(testCase)
E = (-100:4:1800)'; y = zeros(size(E)); y(E > 0) = 3e7 * E(E > 0) .^ -3.2;
[bg, r, A] = qe_powerlaw_tail(E, y, [150 250]);
verifyEqual(testCase, r, 3.2, 'AbsTol', 1e-9);
verifyEqual(testCase, A, 3e7, 'RelTol', 1e-9);
verifyEqual(testCase, bg, y, 'RelTol', 1e-9);
end


function testPeakApexFindsAPeakOnATail(testCase)
E = (0:4:1800)';
y = 1e4 * max(E, 4) .^ -1 + 800 * exp(-(E - 600).^2 / (2 * 100^2));
verifyTrue(testCase, isnan(qe_peak_apex(E, 1e4 * max(E, 4) .^ -1, [250 1000])));   % tail only
tail = qe_powerlaw_tail(E, y, [150 250]);
e = qe_peak_apex(E, y - tail, [250 1000]);
verifyEqual(testCase, e, 600, 'AbsTol', 10);
end


function testPlotHeatmapNormalisesColumns(testCase)
E = (0:4:2000)'; q = (-10:10) * 0.0005;
M = (1 + abs(q)) .* exp(-(E - 600).^2 / 2e4) + 1;
base = fullfile(tempname, 'map');
info = qe_plot_heatmap(E, q, M, base, q_window=[-0.003 0.003]);
verifyTrue(testCase, isfile([base '.png']) && isfile([base '.pdf']));
verifyEqual(testCase, numel(info.q), 13);
verifyEqual(testCase, trapz(info.E, info.map, 1), ones(1, 13), 'RelTol', 1e-12);
verifyLessThan(testCase, info.clim(1), info.clim(2));
rmdir(fileparts(base), 's');
end
