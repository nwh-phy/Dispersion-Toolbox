function tests = test_bisb_slot_regrid
tests = functiontests(localfunctions);
end


function setupOnce(~)
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts'));
end


function testIdentityOnIntegerGrid(testCase)
E = (-100:4:400)'; q = 1:50; X = rand(numel(E), numel(q), 3);
Y = bisb_slot_regrid(E, q, X, zeros(1, 50), 20, q - 20);
verifyEqual(testCase, Y, X, 'AbsTol', 1e-12);
end


function testElasticZeroMovesToZero(testCase)
E = (-200:4:400)'; q = 1:30; z = linspace(-25, 38, 30) + 0.7;
X = exp(-(E - z).^2 / (2 * 12^2));
Y = bisb_slot_regrid(E, q, X, z, 15, q - 15);
m = abs(E) <= 40; c = sum(E(m) .* Y(m, :)) ./ sum(Y(m, :));
verifyLessThan(testCase, max(abs(c)), 0.5);
end


function testEnergyDependentDistortionIsRemoved(testCase)
% Two features at true k0 drawn where the distortion model puts them;
% after regridding each sits at k0 at every energy.
E = (0:20:3000)'; c0 = 120.4; q = 1:250; gamma = 0.55; beta = -0.0044; s = 2.5;
X = zeros(numel(E), numel(q));
for k0 = [-90 60]
    X = X + exp(-(q - (c0 + k0 + (gamma + beta * k0) * E / 1000)).^2 / (2 * s^2));
end
k = -110:100; Y = bisb_slot_regrid(E, q, X, zeros(size(q)), c0, k, gamma, beta);
for k0 = [-90 60]
    w = abs(k - k0) <= 10; c = sum(k(w) .* Y(:, w), 2) ./ sum(Y(:, w), 2);
    verifyLessThan(testCase, max(abs(c - k0)), 0.1);
end
% Without the correction the far feature drifts by more than one channel.
Y0 = bisb_slot_regrid(E, q, X, zeros(size(q)), c0, k);
w = abs(k + 90) <= 10; c = sum(k(w) .* Y0(:, w), 2) ./ sum(Y0(:, w), 2);
verifyGreaterThan(testCase, max(abs(c + 90)), 1);
end


function testTrailingDimensionsMatchSeparateCalls(testCase)
E = (-40:4:200)'; q = 1:40; X = rand(numel(E), numel(q), 2, 2); z = randn(1, 40) * 5;
Y = bisb_slot_regrid(E, q, X, z, 18.3, -15:15, 0.4, -0.003);
for r = 1:4
    Yr = bisb_slot_regrid(E, q, X(:, :, r), z, 18.3, -15:15, 0.4, -0.003);
    verifyEqual(testCase, Y(:, :, r), Yr, 'AbsTol', 1e-12);
end
verifyTrue(testCase, any(isnan(Y(:))), 'Shifted energies outside the support are NaN');
end
