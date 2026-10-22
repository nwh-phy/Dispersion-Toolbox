function tests = test_qe_auto_fit_joint
% qe_auto_fit_joint on synthetic counts with known dispersions (ZLP, 2D kinematic
% prefactor, Poisson noise): tracking on both sides of q = 0, a flagged excursion,
% seed checks and cancelling.
tests = functiontests(localfunctions);
end


function setupOnce(~)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
end


function testTracksTwoBranchesOnBothSidesOfQZero(testCase)
truth = {@(q) 700 + 6e4 * abs(q), @(q) 1700 + 2e4 * abs(q)};
[a3, e, q] = synth_joint_dataset(nq=41, E0=truth, gamma=[400 500], height=[2e3 1.5e3]);
qe = struct('intensity', a3, 'energy_meV', e, 'q_Ainv', q);
out = qe_auto_fit_joint(qe, channel(q, 0.005), n_peaks=2, signal_window=[250 2200]);

fitted = [out.channels.q_Ainv];
verifyEqual(testCase, sort(fitted), q(q ~= 0), 'AbsTol', 1e-12);   % every channel once, q = 0 never
verifyEqual(testCase, out.n_success, 40);
verifyFalse(testCase, out.cancelled);
for k = 1:2
    rows = out.peaks(out.peaks(:, 13) == k, :);
    verifyGreaterThanOrEqual(testCase, size(rows, 1), 36, sprintf('branch %d keeps most channels', k));
    verifyLessThan(testCase, max(abs(rows(:, 2) - truth{k}(rows(:, 1)))), 30, sprintf('branch %d E0', k));
    verifyEqual(testCase, out.branch_assignment.branches{k}, sortrows(rows, 1));
end
verifyTrue(testCase, all(isnan(out.peaks(:, 6:11)), 'all'));      % no intervals
verifyTrue(testCase, all(out.peaks(:, 14) < out.peaks(:, 2)));     % Drude-Lorentz apex below E0
verifyTrue(testCase, all(ismember(out.review_reason, ...
    {'failed', 'collapse', 'boundary', 'zero_amplitude', 'jump', 'poor_fit', 'low_R2', 'wide'})));
end


function testFlagsAnExcursionAndPicksUpTheBranchAfterIt(testCase)
% Two channels where the peak sits 400 meV above the smooth branch: they go to
% review as jumps, and the sweep follows the branch again after them.
smooth = @(q) 700 + 6e4 * abs(q);
excursion = @(q) smooth(q) + 400 * any(abs(q - [0.0075 0.008]) < 1e-9);
[a3, e, q] = synth_joint_dataset(nq=41, E0={excursion});
qe = struct('intensity', a3, 'energy_meV', e, 'q_Ainv', q);
out = qe_auto_fit_joint(qe, channel(q, 0.004), signal_window=[250 2200], q_range=[0.002 0.01]);

verifyEmpty(testCase, out.mirror_index);
verifyEqual(testCase, out.review(:, 1), [0.0075; 0.008], 'AbsTol', 1e-12);
verifyEqual(testCase, out.review_reason, {'jump'; 'jump'});
after = out.peaks(out.peaks(:, 1) > 0.00875, :);
verifyNumElements(testCase, after(:, 1), 3);
verifyLessThan(testCase, max(abs(after(:, 2) - smooth(after(:, 1)))), 30);
end


function testRefusesAnUnusableSeed(testCase)
[a3, e, q] = synth_joint_dataset(nq=21);
qe = struct('intensity', a3, 'energy_meV', e, 'q_Ainv', q);
verifyError(testCase, @() qe_auto_fit_joint(qe, channel(q, 0)), 'qe_auto_fit_joint:Seed');
verifyError(testCase, @() qe_auto_fit_joint(qe, channel(q, 0.003), q_range=[-0.002 0.002]), ...
    'qe_auto_fit_joint:Seed');
% The peak (E0 = 880 meV at q = 0.003) lies below a 1500-2200 meV window.
verifyError(testCase, @() qe_auto_fit_joint(qe, channel(q, 0.003), signal_window=[1500 2200]), ...
    'qe_auto_fit_joint:BadSeed');
end


function testCancelStopsTheRunAndProgressCountsChannels(testCase)
[a3, e, q] = synth_joint_dataset(nq=21);
qe = struct('intensity', a3, 'energy_meV', e, 'q_Ainv', q);
fractions = [];
out = qe_auto_fit_joint(qe, channel(q, 0.003), signal_window=[250 2200], ...
    progress_fn=@record, cancel_fn=@stop);
verifyTrue(testCase, out.cancelled);
verifyNumElements(testCase, out.channels, 3);
verifyEqual(testCase, fractions, (1:3) / 20, 'AbsTol', 1e-12);

    function record(fraction, ~)
        fractions(end + 1) = fraction;
    end

    function yes = stop()
        yes = numel(fractions) >= 3;
    end
end


function k = channel(q, value)
[~, k] = min(abs(q - value));
end
