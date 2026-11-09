function tests = test_qe_transport_frames
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testNativeBinsWithoutDriftReturnTheChannels(testCase)
dq = 0.0005; Q = (-4:4) * dq;
F = reshape(1:3*2*numel(Q), 3, 2, numel(Q));
[frames, err] = qe_transport_frames(F, Q, dq, Q, dq * ones(size(Q)), [0 0]);
verifyClass(testCase, frames, 'single');
verifyEqual(testCase, double(frames), F);
verifyLessThan(testCase, err, 1e-12);
end


function testWholeChannelDriftMovesCountsToTheNeighbour(testCase)
dq = 0.0005; Q = (-4:4) * dq;
F = zeros(1, 2, numel(Q)); F(1, :, 6) = 7;      % counts at Q = +dq in both frames
% frame 2: beam moved by +dq, so physical q = 0 now sits on native channel Q = +dq
[frames, err] = qe_transport_frames(F, Q, dq, 0, dq, [0 dq]);
verifyEqual(testCase, reshape(double(frames), 1, []), [0 7], 'AbsTol', 1e-12);
verifyLessThan(testCase, err, 1e-12);
end


function testHalfChannelDriftAveragesTheTwoNeighbours(testCase)
rng(1);
dq = 0.0005; Q = (-6:6) * dq;
F = rand(4, 2, numel(Q));
centers = (-3:3) * dq;
[frames, err] = qe_transport_frames(F, Q, dq, centers, dq * ones(size(centers)), [0 dq/2]);
for b = 1:numel(centers)
    j = find(abs(Q - centers(b)) < dq/4);
    verifyEqual(testCase, double(frames(:, 1, b)), F(:, 1, j), 'RelTol', 1e-6);
    verifyEqual(testCase, double(frames(:, 2, b)), (F(:, 2, j) + F(:, 2, j+1)) / 2, 'RelTol', 1e-6);
end
verifyLessThan(testCase, err, 1e-12);
end


function testExpectedCountsReportTheWeightError(testCase)
dq = 0.0005; Q = (-2:2) * dq;
F = ones(1, 1, numel(Q));
[~, err] = qe_transport_frames(F, Q, dq, 0, 3 * dq, 0, 2);   % a 3-channel bin declared as 2
verifyEqual(testCase, err, 1, 'AbsTol', 1e-12);
end
