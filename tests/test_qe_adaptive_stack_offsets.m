function tests = test_qe_adaptive_stack_offsets
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testAdaptiveOffsetsPreventCrossingWhenCenterTraceIsTall(testCase)
traces = [
    10   0   0;
     0  10   0;
    10   0   0];
min_gap = 1e-4;

offsets = qe_adaptive_stack_offsets(traces, min_gap);
shifted = traces + reshape(offsets, 1, []);

verifyEqual(testCase, offsets(1), 0, 'AbsTol', 1e-12);
verifyGreaterThan(testCase, offsets(2), 9.999);
verifyGreaterThan(testCase, offsets(3), offsets(2) + 9.999);
verifyGreaterThanOrEqual(testCase, min(shifted(:, 2) - shifted(:, 1)), min_gap * 0.99);
verifyGreaterThanOrEqual(testCase, min(shifted(:, 3) - shifted(:, 2)), min_gap * 0.99);
end


function testAdaptiveOffsetsCheckAllEarlierTracesNotOnlyNeighbor(testCase)
traces = [
    100   0   40;
      0   0   40;
      0 100   40];
min_gap = 0.5;

offsets = qe_adaptive_stack_offsets(traces, min_gap);
shifted = traces + reshape(offsets, 1, []);

for col = 2:size(shifted, 2)
    previous_envelope = max(shifted(:, 1:col-1), [], 2);
    verifyGreaterThanOrEqual(testCase, min(shifted(:, col) - previous_envelope), min_gap * 0.99);
end
end


function testConstantTracesUseRequestedMinimumGap(testCase)
traces = ones(4, 3);
min_gap = 1e-4;

offsets = qe_adaptive_stack_offsets(traces, min_gap);

verifyEqual(testCase, offsets, [0 min_gap 2 * min_gap], 'AbsTol', 1e-12);
end
