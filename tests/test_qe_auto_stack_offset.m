function tests = test_qe_auto_stack_offset
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
end


function testOffsetTracksTraceVariationNotLargeCommonBaseline(testCase)
traces = 1e7 + [
    0   0   0;
    8  10  12;
   16  20  24;
    8  10  12;
    0   0   0];

offset = qe_auto_stack_offset(traces);

verifyGreaterThan(testCase, offset, 10);
verifyLessThan(testCase, offset, 40);
end


function testOffsetUsesSelectedTraceColumnsOnly(testCase)
small_traces = [
    0  0;
    2  4;
    4  8;
    2  4;
    0  0];
unused_huge_trace = [0; 1e6; 2e6; 1e6; 0];

offset_without_huge = qe_auto_stack_offset(small_traces);
offset_with_huge = qe_auto_stack_offset([small_traces unused_huge_trace]);

verifyLessThan(testCase, offset_without_huge, 12);
verifyGreaterThan(testCase, offset_with_huge, offset_without_huge);
end


function testConstantTraceFallsBackToPositiveSpacing(testCase)
offset = qe_auto_stack_offset(ones(5, 3) * 42);

verifyGreaterThan(testCase, offset, 0);
verifyTrue(testCase, isfinite(offset));
end
