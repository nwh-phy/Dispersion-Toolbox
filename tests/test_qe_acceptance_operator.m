function tests = test_qe_acceptance_operator
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.Qm = [0:1e-5:2e-4, 2.25e-4:2.5e-5:2e-3, 2.1e-3:1e-4:1e-2, 1.05e-2:5e-4:4e-2].';
end


function testLinearRadialFunctionIsExact(testCase)
% hat interpolation integrates C(Q) = Q exactly; compare with direct sampling
Qm = testCase.TestData.Qm; dq = 0.0005; h = 0.0005;
ax = [-0.0003; 0.0001]; wx = [2; 1]; ay = [-0.0002; 0.0003]; wy = [0.4; 0.6];
op = qe_acceptance_operator(0.001, ax, wx, [1; 1], ay, wy, Qm, dq=dq, h=h, n_sx=7, n_sy=5);
G = squeeze(op.W(1, 1, :)).' * Qm;
sx = 0.001 + dq * (((1:7) - 0.5) / 7 - 0.5); sy = 2 * h * (((1:5) - 0.5) / 5 - 0.5);
ref = 0;
for n = 1:2, for k = 1:2
    ref = ref + wx(n) * wy(k) * mean(sqrt((sx(:) - ax(n)) .^ 2 + (sy(:).' - ay(k)) .^ 2), 'all');
end, end
ref = ref * dq * 2 * h;
verifyEqual(testCase, G, ref, 'RelTol', 1e-12);
verifyEqual(testCase, sum(op.W(:)), dq * 2 * h * sum(wx) * sum(wy), 'RelTol', 1e-12);
end


function testSlicesSplitTheWeights(testCase)
Qm = testCase.TestData.Qm;
ax = [-0.001; -0.0009; 0.0002]; wx = [1; 2; 3]; rl = [1; 1; 2];
op1 = qe_acceptance_operator(0.0015, ax, wx, rl, 0, 1, Qm);
op0 = qe_acceptance_operator(0.0015, ax, wx, ones(3, 1), 0, 1, Qm);
verifyEqual(testCase, squeeze(sum(op1.W, 2)), squeeze(op0.W), 'AbsTol', 1e-18);
verifyEqual(testCase, sum(op1.W(1, 1, :)) / sum(op1.W(1, 2, :)), 3 / 3, 'RelTol', 1e-12);
end


function testDetailedBalanceAndFiniteZeroLoss(testCase)
intr = qe_lowq_joint_forward('intrinsic');
th = struct('E0sq', [800 790 770 740 720; 1550 1600 1650 1680 1600] .^ 2, 'G', [380 500 580; 700 1300 1750]);
w = -200:4:200; Q = [0 1e-4 1e-3 4e-3].';
C = qe_lowq_joint_forward('C', Q, w, th, intr);
verifyTrue(testCase, all(isfinite(C(:))));
verifyEqual(testCase, C(1, :, :), zeros(1, numel(w), 6), 'AbsTol', 0);  % Q = 0 row
i = find(w == 100); j = find(w == -100); a = C(:, j, :); b = C(:, i, :); nz = b ~= 0;
verifyEqual(testCase, a(nz) ./ b(nz), exp(-100 / intr.kT) * ones(nnz(nz), 1), 'RelTol', 1e-10);
verifyEqual(testCase, a(~nz), zeros(nnz(~nz), 1));
end


function testSliceEqualsIntegratedWhenKernelsAgree(testCase)
[ops, R, opt] = smallSetup(testCase, 'same');
th = refTheta(); intr = qe_lowq_joint_forward('intrinsic');
opt.variant = 'slice'; cs = qe_lowq_joint_forward('eval', qe_lowq_joint_forward('build', ops, R, opt), th, intr);
opt.variant = 'integrated'; ci = qe_lowq_joint_forward('eval', qe_lowq_joint_forward('build', ops, R, opt), th, intr);
opt.variant = 'exit'; cx = qe_lowq_joint_forward('eval', qe_lowq_joint_forward('build', ops, R, opt), th, intr);
verifyEqual(testCase, cs, ci, 'RelTol', 1e-9, 'AbsTol', 1e-12 * max(abs(cs(:))));
verifyEqual(testCase, cs, cx, 'RelTol', 1e-9, 'AbsTol', 1e-12 * max(abs(cs(:))));
end


function testSlicesWithShiftedKernelsMoveTheSignal(testCase)
% two slices whose kernels sit at -8 and +8 meV: the slice variant is the weighted
% sum of the two shifted responses; the integrated variant uses the mixed kernel
[ops, R, opt] = smallSetup(testCase, 'shifted');
th = refTheta(); intr = qe_lowq_joint_forward('intrinsic');
opt.variant = 'slice'; cs = qe_lowq_joint_forward('eval', qe_lowq_joint_forward('build', ops, R, opt), th, intr);
% direct: per slice, sum_m W C, then convolve with that slice's kernel
C = qe_lowq_joint_forward('C', ops{1}.Qm, opt.omega, th, intr); n = numel(opt.omega) + numel(opt.eps) - 1;
ref = zeros(n, 1);
for r = 1:2
    G = squeeze(ops{1}.W(1, r, :)).' * squeeze(C(:, :, 1));
    ref = ref + conv(G(:), R{1}(r, :).');
end
verifyEqual(testCase, cs(:, 1, 1, 1), ref, 'RelTol', 1e-8, 'AbsTol', 1e-10 * max(abs(ref)));
end


function testNarrowLossLineAppearsAtItsEnergy(testCase)
% a narrow, weakly dispersing line at 500 meV is seen at 500 meV plus the kernel offset
[ops, R, opt] = smallSetup(testCase, 'shifted');
th = struct('E0sq', repmat(500 ^ 2, 2, 5), 'G', repmat(4, 2, 3)); intr = qe_lowq_joint_forward('intrinsic');
opt.variant = 'exit'; opt.exit_slice = 2;          % exit kernel sits at +8 meV
c = qe_lowq_joint_forward('eval', qe_lowq_joint_forward('build', ops, R, opt), th, intr);
E = qe_lowq_joint_forward('build', ops, R, opt).E_out(:); y = c(:, 1, 1, 1); m = E > 300;
[~, i] = max(y .* m); verifyEqual(testCase, E(i), 508, 'AbsTol', 2);
end


function [ops, R, opt] = smallSetup(testCase, kind)
Qm = testCase.TestData.Qm;
ax = [-0.0002; -0.0001; 0.0001; 0.0002]; wx = [1; 2; 2; 1]; rl = [1; 1; 2; 2];
op = qe_acceptance_operator(0.001, ax, wx, rl, [-0.0003; 0; 0.0003], [0.25; 0.5; 0.25], Qm);
opt = struct('omega', -400:4:2400, 'eps', -40:4:40, 'slice_weight', {{[3 3]}}, 'exit_slice', 1);
e = opt.eps; g = @(c) exp(-(e - c) .^ 2 / (2 * 8 ^ 2)); k1 = g(0); k2 = g(0);
if strcmp(kind, 'shifted'), k1 = g(-8); k2 = g(8); end
R = {[k1 / sum(k1); k2 / sum(k2)]}; ops = {op};
end


function th = refTheta()
th = struct('E0sq', [800 790 770 740 720; 1550 1600 1650 1680 1600] .^ 2, 'G', [380 500 580; 700 1300 1750]);
end
