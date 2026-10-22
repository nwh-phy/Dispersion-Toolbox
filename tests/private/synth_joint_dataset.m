function [a3, e, q] = synth_joint_dataset(options)
%SYNTH_JOINT_DATASET  Synthetic 10w q-E counts for the ZLP joint-fit tests.
%
%   [a3, e, q] = synth_joint_dataset(E0={@(q) 700 + 6e4 * abs(q)})
%
%   nq channels of width dq around q = 0 (the middle channel). Each channel is a
%   ZLP, strongest at q = 0, plus one Drude-Lorentz peak per entry of E0 (a
%   function of q) with width gamma and the 2D kinematic prefactor, convolved
%   with the ZLP and scaled to height counts at its maximum, on a flat 20 counts.
%   Gaussian-approximated Poisson noise from a fixed seed.

arguments
    options.nq (1,1) double {mustBeInteger, mustBePositive} = 101
    options.dq (1,1) double {mustBePositive} = 0.0005
    options.E0 cell = {@(q) 700 + 6e4 * abs(q)}
    options.gamma (1,:) double = 400
    options.height (1,:) double = 2e3
    options.seed (1,1) double = 20261022
end

e = (-300:4:2400).';
nq = options.nq;
q0 = (nq + 1) / 2;
q = ((1:nq) - q0) * options.dq;
rng(options.seed, 'twister');
a3 = zeros(numel(e), nq);
pos = e > 0;
for j = 1:nq
    zlp = 1e6 * (1 + 0.2 * exp(-(j - q0) ^ 2 / 8)) * (1 + (e / 18) .^ 2 / 1.5) .^ (-1.5);
    kernel = zlp / sum(zlp);
    K = qe_kinematic_prefactor(abs(q(j)), dq=options.dq, form='2d', beam_kV=30);
    mu = zlp;
    for p = 1:numel(options.E0)
        E0 = options.E0{p}(q(j));
        G = options.gamma(p);
        peak = zeros(size(e));
        peak(pos) = e(pos) * G ./ ((e(pos) .^ 2 - E0 ^ 2) .^ 2 + e(pos) .^ 2 * G ^ 2) .* K(e(pos));
        peak = conv(peak, kernel(e >= -120 & e <= 120), 'same');
        mu = mu + peak / max(peak) * options.height(p);
    end
    mu = mu + 20;
    a3(:, j) = mu + sqrt(mu) .* randn(size(mu));
end
end
