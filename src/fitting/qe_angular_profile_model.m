function out = qe_angular_profile_model(mode, x, dq, p)
%QE_ANGULAR_PROFILE_MODEL  Two split-Pearson-VII angular profile of the elastic core.
%
%   d = qe_angular_profile_model('density', Q, dq, p)     counts per channel width at the points Q
%   m = qe_angular_profile_model('channel', q, dq, p)     counts in channels centred at q (width dq)
%   c = qe_angular_profile_model('components', Q, dq, p)  [narrow broad] densities at the points Q
%
%   p = [log A1, mu1/dq, log(A2/A1), (mu2 - mu1)/dq, log(g1m/dq), log(g1p/dq), log m1,
%        log(g2m/dq), log(g2p/dq), log m2]
%   Each component is A (1 + ((x - mu)/g)^2)^(-m) with g = g_minus for x < mu and
%   g_plus otherwise; A1 is the narrow-component density (counts per channel width)
%   at its mode mu1. g is a scale: the half width at half maximum is
%   g sqrt(2^(1/m) - 1). 'channel' averages 21 points per channel, as used for the
%   fits in b1_p2a_inputs (P2a step 1), so fitted parameters can be reused here.
%   Q, q in A^-1, dq the channel width (A^-1).
%
%   See also: b1_p2a_inputs, qe_acceptance_operator

switch mode
    case 'density'
        out = local_density(x, dq, p);
    case 'components'
        out = [local_spv(x, dq, [p(1) p(2) p(5:7)]), local_spv(x, dq, [p(1) + p(3) p(2) + p(4) p(8:10)])];
    case 'channel'
        ns = 21; xs = ((1:ns) - 0.5) / ns - 0.5;
        out = mean(local_density(x(:) + dq * xs, dq, p), 2);
    otherwise
        error('qe_angular_profile_model:Mode', 'Unknown mode %s', mode);
end
end

function d = local_density(Q, dq, p)
d = local_spv(Q, dq, [p(1) p(2) p(5:7)]) + local_spv(Q, dq, [p(1) + p(3) p(2) + p(4) p(8:10)]);
end

function d = local_spv(Q, dq, p)
mu = p(2) * dq; g = exp(p(3)) * dq * ones(size(Q)); g(Q > mu) = exp(p(4)) * dq;
d = exp(p(1)) * (1 + ((Q - mu) ./ g) .^ 2) .^ (-exp(p(5)));
end
