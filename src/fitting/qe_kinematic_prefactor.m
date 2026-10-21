function [Kfun, info] = qe_kinematic_prefactor(member_q, options)
%QE_KINEMATIC_PREFACTOR  Energy-dependent q-EELS prefactor averaged over a q bin.
%
%   Kfun = qe_kinematic_prefactor(member_q, beam_kV=30, form='2d')
%   K = Kfun(E)          % E = |energy loss| in meV, K in A^3 ('2d') or A^2 ('3d')
%
%   Fast electron at normal incidence on a thin film, quasi-static and
%   non-relativistic (Do et al., Nat. Commun. 16, 5801 (2025), SI Eq. S13,
%   after Rodriguez Echarri et al., PRR 2, 023096 (2020)):
%       '2d'  I_kin ~ q / (q^2 + qE^2)^2     (q d << 1)
%       '3d'  I_kin ~ 1 / (q^2 + qE^2)       (bulk-like limit)
%   with qE = omega / v. At 30 kV qE = 1.54e-3 A^-1 per eV, comparable with the
%   smallest measured q, so the factor reshapes the spectrum inside one bin.
%
%   The in-plane momentum q = sqrt(qx^2 + qy^2) is averaged over the bin: qx
%   uniform inside each member channel (width dq), qy uniform within +/- h_perp
%   (spectrometer acceptance across the dispersion axis); both are convolved
%   with a Gaussian probe angular spread sigma_probe (A^-1).
%   sigma_probe_y sets a different spread across the dispersion axis (default:
%   sigma_probe); n_probe / probe_extent set the Gaussian sampling (default 13
%   points over +/-3 sigma), n_sub the samples per channel width and q_round the
%   node merging; raise them to check quadrature convergence.
%
%   See also: qe_zlp_joint_fit

arguments
 member_q (:,1) double
 options.dq (1,1) double {mustBePositive} = 0.0005
 options.beam_kV (1,1) double {mustBePositive} = 30
 options.form char {mustBeMember(options.form,{'2d','3d'})} = '2d'
 options.h_perp (1,1) double {mustBeNonnegative} = 0
 options.sigma_probe (1,1) double {mustBeNonnegative} = 0
 options.n_sub (1,1) double {mustBePositive,mustBeInteger} = 9
 options.sigma_probe_y (1,1) double = NaN      % NaN: same as sigma_probe
 options.n_probe (1,1) double {mustBePositive,mustBeInteger} = 13
 options.probe_extent (1,1) double {mustBePositive} = 3
 options.q_round (1,1) double {mustBePositive} = 1e-6
end

gam = 1 + options.beam_kV * 1e3 / 510998.95;
v = sqrt(1 - 1 / gam ^ 2) * 299792458;
qE_per_meV = 1e-3 / 6.582119569e-16 / v * 1e-10;

ns = options.n_sub;
box = ((1:ns) - 0.5) / ns - 0.5;
sx = options.sigma_probe; sy = options.sigma_probe_y;
if isnan(sy), sy = sx; end
[zx, wzx] = local_probe(sx, options.n_probe, options.probe_extent);
[zy, wzy] = local_probe(sy, options.n_probe, options.probe_extent);
qx = member_q(:) + options.dq * box;
qx = qx(:) + sx * zx;
wx = ones(numel(member_q) * ns, 1) * wzx;
if options.h_perp > 0
    qy = options.h_perp * 2 * box(:) + sy * zy;
    wy = ones(ns, 1) * wzy;
else
    qy = sy * zy;
    wy = wzy;
end
Q = sqrt(qx(:) .^ 2 + qy(:).' .^ 2);
W = wx(:) * wy(:).';
[qn, ~, ic] = unique(round(Q(:) / options.q_round) * options.q_round);
wn = accumarray(ic, W(:));
wn = wn / sum(wn);
keep = qn > 0;
qn = qn(keep); wn = wn(keep) / sum(wn(keep));

form = options.form;
Kfun = @(E) local_eval(E, qn, wn, qE_per_meV, form);
info = struct('form', form, 'beam_kV', options.beam_kV, 'qE_per_meV', qE_per_meV, ...
    'q_nodes', qn, 'q_weights', wn, 'q_mean', sum(wn .* qn), 'h_perp', options.h_perp, ...
    'sigma_probe', options.sigma_probe, 'sigma_probe_y', sy, 'n_sub', ns, 'n_probe', options.n_probe, ...
    'probe_extent', options.probe_extent, 'member_q', member_q(:).');
end


function [z, w] = local_probe(s, n, ext)
% Gaussian probe samples (unnormalised weights); a delta when s = 0.
if s > 0
    z = linspace(-ext, ext, n);
    w = exp(-z .^ 2 / 2);
else
    z = 0; w = 1;
end
end


function K = local_eval(E, qn, wn, qE_per_meV, form)
qE2 = (qE_per_meV * abs(E(:))) .^ 2;
q = qn(:).';
switch form
    case '2d'
        K = (q ./ (q .^ 2 + qE2) .^ 2) * wn(:);
    otherwise
        K = (1 ./ (q .^ 2 + qE2)) * wn(:);
end
K = reshape(K, size(E));
end
