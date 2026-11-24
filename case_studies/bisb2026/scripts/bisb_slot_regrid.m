function [Y, cq] = bisb_slot_regrid(E, q, X, z, c0, k, gamma, beta)
%BISB_SLOT_REGRID Per-channel energy zero and energy-dependent q distortion.
%   Y = bisb_slot_regrid(E, q, X, z, c0, k, gamma, beta)
%   X is energy x detector channel x any trailing dimensions (spatial rows).
%   1. Channel j is shifted in energy so that its elastic zero z(j) (meV)
%      moves to E = 0:  Xe(E, j) = X(E + z(j), j).
%   2. The channel axis is resampled on k, channels from the direct beam c0,
%      with the detector distortion measured in the loss signal
%          c_measured = c0 + k + (gamma + beta*k) * E/1000,   E in meV,
%      a shift gamma (channels/eV) plus a magnification change beta (1/eV).
%   Linear interpolation in both steps; output points that need data outside
%   the measured energy or channel range are NaN. cq returns c_measured
%   (energy x numel(k)).
arguments
    E (:,1) double
    q (1,:) double
    X double
    z (1,:) double
    c0 (1,1) double
    k (1,:) double
    gamma (1,1) double = 0
    beta (1,1) double = 0
end
sz = size(X); n = prod(sz(3:end));
assert(sz(1) == numel(E) && sz(2) == numel(q) && numel(z) == numel(q), ...
    'bisb_slot_regrid:Size', 'X must be numel(E) x numel(q) x ...; z one value per channel');
X = reshape(double(X), sz(1), sz(2), n); Xe = nan(size(X));
for j = 1:numel(q)
    Xe(:, j, :) = reshape(interp1(E, reshape(X(:, j, :), sz(1), n), E + z(j), 'linear', NaN), sz(1), 1, n);
end
cq = c0 + k + (gamma + beta * k) .* E / 1000;
Eg = repmat(E, 1, numel(k)); Y = nan(sz(1), numel(k), n);
for r = 1:n
    G = griddedInterpolant({E, q(:)}, Xe(:, :, r), 'linear', 'none');
    Y(:, :, r) = G(Eg, cq);
end
Y = reshape(Y, [sz(1) numel(k) sz(3:end)]);
end
