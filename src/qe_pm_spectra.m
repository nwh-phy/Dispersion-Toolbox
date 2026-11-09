function [S, found] = qe_pm_spectra(M, q, abs_q, side)
%QE_PM_SPECTRA  Spectra at given |q| from the +q and -q channels of a q-E map.
%   S = qe_pm_spectra(M, q, abs_q) returns E x numel(abs_q): for each target
%   |q| the mean of the channel at +|q| and the channel at -|q| (the single
%   channel at q = 0). M is E x q, q the channel centres in 1/A; a channel
%   matches a target when |q - target| < dq/2 (dq = median channel spacing).
%   side = 'mean' (default), 'plus' or 'minus'. found(i) is false and S(:, i)
%   NaN when a required channel is missing.
%
%   See also: qe_prepare_stacked_spectra
arguments
    M {mustBeNumeric}
    q (1,:) double
    abs_q (1,:) double
    side (1,:) char {mustBeMember(side, {'mean', 'plus', 'minus'})} = 'mean'
end
assert(size(M, 2) == numel(q), 'qe_pm_spectra:Size', 'M must be E x numel(q)');
dq = median(diff(q));
S = nan(size(M, 1), numel(abs_q)); found = false(1, numel(abs_q));
for i = 1:numel(abs_q)
    jp = find(abs(q - abs_q(i)) < dq/2, 1); jm = find(abs(q + abs_q(i)) < dq/2, 1);
    switch side
        case 'mean', ok = ~isempty(jp) && ~isempty(jm); if ok, S(:, i) = (double(M(:, jp)) + double(M(:, jm))) / 2; end
        case 'plus', ok = ~isempty(jp); if ok, S(:, i) = double(M(:, jp)); end
        case 'minus', ok = ~isempty(jm); if ok, S(:, i) = double(M(:, jm)); end
    end
    found(i) = ok;
end
end
