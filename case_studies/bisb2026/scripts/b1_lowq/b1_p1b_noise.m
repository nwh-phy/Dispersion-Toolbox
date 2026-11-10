function [sigma, g, info] = b1_p1b_noise(F, E, ref, holdout, method, npts)
% Noise of a frame-summed spectrum, sigma^2 = g(E) * counts, from the frames
% F (E x T). ref (T x 1) is the beam reference per frame:
%   method 'first_diff' : g from first differences of the reference-corrected
%                         frames (var(diff)/2 / mean), removes slow drifts;
%   method 'common'     : g = var/mean of the reference-corrected frames;
%   method 'channel'    : reference = this channel's own frame totals (the
%                         model used up to P1; ref is ignored).
% The running median (npts points) of g never mixes holdout and non-holdout
% energies, so data inside holdout windows do not change sigma outside.
if nargin < 6, npts = 21; end
F = double(F);
switch method
    case 'channel'
        s = sum(F, 1);
    otherwise
        s = ref(:).';
end
Fc = F .* (mean(s) ./ s);
m = mean(F, 2);
switch method
    case 'first_diff'
        ratio = var(diff(Fc, 1, 2), 0, 2) / 2 ./ max(m, eps);
    otherwise
        ratio = var(Fc, 0, 2) ./ max(m, eps);
end
hm = false(size(E(:)));
for k = 1:size(holdout, 1)
    hm = hm | (E(:) >= holdout(k, 1) & E(:) <= holdout(k, 2));
end
ro = ratio; ro(hm) = NaN; rh = ratio; rh(~hm) = NaN;
g = movmedian(ro, npts, 1, 'omitnan');
gh = movmedian(rh, npts, 1, 'omitnan');
g(hm) = gh(hm);
g(~isfinite(g) | g <= 0) = 1;
sigma = sqrt(g .* max(sum(F, 2), 1));
info = struct('method', method, 'ratio', ratio, 'npts', npts);
end
