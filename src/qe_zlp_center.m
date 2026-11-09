function [c, fw] = qe_zlp_center(y, E)
%QE_ZLP_CENTER  ZLP centre and FWHM of one spectrum.
%   [c, fw] = qe_zlp_center(y, E): centre c (meV) of the log-parabola through
%   the maximum within |E| < 40 meV and its two neighbours; fw = FWHM (meV) by
%   linear interpolation of the half-maximum crossings (NaN if not found).
%   E must be uniformly spaced. This is the ZLP-centre rule of P1b/P2a
%   (b1_p1b_diag, b1_p2a_inputs) and of qe_elastic_line_model, kept in one place.
%
%   See also: qe_elastic_line_model, qe_zlp_integer_align
y = y(:); dE = E(2) - E(1); [~, i0] = max(y .* (abs(E) < 40));
ly = log(max(y(i0-1:i0+1), 1)); d = (ly(3) - ly(1)) / (2 * (2 * ly(2) - ly(1) - ly(3))); c = E(i0) + d * dE;
h = y(i0) / 2; l1 = find(y(1:i0) < h, 1, 'last'); r1 = i0 - 1 + find(y(i0:end) < h, 1, 'first'); fw = NaN;
if ~isempty(l1) && ~isempty(r1)
    fw = interp1(y([r1-1 r1]), E([r1-1 r1]), h) - interp1(y([l1 l1+1]), E([l1 l1+1]), h);
end
end
