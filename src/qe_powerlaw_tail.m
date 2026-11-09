function [bg, r, A] = qe_powerlaw_tail(E, y, window)
%QE_POWERLAW_TAIL  Power-law ZLP tail A*E^-r fitted in a loss window.
%   [bg, r, A] = qe_powerlaw_tail(E, y, [E1 E2]) fits log(y) = log(A) - r*log(E)
%   by least squares on E1 <= E <= E2 (y > 0) and returns the tail bg on the
%   grid E (0 where E <= 0). Subtract it to see a loss peak above the tail.
%   Check y - bg: a window that contains a loss feature (e.g. the ~130 meV
%   phonon of amorphous-carbon contamination) gives a too flat tail and a
%   clearly negative difference above the window.
%
%   See also: fit_loss_function, qe_preprocess
E = E(:); y = double(y(:));
w = E >= window(1) & E <= window(2) & y > 0;
assert(nnz(w) >= 3, 'qe_powerlaw_tail:Window', 'Fewer than three positive points in the fit window');
p = polyfit(log(E(w)), log(y(w)), 1);
r = -p(1); A = exp(p(2));
bg = zeros(size(y)); pos = E > 0; bg(pos) = exp(polyval(p, log(E(pos))));
end
