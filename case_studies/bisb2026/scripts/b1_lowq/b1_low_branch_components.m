function T = b1_low_branch_components(folder)
%B1_LOW_BRANCH_COMPONENTS  Per-bin component table of one b1_low_branch_matlab session folder.
%   T = b1_low_branch_components(folder) reads folder/fit/bin_*.mat and folder/fit_summary.csv
%   (the session folder written by the fit and report stages) and returns one row per bin:
%   count apex of the one-B1 fit and of the low and high components of the two-B1 fit
%   (maximum of the fitted count component in 300-1800 meV, the definition of
%   low_count_apex_meV in b1_low_branch_matlab), the Drude-Lorentz E0 and widths, the low
%   component's share of the B1 counts, the held-out gain of two over one B1, and flagged
%   = low component at a bound, censored, merged with the high one or weak. Nothing is refitted.
%
%   See also: b1_low_branch_matlab
arguments
    folder (1,:) char {mustBeFolder}
end
S = readtable(fullfile(folder, 'fit_summary.csv'));
files = dir(fullfile(folder, 'fit', 'bin_*.mat'));
T = table();
for i = 1:numel(files)
    L = load(fullfile(files(i).folder, files(i).name), 'fits', 'rows');
    r2 = L.rows(2); f1 = L.fits{1}; f2 = L.fits{2};
    s = S(S.bin == r2.bin & S.n_peaks == 2, :);
    flagged = local_bool(s.low_boundary) || local_bool(s.low_count_censored) || ...
        local_bool(s.components_merge) || local_bool(s.low_component_weak);
    T = [T; table(r2.bin, r2.q_Ainv, abs(r2.q_Ainv), r2.N, ...
        local_apex(f1, 1), local_par(f1, 1, 1), local_par(f1, 1, 2), ...
        r2.low_count_apex_meV, local_apex(f2, 2), r2.low_E0_meV, r2.high_E0_meV, ...
        local_par(f2, 1, 2), local_par(f2, 2, 2), r2.low_area_fraction, s.test_gain_fraction, ...
        f1.success && f2.success, flagged, ...
        'VariableNames', {'bin', 'q_Ainv', 'abs_q', 'N', 'one_count_apex_meV', 'one_E0_meV', 'one_width_meV', ...
        'low_count_apex_meV', 'high_count_apex_meV', 'low_E0_meV', 'high_E0_meV', 'low_width_meV', ...
        'high_width_meV', 'low_area_fraction', 'test_gain_fraction', 'success', 'flagged'})]; %#ok<AGROW>
end
T = sortrows(T, {'N', 'q_Ainv'});
end

function a = local_apex(f, k)
a = NaN;
if ~f.success || size(f.peaks, 2) < k, return; end
m = f.energy_meV >= 300 & f.energy_meV <= 1800; e = f.energy_meV(m);
[~, i] = max(f.peaks(m, k)); a = e(i);
end

function v = local_par(f, k, c)
v = NaN;
if f.success && size(f.parameters, 1) >= k, v = f.parameters(k, c); end
end

function b = local_bool(x)
if iscell(x) || isstring(x), b = strcmpi(string(x), "true"); else, b = logical(x); end
b = any(b);
end
