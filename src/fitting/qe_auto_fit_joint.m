function out = qe_auto_fit_joint(qe, seed_index, options)
%QE_AUTO_FIT_JOINT  ZLP joint fit of every q channel, tracked from a seed channel.
%
%   out = qe_auto_fit_joint(qe, seed_index, n_peaks=2, signal_window=[250 1800])
%
%   qe holds unprocessed counts: intensity (energy x q, ZLP included), energy_meV
%   and q_Ainv. Every channel in q_range except q = 0 is fitted once with
%   qe_fit_spectrum_joint. The seed channel gets n_starts starts, the first from
%   initial_E0 when that has n_peaks values. Four sweeps follow, each in order of
%   |q|: outward and inward from the seed on its side, and outward and inward from
%   the channel nearest -q of the seed on the other side, whose own fit starts from
%   the seed. Within a sweep a channel starts from the parameters of the last
%   channel whose peaks the sweep followed (one start).
%
%   Peaks are numbered by E0 within a channel, and peak k is branch k. A peak is
%   flagged 'failed' (no converged fit), 'collapse' (two peaks within one energy
%   step), 'boundary' (a parameter at its bound), 'zero_amplitude', 'jump' (E0 more
%   than max_shift from both the last E0 its sweep followed and the median of the
%   last three; each channel since the last one adds max_shift to that allowance),
%   'poor_fit' (chi2_red above max_chi2_ratio x that of the seed), 'low_R2' (R^2 of
%   the channel below min_R2) or 'wide' (Gamma/E0 above max_gamma_ratio). A sweep
%   follows the peaks without the first six flags. Those six, or chi2_red above
%   max(1.5, 1.3 x that of the last channel it followed in full), make the channel
%   fit again with retry_starts starts, the warm start among them, and the lower
%   cost wins; a converged channel with R^2 below min_R2 is noise-dominated and
%   not fitted again.
%
%   out.peaks (unflagged) and out.review (flagged) have one row per peak in the
%   qe_auto_fit layout, [q, E0, Gamma, R2, A, six NaN (no intervals), height, k,
%   apex]: the model E0 is the dispersion energy, height the peak's maximum in
%   counts, k the branch id, apex the loss-function maximum for comparison.
%   out.review_reason holds the flag of each review row, out.branch_assignment
%   qe_assign_peak_branches_by_windows of out.peaks with one branch per peak
%   number, out.channels one entry per fitted channel. progress_fn(fraction,
%   message) is called after each channel and cancel_fn() before it; when
%   cancel_fn returns true the run stops and out.cancelled is true.
%
%   See also: qe_fit_spectrum_joint, qe_auto_fit, qe_assign_peak_branches_by_windows

arguments
    qe (1,1) struct
    seed_index (1,1) double {mustBeInteger, mustBePositive}
    options.n_peaks (1,1) double {mustBeInteger, mustBePositive} = 1
    options.peak_model char {mustBeMember(options.peak_model, {'lorentz', 'lorentz_symmetric'})} = 'lorentz'
    options.signal_window (1,2) double = [300 1800]
    options.kinematic (1,1) logical = true
    options.beam_kV (1,1) double {mustBePositive} = 30
    options.weights char {mustBeMember(options.weights, {'poisson', 'none'})} = 'poisson'
    options.q_range (1,2) double = [-Inf Inf]
    options.initial_E0 (1,:) double = []
    options.max_shift (1,1) double {mustBePositive} = 80
    options.n_starts (1,1) double {mustBeInteger, mustBePositive} = 8
    options.retry_starts (1,1) double {mustBeInteger, mustBePositive} = 4
    options.max_chi2_ratio (1,1) double {mustBePositive} = 2
    options.min_R2 (1,1) double = 0.3
    options.max_gamma_ratio (1,1) double {mustBePositive} = 2
    options.progress_fn = []
    options.cancel_fn = []
end

E = double(qe.energy_meV(:));
q = double(qe.q_Ainv(:)).';
assert(seed_index <= numel(q), 'qe_auto_fit_joint:Seed', 'The seed channel is outside the q axis.');
dq = median(abs(diff(q)), 'omitnan');
q_range = sort(options.q_range);
fitted = q >= q_range(1) - dq / 2 & q <= q_range(2) + dq / 2 & abs(q) > dq / 2;
q0 = q(seed_index);
if abs(q0) <= dq / 2
    error('qe_auto_fit_joint:Seed', 'The seed channel is q = 0; select a channel with q ~= 0.');
end
if ~fitted(seed_index)
    error('qe_auto_fit_joint:Seed', 'The seed channel q = %.4f is outside the fitted q range.', q0);
end
fit_args = {'n_peaks', options.n_peaks, 'peak_model', options.peak_model, ...
    'signal_window', options.signal_window, 'kinematic', options.kinematic, ...
    'beam_kV', options.beam_kV, 'dq', dq, 'weights', options.weights};

% Sweeps in order of |q|: seed side outward and inward, then the mirror side.
side = byAbsQ(find(fitted & sign(q) == sign(q0)), q);
p = find(side == seed_index);
sweeps = {side(p + 1:end), fliplr(side(1:p - 1))};
other = byAbsQ(find(fitted & sign(q) == -sign(q0)), q);
mirror_index = [];
mirror_sweeps = {};
if ~isempty(other)
    [~, m] = min(abs(q(other) + q0));
    mirror_index = other(m);
    mirror_sweeps = {other(m + 1:end), fliplr(other(1:m - 1))};
end
run = struct('qe', qe, 'E', E, 'q', q, 'fit_args', {fit_args}, 'options', options, ...
    'n_total', nnz(fitted), 'done', 0, 'chi2_limit', Inf);

[channels, seed_track] = fitChannel(run, seed_index, []);
bad = ~cellfun(@isempty, channels.flags);
if any(bad)
    error('qe_auto_fit_joint:BadSeed', ['The fit of the seed channel q = %.4f is not usable (%s). ' ...
        'Select another channel or enter Guesses.'], q0, strjoin(unique(channels.flags(bad)), ', '));
end
run.chi2_limit = options.max_chi2_ratio * channels.chi2_red;
run.done = 1;
report(run, q0);
cancelled = false;
for s = 1:numel(sweeps)
    [channels, run, cancelled] = sweep(run, sweeps{s}, seed_track, channels);
    if cancelled, break; end
end
if ~cancelled && ~isempty(mirror_index)
    [channels, run, cancelled, mirror_track] = sweep(run, mirror_index, seed_track, channels);
    for s = 1:numel(mirror_sweeps)
        if cancelled, break; end
        [channels, run, cancelled] = sweep(run, mirror_sweeps{s}, mirror_track, channels);
    end
end

[peaks, review, reasons] = peakRows(channels, options.n_peaks);
np = options.n_peaks;
specs = struct('name', arrayfun(@(k) sprintf('Branch %d', k), 1:np, 'UniformOutput', false), ...
    'energy_window_meV', {sort(options.signal_window)}, 'enabled', true, 'branch_index', num2cell(1:np));
out = struct();
out.peaks = peaks;
out.review = review;
out.review_reason = reasons;
out.branch_assignment = qe_assign_peak_branches_by_windows(peaks, specs(:), ...
    struct('q_skip_Ainv', 0, 'min_R2', -Inf, 'max_gamma_ratio', Inf, 'score_column', 12));
out.channels = channels;
out.seed_index = seed_index;
out.mirror_index = mirror_index;
out.n_success = nnz([channels.success]);
out.cancelled = cancelled;
out.settings = rmfield(options, {'progress_fn', 'cancel_fn'});
out.settings.dq = dq;
end


function idx = byAbsQ(idx, q)
[~, order] = sort(abs(q(idx)));
idx = idx(order);
end


function [channels, run, cancelled, track] = sweep(run, path, track, channels)
cancelled = false;
for k = path
    if ~isempty(run.options.cancel_fn) && run.options.cancel_fn()
        cancelled = true;
        return
    end
    [channels(end + 1), track] = fitChannel(run, k, track); %#ok<AGROW>
    run.done = run.done + 1;
    report(run, run.q(k));
end
end


function report(run, qk)
if ~isempty(run.options.progress_fn)
    run.options.progress_fn(run.done / run.n_total, ...
        sprintf('ZLP joint auto fit: %d of %d channels (q = %.4f)', run.done, run.n_total, qk));
end
end


function [ch, track] = fitChannel(run, k, track)
% Seed (track empty): multistart. Otherwise one warm start, and a multistart
% that includes it when a peak loses the track or chi2 rises.
o = run.options;
Y = double(run.qe.intensity(:, k));
[retried, improved] = deal(false);
if isempty(track)
    r = fitOnce(run, k, Y, o.n_starts, [], o.initial_E0);
else
    r = fitOnce(run, k, Y, 1, track.u, track.history(:, end).');
    lost = ismember(flagsOf(r, track, run), trackingFlags());
    noisy = r.success && r.R_squared < o.min_R2;   % more starts do not lift a noise-dominated channel
    retried = ~noisy && (any(lost) || r.chi2_red > max(1.5, 1.3 * track.chi2));
    if retried
        again = fitOnce(run, k, Y, o.retry_starts, track.u, track.history(:, end).');
        improved = again.success && (~r.success || again.joint.cost < r.joint.cost);
        if improved, r = again; end
    end
end
flags = flagsOf(r, track, run);
follow = ~ismember(flags, trackingFlags());
if isempty(track)
    if all(follow)
        track = struct('history', [nan(o.n_peaks, 2), r.omega_p(:)], 'skipped', zeros(o.n_peaks, 1), ...
            'u', r.joint.u, 'chi2', r.chi2_red);
    end
else
    track.history(follow, :) = [track.history(follow, 2:end), r.omega_p(follow)];
    track.skipped = (track.skipped + 1) .* ~follow;
    if all(follow)
        track.u = r.joint.u;
        track.chi2 = r.chi2_red;
    end
end
height = nan(1, o.n_peaks);
if r.success
    height = cellfun(@max, r.peak_curves);
end
ch = struct('index', k, 'q_Ainv', run.q(k), 'success', r.success, 'status', char(string(r.status)), ...
    'retried', retried, 'improved', improved, ...
    'E0', r.omega_p(:).', 'gamma', r.gamma(:).', 'amplitude', r.amplitude(:).', ...
    'apex', r.apex_energy_meV(:).', 'height', height, 'R_squared', r.R_squared, ...
    'chi2_red', r.chi2_red, 'flags', {flags(:).'});
end


function r = fitOnce(run, k, Y, n_starts, warm_u, E0)
r = qe_fit_spectrum_joint(run.E, Y, run.q(k), run.fit_args{:}, 'n_starts', n_starts, ...
    'warm_u', warm_u, 'initial_E0', E0);
end


function names = trackingFlags()
% Flags that mean the fit lost the peak; low_R2 and wide only keep it out of the branch.
names = {'failed', 'collapse', 'boundary', 'zero_amplitude', 'jump', 'poor_fit'};
end


function flags = flagsOf(r, track, run)
o = run.options;
flags = repmat({''}, o.n_peaks, 1);
if ~r.success
    flags(:) = {'failed'};
    return
end
if strcmp(r.status, 'collapse_sampling_scale')
    flags(:) = {'collapse'};
    return
end
for p = 1:o.n_peaks
    if any(r.joint.peak_boundary(p, :))
        flags{p} = 'boundary';
    elseif r.joint.zero_amplitude(p)
        flags{p} = 'zero_amplitude';
    elseif ~isempty(track) && jumped(r.omega_p(p), track, p, o.max_shift)
        flags{p} = 'jump';
    elseif r.chi2_red > run.chi2_limit
        flags{p} = 'poor_fit';
    elseif r.R_squared < o.min_R2
        flags{p} = 'low_R2';
    elseif r.gamma_ratio(p) > o.max_gamma_ratio
        flags{p} = 'wide';
    end
end
end


function yes = jumped(E0, track, p, max_shift)
% Far from the last followed E0 and from the median of the last three, so that
% neither one outlying channel nor a steep branch reads as a jump.
allowed = max_shift * (1 + track.skipped(p));
yes = abs(E0 - track.history(p, end)) > allowed && ...
    abs(E0 - median(track.history(p, :), 'omitnan')) > allowed;
end


function [peaks, review, reasons] = peakRows(channels, np)
peaks = zeros(0, 14);
review = zeros(0, 14);
reasons = cell(0, 1);
for c = channels(:).'
    for p = 1:np
        row = [c.q_Ainv, c.E0(p), c.gamma(p), c.R_squared, c.amplitude(p), nan(1, 6), ...
            c.height(p), p, c.apex(p)];
        if isempty(c.flags{p})
            peaks(end + 1, :) = row; %#ok<AGROW>
        elseif isfinite(c.E0(p))
            review(end + 1, :) = row; %#ok<AGROW>
            reasons{end + 1, 1} = c.flags{p}; %#ok<AGROW>
        end
    end
end
[review, order] = sortrows(review, [1 13]);
reasons = reasons(order);
end
