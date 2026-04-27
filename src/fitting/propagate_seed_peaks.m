function results = propagate_seed_peaks(intensity, energy_axis, channel_axis, options)
%PROPAGATE_SEED_PEAKS  Propagate initial peak guesses across q/x channels.
%
%   results = propagate_seed_peaks(intensity, energy_axis, channel_axis, ...)
%
%   Given a "seed" channel with known peak positions, this function fits
%   that channel and then propagates the fitted peak positions as initial
%   guesses to neighboring channels, stepping outward in both directions.
%
%   This solves the key problem in phonon EELS: near the Γ-point, peaks
%   overlap and auto-detection fails. By seeding from a clear channel
%   (typically at moderate |q|) and propagating inward/outward, the
%   algorithm tracks peaks through regions where blind findpeaks would fail.
%
%   Name-value options:
%       seed_guesses    - Initial peak position guesses [meV] (required)
%       seed_idx        - Index of seed channel (default: auto from q_axis)
%       direction       - 'both' | 'positive' | 'negative' (default: 'both')
%       max_shift       - Max peak position shift per step [meV] (default: 80)
%       min_amplitude   - Stop propagating if peak amp < this fraction of seed (default: 0.1)
%       min_R2          - Stop propagating if R² drops below (default: 0.3)
%       peak_model      - Peak model name (default: 'lorentz')
%       pre_subtracted  - True if background already subtracted (default: false)
%       E_min           - Fitting window lower bound [meV] (default: 50)
%       E_max           - Fitting window upper bound [meV] (default: Inf)
%       smooth_width    - Smoothing for peak detection fallback (default: 25)
%       verbose         - Print progress (default: true)
%
%   Output:
%     results struct with fields:
%       peaks       - [N × 5] array: [channel_val, E0, width, R², amplitude]
%       fit_details - cell array of per-channel fit results
%       seed_idx    - index of the seed channel used
%       n_success   - number of successfully fitted channels
%
%   See also: fit_loss_function, peak_models

arguments
    intensity       (:,:) double        % [N_energy × N_channels]
    energy_axis     (:,1) double        % meV
    channel_axis    (:,1) double        % q (Å⁻¹) or x (nm)
    options.seed_guesses  (:,1) double  % Required: initial peak positions [meV]
    options.seed_idx      (1,1) double {mustBePositive, mustBeInteger} = 0
    options.direction     (1,:) char {mustBeMember(options.direction, ...
        {'both', 'positive', 'negative'})} = 'both'
    options.max_shift     (1,1) double {mustBePositive} = 80
    options.min_amplitude (1,1) double = 0.10
    options.min_R2        (1,1) double = 0.30
    options.peak_model    (1,:) char = 'lorentz'
    options.pre_subtracted (1,1) logical = false
    options.E_min         (1,1) double = 50
    options.E_max         (1,1) double = Inf
    options.smooth_width  (1,1) double {mustBePositive} = 25
    options.verbose       (1,1) logical = true
end

n_channels = numel(channel_axis);
n_seeds = numel(options.seed_guesses);

% Auto-detect seed index if not provided: pick the channel closest to
% the median |q| (or median x), avoiding the extremes and Γ-point
if options.seed_idx == 0
    abs_vals = abs(channel_axis);
    target = median(abs_vals(abs_vals > 0));
    [~, options.seed_idx] = min(abs(abs_vals - target));
end

seed_idx = options.seed_idx;

if options.verbose
    fprintf('Seed propagation: %d peaks from channel %d (val=%.4f)\n', ...
        n_seeds, seed_idx, channel_axis(seed_idx));
    fprintf('  Model: %s, max_shift=%.0f meV, direction=%s\n', ...
        options.peak_model, options.max_shift, options.direction);
end

%% Fit seed channel
seed_spectrum = intensity(:, seed_idx);
seed_result = fit_loss_function(energy_axis, seed_spectrum, ...
    'E_min', options.E_min, 'E_max', options.E_max, ...
    'initial_guesses', options.seed_guesses, ...
    'peak_model', options.peak_model, ...
    'pre_subtracted', options.pre_subtracted, ...
    'smooth_width', options.smooth_width, ...
    'max_peaks', n_seeds + 1);  % allow +1 for discovery

if seed_result.R_squared < options.min_R2
    warning('propagate_seed_peaks:PoorSeed', ...
        'Seed fit has low R² (%.3f). Results may be unreliable.', ...
        seed_result.R_squared);
end

% Reference amplitudes for stopping criterion
seed_amps = seed_result.amplitude;

if options.verbose
    fprintf('  Seed fit: %d peaks, R²=%.4f\n', seed_result.n_peaks, seed_result.R_squared);
    for p = 1:seed_result.n_peaks
        fprintf('    Peak %d: E0=%.0f meV, Γ=%.0f meV\n', ...
            p, seed_result.omega_p(p), seed_result.gamma(p));
    end
end

%% Initialize storage
all_peaks = [];  % [channel_val, omega_p, gamma, R², amplitude]
fit_details = cell(n_channels, 1);

% Store seed results
fit_details{seed_idx} = seed_result;
for p = 1:seed_result.n_peaks
    % [channel_val, omega_p, gamma, R², amplitude, branch_id]
    all_peaks(end+1, :) = [channel_axis(seed_idx), ...
        seed_result.omega_p(p), seed_result.gamma(p), ...
        seed_result.R_squared, seed_result.amplitude(p), p]; %#ok<AGROW>
end

%% Propagate in requested directions
directions = {};
if strcmp(options.direction, 'both') || strcmp(options.direction, 'positive')
    directions{end+1} = 'positive';
end
if strcmp(options.direction, 'both') || strcmp(options.direction, 'negative')
    directions{end+1} = 'negative';
end

n_success = 1;  % seed counts

for d = 1:numel(directions)
    dir_name = directions{d};
    if strcmp(dir_name, 'positive')
        step = 1;
        range = (seed_idx + 1):n_channels;
    else
        step = -1;
        range = (seed_idx - 1):-1:1;
    end

    % Start with seed's fitted peaks as the current guesses
    current_guesses = seed_result.omega_p(:);
    current_amps = seed_result.amplitude(:);
    active_peaks = true(numel(current_guesses), 1);  % track which peaks are still alive

    if options.verbose
        fprintf('  Propagating %s (%d channels)...\n', dir_name, numel(range));
    end

    for idx = range
        spectrum = intensity(:, idx);

        % Skip zero/NaN channels
        if all(spectrum == 0) || all(isnan(spectrum))
            continue
        end

        % Use only active peak guesses
        guesses = current_guesses(active_peaks);
        if isempty(guesses)
            if options.verbose
                fprintf('    Channel %d: all peaks died, stopping %s\n', idx, dir_name);
            end
            break
        end

        try
            result = fit_loss_function(energy_axis, spectrum, ...
                'E_min', options.E_min, 'E_max', options.E_max, ...
                'initial_guesses', guesses, ...
                'peak_model', options.peak_model, ...
                'pre_subtracted', options.pre_subtracted, ...
                'smooth_width', options.smooth_width, ...
                'max_peaks', numel(guesses) + 1);

            fit_details{idx} = result;

            % Quality check
            if result.R_squared < options.min_R2
                if options.verbose
                    fprintf('    Channel %d: R²=%.3f < %.3f, stopping %s\n', ...
                        idx, result.R_squared, options.min_R2, dir_name);
                end
                break
            end

            n_success = n_success + 1;

            % Match fitted peaks to current guesses (nearest-energy assignment)
            new_guesses = current_guesses;
            new_amps = current_amps;
            active_idx = find(active_peaks);

            for g = 1:numel(active_idx)
                gi = active_idx(g);
                old_pos = current_guesses(gi);

                % Find nearest fitted peak
                if result.n_peaks > 0
                    [min_dist, best_p] = min(abs(result.omega_p - old_pos));

                    if min_dist <= options.max_shift
                        % Accept: update guess position
                        new_guesses(gi) = result.omega_p(best_p);
                        new_amps(gi) = result.amplitude(best_p);

                        % Store the peak: include branch_id = gi
                        all_peaks(end+1, :) = [channel_axis(idx), ...
                            result.omega_p(best_p), result.gamma(best_p), ...
                            result.R_squared, result.amplitude(best_p), gi]; %#ok<AGROW>
                    else
                        % Shift too large: keep old guess but mark as suspect
                        new_guesses(gi) = old_pos;
                    end

                    % Check amplitude death
                    if new_amps(gi) < options.min_amplitude * seed_amps(min(gi, numel(seed_amps)))
                        active_peaks(gi) = false;
                        if options.verbose
                            fprintf('    Channel %d: peak at %.0f meV died (amp)\n', ...
                                idx, old_pos);
                        end
                    end
                else
                    active_peaks(gi) = false;
                end
            end

            current_guesses = new_guesses;
            current_amps = new_amps;

        catch
            % Fit failed; skip this channel but don't stop propagation
            if options.verbose
                fprintf('    Channel %d: fit exception, skipping\n', idx);
            end
        end
    end
end

if options.verbose
    fprintf('  Done: %d channels fitted, %d total peak observations\n', ...
        n_success, size(all_peaks, 1));
end

%% Build output
results = struct();
results.peaks = all_peaks;
results.fit_details = fit_details;
results.seed_idx = seed_idx;
results.n_success = n_success;
results.peak_model = options.peak_model;
results.channel_axis = channel_axis;
results.energy_axis = energy_axis;

end
