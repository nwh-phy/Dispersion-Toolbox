function fig = qe_loss_map(qe_processed, phys, opts)
%QE_LOSS_MAP  Apply a q-dependent I_kin correction to a q-E EELS map.
%
%   fig = qe_loss_map(qe_processed, phys)
%   fig = qe_loss_map(qe_processed, phys, opts)
%
%   Divides the preprocessed EELS map by I_kin(q) to produce an
%   I_kin-corrected EELS map.
%
%   In experimental mode this is a loss-function-like visualization. It is
%   useful for comparing q-dependent spectral trends, but it should not be
%   over-interpreted as a full reconstruction of Im[-1/ε(ω,q)] away from
%   the fitted plasmon branch.
%
%   Two modes:
%     'simple'       — use I_kin ≈ C·q⁻³ (pure 2D assumption)
%     'experimental' — use experimentally extracted I_kin(q) from
%                      qe_physics_extract (self-consistent, handles
%                      the 2D→3D crossover)
%
%   CRITICAL WARNING:
%     Using 'simple' mode at large q (> ~0.05 Å⁻¹) produces ARTIFACTS:
%     I_kin actually scales as q⁻² in the 3D regime, so dividing by q⁻³
%     overcorrects, causing anomalous intensity increase. This is exactly
%     the artifact the Do et al. paper discusses in Fig. 3a vs 3e.
%
%   Inputs:
%     qe_processed — qeData struct (after preprocessing)
%     phys         — output of qe_physics_extract (or [] for simple mode)
%     opts         — struct with optional fields:
%       .mode       — 'simple' or 'experimental' (default: auto)
%       .E_range    — [E_min, E_max] meV for display (default: full)
%       .q_range    — [q_min, q_max] Å⁻¹ for display (default: full)
%       .cmap       — colormap (default: turbo)
%       .log_scale  — use log color scale (default: true)
%
%   See also: qe_physics_extract, qe_gamma_dashboard

arguments
    qe_processed struct
    phys                = []
    opts         struct = struct()
end

%% Options
if ~isfield(opts, 'mode')
    if isempty(phys)
        opts.mode = 'simple';
    else
        opts.mode = 'experimental';
    end
end
if ~isfield(opts, 'E_range'), opts.E_range = [-Inf, Inf]; end
if ~isfield(opts, 'q_range'), opts.q_range = [-Inf, Inf]; end
if ~isfield(opts, 'cmap'),    opts.cmap = turbo; end
if ~isfield(opts, 'log_scale'), opts.log_scale = true; end

%% Extract data
energy_meV = double(qe_processed.energy_meV(:));
q_Ainv     = double(qe_processed.q_Ainv(:));
intensity  = double(qe_processed.intensity);

% Apply display masks
E_mask = energy_meV >= opts.E_range(1) & energy_meV <= opts.E_range(2);
q_mask = q_Ainv >= opts.q_range(1) & q_Ainv <= opts.q_range(2);

energy_disp = energy_meV(E_mask);
q_disp      = q_Ainv(q_mask);
map         = intensity(E_mask, q_mask);

%% Build I_kin(q) correction vector
switch lower(opts.mode)
    case 'simple'
        % I_kin ≈ C · |q|⁻³
        % We don't know C, so normalize s.t. I_kin(q_ref) = 1
        q_abs = abs(q_disp);
        q_abs(q_abs < eps) = eps;  % avoid division by zero
        I_kin_q = q_abs.^(-3);
        I_kin_q = I_kin_q / max(I_kin_q);  % normalized
        mode_label = 'I_{kin}-corrected map — simple q^{-3} correction';
        
    case 'experimental'
        if isempty(phys) || ~isfield(phys, 'I_kin')
            error('qe_loss_map:noPhysics', ...
                'Experimental mode requires phys struct from qe_physics_extract.');
        end
        
        % Interpolate I_kin onto the qe display grid
        % Deduplicate phys.q (both +q and -q → same |q|)
        [q_uniq, ia] = unique(phys.q);
        Ik_uniq = phys.I_kin(ia);
        
        q_abs_disp = abs(q_disp);
        if numel(q_uniq) >= 2
            % interp1 with 'linear' returns NaN outside the data range
            I_kin_interp = interp1(q_uniq, Ik_uniq, q_abs_disp, 'linear');
            
            % For q slightly outside the fitted range, fill with nearest
            % neighbor (clamped). For q far beyond (>1.5× coverage), keep
            % NaN to prevent bright-edge artifacts (Issue #3).
            q_lo = min(q_uniq);
            q_hi = max(q_uniq);
            nan_mask = isnan(I_kin_interp);
            near_mask = nan_mask & q_abs_disp >= q_lo*0.5 & q_abs_disp <= q_hi*1.5;
            far_mask  = nan_mask & ~near_mask;
            
            % Fill near-range NaN with edge values (nearest-neighbor clamp)
            I_kin_interp(near_mask & q_abs_disp < q_lo) = Ik_uniq(1);
            I_kin_interp(near_mask & q_abs_disp > q_hi) = Ik_uniq(end);
            
            % Far-range stays NaN → loss_map will show as blank (not bright)
            I_kin_interp(far_mask) = NaN;
        else
            % Single point — use constant
            I_kin_interp = ones(size(q_abs_disp)) * Ik_uniq(1);
        end
        I_kin_interp(I_kin_interp <= 0) = eps;
        I_kin_q = I_kin_interp;
        mode_label = sprintf('I_{kin}-corrected map — experimental correction (ρ₀=%.0fÅ)', phys.rho0);
        
    otherwise
        error('qe_loss_map:badMode', 'Mode must be ''simple'' or ''experimental''.');
end

%% Divide map by I_kin(q) column-wise
%   S_corr(ω,q) = S(ω,q) / I_kin(q)
loss_map = map ./ I_kin_q(:)';

% Clip negatives
loss_map = max(loss_map, 0);

%% Create figure
fig = figure('Name', 'I_kin-Corrected Map', ...
    'NumberTitle', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 700]);

if opts.log_scale
    display_map = log10(max(loss_map, eps));
    clabel = 'log_{10}(S / I_{kin})';
else
    display_map = loss_map;
    clabel = 'S / I_{kin} (a.u.)';
end

% Auto-scale: use [5th, 95th] percentile
vals = display_map(:);
vals = vals(isfinite(vals));
if ~isempty(vals)
    clim_lo = prctile(vals, 5);
    clim_hi = prctile(vals, 95);
    if clim_hi <= clim_lo
        clim_hi = clim_lo + max(abs(clim_lo)*0.1, 0.5);  % widen to avoid error
    end
else
    clim_lo = 0; clim_hi = 1;
end

%% Plot
ax1 = subplot(1, 2, 1);
imagesc(ax1, q_disp, energy_disp, display_map);
axis(ax1, 'xy');
colormap(ax1, opts.cmap);
clim(ax1, [clim_lo, clim_hi]);
cb = colorbar(ax1);
cb.Label.String = clabel;
xlabel(ax1, 'q (Å⁻¹)');
ylabel(ax1, 'Energy (meV)');
title(ax1, mode_label, 'Interpreter', 'tex');

% Plot 2: original map for comparison
ax2 = subplot(1, 2, 2);
if opts.log_scale
    orig_display = log10(max(map, eps));
else
    orig_display = map;
end
imagesc(ax2, q_disp, energy_disp, orig_display);
axis(ax2, 'xy');
colormap(ax2, opts.cmap);
if ~isempty(vals)
    orig_vals = orig_display(:);
    orig_vals = orig_vals(isfinite(orig_vals));
    orig_lo = prctile(orig_vals, 5);
    orig_hi = prctile(orig_vals, 95);
    if orig_hi > orig_lo
        clim(ax2, [orig_lo, orig_hi]);
    end
end
cb2 = colorbar(ax2);
cb2.Label.String = clabel;
xlabel(ax2, 'q (Å⁻¹)');
ylabel(ax2, 'Energy (meV)');
title(ax2, 'Raw EELS probability d^2P/d\Omega dE');

% Overlay dispersion if available from phys
if ~isempty(phys) && isfield(phys, 'q') && isfield(phys, 'omega_p')
    hold(ax1, 'on');
    plot(ax1, phys.q, phys.omega_p, 'w.', 'MarkerSize', 8);
    plot(ax1, -phys.q, phys.omega_p, 'w.', 'MarkerSize', 8);
    hold(ax1, 'off');
    
    hold(ax2, 'on');
    plot(ax2, phys.q, phys.omega_p, 'w.', 'MarkerSize', 8);
    plot(ax2, -phys.q, phys.omega_p, 'w.', 'MarkerSize', 8);
    hold(ax2, 'off');
end

sgtitle(fig, 'Momentum-Dependent I_{kin} Correction (Do et al. 2025-inspired)', ...
    'FontWeight', 'bold');
end
