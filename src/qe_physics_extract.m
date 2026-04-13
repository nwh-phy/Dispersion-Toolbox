function phys = qe_physics_extract(branches, opts)
%QE_PHYSICS_EXTRACT  Self-consistent I_kin / ε_inter / Loss Function extraction.
%
%   phys = qe_physics_extract(branches)
%   phys = qe_physics_extract(branches, opts)
%
%   Implements Phase 3 of Do et al. (2025) pipeline — the "magic step":
%   extract material-intrinsic optical response from measured EELS by
%   self-consistently decoupling the EELS prefactor I_kin(q) from the
%   interband screening ε_inter(q) using only experimental data.
%
%   Theory (Do et al. 2025, Eqs. 1-3):
%     d²P/dΩdE = I_kin(q) × Im[-1/ε(q,ω)]           ... (Eq. 1)
%     ω_p(q)   = sqrt(D·q / (2π·ε₀·ε_inter(q)))     ... (Eq. 2)
%     A(q)     = I_kin(q) × 1/ε_inter(q)             ... (Eq. 3)
%
%   Pipeline:
%     1. Fit Drude weight D from low-q ω_p²(q) ∝ q regime
%     2. Compute ε_inter(q) = D·q / (2π·ε₀·ω_p²(q)) from measured ω_p
%     3. Fit Keldysh model: ε_inter(q) = 1 + ρ₀·q → extract ρ₀
%     4. Compute I_kin(q) = A(q) × ε_inter(q)
%     5. Fit piecewise power-law to I_kin in log-log space
%     6. Detect dimensional crossover q* (2D→3D)
%
%   Inputs:
%     branches — cell array of Nx12 branch matrices from qe_auto_fit
%                columns: [q E Γ R² A E_lo E_hi G_lo G_hi A_lo A_hi raw_h]
%     opts     — (optional) struct with fields:
%       .branch_index    — which branch to analyze (default: 1 = plasmon)
%       .q_min_fit       — min |q| for Drude weight fit (default: 0.01 Å⁻¹)
%       .q_max_linear    — max |q| for linear ω_p² regime (default: 0.06 Å⁻¹)
%       .use_raw_A       — use raw peak heights (col 12) vs model A (col 5)
%       .rho0_fallback   — fallback ρ₀ if fit fails (default: 25 Å)
%       .energy_unit     — 'meV' or 'eV' (default: 'meV')
%
%   Output:
%     phys — struct with fields:
%       .q              — |q| array (Å⁻¹)
%       .omega_p        — plasmon energy (meV)
%       .gamma          — damping rate (meV)
%       .A_raw          — raw peak amplitude
%       .Drude_weight   — fitted D (eV·Å⁻¹), from ω_p²∝q
%       .epsilon_inter  — ε_inter(q) from ω_p(q)
%       .rho0           — Keldysh screening length (Å)
%       .I_kin          — extracted EELS prefactor
%       .I_kin_slope_lo — low-q log-log slope (expect ≈ -3 for 2D)
%       .I_kin_slope_hi — high-q log-log slope (expect ≈ -2 for 3D)
%       .q_crossover    — dimensional crossover momentum (Å⁻¹)
%       .loss_function  — L(q) = 1/ε_inter(q)
%       .quality_factor — Q = ω_p / Γ
%
%   See also: qe_gamma_dashboard, qe_loss_map, fit_loss_function

arguments
    branches    cell
    opts        struct = struct()
end

%% Default options
if ~isfield(opts, 'branch_index'),  opts.branch_index = 1;         end
if ~isfield(opts, 'q_min_fit'),     opts.q_min_fit = 0.01;         end
if ~isfield(opts, 'q_max_linear'),  opts.q_max_linear = 0.06;      end
if ~isfield(opts, 'use_raw_A'),     opts.use_raw_A = true;         end
if ~isfield(opts, 'rho0_fallback'), opts.rho0_fallback = 25;       end
if ~isfield(opts, 'energy_unit'),   opts.energy_unit = 'meV';      end

bi = opts.branch_index;
if bi > numel(branches) || isempty(branches{bi})
    error('qe_physics_extract:noBranch', ...
        'Branch %d is empty or does not exist.', bi);
end
br = branches{bi};

%% Extract columns from branch matrix
q_signed = br(:,1);
q_abs    = abs(q_signed);
omega_p  = br(:,2);    % peak energy (meV)
gamma    = br(:,3);    % damping (meV)

% Select amplitude source
if opts.use_raw_A && size(br,2) >= 12 && ~all(isnan(br(:,12)))
    A_raw = br(:,12);
else
    A_raw = br(:,5);
end

% Filter: valid data only
valid = q_abs > 0 & omega_p > 0 & A_raw > 0 & ...
        isfinite(q_abs) & isfinite(omega_p) & isfinite(A_raw);
q_abs   = q_abs(valid);
omega_p = omega_p(valid);
gamma   = gamma(valid);
A_raw   = A_raw(valid);

% Sort by q
[q_abs, sort_idx] = sort(q_abs);
omega_p = omega_p(sort_idx);
gamma   = gamma(sort_idx);
A_raw   = A_raw(sort_idx);

% Convert energy units for internal calculation
if strcmpi(opts.energy_unit, 'meV')
    omega_eV = omega_p * 1e-3;   % meV → eV
else
    omega_eV = omega_p;
end

%% Step 1: Fit Drude weight D from ω_p²(q) ∝ q
%
%   From Eq. 2: ω_p² = D·q / (2π·ε₀·ε_inter(q))
%   In the low-q limit where ε_inter ≈ 1: ω_p² ≈ (D / 2πε₀) · q
%   So D = slope × 2π × ε₀, with slope = dω_p²/dq

linear_mask = q_abs >= opts.q_min_fit & q_abs <= opts.q_max_linear;

if sum(linear_mask) >= 3
    % Fit ω_p² = slope * q + intercept (should have near-zero intercept)
    P_drude = polyfit(q_abs(linear_mask), omega_eV(linear_mask).^2, 1);
    slope_eV2_per_Ainv = P_drude(1);
    
    % D = slope × 2π × ε₀
    % ε₀ = 8.854e-12 F/m = 8.854e-12 C²/(N·m²)
    % But in atomic units: slope is eV²·Å, convert to SI:
    %   slope_SI = slope_eV2_per_Ainv × (1.6e-19)² / (1e-10)  [J²·m / m]
    %   D_SI = slope_SI × 2π × 8.854e-12
    % For practical purposes, we store the dimensionless ratio:
    Drude_weight = slope_eV2_per_Ainv;  % eV²·Å (needs unit context)
    
    fprintf('  Drude weight fit: D = %.3f eV²·Å (from %d points)\n', ...
        Drude_weight, sum(linear_mask));
    fprintf('  Linear fit: ω_p² ≈ %.3f·q + %.4f\n', P_drude(1), P_drude(2));
else
    warning('qe_physics_extract:insufficientLinear', ...
        'Not enough points in linear regime [%.3f, %.3f] Å⁻¹ for Drude weight fit.', ...
        opts.q_min_fit, opts.q_max_linear);
    Drude_weight = NaN;
end

%% Step 2: Compute ε_inter(q) from measured ω_p(q)
%
%   ε_inter(q) = D·q / (2π·ε₀·ω_p²)
%   Using our dimensionless convention: ε_inter = (D/ω_p²) · q
%   where D is the fitted slope from Step 1.

if isfinite(Drude_weight) && Drude_weight > 0
    epsilon_inter = (Drude_weight .* q_abs) ./ (omega_eV.^2);
else
    % Fallback: use Keldysh model directly
    epsilon_inter = 1 + opts.rho0_fallback * q_abs;
    fprintf('  [fallback] Using Keldysh with ρ₀ = %.0f Å\n', opts.rho0_fallback);
end

%% Step 3: Fit Keldysh model ε_inter(q) = 1 + ρ₀·q
%
%   Linear fit: ε_inter - 1 = ρ₀ · q

keldysh_mask = q_abs >= opts.q_min_fit & q_abs <= opts.q_max_linear;
if sum(keldysh_mask) >= 3
    eps_minus_1 = epsilon_inter(keldysh_mask) - 1;
    q_kel = q_abs(keldysh_mask);
    P_kel = polyfit(q_kel, eps_minus_1, 1);
    rho0 = P_kel(1);
    if rho0 < 0
        rho0 = opts.rho0_fallback;
        fprintf('  [warning] ρ₀ fit gave negative value, using fallback %.0f Å\n', rho0);
    else
        fprintf('  Keldysh fit: ε_inter = 1 + %.1f·q (ρ₀ = %.1f Å)\n', rho0, rho0);
    end
else
    rho0 = opts.rho0_fallback;
    fprintf('  [fallback] ρ₀ = %.0f Å (insufficient data for fit)\n', rho0);
end

%% Step 4: Extract I_kin(q) = A(q) × ε_inter(q)
%
%   Do et al. Eq. 3 rearranged:
%     A(q) = I_kin(q) / ε_inter(q)
%   Therefore:
%     I_kin(q) = A(q) × ε_inter(q)
%
%   For a 2D electron gas with qd << 1: I_kin ∝ q⁻³
%   For bulk-like regime: I_kin ∝ q⁻²

I_kin = A_raw .* epsilon_inter;

%% Step 5: Piecewise power-law fit to I_kin
%
%   In log-log space: log(I_kin) = α·log(q) + const
%   α ≈ -3 for 2D, α ≈ -2 for 3D

ik_valid = I_kin > 0 & q_abs > 0;
log_q = log(q_abs(ik_valid));
log_Ik = log(I_kin(ik_valid));

% Global fit
if sum(ik_valid) >= 3
    P_global = polyfit(log_q, log_Ik, 1);
    fprintf('  I_kin global slope: q^{%.2f}\n', P_global(1));
else
    P_global = [NaN, NaN];
end

% Piecewise: split at median q
q_median = median(q_abs(ik_valid));
lo_mask = q_abs(ik_valid) <= q_median;
hi_mask = q_abs(ik_valid) > q_median;

if sum(lo_mask) >= 3
    P_lo = polyfit(log_q(lo_mask), log_Ik(lo_mask), 1);
    I_kin_slope_lo = P_lo(1);
    fprintf('  I_kin low-q slope:  q^{%.2f}  (expect ~ -3 for 2D)\n', I_kin_slope_lo);
else
    I_kin_slope_lo = NaN;
end

if sum(hi_mask) >= 3
    P_hi = polyfit(log_q(hi_mask), log_Ik(hi_mask), 1);
    I_kin_slope_hi = P_hi(1);
    fprintf('  I_kin high-q slope: q^{%.2f}  (expect ~ -2 for 3D)\n', I_kin_slope_hi);
else
    I_kin_slope_hi = NaN;
end

%% Step 6: Dimensional crossover detection
%
%   Find q* where the slope transitions from -3 to -2.
%   Simple approach: intersection of the two power-law fits.

if isfinite(I_kin_slope_lo) && isfinite(I_kin_slope_hi) && ...
   abs(I_kin_slope_lo - I_kin_slope_hi) > 0.1  % slopes are different
    % Intersection: P_lo(1)*log(q) + P_lo(2) = P_hi(1)*log(q) + P_hi(2)
    log_q_cross = (P_hi(2) - P_lo(2)) / (P_lo(1) - P_hi(1));
    q_crossover = exp(log_q_cross);
    fprintf('  Dimensional crossover: q* ≈ %.4f Å⁻¹\n', q_crossover);
else
    q_crossover = NaN;
end

%% Loss function: L(q) = 1/ε_inter(q)
%   This is the q-dependent part of Im[-1/ε₂D(ω,q)] at the plasmon peak.
loss_function = 1 ./ epsilon_inter;

%% Quality factor
quality_factor = omega_p ./ gamma;

%% Pack output
phys = struct();
phys.q              = q_abs;
phys.omega_p        = omega_p;
phys.gamma          = gamma;
phys.A_raw          = A_raw;
phys.Drude_weight   = Drude_weight;
phys.epsilon_inter  = epsilon_inter;
phys.rho0           = rho0;
phys.I_kin          = I_kin;
phys.I_kin_slope_lo = I_kin_slope_lo;
phys.I_kin_slope_hi = I_kin_slope_hi;
phys.q_crossover    = q_crossover;
phys.loss_function  = loss_function;
phys.quality_factor = quality_factor;

% Store fit info for downstream use
phys.Drude_fit      = struct('slope', Drude_weight, ...
                             'q_range', [opts.q_min_fit, opts.q_max_linear]);
phys.Keldysh_fit    = struct('rho0', rho0, ...
                             'q_range', [opts.q_min_fit, opts.q_max_linear]);
phys.I_kin_fit      = struct('slope_lo', I_kin_slope_lo, ...
                             'slope_hi', I_kin_slope_hi, ...
                             'q_crossover', q_crossover);
phys.opts           = opts;

fprintf('  Physics extraction complete: %d data points\n', numel(q_abs));
end
