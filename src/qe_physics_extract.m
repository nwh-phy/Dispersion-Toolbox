function phys = qe_physics_extract(branches, opts)
%QE_PHYSICS_EXTRACT  Extract screening and kinematic prefactor from one branch.
%
%   phys = qe_physics_extract(branches)
%   phys = qe_physics_extract(branches, opts)
%
%   This routine takes one dispersion branch, fits the quasi-2D plasmon
%   dispersion model, and uses the fitted screening factor to separate the
%   measured peak intensity A(q) into:
%
%       A(q) = I_kin(q) / epsilon_inter(q)
%
%   The key point is that epsilon_inter(q) is obtained from a global
%   quasi-2D fit of omega_p(q), rather than from a low-q linearized slope.
%   That avoids the circular "fit D from low-q, then infer epsilon_inter
%   from the same low-q assumption" logic that is too fragile for
%   quantitative analysis.
%
%   Inputs:
%     branches - cell array of Nx12 branch matrices
%     opts     - optional struct:
%       .branch_index      which branch to analyze (default: 1)
%       .q_min_fit         minimum |q| for dispersion fit (default: 0.01 Å^-1)
%       .q_max_linear      maximum |q| used for screening fit (default: 0.06 Å^-1)
%       .use_raw_A         use column 12 measured peak heights when present
%       .rho0_fallback     fallback rho0 if fit fails (default: 25 Å)
%       .energy_unit       'meV' or 'eV' for branch energies (default: 'meV')
%       .epsilon_s         substrate dielectric constant (default: 1)
%       .symmetry_average  average +q/-q before fitting (default: true)
%
%   Output:
%     phys.q                  unique |q| values (Å^-1)
%     phys.omega_p            branch energies in the input energy unit
%     phys.gamma              damping values in the input energy unit
%     phys.A_raw              measured peak intensity
%     phys.Drude_weight       fitted quasi-2D A parameter (meV^2 * Å)
%     phys.epsilon_inter      model screening epsilon_bg + rho0*q
%     phys.epsilon_inter_exp  experimental A*q/omega_p^2 estimate
%     phys.rho0               fitted screening length (Å)
%     phys.I_kin              extracted kinematic prefactor
%     phys.loss_function      screening factor 1/epsilon_inter(q)
%     phys.dispersion_fit     full result from fit_quasi2d_plasmon

arguments
    branches cell
    opts struct = struct()
end

if ~isfield(opts, 'branch_index'),     opts.branch_index = 1;        end
if ~isfield(opts, 'q_min_fit'),        opts.q_min_fit = 0.01;        end
if ~isfield(opts, 'q_max_linear'),     opts.q_max_linear = 0.06;     end
if ~isfield(opts, 'use_raw_A'),        opts.use_raw_A = true;        end
if ~isfield(opts, 'rho0_fallback'),    opts.rho0_fallback = 25;      end
if ~isfield(opts, 'energy_unit'),      opts.energy_unit = 'meV';     end
if ~isfield(opts, 'epsilon_s'),        opts.epsilon_s = 1;           end
if ~isfield(opts, 'symmetry_average'), opts.symmetry_average = true; end

bi = opts.branch_index;
if bi > numel(branches) || isempty(branches{bi})
    error('qe_physics_extract:noBranch', ...
        'Branch %d is empty or does not exist.', bi);
end
br = branches{bi};

q_signed = br(:, 1);
q_abs = abs(q_signed);
omega_in = br(:, 2);
gamma_in = br(:, 3);

if opts.use_raw_A && size(br, 2) >= 12 && ~all(isnan(br(:, 12)))
    A_raw = br(:, 12);
else
    A_raw = br(:, 5);
end

valid = q_abs > 0 & omega_in > 0 & A_raw > 0 & ...
        isfinite(q_abs) & isfinite(omega_in) & isfinite(gamma_in) & isfinite(A_raw);
q_abs = q_abs(valid);
omega_in = omega_in(valid);
gamma_in = gamma_in(valid);
A_raw = A_raw(valid);
q_signed = q_signed(valid);

if isempty(q_abs)
    error('qe_physics_extract:noValidPoints', ...
        'No valid branch points remain after filtering.');
end

[~, q_abs, omega_in, gamma_in, A_raw, n_candidates_removed] = ...
    collapse_by_signed_q(q_signed, q_abs, omega_in, gamma_in, A_raw);

if opts.symmetry_average
    [q_abs, omega_in, gamma_in, A_raw, q_counts] = average_by_abs_q(q_abs, omega_in, gamma_in, A_raw);
else
    [q_abs, sort_idx] = sort(q_abs);
    omega_in = omega_in(sort_idx);
    gamma_in = gamma_in(sort_idx);
    A_raw = A_raw(sort_idx);
    q_counts = ones(size(q_abs));
end

switch lower(opts.energy_unit)
    case 'mev'
        omega_meV = omega_in;
        gamma_out = gamma_in;
        omega_out = omega_in;
    case 'ev'
        omega_meV = omega_in * 1e3;
        gamma_out = gamma_in;
        omega_out = omega_in;
    otherwise
        error('qe_physics_extract:badEnergyUnit', ...
            'energy_unit must be ''meV'' or ''eV''.');
end

fit_mask = q_abs >= opts.q_min_fit & q_abs <= opts.q_max_linear;
if sum(fit_mask) < 3
    fit_mask = q_abs > 0;
end

weights = ones(sum(fit_mask), 1);
dispersion_fit = struct();
eps_bg = (1 + opts.epsilon_s) / 2;

try
    dispersion_fit = fit_quasi2d_plasmon( ...
        q_abs(fit_mask), omega_meV(fit_mask), ...
        confidence=weights, ...
        epsilon_s=opts.epsilon_s, ...
        rho0_init=opts.rho0_fallback);
    Drude_weight = dispersion_fit.A;
    rho0 = dispersion_fit.rho0;
catch ME
    warning('qe_physics_extract:dispersionFitFailed', ...
        'Quasi-2D dispersion fit failed: %s', ME.message);
    rho0 = opts.rho0_fallback;
    epsilon_inter_fallback = eps_bg + rho0 * q_abs;
    Drude_weight = median((omega_meV.^2) .* epsilon_inter_fallback ./ max(q_abs, eps), 'omitnan');
end

if ~isfinite(Drude_weight) || Drude_weight <= 0
    epsilon_inter_fallback = eps_bg + opts.rho0_fallback * q_abs;
    Drude_weight = median((omega_meV.^2) .* epsilon_inter_fallback ./ max(q_abs, eps), 'omitnan');
end
if ~isfinite(rho0) || rho0 <= 0
    rho0 = opts.rho0_fallback;
end

epsilon_inter_model = eps_bg + rho0 * q_abs;
epsilon_inter_exp = (Drude_weight .* q_abs) ./ max(omega_meV.^2, eps);
epsilon_inter_exp = max(epsilon_inter_exp, eps_bg);
epsilon_inter = epsilon_inter_model;

fprintf('  Screening extraction: A = %.3g meV^2*A, rho0 = %.2f A, eps_bg = %.2f\n', ...
    Drude_weight, rho0, eps_bg);
if n_candidates_removed > 0
    fprintf('  Candidate cleanup: removed %d duplicate peak candidates at fixed signed q\n', ...
        n_candidates_removed);
end
if any(q_counts > 1)
    fprintf('  Symmetry averaging: merged %d signed points into %d |q| bins\n', ...
        sum(q_counts), numel(q_counts));
end

I_kin = A_raw .* epsilon_inter;

ik_valid = I_kin > 0 & q_abs > 0;
log_q = log(q_abs(ik_valid));
log_Ik = log(I_kin(ik_valid));

P_lo = [NaN, NaN];
P_hi = [NaN, NaN];
I_kin_slope_lo = NaN;
I_kin_slope_hi = NaN;

if sum(ik_valid) >= 3
    P_global = polyfit(log_q, log_Ik, 1);
    fprintf('  I_kin global slope: q^{%.2f}\n', P_global(1));
end

q_crossover = NaN;
if sum(ik_valid) >= 6
    q_median = median(q_abs(ik_valid));
    lo_mask = q_abs(ik_valid) <= q_median;
    hi_mask = q_abs(ik_valid) > q_median;

    if sum(lo_mask) >= 3
        P_lo = polyfit(log_q(lo_mask), log_Ik(lo_mask), 1);
        I_kin_slope_lo = P_lo(1);
        fprintf('  I_kin low-q slope:  q^{%.2f}  (expect ~ -3 for 2D)\n', I_kin_slope_lo);
    end
    if sum(hi_mask) >= 3
        P_hi = polyfit(log_q(hi_mask), log_Ik(hi_mask), 1);
        I_kin_slope_hi = P_hi(1);
        fprintf('  I_kin high-q slope: q^{%.2f}  (expect ~ -2 for 3D)\n', I_kin_slope_hi);
    end

    if isfinite(I_kin_slope_lo) && isfinite(I_kin_slope_hi) && ...
            abs(I_kin_slope_lo - I_kin_slope_hi) > 0.1
        log_q_cross = (P_hi(2) - P_lo(2)) / (P_lo(1) - P_hi(1));
        q_crossover = exp(log_q_cross);
        fprintf('  Dimensional crossover: q* ~= %.4f A^-1\n', q_crossover);
    end
end

loss_function = 1 ./ epsilon_inter;
quality_factor = omega_out ./ gamma_out;

phys = struct();
phys.q = q_abs;
phys.q_counts = q_counts;
phys.omega_p = omega_out;
phys.gamma = gamma_out;
phys.A_raw = A_raw;
phys.Drude_weight = Drude_weight;
phys.Drude_weight_units = 'meV^2*Angstrom';
phys.epsilon_bg = eps_bg;
phys.epsilon_inter = epsilon_inter;
phys.epsilon_inter_exp = epsilon_inter_exp;
phys.epsilon_inter_model = epsilon_inter_model;
phys.rho0 = rho0;
phys.I_kin = I_kin;
phys.I_kin_slope_lo = I_kin_slope_lo;
phys.I_kin_slope_hi = I_kin_slope_hi;
phys.q_crossover = q_crossover;
phys.loss_function = loss_function;
phys.quality_factor = quality_factor;
phys.Drude_fit = struct('A', Drude_weight, ...
    'units', 'meV^2*Angstrom', ...
    'q_range', [min(q_abs(fit_mask)), max(q_abs(fit_mask))]);
phys.Keldysh_fit = struct('rho0', rho0, ...
    'epsilon_bg', eps_bg, ...
    'q_range', [min(q_abs(fit_mask)), max(q_abs(fit_mask))]);
phys.I_kin_fit = struct('slope_lo', I_kin_slope_lo, ...
    'slope_hi', I_kin_slope_hi, ...
    'q_crossover', q_crossover);
phys.dispersion_fit = dispersion_fit;
phys.opts = opts;

fprintf('  Physics extraction complete: %d q bins\n', numel(q_abs));
end


function [q_unique, omega_mean, gamma_mean, A_mean, counts] = average_by_abs_q(q_abs, omega_in, gamma_in, A_raw)
[q_unique, ~, group_idx] = unique(q_abs, 'sorted');
omega_mean = accumarray(group_idx, omega_in, [], @mean);
gamma_mean = accumarray(group_idx, gamma_in, [], @mean);
A_mean = accumarray(group_idx, A_raw, [], @mean);
counts = accumarray(group_idx, 1, [], @sum);
end


function [q_signed_out, q_abs_out, omega_out, gamma_out, A_out, n_removed] = ...
        collapse_by_signed_q(q_signed, q_abs, omega_in, gamma_in, A_raw)
[q_signed_out, ~, group_idx] = unique(q_signed, 'sorted');
n_groups = numel(q_signed_out);

q_abs_out = zeros(n_groups, 1);
omega_out = zeros(n_groups, 1);
gamma_out = zeros(n_groups, 1);
A_out = zeros(n_groups, 1);
n_removed = 0;

for gi = 1:n_groups
    rows = find(group_idx == gi);
    if isscalar(rows)
        best_idx = rows;
    else
        [~, local_idx] = max(A_raw(rows));
        best_idx = rows(local_idx);
        n_removed = n_removed + numel(rows) - 1;
    end

    q_abs_out(gi) = q_abs(best_idx);
    omega_out(gi) = omega_in(best_idx);
    gamma_out(gi) = gamma_in(best_idx);
    A_out(gi) = A_raw(best_idx);
end
end
