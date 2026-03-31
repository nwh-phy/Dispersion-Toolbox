function fig = qe_gamma_dashboard(all_peaks, branches)
%QE_GAMMA_DASHBOARD  Linewidth, amplitude & quality factor analysis.
%
%   fig = qe_gamma_dashboard(all_peaks, branches)
%
%   Creates a 3×2 figure showing:
%     (1,1) Γ(q)         — linewidth vs momentum
%     (1,2) Γ(E)         — linewidth vs peak energy
%     (2,1) Amplitude(q)  — model amplitude vs momentum
%     (2,2) Amplitude(E)  — model amplitude vs energy
%     (3,1) Raw A(q)      — log-log raw peak height with power-law fit
%     (3,2) Q-factor(q)   — plasmon quality factor E₀/Γ
%
%   Inputs:
%     all_peaks  — Nx12 matrix [q E Γ R² A E_lo E_hi G_lo G_hi A_lo A_hi raw_h]
%     branches   — cell array of per-branch subsets of all_peaks
%
%   See also: interactive_qe_browser, qe_auto_fit

arguments
    all_peaks   double
    branches    cell
end

n_branches = numel(branches);
branch_colors = get_branch_colors();

% CI column checks
has_gamma_ci = size(all_peaks, 2) >= 9;
has_amp_ci = size(all_peaks, 2) >= 11;

fig = figure('Name', 'Linewidth Γ, Amplitude & EELS Prefactor', ...
    'NumberTitle', 'off', 'Color', 'w', ...
    'Position', [120 40 1100 1000]);

% --- Subplot 1: Γ(q) ---
ax1 = subplot(3, 2, 1, 'Parent', fig);
hold(ax1, 'on');
for b = 1:n_branches
    br = branches{b};
    if isempty(br), continue; end
    col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
    lbl = sprintf('B%d: Γ̄=%.0f meV', b, mean(br(:,3)));
    plot_branch_property(ax1, br(:,1), br(:,3), br, col, lbl, ...
        has_gamma_ci, 8, 9);
end
hold(ax1, 'off');
xlabel(ax1, 'q (1/Å)'); ylabel(ax1, 'Γ (meV)');
title(ax1, 'Linewidth Γ vs Momentum');
legend(ax1, 'Location', 'best'); grid(ax1, 'on');

% --- Subplot 2: Γ(E) ---
ax2 = subplot(3, 2, 2, 'Parent', fig);
hold(ax2, 'on');
for b = 1:n_branches
    br = branches{b};
    if isempty(br), continue; end
    col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
    lbl = sprintf('B%d: Γ̄=%.0f meV', b, mean(br(:,3)));
    plot_branch_property(ax2, br(:,2), br(:,3), br, col, lbl, ...
        has_gamma_ci, 8, 9);
end
hold(ax2, 'off');
xlabel(ax2, 'ω_p (meV)'); ylabel(ax2, 'Γ (meV)');
title(ax2, 'Linewidth Γ vs Peak Energy');
legend(ax2, 'Location', 'best'); grid(ax2, 'on');

% --- Subplot 3: Amplitude(q) ---
ax3 = subplot(3, 2, 3, 'Parent', fig);
hold(ax3, 'on');
for b = 1:n_branches
    br = branches{b};
    if isempty(br), continue; end
    col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
    lbl = sprintf('B%d', b);
    plot_branch_property(ax3, br(:,1), br(:,5), br, col, lbl, ...
        has_amp_ci, 10, 11);
end
hold(ax3, 'off');
xlabel(ax3, 'q (1/Å)'); ylabel(ax3, 'Amplitude (a.u.)');
title(ax3, 'Peak Amplitude vs Momentum');
legend(ax3, 'Location', 'best'); grid(ax3, 'on');

% --- Subplot 4: Amplitude(E) ---
ax4 = subplot(3, 2, 4, 'Parent', fig);
hold(ax4, 'on');
for b = 1:n_branches
    br = branches{b};
    if isempty(br), continue; end
    col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
    lbl = sprintf('B%d', b);
    plot_branch_property(ax4, br(:,2), br(:,5), br, col, lbl, ...
        has_amp_ci, 10, 11);
end
hold(ax4, 'off');
xlabel(ax4, 'ω_p (meV)'); ylabel(ax4, 'Amplitude (a.u.)');
title(ax4, 'Peak Amplitude vs Energy');
legend(ax4, 'Location', 'best'); grid(ax4, 'on');

% --- Subplot 5: Raw Peak Height A(q) log-log with power law ---
ax5 = subplot(3, 2, 5, 'Parent', fig);
hold(ax5, 'on');
for b = 1:n_branches
    br = branches{b};
    if isempty(br), continue; end
    col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
    plot_raw_amplitude_loglog(ax5, br, b, col);
end
hold(ax5, 'off');
xlabel(ax5, 'q (Å^{-1})'); ylabel(ax5, 'Peak Height (a.u.)');
title(ax5, 'Raw Peak Height A(q) — log-log (no areanorm)');
legend(ax5, 'Location', 'best', 'FontSize', 7); grid(ax5, 'on');

% --- Subplot 6: Quality factor Q = E₀/Γ vs q ---
ax6 = subplot(3, 2, 6, 'Parent', fig);
hold(ax6, 'on');
for b = 1:n_branches
    br = branches{b};
    if isempty(br), continue; end
    col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
    Q_factor = br(:,2) ./ br(:,3);
    scatter(ax6, br(:,1), Q_factor, 30, col, 'filled', ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', sprintf('B%d: Q̄=%.1f', b, mean(Q_factor)));
end
hold(ax6, 'off');
xlabel(ax6, 'q (Å^{-1})'); ylabel(ax6, 'Q = ω_p / Γ');
title(ax6, 'Plasmon Quality Factor vs q');
legend(ax6, 'Location', 'best'); grid(ax6, 'on');

sgtitle(fig, sprintf('Γ, Amplitude & Quality Factor (%d branches, %d peaks)', ...
    n_branches, size(all_peaks, 1)), 'FontSize', 14);

end


%% ═══════════ Local helpers ═══════════

function colors = get_branch_colors()
    colors = [0.1 0.5 0.9; 0.9 0.3 0.1; 0.2 0.7 0.3; 0.6 0.2 0.8; 0.9 0.6 0.1];
end


function plot_branch_property(ax, x, y, br, col, lbl, has_ci, ci_lo_col, ci_hi_col)
    % Generic branch scatter with optional error bars from CI columns
    if has_ci && size(br, 2) >= ci_hi_col && ~all(isnan(br(:, ci_lo_col)))
        y_err = (br(:, ci_hi_col) - br(:, ci_lo_col)) / 2;
        y_err(isnan(y_err)) = 0;
        errorbar(ax, x, y, y_err, ...
            'o', 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
            'DisplayName', lbl);
    else
        scatter(ax, x, y, 40, col, 'filled', ...
            'MarkerFaceAlpha', 0.7, 'DisplayName', lbl);
    end
end


function plot_raw_amplitude_loglog(ax, br, b, col)
    q_abs = abs(br(:,1));
    if size(br, 2) >= 12 && ~all(isnan(br(:,12)))
        A_raw = br(:,12);
        src_label = 'raw';
    else
        A_raw = br(:,5);
        src_label = 'model';
    end
    valid = q_abs > 0 & A_raw > 0 & isfinite(A_raw);
    if sum(valid) < 2, return; end

    loglog(ax, q_abs(valid), A_raw(valid), 'o', ...
        'Color', col, 'MarkerFaceColor', col, ...
        'MarkerSize', 5, 'DisplayName', sprintf('B%d (%s)', b, src_label));

    log_q = log(q_abs(valid));
    log_A = log(A_raw(valid));
    P = polyfit(log_q, log_A, 1);
    q_fit_line = linspace(min(q_abs(valid)), max(q_abs(valid)), 100);
    loglog(ax, q_fit_line, exp(P(2))*q_fit_line.^P(1), '--', ...
        'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('B%d: q^{%.1f}', b, P(1)));

    % Piecewise power-law (low/high q)
    q_median = median(q_abs(valid));
    low_mask = q_abs(valid) <= q_median;
    high_mask = q_abs(valid) > q_median;
    if sum(low_mask) >= 3 && sum(high_mask) >= 3
        P_lo = polyfit(log_q(low_mask), log_A(low_mask), 1);
        P_hi = polyfit(log_q(high_mask), log_A(high_mask), 1);
        q_lo = linspace(min(q_abs(valid)), q_median, 50);
        q_hi = linspace(q_median, max(q_abs(valid)), 50);
        loglog(ax, q_lo, exp(P_lo(2))*q_lo.^P_lo(1), ':', ...
            'Color', col*0.6, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('lo: q^{%.1f}', P_lo(1)));
        loglog(ax, q_hi, exp(P_hi(2))*q_hi.^P_hi(1), '-.', ...
            'Color', col*0.6+0.4, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('hi: q^{%.1f}', P_hi(1)));
    end
end
