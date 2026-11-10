function [F, C] = b1_p2b_hscan_summary(run_root, h, Q)
% Tables of the slot-acceptance scan of P2b (reports/p2b_hscan_20261017.md) from the runs
% of b1_p2b_hscan_run. The run folder of every h is found the way b1_p2b_fit names it:
% the hash of b1_p2b_slot_config('hscan', h).
%   F  one row per h and variant: chi2 in the window, its difference to the best row,
%      chi2/dof, the 150-250 meV holdout (chi2 per point, mean residual in sigma), the
%      number of negative tails and, for the main variant, the predicted chi2/n of the
%      left-out block.
%   C  main variant: DL apex and E0 of both branches and the width of the high branch
%      (meV, with the statistical error) at the curve points nearest to Q.
% Defaults: the five h of the report and Q = 0.001 ... 0.008 A^-1.
root = b1_simple_root();
if nargin < 1 || isempty(run_root), run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
if nargin < 2 || isempty(h), h = [0.0005 0.001 0.002 0.004 0.0066]; end
if nargin < 3, Q = [0.001 0.002 0.003 0.004 0.005 0.0065 0.008]; end
measures = {'apex_b1', 'E0_b1', 'apex_b2', 'E0_b2', 'G_b2'};

F = table(); C = table();
for k = 1:numel(h)
    tag = b1_sha256(b1_p2b_slot_config('hscan', h(k)));
    d = fullfile(run_root, 'p2b', ['fit_590_' tag(1:8)]);
    assert(isfile(fullfile(d, 'p2b_fits.csv')), 'b1_p2b_hscan_summary:NoRun', 'no P2b run for h = %g in %s', h(k), d);
    fits = readtable(fullfile(d, 'p2b_fits.csv'), 'Delimiter', ',');
    lobo = readtable(fullfile(d, 'p2b_lobo.csv'), 'Delimiter', ',');
    n = height(fits); is_main = strcmp(fits.variant, 'main');
    lobo_pred = nan(n, 1); lobo_pred(is_main) = lobo.chi2n_predicted(1);
    F = [F; table(repmat(h(k), n, 1), repmat({['fit_590_' tag(1:8)]}, n, 1), fits.variant, fits.chi2, fits.chi2dof, ...
        fits.holdout_chi2_per_point, fits.holdout_mean_residual_sigma, fits.negative_tails, lobo_pred, ...
        'VariableNames', {'h', 'run', 'variant', 'chi2', 'chi2dof', 'holdout_chi2_per_point', ...
        'holdout_mean_residual_sigma', 'negative_tails', 'lobo_chi2n_predicted'})]; %#ok<AGROW>
    curves = readtable(fullfile(d, 'p2b_curves_main.csv'));
    [~, at] = min(abs(curves.Q - Q(:).'), [], 1);
    for m = 1:numel(measures)
        C = [C; table(repmat(h(k), numel(Q), 1), repmat(measures(m), numel(Q), 1), Q(:), curves.Q(at(:)), ...
            curves.(measures{m})(at(:)), curves.([measures{m} '_err'])(at(:)), ...
            'VariableNames', {'h', 'measure', 'Q', 'Q_curve', 'value_meV', 'stat_err_meV'})]; %#ok<AGROW>
    end
end
F.dchi2 = F.chi2 - min(F.chi2); F = movevars(F, 'dchi2', 'After', 'chi2');
disp(F);
W = unstack(C(:, {'h', 'measure', 'Q', 'value_meV'}), 'value_meV', 'Q', 'VariableNamingRule', 'preserve');
disp(W);
end
