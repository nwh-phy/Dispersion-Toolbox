function T = b1_layer_compare(run_root)
% Band averages of the 1L / 2L comparison (the numbers of reports/layer_1L_2L_20261017.md;
% the point-by-point figure is b1_layer_fig): 590 (1L), n0 (1L repeat) and 20w (2L), from
% the 2d_h5 adaptive-binning tables in <run_root>/stage5.
% For every measure (DL apex and E0 of both branches, apex of the measured loss function),
% |q| band and side (both signs, +q, -q): mean, standard error and number of the converged
% bins of each dataset, and the three differences of the means. The band edges sit between
% grid points, so no bin falls on an edge.
% Returns the table and prints the rows of both signs.
root = b1_simple_root();
if nargin < 1 || isempty(run_root), run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
ids = {'590_PL2_10w', 'n0_PL2_10w_repeat', 'no_PL2_20w_2film'}; tag = {'590', 'n0', '20w'};
bands = [0.0021 0.0031; 0.0031 0.0051; 0.0051 0.0071; 0.0071 0.0091; 0.0091 0.0125];
measures = {'apex_low', 'apex_high', 'apex_measured_LF', 'E0_low', 'E0_high'};
sides = {'both', '+q', '-q'};

D = cell(1, 3);
for s = 1:3
    D{s} = readtable(fullfile(run_root, 'stage5', sprintf('e0_stage2_2d_h5_%s.csv', ids{s})));
    assert(numel(unique(round(D{s}.q, 5))) == height(D{s}), 'b1_layer_compare:DuplicateQ', 'repeated q in %s', ids{s});
    D{s} = D{s}(strcmp(D{s}.status, 'converged'), :);
end

rows = cell(0, 16);
for m = 1:numel(measures)
    for b = 1:size(bands, 1)
        for sd = 1:numel(sides)
            st = nan(3, 3);                                   % dataset x (mean, sem, n)
            for s = 1:3
                q = D{s}.q; in = abs(q) >= bands(b, 1) & abs(q) < bands(b, 2);
                if sd == 2, in = in & q > 0; elseif sd == 3, in = in & q < 0; end
                v = D{s}.(measures{m})(in); n = numel(v);
                st(s, :) = [mean(v), std(v) / sqrt(n), n];    % std of one value is 0: no sem
                if n < 2, st(s, 2) = NaN; end
            end
            rows(end + 1, :) = [measures(m), {bands(b, 1), bands(b, 2)}, sides(sd), num2cell(reshape(st.', 1, [])), ...
                {st(2, 1) - st(1, 1), st(3, 1) - st(1, 1), st(3, 1) - st(2, 1)}]; %#ok<AGROW>
        end
    end
end
names = [{'measure', 'q_lo', 'q_hi', 'side'}, reshape([strcat('mean_', tag); strcat('sem_', tag); strcat('n_', tag)], 1, []), ...
    {'n0_minus_590', 'w20_minus_590', 'w20_minus_n0'}];
T = cell2table(rows, 'VariableNames', names);
disp(T(strcmp(T.side, 'both'), :));
end
