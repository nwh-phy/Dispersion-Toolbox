function branches_out = symmetrize_qe_branches(branches_in, opts)
%SYMMETRIZE_QE_BRANCHES Weighted +q/-q averaging for accepted branch points.
%   The output keeps the same branch fields and adds/updates .symmetrized.
%   Column 1 of .symmetrized is |q|, not signed q.

if nargin < 2 || isempty(opts)
    opts = thesis_config().symmetry;
end
branches_out = branches_in;
labels = {'low', 'high'};
for li = 1:numel(labels)
    label = labels{li};
    if ~isfield(branches_in, label)
        continue
    end
    accepted = branches_in.(label).accepted;
    branches_out.(label).symmetrized = local_symmetrize_matrix(accepted, opts);
end
end


function sym = local_symmetrize_matrix(points, opts)
if isempty(points)
    sym = zeros(0, 12);
    return
end
if size(points,2) < 12
    points(:, end+1:12) = NaN;
end
q_abs = abs(points(:,1));
tol = opts.q_tolerance_Ainv;
if tol <= 0
    tol = 1e-9;
end
q_key = round(q_abs ./ tol) .* tol;
keys = unique(q_key);
sym = nan(numel(keys), 12);
for i = 1:numel(keys)
    mask = q_key == keys(i);
    group = points(mask, :);
    weights = local_symmetry_weights(group, opts);
    weights = weights ./ max(sum(weights), opts.min_weight);
    row = local_weighted_row(group, weights);
    row(1) = mean(abs(group(:,1)), 'omitnan');
    row(4) = local_weighted_scalar(group(:,4), weights);
    sym(i, :) = row;
end
[~, order] = sort(sym(:,1), 'ascend');
sym = sym(order, :);
end


function row = local_weighted_row(group, weights)
row = nan(1, size(group, 2));
for col = 1:size(group, 2)
    row(col) = local_weighted_scalar(group(:, col), weights);
end
end


function value = local_weighted_scalar(values, weights)
valid = isfinite(values) & isfinite(weights) & weights > 0;
if ~any(valid)
    value = NaN;
else
    value = sum(values(valid) .* weights(valid)) ./ sum(weights(valid));
end
end


function weights = local_symmetry_weights(group, opts)
r2 = group(:,4);
damping = min(1, group(:,2) ./ max(group(:,3), 1));
height = group(:,12);
if all(~isfinite(height) | height <= 0)
    height = group(:,5);
end
height(~isfinite(height) | height <= 0) = 1;
height = height ./ max(height);
weights = max(r2, 0) .* max(damping, 0) .* max(height, opts.min_weight);
weights(~isfinite(weights) | weights <= 0) = opts.min_weight;
end
