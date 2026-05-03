function model_results = fit_thesis_branch_models(branches, cfg)
%FIT_THESIS_BRANCH_MODELS Fit conservative thesis models to symmetrized branches.
%   Low branch defaults to quasi-2D plasmon as an empirical parameterization;
%   high branch defaults to optical_constant to avoid over-interpreting it as
%   a plasmon-like branch.

if nargin < 2 || isempty(cfg)
    cfg = thesis_config();
end

model_results = struct();
model_results.low = local_fit_one(branches.low, cfg.models.low, cfg.models.epsilon_s);
model_results.high = local_fit_one(branches.high, cfg.models.high, cfg.models.epsilon_s);
end


function out = local_fit_one(branch, model_name, epsilon_s)
out = struct();
out.branch = branch.name;
out.model = model_name;
out.success = false;
out.error = '';
out.fit = struct();
points = branch.symmetrized;
if isempty(points)
    out.error = 'empty_branch';
    return
end
valid = isfinite(points(:,1)) & isfinite(points(:,2)) & points(:,1) > 0 & points(:,2) > 0;
points = points(valid, :);
if size(points, 1) < 2
    out.error = 'insufficient_points';
    return
end
weights = points(:,4);
weights(~isfinite(weights) | weights <= 0) = 1;
try
    fit = fit_dispersion_generic(points(:,1), points(:,2), ...
        'model', model_name, ...
        'confidence', weights, ...
        'epsilon_s', epsilon_s);
    out.success = true;
    out.fit = fit;
catch ME
    out.error = ME.message;
end
end
