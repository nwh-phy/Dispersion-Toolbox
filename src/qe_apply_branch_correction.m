function [branches_out, correction] = qe_apply_branch_correction(branches, q_Ainv, energy_meV, varargin)
%QE_APPLY_BRANCH_CORRECTION  Replace or add one manually corrected branch point.
%
%   [branches_out, correction] = qe_apply_branch_correction(branches, q, E)
%   updates a cell array of auto-fit branches with one manually reviewed peak
%   position.  It is intended for the GUI workflow where automatic peak finding
%   is mostly trusted, but isolated q channels need manual correction.
%
%   Name-value options:
%     Branch        'auto' (default) or a positive branch index.  In 'auto'
%                   mode, the nearest same-q point is replaced; if none exists,
%                   the branch with the closest median energy is used.
%     QTolerance    Same-q tolerance in Å^-1.  Default: 1e-6.
%     PreserveExtra Keep columns 3:end from the replaced auto point.  Default:
%                   true.  Added points receive NaN in extra columns.
%
%   Branch rows are expected to have at least [q_Ainv, energy_meV].  Extra
%   columns such as Γ, R², amplitude, confidence intervals, and raw height are
%   preserved when replacing an existing point so downstream GUI panels do not
%   break, but the returned correction struct records the manual override.

p = inputParser();
p.addParameter('Branch', 'auto');
p.addParameter('QTolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('PreserveExtra', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

branch_spec = p.Results.Branch;
q_tol = p.Results.QTolerance;
preserve_extra = logical(p.Results.PreserveExtra);

if nargin < 1 || isempty(branches)
    branches_out = cell(0, 1);
else
    branches_out = branches(:);
end

validateattributes(q_Ainv, {'numeric'}, {'scalar', 'finite'});
validateattributes(energy_meV, {'numeric'}, {'scalar', 'finite'});

[target_branch, preferred_point] = chooseTargetBranch(branches_out, branch_spec, q_Ainv, energy_meV, q_tol);
branches_out = ensureBranchCount(branches_out, target_branch);

br = branches_out{target_branch};
if isempty(br)
    br = zeros(0, 2);
end
if size(br, 2) < 2
    error('qe_apply_branch_correction:InvalidBranch', ...
        'Branch %d must have at least [q, energy] columns.', target_branch);
end

if isempty(preferred_point)
    preferred_point = nearestSameQPoint(br, q_Ainv, energy_meV, q_tol);
end

correction = baseCorrection(branch_spec, target_branch, q_Ainv, energy_meV);

if ~isempty(preferred_point)
    old_row = br(preferred_point, :);
    new_row = old_row;
    new_row(1:2) = [q_Ainv, energy_meV];
    if ~preserve_extra && size(new_row, 2) > 2
        new_row(3:end) = NaN;
    end
    br(preferred_point, :) = new_row;

    correction.action = 'replace';
    correction.point_index = preferred_point;
    correction.old_q_Ainv = old_row(1);
    correction.old_energy_meV = old_row(2);
else
    n_cols = max(size(br, 2), 2);
    new_row = NaN(1, n_cols);
    new_row(1:2) = [q_Ainv, energy_meV];
    br = [br; new_row]; %#ok<AGROW>

    correction.action = 'add';
    correction.old_q_Ainv = NaN;
    correction.old_energy_meV = NaN;
end

if ~isempty(br)
    [br, sort_idx] = sortrows(br, [1 2]);
else
    sort_idx = [];
end
branches_out{target_branch} = br;

if strcmp(correction.action, 'replace') && ~isempty(sort_idx)
    correction.point_index = find(sort_idx == correction.point_index, 1, 'first');
elseif strcmp(correction.action, 'add')
    correction.point_index = find(abs(br(:,1) - q_Ainv) <= max(q_tol, eps) & ...
                                  abs(br(:,2) - energy_meV) <= eps(max(1, abs(energy_meV))), ...
                                  1, 'first');
    if isempty(correction.point_index)
        correction.point_index = size(br, 1);
    end
end
end


function [target_branch, preferred_point] = chooseTargetBranch(branches, branch_spec, q_Ainv, energy_meV, q_tol)
preferred_point = [];

if isnumeric(branch_spec)
    target_branch = round(branch_spec(1));
    if ~isfinite(target_branch) || target_branch < 1
        error('qe_apply_branch_correction:InvalidBranchSpec', ...
            'Branch must be ''auto'' or a positive branch index.');
    end
    return
end

branch_text = lower(strtrim(char(string(branch_spec))));
if startsWith(branch_text, 'branch')
    nums = regexp(branch_text, '\d+', 'match');
    if ~isempty(nums)
        target_branch = str2double(nums{1});
        return
    end
end

if ~strcmp(branch_text, 'auto')
    error('qe_apply_branch_correction:InvalidBranchSpec', ...
        'Branch must be ''auto'' or a positive branch index.');
end

% Prefer an existing point at the same q and closest to the clicked energy.
best = [Inf, Inf, Inf];  % qdiff, ediff, branch
for bi = 1:numel(branches)
    br = branches{bi};
    if isempty(br) || size(br, 2) < 2, continue; end
    qdiff = abs(br(:,1) - q_Ainv);
    same_q = find(qdiff <= q_tol);
    if isempty(same_q), continue; end
    [ediff, local_idx] = min(abs(br(same_q, 2) - energy_meV));
    candidate = [qdiff(same_q(local_idx)), ediff, bi];
    if candidate(1) < best(1) || ...
            (candidate(1) == best(1) && candidate(2) < best(2))
        best = candidate;
        preferred_point = same_q(local_idx);
    end
end

if isfinite(best(3))
    target_branch = best(3);
    return
end

% If no same-q point exists, add to the branch whose energy scale is closest.
best_energy = [Inf, Inf];  % energy diff, branch
for bi = 1:numel(branches)
    br = branches{bi};
    if isempty(br) || size(br, 2) < 2, continue; end
    e_med = median(br(:,2), 'omitnan');
    ediff = abs(e_med - energy_meV);
    if isfinite(ediff) && ediff < best_energy(1)
        best_energy = [ediff, bi];
    end
end

if isfinite(best_energy(2))
    target_branch = best_energy(2);
else
    target_branch = 1;
end
end


function point_idx = nearestSameQPoint(br, q_Ainv, energy_meV, q_tol)
point_idx = [];
if isempty(br)
    return
end
qdiff = abs(br(:,1) - q_Ainv);
same_q = find(qdiff <= q_tol);
if isempty(same_q)
    return
end
[~, local_idx] = min(abs(br(same_q, 2) - energy_meV));
point_idx = same_q(local_idx);
end


function branches = ensureBranchCount(branches, n)
if numel(branches) < n
    branches{n, 1} = [];
end
end


function correction = baseCorrection(branch_spec, target_branch, q_Ainv, energy_meV)
correction = struct();
correction.action = '';
correction.branch_spec = branch_spec;
correction.branch_index = target_branch;
correction.point_index = NaN;
correction.old_q_Ainv = NaN;
correction.old_energy_meV = NaN;
correction.new_q_Ainv = q_Ainv;
correction.new_energy_meV = energy_meV;
correction.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
correction.source = 'manual_correction';
end
