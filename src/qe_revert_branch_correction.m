function branches_out = qe_revert_branch_correction(branches, correction)
%QE_REVERT_BRANCH_CORRECTION  Undo one qe_apply_branch_correction operation.
%
%   branches_out = qe_revert_branch_correction(branches, correction) reverts a
%   single correction struct returned by qe_apply_branch_correction.  Replace
%   corrections restore the old [q, energy] values; add corrections remove the
%   inserted row.  Extra columns are preserved for replace corrections.

if nargin < 1 || isempty(branches)
    branches_out = cell(0, 1);
else
    branches_out = branches(:);
end
if nargin < 2 || isempty(correction)
    return
end

required = {'action', 'branch_index', 'new_q_Ainv', 'new_energy_meV'};
for k = 1:numel(required)
    if ~isfield(correction, required{k})
        error('qe_revert_branch_correction:InvalidCorrection', ...
            'Correction is missing required field "%s".', required{k});
    end
end

bi = correction.branch_index;
if ~isfinite(bi) || bi < 1 || bi > numel(branches_out)
    return
end

br = branches_out{bi};
if isempty(br)
    return
end

row_idx = locateCorrectedRow(br, correction);
if isempty(row_idx)
    return
end

switch lower(char(correction.action))
    case 'replace'
        if ~isfield(correction, 'old_q_Ainv') || ~isfield(correction, 'old_energy_meV')
            error('qe_revert_branch_correction:InvalidReplaceCorrection', ...
                'Replace correction must contain old_q_Ainv and old_energy_meV.');
        end
        br(row_idx, 1:2) = [correction.old_q_Ainv, correction.old_energy_meV];
        branches_out{bi} = sortrows(br, [1 2]);
    case 'add'
        br(row_idx, :) = [];
        branches_out{bi} = br;
        branches_out = trimTrailingEmptyBranches(branches_out);
    otherwise
        error('qe_revert_branch_correction:UnknownAction', ...
            'Unknown correction action: %s.', char(correction.action));
end
end


function row_idx = locateCorrectedRow(br, correction)
row_idx = [];
if isfield(correction, 'point_index') && isfinite(correction.point_index)
    idx = round(correction.point_index);
    if idx >= 1 && idx <= size(br, 1) && rowMatchesCorrection(br(idx, :), correction)
        row_idx = idx;
        return
    end
end

matches = find(abs(br(:,1) - correction.new_q_Ainv) <= toleranceFor(correction.new_q_Ainv) & ...
               abs(br(:,2) - correction.new_energy_meV) <= toleranceFor(correction.new_energy_meV));
if ~isempty(matches)
    row_idx = matches(1);
end
end


function tf = rowMatchesCorrection(row, correction)
tf = abs(row(1) - correction.new_q_Ainv) <= toleranceFor(correction.new_q_Ainv) && ...
     abs(row(2) - correction.new_energy_meV) <= toleranceFor(correction.new_energy_meV);
end


function tol = toleranceFor(value)
tol = max(1e-9, 10 * eps(max(1, abs(value))));
end


function branches = trimTrailingEmptyBranches(branches)
while ~isempty(branches) && isempty(branches{end})
    branches(end) = [];
end
end
