function pts = qe_flatten_branch_points(branches)
%QE_FLATTEN_BRANCH_POINTS  Return all branch [q, E] points as one sorted matrix.
%
%   pts = qe_flatten_branch_points(branches) concatenates the first two columns
%   of every non-empty branch and sorts by q then energy.  It is useful for GUI
%   counters and legacy Save Pts export, where older code stored only a single
%   unbranched manual_points array.

if nargin < 1 || isempty(branches)
    pts = zeros(0, 2);
    return
end

pts = zeros(0, 2);
branches = branches(:);
for bi = 1:numel(branches)
    br = branches{bi};
    if isempty(br)
        continue
    end
    if size(br, 2) < 2
        error('qe_flatten_branch_points:InvalidBranch', ...
            'Branch %d must have at least [q, energy] columns.', bi);
    end
    pts = [pts; br(:, 1:2)]; %#ok<AGROW>
end

if ~isempty(pts)
    pts = sortrows(pts, [1 2]);
end
end
