function indices = qe_selected_branch_indices(n_branches, selection)
%QE_SELECTED_BRANCH_INDICES Convert a GUI branch selection into indices.
%   INDICES = QE_SELECTED_BRANCH_INDICES(N_BRANCHES, SELECTION) returns
%   1:N_BRANCHES for "All" and the selected branch number for "Branch N".

if nargin < 1 || isempty(n_branches)
    n_branches = 0;
end
if nargin < 2 || isempty(selection)
    selection = 'All';
end

n_branches = max(0, floor(double(n_branches)));
selection = strtrim(char(string(selection)));

if n_branches == 0
    indices = [];
    return
end

if isempty(selection) || strcmpi(selection, 'All')
    indices = 1:n_branches;
    return
end

token = regexp(selection, 'Branch\s+(\d+)', 'tokens', 'once');
if isempty(token)
    indices = [];
    return
end

idx = str2double(token{1});
if isfinite(idx) && idx >= 1 && idx <= n_branches && idx == floor(idx)
    indices = idx;
else
    indices = [];
end
end
