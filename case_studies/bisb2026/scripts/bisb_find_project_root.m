function project_root = bisb_find_project_root(start_dir)
%BISB_FIND_PROJECT_ROOT Find the repository root for BiSb case-study scripts.
%   The root is the nearest parent directory containing both startup.m and src/.

if nargin < 1 || strlength(string(start_dir)) == 0
    start_dir = fileparts(mfilename('fullpath'));
end

current_dir = char(string(start_dir));
while true
    if isfile(fullfile(current_dir, 'startup.m')) && ...
            isfolder(fullfile(current_dir, 'src'))
        project_root = current_dir;
        return;
    end

    parent_dir = fileparts(current_dir);
    if isempty(parent_dir) || strcmp(parent_dir, current_dir)
        error('bisb_find_project_root:NotFound', ...
            'Could not find a project root containing startup.m and src/ from %s.', ...
            start_dir);
    end
    current_dir = parent_dir;
end
end
