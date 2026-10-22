function sources = qe_browser_source_identity(dataset)
%QE_BROWSER_SOURCE_IDENTITY Bind the loaded data to its source bytes at import.
% Raw NPY metadata is part of the source identity when it is available.
paths = {char(dataset.source_path)};
[folder, name, ext] = fileparts(paths{1});
if strcmpi(ext, '.npy') && isfile(fullfile(folder, [name '.json']))
    paths{end+1} = fullfile(folder, [name '.json']);
end
sources = struct('path', {}, 'bytes', {}, 'sha256', {});
for k = 1:numel(paths)
    path = char(java.io.File(paths{k}).getCanonicalPath());
    info = dir(path);
    if isempty(info)
        error('qe_browser_source_identity:MissingSource', 'Source file is missing: %s', path);
    end
    sources(k) = struct('path', path, 'bytes', info.bytes, 'sha256', qe_file_sha256(path));
end
end
