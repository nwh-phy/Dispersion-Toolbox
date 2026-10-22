function missing = qe_check_browser_sources(sources)
%QE_CHECK_BROWSER_SOURCES Reject changed sources; allow offline embedded sessions.
missing = {};
for k = 1:numel(sources)
    path = sources(k).path;
    if ~isfile(path)
        missing{end+1} = path; %#ok<AGROW>
        continue
    end
    info = dir(path);
    if info.bytes ~= sources(k).bytes || ~strcmp(qe_file_sha256(path), sources(k).sha256)
        error('qe_check_browser_sources:SourceChanged', ...
            'Source file changed since import: %s. Its contents do not match the recorded source.', path);
    end
end
end
