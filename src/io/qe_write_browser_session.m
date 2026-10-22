function qe_write_browser_session(path, session)
%QE_WRITE_BROWSER_SESSION Validate and atomically replace a versioned MAT session.
qe_validate_browser_session(session);
qe_check_browser_sources(session.sources);
target = char(java.io.File(char(path)).getCanonicalPath());
check_paths = {session.sources.path};
if any(strcmp(target, check_paths))
    error('qe_write_browser_session:SourceOverwrite', 'A session cannot overwrite its source data.');
end
folder = fileparts(target);
temporary = [tempname(folder) '.mat'];
cleanup = onCleanup(@() removeTemporary(temporary));
save(temporary, 'session', '-v7.3');
[ok, message] = movefile(temporary, target, 'f');
if ~ok, error('qe_write_browser_session:SaveFailed', '%s', message); end
end

function removeTemporary(path)
if isfile(path), delete(path); end
end
