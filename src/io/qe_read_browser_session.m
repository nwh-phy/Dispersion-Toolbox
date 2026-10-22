function [session, missing_sources] = qe_read_browser_session(path)
%QE_READ_BROWSER_SESSION Validate the complete session before touching the GUI.
loaded = load(path, 'session');
if ~isfield(loaded, 'session')
    error('qe_validate_browser_session:InvalidSession', 'File does not contain a browser session.');
end
session = loaded.session;
qe_validate_browser_session(session);
missing_sources = qe_check_browser_sources(session.sources);
end
