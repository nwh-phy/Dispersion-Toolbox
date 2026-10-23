function opHistory = qe_read_browser_history(path)
%QE_READ_BROWSER_HISTORY Read and validate GUI history without UI changes.
loaded = load(path, 'opHistory');
if ~isfield(loaded, 'opHistory')
    error('qe_read_browser_history:InvalidHistory', 'File does not contain opHistory.');
end
opHistory = loaded.opHistory;
qe_validate_browser_history(opHistory);
end
