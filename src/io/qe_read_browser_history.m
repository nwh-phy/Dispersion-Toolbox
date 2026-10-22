function opHistory = qe_read_browser_history(path)
%QE_READ_BROWSER_HISTORY Read and structurally validate GUI history without UI changes.
% Fields added after the original snapshot format remain optional. UI value
% constraints are checked by the GUI's transactional restoration.
loaded = load(path, 'opHistory');
if ~isfield(loaded, 'opHistory') || ~iscell(loaded.opHistory) || ...
        (~isempty(loaded.opHistory) && ~isvector(loaded.opHistory))
    error('qe_read_browser_history:InvalidHistory', 'File does not contain a valid opHistory cell vector.');
end
opHistory = loaded.opHistory;
required = {'qStart', 'qEnd', 'qStep', 'dqOverride', 'energyMin', 'energyMax', ...
    'zlpMin', 'zlpMax', 'refQMin', 'refQMax', 'smoothMode', 'smoothWidth', ...
    'sgOrder', 'sgFrameLen', 'yScale', 'autoY', 'yMin', 'yMax', 'offset', ...
    'traceMode', 'baseLambda', 'areaNorm', 'areaNormMin', 'areaNormMax', ...
    'bgSub', 'bgMethod', 'denoise', 'denoiseMethod', 'denoiseSigma', ...
    'deconv', 'deconvIter', 'viewMode'};
groups = {{'waterfallMin', 'waterfallMax', 'waterfallNorm', 'waterfallNormMin', ...
    'waterfallNormMax', 'waterfallGain', 'waterfallResidual'}, ...
    {'bgWin1Lo', 'bgWin1Hi', 'bgDual', 'bgWin2Lo', 'bgWin2Hi', 'bgIter'}, ...
    {'branch1Min', 'branch1Max', 'branch2Min', 'branch2Max', 'branch3Min', 'branch3Max'}};
for k = 1:numel(opHistory)
    entry = opHistory{k};
    if ~isstruct(entry) || ~isscalar(entry) || ...
            ~all(isfield(entry, {'snapshot', 'label', 'timestamp'}))
        error('qe_read_browser_history:InvalidHistory', 'History entry %d is invalid.', k);
    end
    snap = entry.snapshot;
    if ~isstruct(snap) || ~isscalar(snap) || ~all(isfield(snap, required))
        error('qe_read_browser_history:InvalidSnapshot', ...
            'History entry %d has an incomplete parameter snapshot.', k);
    end
    for g = 1:numel(groups)
        fields = groups{g};
        if isfield(snap, fields{1}) && ~all(isfield(snap, fields))
            error('qe_read_browser_history:InvalidSnapshot', ...
                'History entry %d has an incomplete %s group.', k, fields{1});
        end
    end
end
end
