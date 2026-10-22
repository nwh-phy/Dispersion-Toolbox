function settings = qe_browser_fit_settings(snapshot)
%QE_BROWSER_FIT_SETTINGS The controls that determine a single-spectrum fit.
if strcmp(snapshot.fitMethod, 'zlp_joint')
    names = {'fitMethod', 'energyMin', 'energyMax', 'maxPeaks', 'peakModel', ...
        'kinematicFactor', 'beamKV', 'poissonWeights', 'dqOverride', 'guessText'};
else
    % Keep processing controls, excluding display and branch/export settings.
    names = fieldnames(snapshot);
    display = {'qStart', 'qEnd', 'qStep', 'yScale', 'autoY', 'yMin', 'yMax', ...
        'offset', 'stackOffsetMode', 'traceMode', 'maxShift', 'dispModel', ...
        'fitBranchTarget', 'correctionBranchTarget', 'exportMode', 'exportRatio', ...
        'selectedQIndex', 'selectedQ_Ainv', 'activeStage'};
    names = names(~ismember(names, display) & ~startsWith(names, 'waterfall') & ~startsWith(names, 'branch'));
end
settings = struct();
for k = 1:numel(names)
    settings.(names{k}) = snapshot.(names{k});
end
end
