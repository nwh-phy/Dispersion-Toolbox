function qe_validate_browser_session(session)
%QE_VALIDATE_BROWSER_SESSION Validate version 1 before any live GUI changes.
check(isstruct(session) && isscalar(session), 'Session must be a scalar struct.');
require(session, {'format', 'version', 'saved_at', 'sources', 'data', ...
    'snapshot', 'history', 'analysis', 'pending_fit', 'display_fit', 'data_hashes'});
check(isequal(session.format, 'qe_browser_session') && isequal(session.version, 1), ...
    'Unsupported session format or version.');
check(isText(session.saved_at), 'Invalid save timestamp.');
check(isstruct(session.sources) && ~isempty(session.sources) && isvector(session.sources), 'Missing source identity.');
require(session.sources, {'path', 'bytes', 'sha256'});
for k = 1:numel(session.sources)
    s = session.sources(k);
    check(isText(s.path) && isscalar(s.bytes) && isnumeric(s.bytes) && ...
        isfinite(s.bytes) && s.bytes >= 0 && isHash(s.sha256), 'Invalid source identity.');
end
d = session.data;
require(d, {'dataset', 'eq3dQE', 'physicalQE', 'activeStage'});
require(d.dataset, {'source_path', 'source_kind', 'label', 'dq_Ainv', 'qe'});
validateQE(d.physicalQE);
validateQE(d.eq3dQE);
validateQE(d.dataset.qe);
check(isequaln(d.dataset.qe, d.eq3dQE), 'Dataset and eq3D stage disagree.');
check(isText(d.dataset.source_path) && isText(d.dataset.source_kind) && isText(d.dataset.label), 'Invalid dataset metadata.');
check(isequal(d.dataset.dq_Ainv, d.physicalQE.dq_Ainv) && strcmp(d.activeStage, 'eq3d'), ...
    'Session supports the projected eq3D stage and a consistent dq.');
require(session.data_hashes, {'physical', 'eq3d'});
check(strcmp(qe_data_fingerprint(d.physicalQE), session.data_hashes.physical) && ...
    strcmp(qe_data_fingerprint(d.eq3dQE), session.data_hashes.eq3d), 'Embedded counts or calibrated axes do not match their fingerprint.');
qe_validate_browser_history(session.history);
qe_validate_browser_history({struct('label', 'Session', 'timestamp', '', 'snapshot', session.snapshot)});
snap = session.snapshot;
require(snap, {'selectedQIndex', 'selectedQ_Ainv', 'fitMethod', 'peakModel', 'guessText'});
nq = numel(d.physicalQE.q_Ainv);
check(isIndex(snap.selectedQIndex, nq) && ...
    isequal(snap.selectedQ_Ainv, d.physicalQE.q_Ainv(snap.selectedQIndex)) && ...
    isequal(snap.dqOverride, d.physicalQE.dq_Ainv), 'Snapshot selection or calibration disagrees with the data.');
a = session.analysis;
require(a, {'manual_points', 'manual_branches', 'branchCorrectionLog', ...
    'fitResults', 'autoFitResults', 'acceptedFits'});
validatePoints(a.manual_points);
check(iscell(a.manual_branches) && iscell(a.fitResults) && iscell(a.acceptedFits), 'Invalid branch or fit archive.');
for k = 1:numel(a.manual_branches), validatePoints(a.manual_branches{k}); end
if ~isempty(a.manual_branches)
    check(isequaln(a.manual_points, qe_flatten_branch_points(a.manual_branches)), 'Flattened points disagree with branches.');
end
check(isempty(a.branchCorrectionLog) || isstruct(a.branchCorrectionLog), 'Invalid corrections.');
check(isempty(a.autoFitResults) || isstruct(a.autoFitResults), 'Invalid auto-fit results.');
for k = 1:numel(a.acceptedFits)
    record = a.acceptedFits{k};
    validateFit(record);
    require(record, {'accepted_at', 'accepted_peak_indices', 'accepted_energy_meV', 'target'});
    check(isequal(record.accepted_peak_indices, find(record.keep)) && ...
        isequal(record.accepted_energy_meV, fitEnergies(record.result, record.keep)) && ...
        isText(record.accepted_at) && isText(record.target), 'Invalid accepted-fit audit record.');
end
pf = session.pending_fit;
if ~isempty(pf)
    check(isequaln(pf, session.display_fit), 'Pending fit and displayed fit disagree.');
end
fits = {session.pending_fit, session.display_fit};
for k = 1:numel(fits)
    pf = fits{k};
    if isempty(pf), continue; end
    validateFit(pf);
    check(isIndex(pf.q_index, nq) && pf.q_index == snap.selectedQIndex && ...
        isequal(pf.q_Ainv, d.physicalQE.q_Ainv(pf.q_index)), 'Pending fit belongs to another channel.');
    check(isequaln(pf.provenance.settings, qe_browser_fit_settings(snap)), 'Pending fit settings are stale.');
    check(strcmp(pf.provenance.data_sha256, session.data_hashes.physical), 'Pending fit belongs to different counts or calibration.');
end
validatePlain(session);
end

function validateQE(qe)
require(qe, {'intensity', 'energy_meV', 'q_Ainv', 'q_abs_Ainv', 'q_channel', 'q_zero_index', 'dq_Ainv'});
check(isnumeric(qe.intensity) && isreal(qe.intensity) && ismatrix(qe.intensity) && ~isempty(qe.intensity), 'Invalid q-E counts.');
[ne, nq] = size(qe.intensity);
check(isnumeric(qe.energy_meV) && numel(qe.energy_meV) == ne && all(isfinite(qe.energy_meV(:))) && ...
    all(diff(qe.energy_meV(:)) > 0), 'Invalid energy axis.');
check(isnumeric(qe.q_Ainv) && isrow(qe.q_Ainv) && numel(qe.q_Ainv) == nq && ...
    isnumeric(qe.q_channel) && numel(qe.q_channel) == nq && all(isfinite(qe.q_channel(:))) && ...
    isIndex(qe.q_zero_index, nq) && isscalar(qe.dq_Ainv) && isfinite(qe.dq_Ainv) && qe.dq_Ainv > 0, 'Invalid q axis.');
check(isequal(qe.q_Ainv, ((1:nq) - qe.q_zero_index) * qe.dq_Ainv) && ...
    isequal(qe.q_abs_Ainv, abs(qe.q_Ainv)), 'Inconsistent momentum calibration.');
end

function validatePoints(points)
check(isnumeric(points) && isreal(points) && (isempty(points) || ...
    (ismatrix(points) && size(points, 2) >= 2 && all(isfinite(points(:, 1:2)), 'all'))), 'Invalid curated points.');
end

function validateFit(pf)
require(pf, {'q_Ainv', 'q_index', 'abs_q', 'result', 'keep', 'display_keep', ...
    'fit_input', 'display_scale', 'provenance'});
r = pf.result;
require(r, {'n_peaks', 'omega_p', 'gamma', 'energy_fit', 'curve_fit', 'peak_curves', 'R_squared'});
check(isIndex(r.n_peaks, inf) && numel(r.omega_p) == r.n_peaks && ...
    numel(r.gamma) == r.n_peaks && iscell(r.peak_curves) && numel(r.peak_curves) == r.n_peaks && ...
    numel(pf.keep) == r.n_peaks && numel(pf.display_keep) == r.n_peaks && ...
    islogical(pf.keep) && islogical(pf.display_keep), 'Invalid fit peaks.');
check(isnumeric(r.energy_fit) && isvector(r.energy_fit) && numel(r.energy_fit) == numel(r.curve_fit) && ...
    all(cellfun(@numel, r.peak_curves) == numel(r.energy_fit)), 'Invalid saved fit curves.');
check(isnumeric(pf.display_scale) && isscalar(pf.display_scale) && ...
    isfinite(pf.display_scale) && pf.display_scale > 0, 'Invalid display scaling.');
require(pf.provenance, {'settings', 'data_sha256', 'snapshot', 'sources', 'signal_window_meV', ...
    'energy_kind', 'input_kind', 'fitted_at'});
check(isHash(pf.provenance.data_sha256) && numel(pf.provenance.energy_kind) == r.n_peaks && ...
    isequaln(pf.provenance.settings, qe_browser_fit_settings(pf.provenance.snapshot)), 'Invalid fit provenance.');
energies = fitEnergies(r, true(r.n_peaks, 1));
expected_display = energies >= pf.provenance.signal_window_meV(1);
expected_keep = expected_display;
if isfield(r, 'peak_valid'), expected_keep = expected_keep & r.peak_valid(:); end
check(isequal(pf.display_keep, expected_display) && isequal(pf.keep, expected_keep), 'Saved fit acceptance mask disagrees with peak validity.');
if ~isempty(pf.fit_input)
    require(pf.fit_input, {'E', 'Y'});
    check(isnumeric(pf.fit_input.E) && numel(pf.fit_input.E) == numel(pf.fit_input.Y), 'Invalid saved fit input.');
end
end

function energy = fitEnergies(result, keep)
energy = result.omega_p(:);
if isfield(result, 'apex_energy_meV')
    apex = result.apex_energy_meV(:);
    energy(isfinite(apex)) = apex(isfinite(apex));
end
energy = energy(keep);
end

function validatePlain(value)
% A session contains values, never GUI handles, callbacks, or mutable objects.
if isstruct(value)
    names = fieldnames(value);
    for k = 1:numel(value)
        for n = 1:numel(names), validatePlain(value(k).(names{n})); end
    end
elseif iscell(value)
    for k = 1:numel(value), validatePlain(value{k}); end
elseif istable(value)
    % Branch assignment summaries/rejections are MATLAB value tables.
    names = value.Properties.VariableNames;
    for k = 1:numel(names), validatePlain(value.(names{k})); end
    validatePlain(value.Properties.UserData);
else
    check(isnumeric(value) || islogical(value) || ischar(value) || isstring(value), 'Session contains a nonportable object.');
end
end

function require(value, names)
check(isstruct(value) && ~isempty(value) && all(isfield(value, names)), 'Session is missing required fields.');
end

function yes = isIndex(value, n)
yes = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 1 && value <= n && value == fix(value);
end

function yes = isHash(value)
yes = ischar(value) && isrow(value) && ~isempty(regexp(value, '^[a-f0-9]{64}$', 'once'));
end

function yes = isText(value)
yes = (ischar(value) && (isrow(value) || isempty(value))) || (isstring(value) && isscalar(value) && ~ismissing(value));
end

function check(condition, message)
if ~condition, error('qe_validate_browser_session:InvalidSession', '%s', message); end
end
