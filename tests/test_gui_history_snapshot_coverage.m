function tests = test_gui_history_snapshot_coverage
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.project_root = project_root;
testCase.TestData.gui_src = fileread(fullfile(project_root, 'src', 'interactive_qe_browser.m'));
testCase.TestData.ui_src = fileread(fullfile(project_root, 'src', 'ui', 'qe_browser_ui.m'));
end

function testHistorySaveFinalizesLiveSnapshot(testCase)
src = testCase.TestData.gui_src;
verifyHasText(testCase, src, 'local_ensure_latest_history_snapshot');
verifyTrue(testCase, ~isempty(regexp(src, ...
    'function\s+local_on_save_history[\s\S]*local_ensure_latest_history_snapshot\("Saved current GUI state"\)', ...
    'once')), 'Save History should append a final live snapshot before writing the MAT file.');
end

function testSnapshotCoversAutoFitAndCuratedControls(testCase)
src = testCase.TestData.gui_src;
fields = {
    'prominence', ...
    'autoFitSmoothWidth', ...
    'maxPeaks', ...
    'peakModel', ...
    'bootstrapCiSamples', ...
    'maxShift', ...
    'guessText', ...
    'dispModel', ...
    'correctionBranchTarget', ...
    'exportRatio', ...
    'selectedQIndex', ...
    'selectedQ_Ainv'};
for k = 1:numel(fields)
    verifyHasText(testCase, src, ['snap.' fields{k}], ...
        sprintf('Snapshot should capture and persist %s.', fields{k}));
end
end

function testFitControlsHaveHistoryCallbacks(testCase)
ui_src = testCase.TestData.ui_src;
verifyHasText(testCase, ui_src, 'cb.on_fit_param_change');
verifyGreaterThanOrEqual(testCase, countOccurrences(ui_src, 'cb.on_fit_param_change'), 6);
verifyHasText(testCase, ui_src, 'cb.on_correction_target_change');
end

function verifyHasText(testCase, text, pattern, message)
if nargin < 4
    message = sprintf('Expected source to contain: %s', pattern);
end
verifyTrue(testCase, contains(string(text), string(pattern)), message);
end

function n = countOccurrences(text, pattern)
n = numel(strfind(text, pattern));
end
