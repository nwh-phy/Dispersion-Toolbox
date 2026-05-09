function tests = test_interactive_qe_browser_fast_selection
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.project_root = project_root;
end


function testMapClickUsesFastSelectedQRefreshInsteadOfFullRedraw(testCase)
src_path = fullfile(testCase.TestData.project_root, 'src', 'interactive_qe_browser.m');
src = fileread(src_path);

body = extractLocalFunctionBody(src, 'local_on_map_clicked', 'local_on_single_axes_clicked');

verifyTrue(testCase, contains(body, 'local_update_selected_q_views'));
verifyFalse(testCase, contains(body, 'local_update_all_views'));
end


function testFitSpectrumOverlayUsesSelectedPeakModelName(testCase)
src_path = fullfile(testCase.TestData.project_root, 'src', 'interactive_qe_browser.m');
src = fileread(src_path);

body = extractLocalFunctionBody(src, 'local_on_fit_spectrum', 'local_on_accept_fit');

verifyTrue(testCase, contains(body, 'local_fit_model_label(result)'));
verifyTrue(testCase, contains(body, 'local_fit_overlay_tag()'));
verifyFalse(testCase, contains(body, '''DisplayName'', ''Lorentz fit'''));
end


function testWaterfallUsesAutoOffsetFromSelectedStackTraces(testCase)
src_path = fullfile(testCase.TestData.project_root, 'src', 'interactive_qe_browser.m');
src = fileread(src_path);

body = extractLocalFunctionBody(src, 'local_plot_waterfall', 'local_update_heatmap_branch_overlays');

verifyTrue(testCase, contains(body, '''offset'', 0'));
verifyTrue(testCase, contains(body, 'local_effective_stack_offset(stack.traces)'));
verifyTrue(testCase, contains(body, 'local_apply_stack_offset(stack, stack_offset)'));
verifyTrue(testCase, contains(body, 'qe_adaptive_stack_offsets(stack.traces, stack_offset)'));
end


function testManualOffsetChangeDisablesAutoStackOffset(testCase)
src_path = fullfile(testCase.TestData.project_root, 'src', 'interactive_qe_browser.m');
src = fileread(src_path);

body = extractLocalFunctionBody(src, 'local_on_display_change', 'local_on_norm_method_changed');

verifyTrue(testCase, contains(body, 'local_is_offset_control(src)'));
verifyTrue(testCase, contains(body, 'state.stackOffsetMode = "manual"'));
end


function body = extractLocalFunctionBody(src, start_name, next_name)
start_pat = sprintf('function %s', start_name);
next_pat = sprintf('function %s', next_name);
start_idx = strfind(src, start_pat);
next_idx = strfind(src, next_pat);

assert(~isempty(start_idx), 'Missing local function %s', start_name);
assert(~isempty(next_idx), 'Missing local function %s', next_name);

start_idx = start_idx(1);
next_idx = next_idx(find(next_idx > start_idx, 1));
assert(~isempty(next_idx), 'Missing following local function %s', next_name);

body = src(start_idx:next_idx-1);
end
