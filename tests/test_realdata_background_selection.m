function tests = test_realdata_background_selection
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project_root, 'startup.m'));
testCase.TestData.project_root = project_root;
testCase.TestData.dataset_path = fullfile(project_root, ...
    '20260120 Bi', '590 PL2 10w 0.004 10sx300', 'eq3D.mat');
end


function testDualWindowAutoPrefersPearsonVIIAtQZeroWhenItAvoidsOversubtraction(testCase)
[qe_focus, opts] = loadRepresentativeLowQDataset(testCase);

[~, bg_diag] = qe_preprocess(qe_focus, opts);
q0_diag = bg_diag(qe_focus.q_zero_index);

pearson_idx = find(strcmp(q0_diag.candidate_methods, 'Pearson'), 1);
pearsonvii_idx = find(strcmp(q0_diag.candidate_methods, 'PearsonVII'), 1);
verifyNotEmpty(testCase, pearson_idx);
verifyNotEmpty(testCase, pearsonvii_idx);

pearson = q0_diag.candidate_details(pearson_idx);
pearsonvii = q0_diag.candidate_details(pearsonvii_idx);

verifyLessThan(testCase, pearsonvii.neg_area_fraction, pearson.neg_area_fraction);
verifyLessThan(testCase, pearsonvii.neg_fraction, pearson.neg_fraction);
verifyEqual(testCase, q0_diag.selected_method, 'PearsonVII');
end


function testDualWindowAutoCanLockCentralNeighborhoodToGroupedWinner(testCase)
[qe_focus, opts] = loadRepresentativeLowQDataset(testCase);
opts.bg_auto_group_qmax = 0.01;

[~, bg_diag] = qe_preprocess(qe_focus, opts);
center_idx = find(abs(qe_focus.q_Ainv) <= opts.bg_auto_group_qmax);
verifyNumElements(testCase, center_idx, 5);
selected = string({bg_diag(center_idx).selected_method});
verifyEqual(testCase, unique(selected), "PearsonVII");
end


function testShoulderGroupedAutoCanLockPositiveShoulderBand(testCase)
[qe_focus, opts] = load20wLowQDataset(testCase);
base_opts = opts;
base_opts.bg_auto_group_qmax = 0.01;
group_opts = base_opts;
group_opts.bg_auto_group_ranges = [0.02, 0.05];

[~, bg_diag_base] = qe_preprocess(qe_focus, base_opts);
[~, bg_diag_grouped] = qe_preprocess(qe_focus, group_opts);
pos_idx = find(qe_focus.q_Ainv >= 0.02 & qe_focus.q_Ainv <= 0.05);
verifyGreaterThan(testCase, numel(pos_idx), 1);
verifyEqual(testCase, unique(string({bg_diag_grouped(pos_idx).selected_method})), "PearsonVII");
verifyLessThan(testCase, mean([bg_diag_grouped(pos_idx).neg_peak_fraction], 'omitnan'), ...
    mean([bg_diag_base(pos_idx).neg_peak_fraction], 'omitnan'));
verifyLessThan(testCase, mean([bg_diag_grouped(pos_idx).neg_area_fraction], 'omitnan'), ...
    mean([bg_diag_base(pos_idx).neg_area_fraction], 'omitnan'));
end


function [qe_focus, opts] = loadRepresentativeLowQDataset(testCase)
dataset_path = testCase.TestData.dataset_path;
assumeTrue(testCase, exist(dataset_path, 'file') == 2, ...
    sprintf('Representative dataset missing: %s', dataset_path));

dataset = load_qe_dataset(dataset_path);
qe = dataset.qe;
focus_idx = find(abs(qe.q_Ainv) <= 0.15);
verifyNotEmpty(testCase, focus_idx);

qe_focus = qe;
qe_focus.intensity = qe.intensity(:, focus_idx);
qe_focus.q_channel = qe.q_channel(focus_idx);
qe_focus.q_Ainv = qe.q_Ainv(focus_idx);
qe_focus.q_abs_Ainv = qe.q_abs_Ainv(focus_idx);
[~, local_zero] = min(abs(qe_focus.q_Ainv));
qe_focus.q_zero_index = local_zero;

opts = struct();
opts.do_despike = false;
opts.do_normalize = false;
opts.do_denoise = false;
opts.do_bg_sub = true;
opts.bg_method = 'Auto';
opts.bg_candidate_methods = {'Power', 'Exp2', 'ExpPoly3', 'Pearson', 'PearsonVII'};
opts.bg_win_lo = [50, 300];
opts.bg_win_hi = [2800, 3400];
opts.bg_iterative = false;
opts.do_deconv = false;
end


function [qe_focus, opts] = load20wLowQDataset(testCase)
dataset_path = fullfile(testCase.TestData.project_root, '20260120 Bi', ...
    'no pl2 20w 0.004 10sx300 2film', 'eq3D.mat');
assumeTrue(testCase, exist(dataset_path, 'file') == 2, ...
    sprintf('Representative dataset missing: %s', dataset_path));

dataset = load_qe_dataset(dataset_path);
qe = dataset.qe;
focus_idx = find(abs(qe.q_Ainv) <= 0.15);
verifyNotEmpty(testCase, focus_idx);

qe_focus = qe;
qe_focus.intensity = qe.intensity(:, focus_idx);
qe_focus.q_channel = qe.q_channel(focus_idx);
qe_focus.q_Ainv = qe.q_Ainv(focus_idx);
qe_focus.q_abs_Ainv = qe.q_abs_Ainv(focus_idx);
[~, local_zero] = min(abs(qe_focus.q_Ainv));
qe_focus.q_zero_index = local_zero;

opts = struct();
opts.do_despike = false;
opts.do_normalize = false;
opts.do_denoise = false;
opts.do_bg_sub = true;
opts.bg_method = 'Auto';
opts.bg_candidate_methods = {'Power', 'Exp2', 'ExpPoly3', 'Pearson', 'PearsonVII'};
opts.bg_win_lo = [50, 300];
opts.bg_win_hi = [2800, 3400];
opts.bg_iterative = false;
opts.do_deconv = false;
end
