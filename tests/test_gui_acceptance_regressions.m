function tests = test_gui_acceptance_regressions
% Regressions found by native GUI acceptance: derived views and dark exports.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root, 'startup.m'));
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
data_dir = fullfile(folder, 'acceptance 10w');
mkdir(data_dir);
[a3, e] = synth_joint_dataset();
save(fullfile(data_dir, 'eq3D.mat'), 'a3', 'e');
testCase.TestData.data_dir = data_dir;
end


function testFreshComparisonMatchesControlsAndSession(testCase)
[app, h] = interactive_qe_browser(testCase.TestData.data_dir);
testCase.addTeardown(@delete, app);
verifyComparisonRoundTrip(testCase, h);
end


function testRecalibratedComparisonMatchesControlsAndSession(testCase)
[app, h] = interactive_qe_browser(testCase.TestData.data_dir);
testCase.addTeardown(@delete, app);
h.ui.DqOverrideField.Value = 0.00025;
h.ui.DqOverrideField.ValueChangedFcn(h.ui.DqOverrideField, []);
verifyComparisonRoundTrip(testCase, h);
end


function verifyComparisonRoundTrip(testCase, h)
% Do not press Build Views: the first displayed map must already be current.
before = h.state();
verifyEqual(testCase, before.comparisonQE.reference_abs_q_min_Ainv, h.ui.RefAbsQMinField.Value);
verifyEqual(testCase, before.comparisonQE.reference_abs_q_max_Ainv, h.ui.RefAbsQMaxField.Value);
original_image = findobj(h.ui.ComparisonAxes, 'Type', 'image');
assertNumElements(testCase, original_image, 1);
cdata = original_image.CData;
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
p = fullfile(folder, 'session.mat');
h.save_session(p);
[app, restored] = interactive_qe_browser(p);
testCase.addTeardown(@delete, app);
after = restored.state();
verifyTrue(testCase, isequaln(after.comparisonQE, before.comparisonQE), ...
    'The derived comparison must not change when the session is reopened.');
restored_image = findobj(restored.ui.ComparisonAxes, 'Type', 'image');
assertNumElements(testCase, restored_image, 1);
verifyTrue(testCase, isequaln(restored_image.CData, cdata), ...
    'The initially displayed and restored comparison pixels must match.');
verifyEqual(testCase, after.physicalQE, before.physicalQE);
end


function testDarkDefaultsExportLegibleWhitePanelsWithoutChangingSource(testCase)
% Explicit dark root defaults make this reproduce on light desktops too.
props = {'defaultAxesColor', 'defaultAxesXColor', 'defaultAxesYColor', ...
    'defaultTextColor', 'defaultLegendColor', 'defaultLegendTextColor'};
old = get(groot, props);
testCase.addTeardown(@set, groot, props, old);
set(groot, props, {[0.07 0.07 0.07], [0.85 0.85 0.85], [0.85 0.85 0.85], ...
    [0.85 0.85 0.85], [0.07 0.07 0.07], [0.85 0.85 0.85]});
app = uifigure('Visible', 'off');
testCase.addTeardown(@delete, app);
ax = uiaxes(app);
line = plot(ax, [0 0.5 1], [0.1 0.9 0.2], 'Color', [0.1 0.4 0.8], ...
    'DisplayName', 'Measured spectrum');
xlim(ax, [0 1]); ylim(ax, [0 1]);
title(ax, 'Publication spectrum'); xlabel(ax, 'Energy'); ylabel(ax, 'Intensity');
lg = legend(ax, 'show', 'Location', 'south');
before = struct('axesColor', ax.Color, 'xColor', ax.XColor, 'yColor', ax.YColor, ...
    'titleColor', ax.Title.Color, 'legendColor', lg.Color, 'legendTextColor', lg.TextColor, ...
    'lineColor', line.Color, 'x', line.XData, 'y', line.YData);
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
paths = qe_export_browser_panel(ax, fullfile(folder, 'dark_source'), ...
    aspect_ratio=4/3, resolution=100, write_pdf=true);
pixels = imread(paths.png);
[height, width, ~] = size(pixels);
interior = pixels(round(height*0.15):round(height*0.8), round(width*0.2):round(width*0.9), 1:3);
verifyGreaterThan(testCase, mean(all(interior > 245, 3), 'all'), 0.8, ...
    'The plotting area and legend must have a white background.');
regions = {[0.025 0.115 0.14 0.94], [0.90 0.985 0.14 0.94], [0.15 0.8 0.01 0.11]};
for k = 1:numel(regions)
    b = regions{k};
    crop = pixels(max(1,round(height*b(1))):round(height*b(2)), ...
        max(1,round(width*b(3))):round(width*b(4)), 1:3);
    verifyGreaterThan(testCase, nnz(all(crop < 100, 3)), 15, ...
        'Title and axis labels must remain dark on the white paper.');
end
verifyEqual(testCase, width/height, 4/3, 'AbsTol', 0.01);
verifyTrue(testCase, isfile(paths.pdf));
after = struct('axesColor', ax.Color, 'xColor', ax.XColor, 'yColor', ax.YColor, ...
    'titleColor', ax.Title.Color, 'legendColor', lg.Color, 'legendTextColor', lg.TextColor, ...
    'lineColor', line.Color, 'x', line.XData, 'y', line.YData);
verifyEqual(testCase, after, before);
end
