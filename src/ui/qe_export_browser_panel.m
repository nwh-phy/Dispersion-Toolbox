function paths = qe_export_browser_panel(source, prefix, options)
%QE_EXPORT_BROWSER_PANEL Export UI axes through ordinary graphics axes.
% Copying the UIAxes itself into an invisible ordinary figure can export a
% blank page. Copy its numeric graphics children into ordinary axes instead.
arguments
    source (1,1) matlab.graphics.axis.Axes
    prefix {mustBeTextScalar}
    options.aspect_ratio (1,1) double {mustBePositive} = 4/3
    options.resolution (1,1) double {mustBePositive} = 300
    options.write_pdf (1,1) logical = false
end
fig = figure('Visible', 'off', 'Color', 'white', 'Units', 'pixels', ...
    'Position', [100 100 round(600 * options.aspect_ratio) 600]);
cleanup = onCleanup(@() close(fig));
ax = axes(fig, 'Units', 'normalized', 'Position', [0.14 0.14 0.80 0.74]);
% Include peak markers with HandleVisibility='off', but recreate labels and
% leave the interactive toolbar in the browser.
children = allchild(source);
excluded = gobjects(1, 4);
excluded(1) = source.Title; excluded(2) = source.XLabel;
excluded(3) = source.YLabel; excluded(4) = source.ZLabel;
if isprop(source, 'Subtitle'), excluded(end+1) = source.Subtitle; end
if ~isempty(source.Toolbar), excluded = [excluded source.Toolbar]; end
children = children(~ismember(children, excluded));
copyobj(children, ax);
properties = {'XLim', 'YLim', 'XScale', 'YScale', 'XDir', 'YDir', 'CLim', ...
    'FontName', 'FontSize', 'XGrid', 'YGrid', 'Box', 'TickLabelInterpreter'};
for k = 1:numel(properties), ax.(properties{k}) = source.(properties{k}); end
colormap(ax, colormap(source));
xlabel(ax, source.XLabel.String, 'Interpreter', source.XLabel.Interpreter);
ylabel(ax, source.YLabel.String, 'Interpreter', source.YLabel.Interpreter);
title(ax, wrapTitle(source.Title.String, 70), 'Interpreter', source.Title.Interpreter);
if ~isempty(source.Legend) && isgraphics(source.Legend)
    labels = source.Legend.String;
    handles = gobjects(0); names = {};
    for k = 1:numel(labels)
        match = findobj(ax, '-property', 'DisplayName', 'DisplayName', labels{k});
        if ~isempty(match)
            handles(end+1) = match(1); names{end+1} = labels{k}; %#ok<AGROW>
        end
    end
    if ~isempty(handles)
        legend(ax, handles, names, 'Location', source.Legend.Location, ...
            'Interpreter', source.Legend.Interpreter, 'FontSize', source.Legend.FontSize);
    end
end
drawnow;
% Reserve actual title height so multiline provenance titles do not touch or
% cross the top of the fixed canvas.
ax.Title.Units = 'pixels';
extent = ax.Title.Extent;
position = getpixelposition(ax);
overflow = max(0, position(2) + extent(2) + extent(4) - fig.Position(4) + 16);
if overflow > 0
    layout = ax.Position;
    layout(4) = max(0.2, layout(4) - overflow / fig.Position(4));
    ax.Position = layout;
    drawnow;
end
paths = struct('png', [char(prefix) '.png'], 'pdf', '');
% A fixed paper canvas preserves the requested ratio on R2024b as well;
% exportgraphics Width/Height/Padding are available only in newer releases.
inches = fig.Position(3:4) / 96;
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 inches], 'PaperSize', inches);
print(fig, paths.png, '-dpng', sprintf('-r%d', round(options.resolution)));
if options.write_pdf
    paths.pdf = [char(prefix) '.pdf'];
    print(fig, paths.pdf, '-dpdf', '-painters', sprintf('-r%d', round(options.resolution)));
end
end

function lines = wrapTitle(value, width)
words = split(strjoin(string(value), ' '));
lines = {''};
for k = 1:numel(words)
    word = char(words(k));
    if numel(lines{end}) + numel(word) + 1 > width && ~isempty(lines{end})
        lines{end+1} = word; %#ok<AGROW>
    elseif isempty(lines{end})
        lines{end} = word;
    else
        lines{end} = [lines{end} ' ' word];
    end
end
end
