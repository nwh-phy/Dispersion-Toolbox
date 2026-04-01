%% startup.m — q-EELS Plasmon Dispersion Toolbox 路径初始化
%  在 MATLAB 中 cd 到本项目根目录后运行 startup 即可自动加载所有模块。
%  也可以将本文件路径添加到 MATLAB 的 userpath/startup.m 中实现自动加载。

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'src')));
fprintf('  q-EELS Toolbox loaded (root: %s)\n', project_root);

% --- Load external toolbox paths from user-local config ---
% Users: copy toolbox_config.m.template → toolbox_config.m and set paths.
try
    cfg = toolbox_config();
    if isfield(cfg, 'nion_toolbox_root') && ~isempty(cfg.nion_toolbox_root)
        nion_root = cfg.nion_toolbox_root;

        % 3D toolbox Process/ (background subtraction, etc.)
        p3d = fullfile(nion_root, '3D-EELS TOOLBOX', 'Process');
        if isfolder(p3d), addpath(p3d); end

        % BM3D denoising (separate directory)
        bm3d_dir = fullfile(nion_root, '3D-EELS TOOLBOX', 'BM3D');
        if isfolder(bm3d_dir), addpath(genpath(bm3d_dir)); end

        fprintf('  Nion EELS Toolbox loaded (root: %s)\n', nion_root);
    end
catch
    % toolbox_config.m not found — BM3D will be unavailable
end
