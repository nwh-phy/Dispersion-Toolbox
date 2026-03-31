%% startup.m — q-EELS Plasmon Dispersion Toolbox 路径初始化
%  在 MATLAB 中 cd 到本项目根目录后运行 startup 即可自动加载所有模块。
%  也可以将本文件路径添加到 MATLAB 的 userpath/startup.m 中实现自动加载。

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'src')));
fprintf('  q-EELS Toolbox loaded (root: %s)\n', project_root);
