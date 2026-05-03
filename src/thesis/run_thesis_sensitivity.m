function outputs = run_thesis_sensitivity(cfg)
%RUN_THESIS_SENSITIVITY Run the small thesis sensitivity audit matrix.
%   This is intentionally a compact matrix, not an exhaustive parameter scan.

if nargin < 1 || isempty(cfg)
    cfg = thesis_config();
end
variants = thesis_sensitivity_configs(cfg);
outputs = repmat(struct('name', '', 'output', []), 1, numel(variants));
summary_tables = cell(1, numel(variants));
for vi = 1:numel(variants)
    run_opts = struct();
    run_opts.output_dir = fullfile(cfg.output.sensitivity_dir, variants(vi).name);
    run_opts.write_outputs = true;
    run_opts.verbose = true;
    outputs(vi).name = variants(vi).name;
    outputs(vi).output = run_thesis_pipeline(variants(vi).cfg, run_opts);
    summary_tables{vi} = summarize_thesis_results(outputs(vi).output.results, variants(vi).name);
end

if ~exist(cfg.output.sensitivity_dir, 'dir')
    mkdir(cfg.output.sensitivity_dir);
end
summary = vertcat(summary_tables{:});
writetable(summary, fullfile(cfg.output.sensitivity_dir, 'sensitivity_summary.csv'));
end
