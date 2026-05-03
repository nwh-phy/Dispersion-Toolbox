function tbl = summarize_thesis_results(results, variant_name)
%SUMMARIZE_THESIS_RESULTS Compact table for thesis pipeline and sensitivity runs.

if nargin < 2 || isempty(variant_name)
    variant_name = "baseline";
end
variant_name = string(variant_name);

n = numel(results);
variant = strings(n, 1);
session = strings(n, 1);
success = false(n, 1);
low_n = zeros(n, 1);
high_n = zeros(n, 1);
rejected_n = zeros(n, 1);
raw_peak_n = zeros(n, 1);
fit_success_n = zeros(n, 1);
low_rho0_A = nan(n, 1);
low_Eflat_meV = nan(n, 1);
low_R2 = nan(n, 1);
low_RMSE_meV = nan(n, 1);
high_constant_meV = nan(n, 1);
high_R2 = nan(n, 1);
high_RMSE_meV = nan(n, 1);
error = strings(n, 1);

for i = 1:n
    variant(i) = variant_name;
    session(i) = string(results(i).session_name);
    success(i) = results(i).success;
    low_n(i) = local_get_scalar(results(i).branch_summary, 'low_n', 0);
    high_n(i) = local_get_scalar(results(i).branch_summary, 'high_n', 0);
    rejected_n(i) = local_get_scalar(results(i).branch_summary, 'rejected_n', 0);
    raw_peak_n(i) = local_get_scalar(results(i), 'n_raw_peaks', 0);
    fit_success_n(i) = local_get_scalar(results(i), 'n_fit_success', 0);
    error(i) = string(local_get_scalar(results(i), 'error', ''));

    if isfield(results(i), 'models') && isfield(results(i).models, 'low') && ...
            isfield(results(i).models.low, 'fit')
        low_fit = results(i).models.low.fit;
        low_rho0_A(i) = local_get_scalar(low_fit, 'rho0', NaN);
        low_Eflat_meV(i) = local_get_scalar(low_fit, 'E_flat_meV', NaN);
        low_R2(i) = local_get_scalar(low_fit, 'R_squared', NaN);
        low_RMSE_meV(i) = local_get_scalar(low_fit, 'RMSE_meV', NaN);
    end

    if isfield(results(i), 'models') && isfield(results(i).models, 'high') && ...
            isfield(results(i).models.high, 'fit')
        high_fit = results(i).models.high.fit;
        if isfield(high_fit, 'params') && ~isempty(high_fit.params)
            high_constant_meV(i) = high_fit.params(1);
        end
        high_R2(i) = local_get_scalar(high_fit, 'R_squared', NaN);
        high_RMSE_meV(i) = local_get_scalar(high_fit, 'RMSE_meV', NaN);
    end
end

tbl = table(variant, session, success, low_n, high_n, rejected_n, ...
    raw_peak_n, fit_success_n, low_rho0_A, low_Eflat_meV, low_R2, ...
    low_RMSE_meV, high_constant_meV, high_R2, high_RMSE_meV, error);
end


function value = local_get_scalar(s, field_name, default_value)
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    value = s.(field_name);
else
    value = default_value;
end
if isnumeric(value) || islogical(value)
    value = value(1);
elseif isstring(value)
    value = value(1);
elseif ischar(value)
    value = char(value);
else
    value = default_value;
end
end
