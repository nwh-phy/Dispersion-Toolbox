function [energy, kind] = qe_fit_point_energy(result, kind)
%QE_FIT_POINT_ENERGY  Energy of each fitted peak that goes into the dispersion.
%
%   [energy, kind] = qe_fit_point_energy(result) takes the model E0 ('model_E0')
%   of every peak of a ZLP joint fit (result.method = 'zlp_joint'), whose E0 is
%   the pole the dispersion models describe. For any other fit it takes the
%   loss-function apex where that is finite ('loss_function_apex'), else E0, as
%   the background + peak fit always has.
%
%   energy = qe_fit_point_energy(result, kind) uses kinds recorded earlier (a
%   cell array, one per peak), so a saved fit keeps the energy it was made with.
%
%   See also: qe_fit_spectrum_joint, fit_loss_function

E0 = double(result.omega_p(:));
apex = nan(size(E0));
if isfield(result, 'apex_energy_meV') && numel(result.apex_energy_meV) == numel(E0)
    apex = double(result.apex_energy_meV(:));
end
if nargin < 2
    kind = repmat({'model_E0'}, numel(E0), 1);
    if ~(isfield(result, 'method') && strcmp(result.method, 'zlp_joint'))
        kind(isfinite(apex)) = {'loss_function_apex'};
    end
end
kind = cellstr(kind(:));
energy = E0;
use_apex = strcmp(kind, 'loss_function_apex');
energy(use_apex) = apex(use_apex);
end
