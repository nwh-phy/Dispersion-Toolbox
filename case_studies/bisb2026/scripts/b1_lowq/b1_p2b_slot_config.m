function [cfg, ocfg] = b1_p2b_slot_config(kind, h)
% Configurations of the P2b reruns with the slot acceptance across the dispersion
% direction, as modified copies of b1_p2b_config and b1_p2a_operator_config
% (pass both to b1_p2b_fit). The version tags change, so every rerun gets its own
% output folder.
%   'h66'    h = 0.0066 A^-1 instead of the assumed 0.0005: half width of the 0.125 mm
%            slot at 10w, from cali.txt (slot 132.8, beam 56.8585 per 1 mrad) and the 16:1
%            slot aspect (user confirmed the 0.125 mm slot, 2026-10-17). n_sy 9 -> 53
%            exit-box samples across the slot, spacing 2.5e-4 A^-1 (checked by
%            b1_p2b_h66_conv against 105). Variants: the main model plus the two whose
%            verdict the slot acceptance can change (y distribution, broad component on B1).
%   'hscan'  h given (A^-1). n_sy with spacing <= 2.5e-4 A^-1, at least the original 9.
%            E0 nodes at 0.0005 0.001 0.002 0.004 0.008, width and residue nodes at
%            0.0005 0.002 0.008: with a wide slot the channels collect Q up to
%            sqrt(0.005^2 + h^2) ~ 0.008, so the top node moves from 0.005 to 0.008
%            (b1_p2b_fit fixes the node counts). Variants: main and y: Gaussian.
% Until 2026-10-21 these were copies of the two config files in h66/ and hscan/, put
% first on the path; the configurations here reproduce those runs field by field.
cfg = b1_p2b_config(); ocfg = b1_p2a_operator_config();
switch kind
    case 'h66'
        h = 0.0066;
        cfg.version = [cfg.version '+h66'];
        keep = {'main', 'y: Gaussian', 'broad component off for B1'};
        ocfg.version = [ocfg.version '+h66'];
        ocfg.n_sy = 53;
    case 'hscan'
        assert(nargin == 2 && isfinite(h) && h > 0, 'b1_p2b_slot_config:h', 'hscan needs h > 0 (A^-1)');
        cfg.version = sprintf('%s+hscan-h%.4f-top0.008', cfg.version, h);
        cfg.nodes = struct('E0', [0.0005 0.001 0.002 0.004 0.008], 'G', [0.0005 0.002 0.008], 'A', [0.0005 0.002 0.008]);
        cfg.Q_curve = exp(linspace(log(0.0005), log(0.008), 49));
        keep = {'main', 'y: Gaussian'};
        ocfg.version = sprintf('%s+hscan-h%.4f', ocfg.version, h);
        ocfg.n_sy = max(9, 2 * ceil(h / 2.5e-4) + 1);
    otherwise
        error('b1_p2b_slot_config:kind', 'kind must be ''h66'' or ''hscan''');
end
sel = ismember({cfg.variants.name}, keep);
assert(nnz(sel) == numel(keep), 'b1_p2b_slot_config:variants', 'variant names changed in b1_p2b_config');
cfg.variants = cfg.variants(sel);
ocfg.h = h;
end
