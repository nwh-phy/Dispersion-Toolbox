function cfg = b1_p2a_response_config()
% Configuration of P2a step 2.1 (slice response); its SHA-256 names the output folder.
cfg.version = 'P2a-step2.1b-2026-10-15';   % light-tailed line model, support 60 meV
cfg.fit_window = [-60 20];     % elastic-line fit window around the centre (meV), core constrained (KA-like)
cfg.support = 60;              % kernel support |E - c| <= support (meV)
cfg.support_check = [40 100];  % supports for the convergence check
cfg.capture_w = [8 16];        % core half-widths of P2a step 1 (capture-fraction check)
cfg.blockings = [6 12];
cfg.tail_report = 40;          % model mass beyond |E - c| > 40 meV is reported
end
