function cfg = b1_p2a_inputs_config()
% Configuration of P2a step 1 (input checks, 590); its SHA-256 names the
% output folder p2a/inputs_<session>_<hash8>.
cfg.version = 'P2a-step1b-2026-10-14';     % step 1 + fixes after review round 4, count cube for step 2
cfg.session = '590_PL2_10w';
cfg.q_max = 0.0065;                   % channels read, both signs, incl. q = 0
cfg.core_w = [8 16];                  % core half-widths |E - c(q) - dE_t| <= w (meV), per-frame descriptors
cfg.window_list = [8 12 16 32];       % summed-profile windows for the sigma_probe comparison
cfg.old_window = 12;                  % b1_adaptq_build: unaligned |E| <= 12 meV
cfg.moment_q_max = 0.001;             % moments over the fixed central five channels
cfg.gauss_q_max = 0.001;              % channel-integrated Gaussian, central five channels
cfg.centre_q_max = 0.001;             % channels whose energy centres give dE_t (count-weighted)
cfg.wing_centre_q = [0.002 0.004];    % |q| range of the wing channels (energy shift free of the bending)
cfg.tie_q_max = 0.0015;               % channels for the bending test c_t(q) against dq_t
cfg.tail_q_min = 0.0015;              % model-free tail fraction: counts at |q| >= tail_q_min
cfg.n_blocks = 6;
cfg.n_blocks_check = 12;              % second blocking for the time-convergence check (count cube only)
cfg.attrib_q_max = 0.0055;            % Delta-I attribution channels (neighbours up to q_max)
cfg.bands = [40 150; 150 250; 250 400; 400 1800];
cfg.key_q = [0.0005 0.001 0.002 0.005];
cfg.shape_edges = [40:10:300, 325:25:600, 675:75:1800];
cfg.shape_regions = [40 150; 150 250; 250 400; 400 1800];
cfg.holdout = [150 250];              % excluded from the common-beam reference (p1b/common_ref.mat)
end
