function cfg = bisb_slot_config()
%BISB_SLOT_CONFIG Wide-angle Bi scan Bi1117_slot.
%   slot-v4 keeps the v3 alignment (one ZLP shift per scan position from the
%   direct beam) and adds the detector-geometry corrections of
%   bisb_slot_geometry, configured in cfg.geometry. v3 results keep their own
%   cfg inside analysis.mat.
cfg.version = 'slot-v4';
cfg.dataset = 'Bi1117_slot';
cfg.channels = 65:388; % illuminated slot interior, inspected in raw elastic map
cfg.bin_width = 8;
cfg.background_windows = [160 220;180 240;200 260]; % below the observed peak upturn
cfg.peak_window = [380 1200];
cfg.bootstrap = 200;
cfg.seed = 1117;
cfg.geometry = struct( ...
    'channel_range', [84 385], ...            % slot edges and the clipped outer-left spot (ch 73) excluded
    'zlp_window_meV', 60, 'zlp_smooth', 5, ...
    'streak_fit_k', -110:10:-30, ...          % gain-side ridge of the diagonal streak
    'streak_search_slope', [1.3 3.2], ...     % search band |E| = (1.3..3.2) meV x |k|
    'streak_core_k', 12, 'streak_mask_meV', 30, 'streak_exclusion_meV', 45, ...
    'drift_bands', [300 450; 450 600; 600 800; 800 1000; 1000 1300; 1300 1700], ...
    'bin_width', 7, ...                       % odd: one bin centred on the direct beam
    'shape_window', [350 1600], ...           % above the streak on the +q loss side
    'replica_window', [350 2200], 'shift_grid', -60:1:60, ...
    'sym_slope_tol', 10, ...                  % % per eV allowed in ln(I(+q)/I(-q)) over shape_window
    'sg_frames', [21 41], ...
    'dq_Ainv', 0.0140237499, ...
    'dq_source', 'provisional ARM200 calibration, Bi spacing 0.33 nm (bi_slot_arm_calibration_20261117.md)', ...
    'mos2_a_A', 3.160);
end
