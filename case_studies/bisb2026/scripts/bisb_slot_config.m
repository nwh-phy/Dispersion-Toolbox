function cfg = bisb_slot_config()
%BISB_SLOT_CONFIG Wide-angle Bi scan; no reciprocal-space calibration assumed.
cfg.version = 'slot-v3';
cfg.dataset = 'Bi1117_slot';
cfg.channels = 65:388; % illuminated slot interior, inspected in raw elastic map
cfg.bin_width = 8;
cfg.background_windows = [160 220;180 240;200 260]; % below the observed peak upturn
cfg.peak_window = [380 1200];
cfg.bootstrap = 200;
cfg.seed = 1117;
end
