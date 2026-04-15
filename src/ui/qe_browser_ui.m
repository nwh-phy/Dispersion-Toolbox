function ui_handles = qe_browser_ui(cb)
% QE_BROWSER_UI  Build the Nion-style q‑E Browser GUI layout.
%
%   ui_handles = qe_browser_ui(cb)
%
%   cb  — struct of callback function handles:
%         .on_load_data, .on_build_views, .on_undo, .on_reset_stage,
%         .on_display_change, .on_normalization_change, .on_dq_override_changed,
%         .on_qnorm_changed, .on_norm_method_changed, .on_auto_y_changed,
%         .on_pick_peaks, .on_undo_pt, .on_clear_pts, .on_save_pts,
%         .on_load_pts, .on_split_branches, .on_fit_model, .on_export,
%         .on_auto_fit_dispersion, .on_fit_spectrum, .on_accept_fit,
%         .on_show_gamma, .on_pick_guesses, .on_reassign_points,
%         .on_fit_dispersion, .on_export_dispersion,
%         .on_history_select, .on_save_history, .on_load_history,
%         .on_clear_history
%
%   Returns a struct of all UI handle fields (Figure, axes, buttons, etc).
%
%   This function is purely declarative — it creates widgets and wires
%   callbacks but does NOT read or modify application state.

    % ═══════════ MAIN FIGURE ═══════════
    fig = uifigure( ...
        "Name", "Nion-style q-E Browser", ...
        "Position", [40 40 1780 980], ...
        "Color", [0.97 0.975 0.985]);

    main_grid = uigridlayout(fig, [3 3]);
    main_grid.RowHeight = {330, "1x", "1x"};
    main_grid.ColumnWidth = {"1x", "1x", 280};
    main_grid.Padding = [10 10 10 10];
    main_grid.RowSpacing = 10;
    main_grid.ColumnSpacing = 10;

    % ═══════════ HISTORY PANEL (right side) ═══════════
    history_panel = uipanel(main_grid, "Title", "Operation History");
    history_panel.Layout.Row = [1 3];
    history_panel.Layout.Column = 3;
    history_grid = uigridlayout(history_panel, [3 1]);
    history_grid.RowHeight = {"1x", 28, 28};
    history_grid.Padding = [4 4 4 4];
    history_grid.RowSpacing = 4;

    history_list = uilistbox(history_grid, ...
        "Items", {}, ...
        "ValueChangedFcn", cb.on_history_select);
    history_list.Layout.Row = 1;

    % Save / Load row
    history_btn_grid = uigridlayout(history_grid, [1 2]);
    history_btn_grid.Layout.Row = 2;
    history_btn_grid.ColumnWidth = {"1x", "1x"};
    history_btn_grid.Padding = [0 0 0 0];
    history_btn_grid.ColumnSpacing = 4;

    history_save_btn = uibutton(history_btn_grid, ...
        "Text", "💾 Save", ...
        "ButtonPushedFcn", cb.on_save_history);
    history_save_btn.Layout.Column = 1;

    history_load_btn = uibutton(history_btn_grid, ...
        "Text", "📂 Load", ...
        "ButtonPushedFcn", cb.on_load_history);
    history_load_btn.Layout.Column = 2;

    history_clear_btn = uibutton(history_grid, ...
        "Text", "Clear History", ...
        "ButtonPushedFcn", cb.on_clear_history);
    history_clear_btn.Layout.Row = 3;

    % ═══════════ CONTROL PANEL ═══════════
    control_panel = uipanel(main_grid, "Title", "Controls");
    control_panel.Layout.Row = 1;
    control_panel.Layout.Column = [1 2];

    control_grid = uigridlayout(control_panel, [10 16]);
    control_grid.RowHeight = {24, 24, 24, 24, 24, 24, 24, 24, 24, 24};
    control_grid.ColumnWidth = {80, 80, 62, 100, 62, 80, 62, 80, 62, 80, 80, 80, 62, 80, 80, "1x"};
    control_grid.Padding = [6 4 6 4];
    control_grid.RowSpacing = 5;
    control_grid.ColumnSpacing = 6;

    % ═══════════ ROW 1: Data Loading & Pipeline ═══════════
    load_eq3d = uibutton(control_grid, ...
        "Text", "Load Data", ...
        "ButtonPushedFcn", cb.on_load_data);
    load_eq3d.Layout.Row = 1;
    load_eq3d.Layout.Column = 1;

    build_button = uibutton(control_grid, ...
        "Text", "Build Views", ...
        "ButtonPushedFcn", cb.on_build_views);
    build_button.Layout.Row = 1;
    build_button.Layout.Column = 2;

    undo_button = uibutton(control_grid, ...
        "Text", "↩ Undo", ...
        "Enable", "off", ...
        "ButtonPushedFcn", cb.on_undo);
    undo_button.Layout.Row = 1;
    undo_button.Layout.Column = 3;

    reset_button = uibutton(control_grid, ...
        "Text", "Reset To", ...
        "ButtonPushedFcn", cb.on_reset_stage);
    reset_button.Layout.Row = 1;
    reset_button.Layout.Column = 4;

    reset_stage = uidropdown(control_grid, ...
        "Items", {'raw'}, ...
        "Value", 'raw');
    reset_stage.Layout.Row = 1;
    reset_stage.Layout.Column = 5;

    view_mode = uidropdown(control_grid, ...
        "Items", {'Physical', 'Normalized'}, ...
        "Value", 'Physical', ...
        "ValueChangedFcn", cb.on_display_change);
    view_mode.Layout.Row = 1;
    view_mode.Layout.Column = 6;

    stage_label = uilabel(control_grid, ...
        "Text", "Stage: none", ...
        "HorizontalAlignment", "left");
    stage_label.Layout.Row = 1;
    stage_label.Layout.Column = [7 8];

    source_label = uilabel(control_grid, ...
        "Text", "No data loaded", ...
        "HorizontalAlignment", "left");
    source_label.Layout.Row = 1;
    source_label.Layout.Column = [9 16];

    % ═══════════ ROW 2: Q Axis ═══════════
    add_label(control_grid, 2, 1, "dq (1/A)");
    dq_display = uieditfield(control_grid, "text", ...
        "Editable", "off", ...
        "Value", "");
    dq_display.Layout.Row = 2;
    dq_display.Layout.Column = 2;

    add_label(control_grid, 2, 3, "Ovr dq");
    dq_override = uieditfield(control_grid, "numeric", ...
        "Value", 0.005, ...
        "Limits", [eps inf], ...
        "LowerLimitInclusive", false, ...
        "ValueDisplayFormat", "%.5f", ...
        "ValueChangedFcn", cb.on_dq_override_changed);
    dq_override.Layout.Row = 2;
    dq_override.Layout.Column = 4;

    add_label(control_grid, 2, 5, "q start");
    q_start = uieditfield(control_grid, "numeric", ...
        "Value", 0, ...
        "ValueDisplayFormat", "%.4f", ...
        "ValueChangedFcn", cb.on_display_change);
    q_start.Layout.Row = 2;
    q_start.Layout.Column = 6;

    add_label(control_grid, 2, 7, "q end");
    q_end = uieditfield(control_grid, "numeric", ...
        "Value", 0, ...
        "ValueDisplayFormat", "%.4f", ...
        "ValueChangedFcn", cb.on_display_change);
    q_end.Layout.Row = 2;
    q_end.Layout.Column = 8;

    add_label(control_grid, 2, 9, "q step");
    q_step = uieditfield(control_grid, "numeric", ...
        "Value", 0.005, ...
        "Limits", [eps inf], ...
        "LowerLimitInclusive", false, ...
        "ValueDisplayFormat", "%.4f", ...
        "ValueChangedFcn", cb.on_display_change);
    q_step.Layout.Row = 2;
    q_step.Layout.Column = 10;

    q_norm = uicheckbox(control_grid, ...
        "Text", "Exp q equal.", ...
        "Value", false, ...
        "ValueChangedFcn", cb.on_qnorm_changed);
    q_norm.Layout.Row = 2;
    q_norm.Layout.Column = [11 12];

    % ═══════════ ROW 3: Energy & Y Axis ═══════════
    add_label(control_grid, 3, 1, "E min");
    energy_min = uieditfield(control_grid, "numeric", ...
        "Value", 0, ...
        "ValueDisplayFormat", "%.1f", ...
        "ValueChangedFcn", cb.on_display_change);
    energy_min.Layout.Row = 3;
    energy_min.Layout.Column = 2;

    add_label(control_grid, 3, 3, "E max");
    energy_max = uieditfield(control_grid, "numeric", ...
        "Value", 0, ...
        "ValueDisplayFormat", "%.1f", ...
        "ValueChangedFcn", cb.on_display_change);
    energy_max.Layout.Row = 3;
    energy_max.Layout.Column = 4;

    auto_y = uicheckbox(control_grid, ...
        "Text", "Auto Y", ...
        "Value", true, ...
        "ValueChangedFcn", cb.on_auto_y_changed);
    auto_y.Layout.Row = 3;
    auto_y.Layout.Column = 5;

    add_label(control_grid, 3, 6, "Y min");
    y_min = uieditfield(control_grid, "numeric", ...
        "Value", 0, ...
        "ValueDisplayFormat", "%.4g", ...
        "ValueChangedFcn", cb.on_display_change);
    y_min.Layout.Row = 3;
    y_min.Layout.Column = 7;

    add_label(control_grid, 3, 8, "Y max");
    y_max = uieditfield(control_grid, "numeric", ...
        "Value", 1, ...
        "ValueDisplayFormat", "%.4g", ...
        "ValueChangedFcn", cb.on_display_change);
    y_max.Layout.Row = 3;
    y_max.Layout.Column = 9;

    add_label(control_grid, 3, 10, "Y scale");
    y_scale = uidropdown(control_grid, ...
        "Items", {'linear', 'log'}, ...
        "Value", 'linear', ...
        "ValueChangedFcn", cb.on_display_change);
    y_scale.Layout.Row = 3;
    y_scale.Layout.Column = 11;

    % ═══════════ ROW 4: Display & Trace ═══════════
    add_label(control_grid, 4, 1, "Trace");
    trace_mode = uidropdown(control_grid, ...
        "Items", {'display', 'background-subtracted', 'curvature'}, ...
        "Value", 'display', ...
        "ValueChangedFcn", cb.on_display_change);
    trace_mode.Layout.Row = 4;
    trace_mode.Layout.Column = 2;

    add_label(control_grid, 4, 3, "Offset");
    offset = uieditfield(control_grid, "numeric", ...
        "Value", 0.25, ...
        "ValueDisplayFormat", "%.4g", ...
        "ValueChangedFcn", cb.on_display_change);
    offset.Layout.Row = 4;
    offset.Layout.Column = 4;

    add_label(control_grid, 4, 5, "Base λ");
    baseline_lambda = uieditfield(control_grid, "numeric", ...
        "Value", 1e6, ...
        "Limits", [1e2 inf], ...
        "ValueDisplayFormat", "%.3g", ...
        "ValueChangedFcn", cb.on_display_change);
    baseline_lambda.Layout.Row = 4;
    baseline_lambda.Layout.Column = 6;

    add_label(control_grid, 4, 7, "Smooth");
    smoothing_mode = uidropdown(control_grid, ...
        "Items", {'gaussian', 'off', 'sgolay'}, ...
        "Value", 'gaussian', ...
        "ValueChangedFcn", cb.on_display_change);
    smoothing_mode.Layout.Row = 4;
    smoothing_mode.Layout.Column = 8;

    add_label(control_grid, 4, 9, "Width");
    smoothing_width = uieditfield(control_grid, "numeric", ...
        "Value", 3, ...
        "Limits", [1 inf], ...
        "RoundFractionalValues", "on", ...
        "ValueDisplayFormat", "%.0f", ...
        "ValueChangedFcn", cb.on_display_change);
    smoothing_width.Layout.Row = 4;
    smoothing_width.Layout.Column = 10;

    add_label(control_grid, 4, 11, "SG ord");
    sg_order = uieditfield(control_grid, "numeric", ...
        "Value", 3, ...
        "Limits", [1 7], ...
        "ValueDisplayFormat", "%.0f", ...
        "ValueChangedFcn", cb.on_display_change);
    sg_order.Layout.Row = 4;
    sg_order.Layout.Column = 12;

    add_label(control_grid, 4, 13, "SG len");
    sg_framelen = uieditfield(control_grid, "numeric", ...
        "Value", 15, ...
        "Limits", [5 201], ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "Frame length (must be odd)", ...
        "ValueChangedFcn", cb.on_display_change);
    sg_framelen.Layout.Row = 4;
    sg_framelen.Layout.Column = 14;

    % ═══════════ ROW 5: Normalization & BG ═══════════
    add_label(control_grid, 5, 1, "ZLP min");
    zlp_min = uieditfield(control_grid, "numeric", ...
        "Value", -20, ...
        "ValueDisplayFormat", "%.1f", ...
        "ValueChangedFcn", cb.on_normalization_change);
    zlp_min.Layout.Row = 5;
    zlp_min.Layout.Column = 2;

    add_label(control_grid, 5, 3, "ZLP max");
    zlp_max = uieditfield(control_grid, "numeric", ...
        "Value", 20, ...
        "ValueDisplayFormat", "%.1f", ...
        "ValueChangedFcn", cb.on_normalization_change);
    zlp_max.Layout.Row = 5;
    zlp_max.Layout.Column = 4;

    add_label(control_grid, 5, 5, "Ref|q|min");
    ref_abs_q_min = uieditfield(control_grid, "numeric", ...
        "Value", 0.010, ...
        "Limits", [eps inf], ...
        "LowerLimitInclusive", false, ...
        "ValueDisplayFormat", "%.5f", ...
        "ValueChangedFcn", cb.on_normalization_change);
    ref_abs_q_min.Layout.Row = 5;
    ref_abs_q_min.Layout.Column = 6;

    add_label(control_grid, 5, 7, "Ref|q|max");
    ref_abs_q_max = uieditfield(control_grid, "numeric", ...
        "Value", 0.020, ...
        "Limits", [eps inf], ...
        "LowerLimitInclusive", false, ...
        "ValueDisplayFormat", "%.5f", ...
        "ValueChangedFcn", cb.on_normalization_change);
    ref_abs_q_max.Layout.Row = 5;
    ref_abs_q_max.Layout.Column = 8;

    area_norm = uicheckbox(control_grid, ...
        "Text", "Normalize", ...
        "Value", false, ...
        "ValueChangedFcn", cb.on_display_change);
    area_norm.Layout.Row = 5;
    area_norm.Layout.Column = 9;

    norm_method = uidropdown(control_grid, ...
        "Items", {'Area', 'ZLP Peak'}, ...
        "Value", 'ZLP Peak', ...
        "Tooltip", "Area: divide by integrated spectral weight (destroys A(q)). ZLP Peak: divide by ZLP peak height (preserves A(q)).", ...
        "ValueChangedFcn", cb.on_norm_method_changed);
    norm_method.Layout.Row = 5;
    norm_method.Layout.Column = 10;

    area_norm_min = uieditfield(control_grid, "numeric", ...
        "Value", -50, ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "ZLP Peak: ZLP window low bound (meV). Area: integration window low bound.", ...
        "ValueChangedFcn", cb.on_display_change);
    area_norm_min.Layout.Row = 5;
    area_norm_min.Layout.Column = 11;

    norm_hi_label = uilabel(control_grid, "Text", "~", "HorizontalAlignment", "center");
    norm_hi_label.Layout.Row = 5;
    norm_hi_label.Layout.Column = 12;
    area_norm_max = uieditfield(control_grid, "numeric", ...
        "Value", 50, ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "ZLP Peak: ZLP window high bound (meV). Area: integration window high bound.", ...
        "ValueChangedFcn", cb.on_display_change);
    area_norm_max.Layout.Row = 5;
    area_norm_max.Layout.Column = 13;

    bg_sub = uicheckbox(control_grid, ...
        "Text", "QE BG", ...
        "Value", false, ...
        "Tooltip", "Quasi-elastic background removal: subtract ZLP tail before peak fitting", ...
        "ValueChangedFcn", cb.on_display_change);
    bg_sub.Layout.Row = 5;
    bg_sub.Layout.Column = 14;

    bg_method = uidropdown(control_grid, ...
        "Items", {'Auto', 'Power', 'ExpPoly3', 'Pearson', 'PearsonVII', 'Power2', 'Exp2'}, ...
        "Value", 'Auto', ...
        "ValueChangedFcn", cb.on_display_change);
    bg_method.Layout.Row = 5;
    bg_method.Layout.Column = 15;

    % ═══════════ ROW 6: Denoise & Deconv ═══════════
    denoise_cb = uicheckbox(control_grid, ...
        "Text", "Denoise", ...
        "Value", false, ...
        "ValueChangedFcn", cb.on_display_change);
    denoise_cb.Layout.Row = 6;
    denoise_cb.Layout.Column = 1;

    denoise_method = uidropdown(control_grid, ...
        "Items", {'Wiener2D', 'SavGol', 'BM3D'}, ...
        "Value", 'Wiener2D', ...
        "ValueChangedFcn", cb.on_display_change);
    denoise_method.Layout.Row = 6;
    denoise_method.Layout.Column = 2;

    add_label(control_grid, 6, 3, "Noise σ");
    denoise_sigma = uieditfield(control_grid, "numeric", ...
        "Value", 0, ...
        "Limits", [0 inf], ...
        "ValueDisplayFormat", "%.1f", ...
        "Tooltip", "0 = auto-estimate from data", ...
        "ValueChangedFcn", cb.on_display_change);
    denoise_sigma.Layout.Row = 6;
    denoise_sigma.Layout.Column = 4;

    deconv_cb = uicheckbox(control_grid, ...
        "Text", "Deconv", ...
        "Value", false, ...
        "ValueChangedFcn", cb.on_display_change);
    deconv_cb.Layout.Row = 6;
    deconv_cb.Layout.Column = 5;

    add_label(control_grid, 6, 6, "LR iter");
    deconv_iter = uieditfield(control_grid, "numeric", ...
        "Value", 5, ...
        "Limits", [1 100], ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "Lucy-Richardson iterations", ...
        "ValueChangedFcn", cb.on_display_change);
    deconv_iter.Layout.Row = 6;
    deconv_iter.Layout.Column = 7;

    % --- BG Window & Advanced Controls (Row 6, columns 8-16) ---
    add_label(control_grid, 6, 8, "BG W1");
    bg_win1_lo = uieditfield(control_grid, "numeric", ...
        "Value", 50, ...
        "Limits", [0 inf], ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "Primary fit window low (meV)", ...
        "ValueChangedFcn", cb.on_display_change);
    bg_win1_lo.Layout.Row = 6;
    bg_win1_lo.Layout.Column = 9;

    add_label(control_grid, 6, 10, "~");
    bg_win1_hi = uieditfield(control_grid, "numeric", ...
        "Value", 300, ...
        "Limits", [0 inf], ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "Primary fit window high (meV)", ...
        "ValueChangedFcn", cb.on_display_change);
    bg_win1_hi.Layout.Row = 6;
    bg_win1_hi.Layout.Column = 11;

    bg_dual_cb = uicheckbox(control_grid, ...
        "Text", "Dual", ...
        "Value", false, ...
        "Tooltip", "Enable dual-window fitting (anchor before AND after loss feature)", ...
        "ValueChangedFcn", cb.on_display_change);
    bg_dual_cb.Layout.Row = 6;
    bg_dual_cb.Layout.Column = 12;

    bg_win2_lo = uieditfield(control_grid, "numeric", ...
        "Value", 3000, ...
        "Limits", [0 inf], ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "Secondary fit window low (meV)", ...
        "ValueChangedFcn", cb.on_display_change);
    bg_win2_lo.Layout.Row = 6;
    bg_win2_lo.Layout.Column = 13;

    add_label(control_grid, 6, 14, "~");
    bg_win2_hi = uieditfield(control_grid, "numeric", ...
        "Value", 3500, ...
        "Limits", [0 inf], ...
        "ValueDisplayFormat", "%.0f", ...
        "Tooltip", "Secondary fit window high (meV)", ...
        "ValueChangedFcn", cb.on_display_change);
    bg_win2_hi.Layout.Row = 6;
    bg_win2_hi.Layout.Column = 15;

    bg_iter_cb = uicheckbox(control_grid, ...
        "Text", "Iter", ...
        "Value", false, ...
        "Tooltip", "Iterative refit: detect negative signal → reweight → refit", ...
        "ValueChangedFcn", cb.on_display_change);
    bg_iter_cb.Layout.Row = 6;
    bg_iter_cb.Layout.Column = 16;

    % ═══════════ ROW 7: Manual Peak Picking ═══════════
    pick_pts = uibutton(control_grid, ...
        "Text", "Pick Peaks", ...
        "ButtonPushedFcn", cb.on_pick_peaks);
    pick_pts.Layout.Row = 7;
    pick_pts.Layout.Column = 1;

    undo_pt = uibutton(control_grid, ...
        "Text", "Undo Pt", ...
        "ButtonPushedFcn", cb.on_undo_pt);
    undo_pt.Layout.Row = 7;
    undo_pt.Layout.Column = 2;

    clear_pts = uibutton(control_grid, ...
        "Text", "Clear Pts", ...
        "ButtonPushedFcn", cb.on_clear_pts);
    clear_pts.Layout.Row = 7;
    clear_pts.Layout.Column = 3;

    save_pts = uibutton(control_grid, ...
        "Text", "Save Pts", ...
        "ButtonPushedFcn", cb.on_save_pts);
    save_pts.Layout.Row = 7;
    save_pts.Layout.Column = 4;

    load_pts = uibutton(control_grid, ...
        "Text", "Load Pts", ...
        "ButtonPushedFcn", cb.on_load_pts);
    load_pts.Layout.Row = 7;
    load_pts.Layout.Column = 5;

    split_btn = uibutton(control_grid, ...
        "Text", "Split 3", ...
        "ButtonPushedFcn", cb.on_split_branches);
    split_btn.Layout.Row = 7;
    split_btn.Layout.Column = 6;

    fit_model_btn = uibutton(control_grid, ...
        "Text", "Fit Model", ...
        "ButtonPushedFcn", cb.on_fit_model);
    fit_model_btn.Layout.Row = 7;
    fit_model_btn.Layout.Column = 7;

    export_btn = uibutton(control_grid, ...
        "Text", "📷 Export", ...
        "ButtonPushedFcn", cb.on_export);
    export_btn.Layout.Row = 7;
    export_btn.Layout.Column = 8;

    auto_fit_btn = uibutton(control_grid, ...
        "Text", "Auto Fit ω(q)", ...
        "ButtonPushedFcn", cb.on_auto_fit_dispersion);
    auto_fit_btn.Layout.Row = 7;
    auto_fit_btn.Layout.Column = 9;

    pts_label = uilabel(control_grid, ...
        "Text", "0 pts", ...
        "HorizontalAlignment", "left");
    pts_label.Layout.Row = 7;
    pts_label.Layout.Column = 10;

    export_ratio = uidropdown(control_grid, ...
        "Items", {'1:1', '4:3', '3:2', '16:9', '2:1'}, ...
        "Value", '4:3');
    export_ratio.Layout.Row = 7;
    export_ratio.Layout.Column = 11;

    info_label = uilabel(control_grid, ...
        "Text", "", ...
        "HorizontalAlignment", "left", ...
        "Visible", "off");
    info_label.Layout.Row = 7;
    info_label.Layout.Column = [12 16];

    % ═══════════ ROW 8: Lorentz Fitting Parameters ═══════════
    prom_lbl = uilabel(control_grid, "Text", "Prominence:"); %#ok<NASGU>
    prom_lbl.Layout.Row = 8; prom_lbl.Layout.Column = 1;

    prom_field = uispinner(control_grid, ...
        "Value", 0.15, "Step", 0.05, ...
        "Limits", [0.01 1], "ValueDisplayFormat", "%.2f");
    prom_field.Layout.Row = 8; prom_field.Layout.Column = 2;

    smooth_lbl = uilabel(control_grid, "Text", "Smooth W:"); %#ok<NASGU>
    smooth_lbl.Layout.Row = 8; smooth_lbl.Layout.Column = 3;

    smooth_field = uispinner(control_grid, ...
        "Value", 25, "Step", 5, ...
        "Limits", [1 100], "RoundFractionalValues", "on");
    smooth_field.Layout.Row = 8; smooth_field.Layout.Column = 4;

    maxpk_lbl = uilabel(control_grid, "Text", "Max Peaks:"); %#ok<NASGU>
    maxpk_lbl.Layout.Row = 8; maxpk_lbl.Layout.Column = 5;

    maxpk_field = uispinner(control_grid, ...
        "Value", 3, "Step", 1, ...
        "Limits", [1 10], "RoundFractionalValues", "on");
    maxpk_field.Layout.Row = 8; maxpk_field.Layout.Column = 6;

    fit_spec_btn = uibutton(control_grid, ...
        "Text", "Fit Spectrum", ...
        "ButtonPushedFcn", cb.on_fit_spectrum);
    fit_spec_btn.Layout.Row = 8; fit_spec_btn.Layout.Column = 7;

    accept_fit_btn = uibutton(control_grid, ...
        "Text", "Accept Fit", ...
        "Enable", "off", ...
        "ButtonPushedFcn", cb.on_accept_fit);
    accept_fit_btn.Layout.Row = 8; accept_fit_btn.Layout.Column = 8;

    fit_info_label = uilabel(control_grid, ...
        "Text", "", ...
        "HorizontalAlignment", "left");
    fit_info_label.Layout.Row = 8;
    fit_info_label.Layout.Column = [14 16];

    show_gamma_btn = uibutton(control_grid, ...
        "Text", "Γ & Amp", ...
        "ButtonPushedFcn", cb.on_show_gamma);
    show_gamma_btn.Layout.Row = 8; show_gamma_btn.Layout.Column = 13;

    guess_lbl = uilabel(control_grid, "Text", "Guesses:"); %#ok<NASGU>
    guess_lbl.Layout.Row = 8; guess_lbl.Layout.Column = 9;

    guess_field = uieditfield(control_grid, 'text', ...
        "Value", "", ...
        "Placeholder", "e.g. 800,2500,3200");
    guess_field.Layout.Row = 8;
    guess_field.Layout.Column = [10 12];

    % ═══════════ ROW 9: Peak Model & Seed Propagation ═══════════
    pk_model_lbl = uilabel(control_grid, "Text", "Peak Model:"); %#ok<NASGU>
    pk_model_lbl.Layout.Row = 9; pk_model_lbl.Layout.Column = 1;

    pk_model_dropdown = uidropdown(control_grid, ...
        "Items", {'Lorentz', 'Gaussian', 'Voigt', 'Damped HO'}, ...
        "ItemsData", {'lorentz', 'gaussian', 'voigt', 'damped_ho'}, ...
        "Value", 'lorentz');
    pk_model_dropdown.Layout.Row = 9; pk_model_dropdown.Layout.Column = [2 3];

    shift_lbl = uilabel(control_grid, "Text", "Max Shift:"); %#ok<NASGU>
    shift_lbl.Layout.Row = 9; shift_lbl.Layout.Column = 4;

    max_shift_field = uispinner(control_grid, ...
        "Value", 80, "Step", 10, ...
        "Limits", [10 500], "RoundFractionalValues", "on");
    max_shift_field.Layout.Row = 9; max_shift_field.Layout.Column = 5;

    pick_guesses_btn = uibutton(control_grid, ...
        "Text", "Pick Guesses", ...
        "ButtonPushedFcn", cb.on_pick_guesses);
    pick_guesses_btn.Layout.Row = 9; pick_guesses_btn.Layout.Column = [6 7];

    reassign_btn = uibutton(control_grid, ...
        "Text", "Reassign Pts", ...
        "ButtonPushedFcn", cb.on_reassign_points);
    reassign_btn.Layout.Row = 9; reassign_btn.Layout.Column = [8 9];

    seed_info_label = uilabel(control_grid, ...
        "Text", "", ...
        "HorizontalAlignment", "left");
    seed_info_label.Layout.Row = 9;
    seed_info_label.Layout.Column = [10 16];

    % ═══════════ ROW 10: Dispersion Model Selector ═══════════
    disp_model_lbl = uilabel(control_grid, "Text", "Disp Model:"); %#ok<NASGU>
    disp_model_lbl.Layout.Row = 10; disp_model_lbl.Layout.Column = 1;

    disp_model_dropdown = uidropdown(control_grid, ...
        "Items", {'Quasi-2D Plasmon', 'Acoustic Linear', 'Optical Constant', 'Optical Quadratic'}, ...
        "ItemsData", {'quasi2d_plasmon', 'acoustic_linear', 'optical_constant', 'optical_quadratic'}, ...
        "Value", 'quasi2d_plasmon');
    disp_model_dropdown.Layout.Row = 10; disp_model_dropdown.Layout.Column = [2 4];

    fit_disp_btn = uibutton(control_grid, ...
        "Text", "Fit Dispersion", ...
        "ButtonPushedFcn", cb.on_fit_dispersion);
    fit_disp_btn.Layout.Row = 10; fit_disp_btn.Layout.Column = [5 6];

    export_disp_btn = uibutton(control_grid, ...
        "Text", "Export Disp.", ...
        "ButtonPushedFcn", cb.on_export_dispersion);
    export_disp_btn.Layout.Row = 10; export_disp_btn.Layout.Column = [7 8];

    loss_map_btn = uibutton(control_grid, ...
        "Text", "🌊 I_kin Map", ...
        "Tooltip", "Generate an I_kin-corrected q-E map (loss-function-like in experimental mode)", ...
        "ButtonPushedFcn", cb.on_show_loss_map);
    loss_map_btn.Layout.Row = 10; loss_map_btn.Layout.Column = [9 10];

    export_data_btn = uibutton(control_grid, ...
        "Text", "📦 Export Data", ...
        "Tooltip", "Export preprocessed q-E data (after denoise + BG removal) as .mat", ...
        "ButtonPushedFcn", cb.on_export_data);
    export_data_btn.Layout.Row = 10; export_data_btn.Layout.Column = [11 12];

    disp_info_label = uilabel(control_grid, ...
        "Text", "", ...
        "HorizontalAlignment", "left");
    disp_info_label.Layout.Row = 10;
    disp_info_label.Layout.Column = [13 16];

    % ═══════════ AXES ═══════════
    qe_axes = uiaxes(main_grid);
    qe_axes.Layout.Row = 2;
    qe_axes.Layout.Column = 1;
    title(qe_axes, "Physical q-E Map");
    xlabel(qe_axes, "q (1/A)");
    ylabel(qe_axes, "Energy relative to ZLP (meV)");

    comparison_axes = uiaxes(main_grid);
    comparison_axes.Layout.Row = 3;
    comparison_axes.Layout.Column = 1;
    title(comparison_axes, "Normalized Off-axis Component");
    xlabel(comparison_axes, "q (1/A)");
    ylabel(comparison_axes, "Energy relative to ZLP (meV)");

    single_axes = uiaxes(main_grid);
    single_axes.Layout.Row = 2;
    single_axes.Layout.Column = 2;
    title(single_axes, "Spectrum");
    xlabel(single_axes, "Energy relative to ZLP (meV)");
    ylabel(single_axes, "Intensity");

    lower_right_tabs = uitabgroup(main_grid);
    lower_right_tabs.Layout.Row = 3;
    lower_right_tabs.Layout.Column = 2;

    dispersion_tab = uitab(lower_right_tabs, "Title", "Dispersion");
    dispersion_grid = uigridlayout(dispersion_tab, [1 1]);
    dispersion_grid.Padding = [0 0 0 0];
    dispersion_grid.RowSpacing = 0;
    dispersion_grid.ColumnSpacing = 0;
    dispersion_axes = uiaxes(dispersion_grid);
    dispersion_axes.Layout.Row = 1;
    dispersion_axes.Layout.Column = 1;
    title(dispersion_axes, "Dispersion");
    xlabel(dispersion_axes, "q (1/A)");
    ylabel(dispersion_axes, "Peak energy (meV)");

    waterfall_tab = uitab(lower_right_tabs, "Title", "Stacked Spectra");
    waterfall_grid = uigridlayout(waterfall_tab, [1 1]);
    waterfall_grid.Padding = [0 0 0 0];
    waterfall_grid.RowSpacing = 0;
    waterfall_grid.ColumnSpacing = 0;
    waterfall_axes = uiaxes(waterfall_grid);
    waterfall_axes.Layout.Row = 1;
    waterfall_axes.Layout.Column = 1;
    title(waterfall_axes, "Stacked Spectra");
    xlabel(waterfall_axes, "Energy relative to ZLP (meV)");
    ylabel(waterfall_axes, "Offset intensity");

    attach_toolbar(qe_axes);
    attach_toolbar(comparison_axes);
    attach_toolbar(single_axes);
    attach_toolbar(dispersion_axes);
    attach_toolbar(waterfall_axes);

    % ═══════════ PACK HANDLES ═══════════
    ui_handles = struct();
    ui_handles.Figure = fig;
    ui_handles.LoadEq3DButton = load_eq3d;
    ui_handles.BuildViewsButton = build_button;
    ui_handles.UndoButton = undo_button;
    ui_handles.ResetButton = reset_button;
    ui_handles.ResetStageDropdown = reset_stage;
    ui_handles.ViewModeDropdown = view_mode;
    ui_handles.SourceLabel = source_label;
    ui_handles.DqDisplay = dq_display;
    ui_handles.DqOverrideField = dq_override;
    ui_handles.QNormCheckbox = q_norm;
    ui_handles.ZlpMinField = zlp_min;
    ui_handles.ZlpMaxField = zlp_max;
    ui_handles.RefAbsQMinField = ref_abs_q_min;
    ui_handles.QStartField = q_start;
    ui_handles.QEndField = q_end;
    ui_handles.QStepField = q_step;
    ui_handles.EnergyMinField = energy_min;
    ui_handles.EnergyMaxField = energy_max;
    ui_handles.StageLabel = stage_label;
    ui_handles.SmoothingModeDropdown = smoothing_mode;
    ui_handles.SmoothingWidthField = smoothing_width;
    ui_handles.YScaleDropdown = y_scale;
    ui_handles.AutoYCheckbox = auto_y;
    ui_handles.YMinField = y_min;
    ui_handles.YMaxField = y_max;
    ui_handles.OffsetField = offset;
    ui_handles.TraceModeDropdown = trace_mode;

    ui_handles.BaselineLambdaField = baseline_lambda;
    ui_handles.InfoLabel = info_label;
    ui_handles.RefAbsQMaxField = ref_abs_q_max;
    ui_handles.AreaNormCheckbox = area_norm;
    ui_handles.NormMethodDropdown = norm_method;
    ui_handles.NormHiLabel = norm_hi_label;
    ui_handles.AreaNormMinField = area_norm_min;
    ui_handles.AreaNormMaxField = area_norm_max;
    ui_handles.BgSubCheckbox = bg_sub;
    ui_handles.BgMethodDropdown = bg_method;
    ui_handles.BgWin1LoField = bg_win1_lo;
    ui_handles.BgWin1HiField = bg_win1_hi;
    ui_handles.BgDualCheckbox = bg_dual_cb;
    ui_handles.BgWin2LoField = bg_win2_lo;
    ui_handles.BgWin2HiField = bg_win2_hi;
    ui_handles.BgIterCheckbox = bg_iter_cb;
    ui_handles.QEAxes = qe_axes;
    ui_handles.ComparisonAxes = comparison_axes;
    ui_handles.SingleAxes = single_axes;
    ui_handles.LowerRightTabs = lower_right_tabs;
    ui_handles.DispersionTab = dispersion_tab;
    ui_handles.WaterfallTab = waterfall_tab;
    ui_handles.WaterfallAxes = waterfall_axes;
    ui_handles.DispersionAxes = dispersion_axes;
    ui_handles.ZLPAxes = dispersion_axes;
    ui_handles.PickPtsButton = pick_pts;
    ui_handles.UndoPtButton = undo_pt;
    ui_handles.ClearPtsButton = clear_pts;
    ui_handles.SavePtsButton = save_pts;
    ui_handles.LoadPtsButton = load_pts;
    ui_handles.PtsLabel = pts_label;
    ui_handles.SplitButton = split_btn;
    ui_handles.FitModelButton = fit_model_btn;
    ui_handles.AutoFitButton = auto_fit_btn;
    ui_handles.PromField = prom_field;
    ui_handles.SmoothField = smooth_field;
    ui_handles.MaxPeaksField = maxpk_field;
    ui_handles.FitSpecButton = fit_spec_btn;
    ui_handles.AcceptFitButton = accept_fit_btn;
    ui_handles.FitInfoLabel = fit_info_label;
    ui_handles.ShowGammaButton = show_gamma_btn;
    ui_handles.GuessField = guess_field;
    ui_handles.ExportButton = export_btn;
    ui_handles.ExportRatioDropdown = export_ratio;
    ui_handles.HistoryListBox = history_list;
    ui_handles.HistoryClearButton = history_clear_btn;
    ui_handles.HistorySaveButton = history_save_btn;
    ui_handles.HistoryLoadButton = history_load_btn;
    ui_handles.PeakModelDropdown = pk_model_dropdown;
    ui_handles.MaxShiftField = max_shift_field;
    ui_handles.SeedInfoLabel = seed_info_label;
    ui_handles.PickGuessesButton = pick_guesses_btn;
    ui_handles.DeconvCheckbox = deconv_cb;
    ui_handles.DeconvIterField = deconv_iter;
    ui_handles.DenoiseCheckbox = denoise_cb;
    ui_handles.DenoiseMethodDropdown = denoise_method;
    ui_handles.DenoiseSigmaField = denoise_sigma;
    ui_handles.SGOrderField = sg_order;
    ui_handles.SGFrameLenField = sg_framelen;
    ui_handles.DispModelDropdown = disp_model_dropdown;
    ui_handles.FitDispButton = fit_disp_btn;
    ui_handles.ExportDispButton = export_disp_btn;
    ui_handles.LossMapButton = loss_map_btn;
    ui_handles.ExportDataButton = export_data_btn;
    ui_handles.DispInfoLabel = disp_info_label;
end

% --- Helper: add a label at (row, col) ---
function add_label(parent, row, column, text_value)
    label = uilabel(parent, "Text", text_value, "HorizontalAlignment", "left");
    label.Layout.Row = row;
    label.Layout.Column = column;
end

% --- Helper: attach standard axis toolbar ---
function attach_toolbar(ax)
    axtoolbar(ax, {'zoomin', 'zoomout', 'pan', 'datacursor', 'restoreview'});
end
