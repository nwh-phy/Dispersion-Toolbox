function app = interactive_qe_browser(startPath, historyFile)
arguments
    startPath {mustBeTextScalar} = ""
    historyFile {mustBeTextScalar} = ""
end

state = local_initial_state();

% Add toolbox background subtraction to MATLAB path
toolbox_process_dir = fullfile(getenv('USERPROFILE'), 'Desktop', ...
    'Nion-EELS数据处理toolbox', '3D-EELS TOOLBOX', 'Process');
if isfolder(toolbox_process_dir)
    addpath(toolbox_process_dir);
end

ui = local_build_ui();

local_sync_manual_y_state();

if strlength(string(startPath)) > 0
    local_load_start_path(string(startPath));
end

% Auto-load history if provided
if strlength(string(historyFile)) > 0 && isfile(historyFile)
    loaded = load(char(historyFile));
    if isfield(loaded, 'opHistory') && ~isempty(loaded.opHistory)
        state.opHistory = loaded.opHistory;
        local_refresh_history_list();
        % Restore the last snapshot
        last_snap = state.opHistory{end}.snapshot;
        try
            local_restore_param_snapshot(last_snap);
            local_rebuild_comparison(false);
            local_update_all_views();
        catch
        end
        ui.InfoLabel.Text = sprintf('History restored (%d entries)', numel(state.opHistory));
        ui.InfoLabel.Visible = "on";
    end
end

if nargout > 0
    app = ui.Figure;
end


    function state_out = local_initial_state()
        state_out = struct();
        state_out.dataset = [];
        state_out.raw4d = [];
        state_out.cropped4d = [];
        state_out.aligned4d = [];
        state_out.binned4d = [];
        state_out.eq3dQE = [];
        state_out.physicalQE = [];
        state_out.comparisonQE = [];
        state_out.cropSpec = [];
        state_out.activeStage = "none";
        state_out.stageHistory = {};  % stack of previous activeStage strings
        state_out.paramHistory = {};  % kept for compat, not used
        state_out.opHistory = {};     % cell array of operation history entries
        state_out.selectedQIndex = 1;
        state_out.manualPickArmed = false;
        state_out.guessPickArmed = false;
        state_out.manual_points = [];  % Nx2 [abs_q, energy]
        state_out.manual_branches = {};  % cell array of Nx2, one per branch
        state_out.fitResults = {};       % cell array of fit_quasi2d_plasmon results
        state_out.autoFitResults = [];   % struct array from auto Drude-Lorentz fit
        state_out.pendingFit = [];       % pending single-spectrum Lorentz fit result
        state_out.qeImage = gobjects(1);
        state_out.selectionMarker = gobjects(1);
    end


    function ui_handles = local_build_ui()
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
            "ValueChangedFcn", @local_on_history_select);
        history_list.Layout.Row = 1;

        % Save / Load row
        history_btn_grid = uigridlayout(history_grid, [1 2]);
        history_btn_grid.Layout.Row = 2;
        history_btn_grid.ColumnWidth = {"1x", "1x"};
        history_btn_grid.Padding = [0 0 0 0];
        history_btn_grid.ColumnSpacing = 4;

        history_save_btn = uibutton(history_btn_grid, ...
            "Text", "💾 Save", ...
            "ButtonPushedFcn", @local_on_save_history);
        history_save_btn.Layout.Column = 1;

        history_load_btn = uibutton(history_btn_grid, ...
            "Text", "📂 Load", ...
            "ButtonPushedFcn", @local_on_load_history);
        history_load_btn.Layout.Column = 2;

        history_clear_btn = uibutton(history_grid, ...
            "Text", "Clear History", ...
            "ButtonPushedFcn", @local_on_clear_history);
        history_clear_btn.Layout.Row = 3;

        control_panel = uipanel(main_grid, "Title", "Controls");
        control_panel.Layout.Row = 1;
        control_panel.Layout.Column = [1 2];
        % keep col 3 for history panel

        control_grid = uigridlayout(control_panel, [10 16]);
        control_grid.RowHeight = {24, 24, 24, 24, 24, 24, 24, 24, 24, 24};
        control_grid.ColumnWidth = {80, 80, 62, 100, 62, 80, 62, 80, 62, 80, 80, 80, 62, 80, 80, "1x"};
        control_grid.Padding = [6 4 6 4];
        control_grid.RowSpacing = 5;
        control_grid.ColumnSpacing = 6;

        % ═══════════ ROW 1: Data Loading & Pipeline ═══════════
        load_eq3d = uibutton(control_grid, ...
            "Text", "Load Data", ...
            "ButtonPushedFcn", @local_on_load_data);
        load_eq3d.Layout.Row = 1;
        load_eq3d.Layout.Column = 1;

        build_button = uibutton(control_grid, ...
            "Text", "Build Views", ...
            "ButtonPushedFcn", @local_on_build_views);
        build_button.Layout.Row = 1;
        build_button.Layout.Column = 2;

        undo_button = uibutton(control_grid, ...
            "Text", "↩ Undo", ...
            "Enable", "off", ...
            "ButtonPushedFcn", @local_on_undo);
        undo_button.Layout.Row = 1;
        undo_button.Layout.Column = 3;

        reset_button = uibutton(control_grid, ...
            "Text", "Reset To", ...
            "ButtonPushedFcn", @local_on_reset_stage);
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
            "ValueChangedFcn", @local_on_display_change);
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
        local_add_label(control_grid, 2, 1, "dq (1/A)");
        dq_display = uieditfield(control_grid, "text", ...
            "Editable", "off", ...
            "Value", "");
        dq_display.Layout.Row = 2;
        dq_display.Layout.Column = 2;

        local_add_label(control_grid, 2, 3, "Ovr dq");
        dq_override = uieditfield(control_grid, "numeric", ...
            "Value", 0.005, ...
            "Limits", [eps inf], ...
            "LowerLimitInclusive", false, ...
            "ValueDisplayFormat", "%.5f", ...
            "ValueChangedFcn", @local_on_dq_override_changed);
        dq_override.Layout.Row = 2;
        dq_override.Layout.Column = 4;

        local_add_label(control_grid, 2, 5, "q start");
        q_start = uieditfield(control_grid, "numeric", ...
            "Value", 0, ...
            "ValueDisplayFormat", "%.4f", ...
            "ValueChangedFcn", @local_on_display_change);
        q_start.Layout.Row = 2;
        q_start.Layout.Column = 6;

        local_add_label(control_grid, 2, 7, "q end");
        q_end = uieditfield(control_grid, "numeric", ...
            "Value", 0, ...
            "ValueDisplayFormat", "%.4f", ...
            "ValueChangedFcn", @local_on_display_change);
        q_end.Layout.Row = 2;
        q_end.Layout.Column = 8;

        local_add_label(control_grid, 2, 9, "q step");
        q_step = uieditfield(control_grid, "numeric", ...
            "Value", 0.005, ...
            "Limits", [eps inf], ...
            "LowerLimitInclusive", false, ...
            "ValueDisplayFormat", "%.4f", ...
            "ValueChangedFcn", @local_on_display_change);
        q_step.Layout.Row = 2;
        q_step.Layout.Column = 10;

        q_norm = uicheckbox(control_grid, ...
            "Text", "Exp q equal.", ...
            "Value", false, ...
            "ValueChangedFcn", @local_on_qnorm_changed);
        q_norm.Layout.Row = 2;
        q_norm.Layout.Column = [11 12];

        % ═══════════ ROW 3: Energy & Y Axis ═══════════
        local_add_label(control_grid, 3, 1, "E min");
        energy_min = uieditfield(control_grid, "numeric", ...
            "Value", 0, ...
            "ValueDisplayFormat", "%.1f", ...
            "ValueChangedFcn", @local_on_display_change);
        energy_min.Layout.Row = 3;
        energy_min.Layout.Column = 2;

        local_add_label(control_grid, 3, 3, "E max");
        energy_max = uieditfield(control_grid, "numeric", ...
            "Value", 0, ...
            "ValueDisplayFormat", "%.1f", ...
            "ValueChangedFcn", @local_on_display_change);
        energy_max.Layout.Row = 3;
        energy_max.Layout.Column = 4;

        auto_y = uicheckbox(control_grid, ...
            "Text", "Auto Y", ...
            "Value", true, ...
            "ValueChangedFcn", @local_on_auto_y_changed);
        auto_y.Layout.Row = 3;
        auto_y.Layout.Column = 5;

        local_add_label(control_grid, 3, 6, "Y min");
        y_min = uieditfield(control_grid, "numeric", ...
            "Value", 0, ...
            "ValueDisplayFormat", "%.4g", ...
            "ValueChangedFcn", @local_on_display_change);
        y_min.Layout.Row = 3;
        y_min.Layout.Column = 7;

        local_add_label(control_grid, 3, 8, "Y max");
        y_max = uieditfield(control_grid, "numeric", ...
            "Value", 1, ...
            "ValueDisplayFormat", "%.4g", ...
            "ValueChangedFcn", @local_on_display_change);
        y_max.Layout.Row = 3;
        y_max.Layout.Column = 9;

        local_add_label(control_grid, 3, 10, "Y scale");
        y_scale = uidropdown(control_grid, ...
            "Items", {'linear', 'log'}, ...
            "Value", 'linear', ...
            "ValueChangedFcn", @local_on_display_change);
        y_scale.Layout.Row = 3;
        y_scale.Layout.Column = 11;

        % ═══════════ ROW 4: Display & Trace ═══════════
        local_add_label(control_grid, 4, 1, "Trace");
        trace_mode = uidropdown(control_grid, ...
            "Items", {'display', 'background-subtracted', 'curvature'}, ...
            "Value", 'display', ...
            "ValueChangedFcn", @local_on_display_change);
        trace_mode.Layout.Row = 4;
        trace_mode.Layout.Column = 2;

        local_add_label(control_grid, 4, 3, "Offset");
        offset = uieditfield(control_grid, "numeric", ...
            "Value", 0.25, ...
            "ValueDisplayFormat", "%.4g", ...
            "ValueChangedFcn", @local_on_display_change);
        offset.Layout.Row = 4;
        offset.Layout.Column = 4;

        local_add_label(control_grid, 4, 5, "Base λ");
        baseline_lambda = uieditfield(control_grid, "numeric", ...
            "Value", 1e6, ...
            "Limits", [1e2 inf], ...
            "ValueDisplayFormat", "%.3g", ...
            "ValueChangedFcn", @local_on_display_change);
        baseline_lambda.Layout.Row = 4;
        baseline_lambda.Layout.Column = 6;

        local_add_label(control_grid, 4, 7, "Smooth");
        smoothing_mode = uidropdown(control_grid, ...
            "Items", {'gaussian', 'off', 'sgolay'}, ...
            "Value", 'gaussian', ...
            "ValueChangedFcn", @local_on_display_change);
        smoothing_mode.Layout.Row = 4;
        smoothing_mode.Layout.Column = 8;

        local_add_label(control_grid, 4, 9, "Width");
        smoothing_width = uieditfield(control_grid, "numeric", ...
            "Value", 3, ...
            "Limits", [1 inf], ...
            "RoundFractionalValues", "on", ...
            "ValueDisplayFormat", "%.0f", ...
            "ValueChangedFcn", @local_on_display_change);
        smoothing_width.Layout.Row = 4;
        smoothing_width.Layout.Column = 10;

        local_add_label(control_grid, 4, 11, "SG ord");
        sg_order = uieditfield(control_grid, "numeric", ...
            "Value", 3, ...
            "Limits", [1 7], ...
            "ValueDisplayFormat", "%.0f", ...
            "ValueChangedFcn", @local_on_display_change);
        sg_order.Layout.Row = 4;
        sg_order.Layout.Column = 12;

        local_add_label(control_grid, 4, 13, "SG len");
        sg_framelen = uieditfield(control_grid, "numeric", ...
            "Value", 15, ...
            "Limits", [5 201], ...
            "ValueDisplayFormat", "%.0f", ...
            "Tooltip", "Frame length (must be odd)", ...
            "ValueChangedFcn", @local_on_display_change);
        sg_framelen.Layout.Row = 4;
        sg_framelen.Layout.Column = 14;

        % ═══════════ ROW 5: Normalization & BG ═══════════
        local_add_label(control_grid, 5, 1, "ZLP min");
        zlp_min = uieditfield(control_grid, "numeric", ...
            "Value", -20, ...
            "ValueDisplayFormat", "%.1f", ...
            "ValueChangedFcn", @local_on_normalization_change);
        zlp_min.Layout.Row = 5;
        zlp_min.Layout.Column = 2;

        local_add_label(control_grid, 5, 3, "ZLP max");
        zlp_max = uieditfield(control_grid, "numeric", ...
            "Value", 20, ...
            "ValueDisplayFormat", "%.1f", ...
            "ValueChangedFcn", @local_on_normalization_change);
        zlp_max.Layout.Row = 5;
        zlp_max.Layout.Column = 4;

        local_add_label(control_grid, 5, 5, "Ref|q|min");
        ref_abs_q_min = uieditfield(control_grid, "numeric", ...
            "Value", 0.010, ...
            "Limits", [eps inf], ...
            "LowerLimitInclusive", false, ...
            "ValueDisplayFormat", "%.5f", ...
            "ValueChangedFcn", @local_on_normalization_change);
        ref_abs_q_min.Layout.Row = 5;
        ref_abs_q_min.Layout.Column = 6;

        local_add_label(control_grid, 5, 7, "Ref|q|max");
        ref_abs_q_max = uieditfield(control_grid, "numeric", ...
            "Value", 0.020, ...
            "Limits", [eps inf], ...
            "LowerLimitInclusive", false, ...
            "ValueDisplayFormat", "%.5f", ...
            "ValueChangedFcn", @local_on_normalization_change);
        ref_abs_q_max.Layout.Row = 5;
        ref_abs_q_max.Layout.Column = 8;

        area_norm = uicheckbox(control_grid, ...
            "Text", "Normalize", ...
            "Value", false, ...
            "ValueChangedFcn", @local_on_display_change);
        area_norm.Layout.Row = 5;
        area_norm.Layout.Column = 9;

        norm_method = uidropdown(control_grid, ...
            "Items", {'Area', 'ZLP Peak'}, ...
            "Value", 'ZLP Peak', ...
            "Tooltip", "Area: divide by integrated spectral weight (destroys A(q)). ZLP Peak: divide by ZLP peak height (preserves A(q)).", ...
            "ValueChangedFcn", @local_on_norm_method_changed);
        norm_method.Layout.Row = 5;
        norm_method.Layout.Column = 10;

        area_norm_min = uieditfield(control_grid, "numeric", ...
            "Value", -50, ...
            "ValueDisplayFormat", "%.0f", ...
            "Tooltip", "ZLP Peak: ZLP window low bound (meV). Area: integration window low bound.", ...
            "ValueChangedFcn", @local_on_display_change);
        area_norm_min.Layout.Row = 5;
        area_norm_min.Layout.Column = 11;

        norm_hi_label = uilabel(control_grid, "Text", "~", "HorizontalAlignment", "center");
        norm_hi_label.Layout.Row = 5;
        norm_hi_label.Layout.Column = 12;
        area_norm_max = uieditfield(control_grid, "numeric", ...
            "Value", 50, ...
            "ValueDisplayFormat", "%.0f", ...
            "Tooltip", "ZLP Peak: ZLP window high bound (meV). Area: integration window high bound.", ...
            "ValueChangedFcn", @local_on_display_change);
        area_norm_max.Layout.Row = 5;
        area_norm_max.Layout.Column = 13;

        bg_sub = uicheckbox(control_grid, ...
            "Text", "BG Sub", ...
            "Value", false, ...
            "ValueChangedFcn", @local_on_display_change);
        bg_sub.Layout.Row = 5;
        bg_sub.Layout.Column = 14;

        bg_method = uidropdown(control_grid, ...
            "Items", {'Power', 'ExpPoly3', 'Pearson'}, ...
            "Value", 'Power', ...
            "ValueChangedFcn", @local_on_display_change);
        bg_method.Layout.Row = 5;
        bg_method.Layout.Column = 15;

        % ═══════════ ROW 6: Denoise & Deconv ═══════════
        denoise_cb = uicheckbox(control_grid, ...
            "Text", "Denoise", ...
            "Value", false, ...
            "ValueChangedFcn", @local_on_display_change);
        denoise_cb.Layout.Row = 6;
        denoise_cb.Layout.Column = 1;

        denoise_method = uidropdown(control_grid, ...
            "Items", {'Wiener2D', 'SavGol'}, ...
            "Value", 'Wiener2D', ...
            "ValueChangedFcn", @local_on_display_change);
        denoise_method.Layout.Row = 6;
        denoise_method.Layout.Column = 2;

        local_add_label(control_grid, 6, 3, "Noise σ");
        denoise_sigma = uieditfield(control_grid, "numeric", ...
            "Value", 0, ...
            "Limits", [0 inf], ...
            "ValueDisplayFormat", "%.1f", ...
            "Tooltip", "0 = auto-estimate from data", ...
            "ValueChangedFcn", @local_on_display_change);
        denoise_sigma.Layout.Row = 6;
        denoise_sigma.Layout.Column = 4;

        deconv_cb = uicheckbox(control_grid, ...
            "Text", "Deconv", ...
            "Value", false, ...
            "ValueChangedFcn", @local_on_display_change);
        deconv_cb.Layout.Row = 6;
        deconv_cb.Layout.Column = 5;

        local_add_label(control_grid, 6, 6, "LR iter");
        deconv_iter = uieditfield(control_grid, "numeric", ...
            "Value", 5, ...
            "Limits", [1 100], ...
            "ValueDisplayFormat", "%.0f", ...
            "ValueChangedFcn", @local_on_display_change);
        deconv_iter.Layout.Row = 6;
        deconv_iter.Layout.Column = 7;


        % ═══════════ ROW 7: Manual Peak Picking ═══════════
        pick_pts = uibutton(control_grid, ...
            "Text", "Pick Peaks", ...
            "ButtonPushedFcn", @local_on_pick_peaks);
        pick_pts.Layout.Row = 7;
        pick_pts.Layout.Column = 1;

        undo_pt = uibutton(control_grid, ...
            "Text", "Undo Pt", ...
            "ButtonPushedFcn", @local_on_undo_pt);
        undo_pt.Layout.Row = 7;
        undo_pt.Layout.Column = 2;

        clear_pts = uibutton(control_grid, ...
            "Text", "Clear Pts", ...
            "ButtonPushedFcn", @local_on_clear_pts);
        clear_pts.Layout.Row = 7;
        clear_pts.Layout.Column = 3;

        save_pts = uibutton(control_grid, ...
            "Text", "Save Pts", ...
            "ButtonPushedFcn", @local_on_save_pts);
        save_pts.Layout.Row = 7;
        save_pts.Layout.Column = 4;

        load_pts = uibutton(control_grid, ...
            "Text", "Load Pts", ...
            "ButtonPushedFcn", @local_on_load_pts);
        load_pts.Layout.Row = 7;
        load_pts.Layout.Column = 5;

        split_btn = uibutton(control_grid, ...
            "Text", "Split 3", ...
            "ButtonPushedFcn", @local_on_split_branches);
        split_btn.Layout.Row = 7;
        split_btn.Layout.Column = 6;

        fit_model_btn = uibutton(control_grid, ...
            "Text", "Fit Model", ...
            "ButtonPushedFcn", @local_on_fit_model);
        fit_model_btn.Layout.Row = 7;
        fit_model_btn.Layout.Column = 7;

        export_btn = uibutton(control_grid, ...
            "Text", "📷 Export", ...
            "ButtonPushedFcn", @local_on_export);
        export_btn.Layout.Row = 7;
        export_btn.Layout.Column = 8;

        auto_fit_btn = uibutton(control_grid, ...
            "Text", "Auto Fit ω(q)", ...
            "ButtonPushedFcn", @local_on_auto_fit_dispersion);
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
        prom_lbl = uilabel(control_grid, "Text", "Prominence:");
        prom_lbl.Layout.Row = 8; prom_lbl.Layout.Column = 1;

        prom_field = uispinner(control_grid, ...
            "Value", 0.15, "Step", 0.05, ...
            "Limits", [0.01 1], "ValueDisplayFormat", "%.2f");
        prom_field.Layout.Row = 8; prom_field.Layout.Column = 2;

        smooth_lbl = uilabel(control_grid, "Text", "Smooth W:");
        smooth_lbl.Layout.Row = 8; smooth_lbl.Layout.Column = 3;

        smooth_field = uispinner(control_grid, ...
            "Value", 25, "Step", 5, ...
            "Limits", [1 100], "RoundFractionalValues", "on");
        smooth_field.Layout.Row = 8; smooth_field.Layout.Column = 4;

        maxpk_lbl = uilabel(control_grid, "Text", "Max Peaks:");
        maxpk_lbl.Layout.Row = 8; maxpk_lbl.Layout.Column = 5;

        maxpk_field = uispinner(control_grid, ...
            "Value", 3, "Step", 1, ...
            "Limits", [1 10], "RoundFractionalValues", "on");
        maxpk_field.Layout.Row = 8; maxpk_field.Layout.Column = 6;

        fit_spec_btn = uibutton(control_grid, ...
            "Text", "Fit Spectrum", ...
            "ButtonPushedFcn", @local_on_fit_spectrum);
        fit_spec_btn.Layout.Row = 8; fit_spec_btn.Layout.Column = 7;

        accept_fit_btn = uibutton(control_grid, ...
            "Text", "Accept Fit", ...
            "Enable", "off", ...
            "ButtonPushedFcn", @local_on_accept_fit);
        accept_fit_btn.Layout.Row = 8; accept_fit_btn.Layout.Column = 8;

        fit_info_label = uilabel(control_grid, ...
            "Text", "", ...
            "HorizontalAlignment", "left");
        fit_info_label.Layout.Row = 8;
        fit_info_label.Layout.Column = [14 16];

        show_gamma_btn = uibutton(control_grid, ...
            "Text", "Γ & Amp", ...
            "ButtonPushedFcn", @local_on_show_gamma);
        show_gamma_btn.Layout.Row = 8; show_gamma_btn.Layout.Column = 13;

        guess_lbl = uilabel(control_grid, "Text", "Guesses:");
        guess_lbl.Layout.Row = 8; guess_lbl.Layout.Column = 9;

        guess_field = uieditfield(control_grid, 'text', ...
            "Value", "", ...
            "Placeholder", "e.g. 800,2500,3200");
        guess_field.Layout.Row = 8;
        guess_field.Layout.Column = [10 12];

        % ═══════════ ROW 9: Peak Model & Seed Propagation ═══════════
        pk_model_lbl = uilabel(control_grid, "Text", "Peak Model:");
        pk_model_lbl.Layout.Row = 9; pk_model_lbl.Layout.Column = 1;

        pk_model_dropdown = uidropdown(control_grid, ...
            "Items", {'Lorentz', 'Gaussian', 'Voigt', 'Damped HO'}, ...
            "ItemsData", {'lorentz', 'gaussian', 'voigt', 'damped_ho'}, ...
            "Value", 'lorentz');
        pk_model_dropdown.Layout.Row = 9; pk_model_dropdown.Layout.Column = [2 3];

        shift_lbl = uilabel(control_grid, "Text", "Max Shift:");
        shift_lbl.Layout.Row = 9; shift_lbl.Layout.Column = 4;

        max_shift_field = uispinner(control_grid, ...
            "Value", 80, "Step", 10, ...
            "Limits", [10 500], "RoundFractionalValues", "on");
        max_shift_field.Layout.Row = 9; max_shift_field.Layout.Column = 5;

        pick_guesses_btn = uibutton(control_grid, ...
            "Text", "Pick Guesses", ...
            "ButtonPushedFcn", @local_on_pick_guesses);
        pick_guesses_btn.Layout.Row = 9; pick_guesses_btn.Layout.Column = [6 7];

        reassign_btn = uibutton(control_grid, ...
            "Text", "Reassign Pts", ...
            "ButtonPushedFcn", @local_on_reassign_points);
        reassign_btn.Layout.Row = 9; reassign_btn.Layout.Column = [8 9];

        seed_info_label = uilabel(control_grid, ...
            "Text", "", ...
            "HorizontalAlignment", "left");
        seed_info_label.Layout.Row = 9;
        seed_info_label.Layout.Column = [10 16];

        % ═══════════ ROW 10: Dispersion Model Selector ═══════════
        disp_model_lbl = uilabel(control_grid, "Text", "Disp Model:");
        disp_model_lbl.Layout.Row = 10; disp_model_lbl.Layout.Column = 1;

        disp_model_dropdown = uidropdown(control_grid, ...
            "Items", {'Quasi-2D Plasmon', 'Acoustic Linear', 'Optical Constant', 'Optical Quadratic'}, ...
            "ItemsData", {'quasi2d_plasmon', 'acoustic_linear', 'optical_constant', 'optical_quadratic'}, ...
            "Value", 'quasi2d_plasmon');
        disp_model_dropdown.Layout.Row = 10; disp_model_dropdown.Layout.Column = [2 4];

        fit_disp_btn = uibutton(control_grid, ...
            "Text", "Fit Dispersion", ...
            "ButtonPushedFcn", @local_on_fit_dispersion);
        fit_disp_btn.Layout.Row = 10; fit_disp_btn.Layout.Column = [5 6];

        export_disp_btn = uibutton(control_grid, ...
            "Text", "Export Disp.", ...
            "ButtonPushedFcn", @local_on_export_dispersion);
        export_disp_btn.Layout.Row = 10; export_disp_btn.Layout.Column = [7 8];

        disp_info_label = uilabel(control_grid, ...
            "Text", "", ...
            "HorizontalAlignment", "left");
        disp_info_label.Layout.Row = 10;
        disp_info_label.Layout.Column = [9 16];

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

        dispersion_axes = uiaxes(main_grid);
        dispersion_axes.Layout.Row = 3;
        dispersion_axes.Layout.Column = 2;
        title(dispersion_axes, "Dispersion");
        xlabel(dispersion_axes, "q (1/A)");
        ylabel(dispersion_axes, "Peak energy (meV)");

        local_attach_toolbar(qe_axes);
        local_attach_toolbar(comparison_axes);
        local_attach_toolbar(single_axes);
        local_attach_toolbar(dispersion_axes);

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
        ui_handles.QEAxes = qe_axes;
        ui_handles.ComparisonAxes = comparison_axes;
        ui_handles.SingleAxes = single_axes;
        ui_handles.WaterfallAxes = gobjects(1);
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
        ui_handles.DispInfoLabel = disp_info_label;
    end


    function local_add_label(parent, row, column, text_value)
        label = uilabel(parent, "Text", text_value, "HorizontalAlignment", "left");
        label.Layout.Row = row;
        label.Layout.Column = column;
    end


    function local_attach_toolbar(ax)
        axtoolbar(ax, {'zoomin', 'zoomout', 'pan', 'datacursor', 'restoreview'});
    end


    function local_load_start_path(path_to_load)
        local_load_dataset(path_to_load, true);
    end


    function local_on_load_data(~, ~)
        start_folder = local_get_dialog_folder();
        [file_name, folder_name] = uigetfile( ...
            {"*.mat;*.npy", "EELS Data (*.mat, *.npy)"; ...
             "*.mat", "MAT files (*.mat)"; ...
             "*.npy", "Raw Nion NPY (*.npy)"; ...
             "*.*", "All files"}, ...
            "Select EELS Data File", ...
            start_folder);
        if isequal(file_name, 0)
            return
        end
        full_path = fullfile(folder_name, file_name);
        [~, ~, ext] = fileparts(file_name);

        % Determine if this is raw data needing preprocessing
        is_raw = false;
        if strcmpi(ext, '.npy')
            is_raw = true;
        elseif strcmpi(ext, '.mat')
            % Check if it's 4D raw or eq3D
            try
                info = whos('-file', full_path);
                names = {info.name};
                classes = {info.class};
                % EELS4D object → raw
                if any(strcmp(classes, 'EELS4D'))
                    is_raw = true;
                % Has a3+e with 2D a3 → eq3D
                elseif ismember('a3', names) && ismember('e', names)
                    a3_idx = find(strcmp(names, 'a3'), 1);
                    if numel(info(a3_idx).size) >= 3
                        is_raw = true;  % 3D/4D a3
                    end
                % Any 3D+ array > 1 MB → raw
                else
                    for ii = 1:numel(info)
                        if numel(info(ii).size) >= 3 && info(ii).bytes > 1e6
                            is_raw = true;
                            break
                        end
                    end
                end
            catch
            end
        end

        if is_raw
            % Ask for q crop range
            answer = inputdlg( ...
                {'Q crop start (pixel, leave 0 for auto):', ...
                 'Q crop end (pixel, leave 0 for auto):'}, ...
                'Raw Import — Q Crop Range', [1 45], {'0', '0'});
            if isempty(answer)
                return
            end
            q_lo = str2double(answer{1});
            q_hi = str2double(answer{2});
            if ~isfinite(q_lo) || q_lo <= 0; q_lo = NaN; end
            if ~isfinite(q_hi) || q_hi <= 0; q_hi = NaN; end

            ui.InfoLabel.Text = "Loading raw data... (this may take a few seconds)";
            ui.InfoLabel.Visible = "on";
            drawnow;

            try
                lower_path = lower(char(full_path));
                if contains(lower_path, "20w")
                    ui.DqOverrideField.Value = 0.0025;
                elseif contains(lower_path, "10w")
                    ui.DqOverrideField.Value = 0.005;
                end
                dataset = load_qe_dataset(full_path, ui.DqOverrideField.Value, ...
                    q_crop=[q_lo q_hi]);
                state.dataset = dataset;
                state.selectedQIndex = 1;
                state.eq3dQE = [];
                state.physicalQE = [];
                state.comparisonQE = [];
                state.eq3dQE = dataset.qe;
                state.physicalQE = dataset.qe;
                state.activeStage = "eq3d";
                local_rebuild_comparison(false);
                ui.ViewModeDropdown.Value = 'Physical';
                local_sync_controls_from_qe(local_get_reference_qe());
                [~, fname] = fileparts(char(full_path));
                local_log_operation(sprintf('Loaded raw: %s (dq=%.4f)', fname, dataset.dq_Ainv));
                local_update_all_views();
            catch ME
                local_show_error(ME, false);
            end
        else
            % .mat path — original behavior
            local_load_dataset(full_path, false);
        end
    end


    function local_load_dataset(path_to_load, throw_on_error)
        try
            % Auto-detect dq from path and update the override field
            lower_path = lower(char(path_to_load));
            if contains(lower_path, "20w")
                ui.DqOverrideField.Value = 0.0025;
            elseif contains(lower_path, "10w")
                ui.DqOverrideField.Value = 0.005;
            end
            dataset = load_qe_dataset(path_to_load, ui.DqOverrideField.Value);
            state.dataset = dataset;
            state.selectedQIndex = 1;
            state.eq3dQE = [];
            state.physicalQE = [];
            state.comparisonQE = [];

            state.eq3dQE = dataset.qe;
            state.physicalQE = dataset.qe;
            state.activeStage = "eq3d";
            local_rebuild_comparison(false);
            ui.ViewModeDropdown.Value = 'Physical';

            local_sync_controls_from_qe(local_get_reference_qe());

            % Log to operation history
            [~, fname] = fileparts(char(path_to_load));
            local_log_operation(sprintf('Loaded: %s (dq=%.4f)', fname, dataset.dq_Ainv));
            
            local_update_all_views();
        catch ME
            local_show_error(ME, throw_on_error);
        end
    end


    function folder_name = local_get_dialog_folder()
        if ~isempty(state.dataset) && isfield(state.dataset, "source_path") ...
                && strlength(string(state.dataset.source_path)) > 0
            folder_name = fileparts(char(state.dataset.source_path));
            if isfolder(folder_name)
                return
            end
        end
        folder_name = pwd;
    end



    function local_on_build_views(~, ~)
        try
            local_rebuild_comparison(true);
            local_update_all_views();
        catch ME
            local_show_error(ME, false);
        end
    end


    function local_on_reset_stage(~, ~)
        target_stage = string(ui.ResetStageDropdown.Value);
        try
            switch target_stage
                case "eq3d"
                    if isempty(state.eq3dQE)
                        error("interactive_qe_browser:MissingEq3DStage", ...
                            "eq3D stage is not available.");
                    end
                    state.physicalQE = state.eq3dQE;
                    state.comparisonQE = [];
                    state.activeStage = "eq3d";

                case "raw"
                    if isempty(state.raw4d)
                        error("interactive_qe_browser:MissingRawStage", ...
                            "Raw stage is not available.");
                    end
                    state.cropSpec = [];
                    state.cropped4d = [];
                    state.aligned4d = [];
                    state.binned4d = [];
                    state.physicalQE = project_4d_to_qe_struct( ...
                        state.raw4d, ...
                        local_current_dq(), ...
                        source_path=string(state.dataset.source_path), ...
                        source_kind="raw", ...
                        label=string(state.dataset.label), ...
                        view_kind="physical", ...
                        stage_name="raw");
                    state.comparisonQE = [];
                    state.activeStage = "raw";

                case "crop"
                    if isempty(state.cropped4d)
                        error("interactive_qe_browser:MissingCropStage", ...
                            "Crop stage is not available.");
                    end
                    state.aligned4d = [];
                    state.binned4d = [];
                    state.physicalQE = project_4d_to_qe_struct( ...
                        state.cropped4d, ...
                        local_current_dq(), ...
                        source_path=string(state.dataset.source_path), ...
                        source_kind="raw", ...
                        label=string(state.dataset.label), ...
                        view_kind="physical", ...
                        stage_name="crop");
                    state.comparisonQE = [];
                    state.activeStage = "crop";

                case "align"
                    if isempty(state.aligned4d)
                        error("interactive_qe_browser:MissingAlignStage", ...
                            "Align stage is not available.");
                    end
                    state.binned4d = [];
                    state.physicalQE = project_4d_to_qe_struct( ...
                        state.aligned4d, ...
                        local_current_dq(), ...
                        source_path=string(state.dataset.source_path), ...
                        source_kind="raw", ...
                        label=string(state.dataset.label), ...
                        view_kind="physical", ...
                        stage_name="align");
                    state.comparisonQE = [];
                    state.activeStage = "align";

                case "bin"
                    if isempty(state.binned4d)
                        error("interactive_qe_browser:MissingBinStage", ...
                            "Bin stage is not available.");
                    end
                    state.physicalQE = project_4d_to_qe_struct( ...
                        state.binned4d, ...
                        local_current_dq(), ...
                        source_path=string(state.dataset.source_path), ...
                        source_kind="raw", ...
                        label=string(state.dataset.label), ...
                        view_kind="physical", ...
                        stage_name="bin");
                    state.comparisonQE = [];
                    state.activeStage = "bin";

                otherwise
                    error("interactive_qe_browser:UnknownResetStage", ...
                        "Unknown reset stage '%s'.", target_stage);
            end

            ui.ViewModeDropdown.Value = 'Physical';
            state.selectedQIndex = 1;
            
            local_sync_controls_from_qe(local_get_reference_qe());
            
            local_update_all_views();
        catch ME
            local_show_error(ME, false);
        end
    end


    function local_push_stage_history()
        % Push current stage onto history stack before a destructive operation
        state.stageHistory{end+1} = char(state.activeStage);
    end


    function local_log_operation(label)
        % Log an operation to the history with a full UI snapshot
        snap = local_capture_snapshot();

        entry = struct();
        entry.label = char(label);
        entry.timestamp = char(datetime('now', 'Format', 'HH:mm:ss'));
        entry.snapshot = snap;
        entry.stage = char(state.activeStage);

        state.opHistory{end+1} = entry;

        % Cap at 50 entries
        if numel(state.opHistory) > 50
            state.opHistory(1) = [];
        end

        local_refresh_history_list();
    end


    function local_push_param_snapshot()
        % Backwards-compatible wrapper — logs as "Param change"
        local_log_operation("Param change");
    end


    function snap = local_capture_snapshot()
        % Capture all UI parameter values into a snapshot struct
        snap = struct();
        snap.qStart = ui.QStartField.Value;
        snap.qEnd = ui.QEndField.Value;
        snap.qStep = ui.QStepField.Value;
        snap.dqOverride = ui.DqOverrideField.Value;
        snap.energyMin = ui.EnergyMinField.Value;
        snap.energyMax = ui.EnergyMaxField.Value;
        snap.zlpMin = ui.ZlpMinField.Value;
        snap.zlpMax = ui.ZlpMaxField.Value;
        snap.refQMin = ui.RefAbsQMinField.Value;
        snap.refQMax = ui.RefAbsQMaxField.Value;
        snap.smoothMode = char(ui.SmoothingModeDropdown.Value);
        snap.smoothWidth = ui.SmoothingWidthField.Value;
        snap.sgOrder = ui.SGOrderField.Value;
        snap.sgFrameLen = ui.SGFrameLenField.Value;
        snap.yScale = char(ui.YScaleDropdown.Value);
        snap.autoY = ui.AutoYCheckbox.Value;
        snap.yMin = ui.YMinField.Value;
        snap.yMax = ui.YMaxField.Value;
        snap.offset = ui.OffsetField.Value;
        snap.traceMode = char(ui.TraceModeDropdown.Value);
        snap.baseLambda = ui.BaselineLambdaField.Value;
        snap.areaNorm = ui.AreaNormCheckbox.Value;
        snap.areaNormMin = ui.AreaNormMinField.Value;
        snap.areaNormMax = ui.AreaNormMaxField.Value;
        snap.bgSub = ui.BgSubCheckbox.Value;
        snap.bgMethod = char(ui.BgMethodDropdown.Value);
        snap.denoise = ui.DenoiseCheckbox.Value;
        snap.denoiseMethod = char(ui.DenoiseMethodDropdown.Value);
        snap.denoiseSigma = ui.DenoiseSigmaField.Value;
        snap.deconv = ui.DeconvCheckbox.Value;
        snap.deconvIter = ui.DeconvIterField.Value;

        snap.viewMode = char(ui.ViewModeDropdown.Value);
        snap.activeStage = char(state.activeStage);
    end


    function local_restore_param_snapshot(snap)
        % Restore all UI controls from a snapshot struct
        ui.QStartField.Value = snap.qStart;
        ui.QEndField.Value = snap.qEnd;
        ui.QStepField.Value = snap.qStep;
        ui.DqOverrideField.Value = snap.dqOverride;
        ui.EnergyMinField.Value = snap.energyMin;
        ui.EnergyMaxField.Value = snap.energyMax;
        ui.ZlpMinField.Value = snap.zlpMin;
        ui.ZlpMaxField.Value = snap.zlpMax;
        ui.RefAbsQMinField.Value = snap.refQMin;
        ui.RefAbsQMaxField.Value = snap.refQMax;
        ui.SmoothingModeDropdown.Value = snap.smoothMode;
        ui.SmoothingWidthField.Value = snap.smoothWidth;
        ui.SGOrderField.Value = snap.sgOrder;
        ui.SGFrameLenField.Value = snap.sgFrameLen;
        ui.YScaleDropdown.Value = snap.yScale;
        ui.AutoYCheckbox.Value = snap.autoY;
        ui.YMinField.Value = snap.yMin;
        ui.YMaxField.Value = snap.yMax;
        ui.OffsetField.Value = snap.offset;
        ui.TraceModeDropdown.Value = snap.traceMode;
        ui.BaselineLambdaField.Value = snap.baseLambda;
        ui.AreaNormCheckbox.Value = snap.areaNorm;
        ui.AreaNormMinField.Value = snap.areaNormMin;
        ui.AreaNormMaxField.Value = snap.areaNormMax;
        ui.BgSubCheckbox.Value = snap.bgSub;
        ui.BgMethodDropdown.Value = snap.bgMethod;
        ui.DenoiseCheckbox.Value = snap.denoise;
        ui.DenoiseMethodDropdown.Value = snap.denoiseMethod;
        ui.DenoiseSigmaField.Value = snap.denoiseSigma;
        ui.DeconvCheckbox.Value = snap.deconv;
        ui.DeconvIterField.Value = snap.deconvIter;

        ui.ViewModeDropdown.Value = snap.viewMode;
    end


    function local_refresh_history_list()
        % Update the history listbox items from opHistory
        n = numel(state.opHistory);
        items = cell(1, n);
        itemData = cell(1, n);
        for k = 1:n
            e = state.opHistory{k};
            items{k} = sprintf('%d. [%s] %s', k, e.timestamp, e.label);
            itemData{k} = k;
        end
        ui.HistoryListBox.Items = items;
        ui.HistoryListBox.ItemsData = itemData;
        if n > 0
            ui.UndoButton.Enable = "on";
            % Scroll to bottom (select last item briefly)
            ui.HistoryListBox.Value = n;
        else
            ui.UndoButton.Enable = "off";
        end
    end


    function local_on_history_select(~, evt)
        % Roll back to the selected history entry
        if isempty(evt.Value) || isempty(state.opHistory)
            return
        end
        sel_idx = evt.Value;
        if sel_idx < 1 || sel_idx > numel(state.opHistory)
            return
        end

        entry = state.opHistory{sel_idx};

        % Trim all entries AFTER the selected one
        state.opHistory = state.opHistory(1:sel_idx);

        try
            local_restore_param_snapshot(entry.snapshot);
            local_rebuild_comparison(false);
            local_refresh_history_list();
            local_update_all_views();
        catch ME
            local_show_error(ME, false);
        end
    end


    function local_on_clear_history(~, ~)
        state.opHistory = {};
        local_refresh_history_list();
    end


    function local_on_save_history(~, ~)
        % Save opHistory to a .mat file in the data folder
        if isempty(state.opHistory)
            ui.InfoLabel.Text = "No history to save";
            ui.InfoLabel.Visible = "on";
            return
        end

        % Default save location: same folder as loaded data
        default_dir = pwd;
        if ~isempty(state.dataset) && isfield(state.dataset, 'source_path') ...
                && strlength(string(state.dataset.source_path)) > 0
            default_dir = fileparts(char(state.dataset.source_path));
        end
        default_file = fullfile(default_dir, 'op_history.mat');

        [fname, fpath] = uiputfile('*.mat', 'Save Operation History', default_file);
        if isequal(fname, 0)
            return
        end

        save_path = fullfile(fpath, fname);
        opHistory = state.opHistory; %#ok<NASGU>
        save(save_path, 'opHistory');
        ui.InfoLabel.Text = sprintf('History saved → %s', fname);
        ui.InfoLabel.Visible = "on";
        fprintf('  History saved to: %s\n', save_path);
    end


    function local_on_load_history(~, ~)
        % Load opHistory from a .mat file and restore the last state
        default_dir = pwd;
        if ~isempty(state.dataset) && isfield(state.dataset, 'source_path') ...
                && strlength(string(state.dataset.source_path)) > 0
            default_dir = fileparts(char(state.dataset.source_path));
        end

        [fname, fpath] = uigetfile('*.mat', 'Load Operation History', ...
            fullfile(default_dir, 'op_history.mat'));
        if isequal(fname, 0)
            return
        end

        load_path = fullfile(fpath, fname);
        loaded = load(load_path);
        if ~isfield(loaded, 'opHistory')
            ui.InfoLabel.Text = "File does not contain opHistory";
            ui.InfoLabel.Visible = "on";
            return
        end

        state.opHistory = loaded.opHistory;
        local_refresh_history_list();

        % Restore the last snapshot if available
        if ~isempty(state.opHistory)
            last_entry = state.opHistory{end};
            try
                local_restore_param_snapshot(last_entry.snapshot);
                local_rebuild_comparison(false);
                local_update_all_views();
            catch ME
                fprintf('  Warning: could not fully restore state: %s\n', ME.message);
            end
        end

        ui.InfoLabel.Text = sprintf('History loaded ← %s (%d entries)', fname, numel(state.opHistory));
        ui.InfoLabel.Visible = "on";
        fprintf('  History loaded from: %s (%d entries)\n', load_path, numel(state.opHistory));
    end



    function local_on_undo(~, ~)
        % Quick undo: roll back to the entry BEFORE the last one
        if numel(state.opHistory) < 2
            return
        end

        target_idx = numel(state.opHistory) - 1;
        entry = state.opHistory{target_idx};
        state.opHistory = state.opHistory(1:target_idx);

        try
            local_restore_param_snapshot(entry.snapshot);
            local_rebuild_comparison(false);
            local_refresh_history_list();
            local_update_all_views();
        catch ME
            local_show_error(ME, false);
        end
    end


    function label = local_describe_changes()
        % Compare current UI values against last snapshot to generate
        % a specific human-readable label describing what changed
        current = local_capture_snapshot();

        if isempty(state.opHistory)
            label = "Initial state";
            return
        end

        prev = state.opHistory{end}.snapshot;

        % Field → short display name mapping
        names = struct( ...
            'qStart', 'q start', ...
            'qEnd', 'q end', ...
            'qStep', 'q step', ...
            'dqOverride', 'dq', ...
            'energyMin', 'E min', ...
            'energyMax', 'E max', ...
            'zlpMin', 'ZLP min', ...
            'zlpMax', 'ZLP max', ...
            'refQMin', 'Ref q min', ...
            'refQMax', 'Ref q max', ...
            'smoothMode', 'Smooth', ...
            'smoothWidth', 'Smooth W', ...
            'sgOrder', 'SG order', ...
            'sgFrameLen', 'SG frame', ...
            'yScale', 'Y scale', ...
            'autoY', 'Auto-Y', ...
            'yMin', 'Y min', ...
            'yMax', 'Y max', ...
            'offset', 'Offset', ...
            'traceMode', 'Trace', ...
            'baseLambda', 'Baseline λ', ...
            'areaNorm', 'Area norm', ...
            'areaNormMin', 'Norm min', ...
            'areaNormMax', 'Norm max', ...
            'bgSub', 'BG sub', ...
            'bgMethod', 'BG method', ...
            'denoise', 'Denoise', ...
            'denoiseMethod', 'Denoise method', ...
            'denoiseSigma', 'Denoise σ', ...
            'deconv', 'Deconv', ...
            'deconvIter', 'Deconv iter', ...
            'searchHW', 'Search HW', ...
            'minConf', 'Min conf', ...
            'viewMode', 'View');

        fields = fieldnames(names);
        parts = {};
        for k = 1:numel(fields)
            f = fields{k};
            if ~isfield(prev, f) || ~isfield(current, f)
                continue
            end
            old_val = prev.(f);
            new_val = current.(f);
            if isequal(old_val, new_val)
                continue
            end
            short = names.(f);
            % Format based on type
            if islogical(old_val) || (isnumeric(old_val) && (old_val==0||old_val==1) && (new_val==0||new_val==1) && contains(f, {'denoise','deconv','bgSub','areaNorm','autoY'}))
                if new_val
                    parts{end+1} = sprintf('%s ON', short); %#ok<AGROW>
                else
                    parts{end+1} = sprintf('%s OFF', short); %#ok<AGROW>
                end
            elseif ischar(old_val) || isstring(old_val)
                parts{end+1} = sprintf('%s: %s', short, char(new_val)); %#ok<AGROW>
            else
                parts{end+1} = sprintf('%s: %g→%g', short, old_val, new_val); %#ok<AGROW>
            end
        end

        if isempty(parts)
            label = "No change";
        elseif numel(parts) <= 3
            label = strjoin(parts, ', ');
        else
            label = sprintf('%s (+%d more)', strjoin(parts(1:2), ', '), numel(parts)-2);
        end
    end


    function local_on_dq_override_changed(~, ~)
        if isempty(state.dataset)
            return
        end

        dq_Ainv = local_current_dq();
        state.dataset.dq_Ainv = dq_Ainv;
        state.eq3dQE = local_reaxis_qe(state.eq3dQE, dq_Ainv);
        state.physicalQE = local_reaxis_qe(state.physicalQE, dq_Ainv);
        state.comparisonQE = local_reaxis_qe(state.comparisonQE, dq_Ainv);

        local_sync_controls_from_qe(local_get_reference_qe());
        
        local_update_all_views();
    end


    function local_on_qnorm_changed(~, ~)
        if ~local_has_raw_source()
            ui.QNormCheckbox.Value = false;
            return
        end
        try
            local_log_operation(local_describe_changes());
            local_rebuild_comparison(true);
            local_update_all_views();
        catch ME
            local_show_error(ME, false);
        end
    end


    function local_on_normalization_change(~, ~)
        try
            local_log_operation(local_describe_changes());
            if ~isempty(state.dataset)
                local_rebuild_comparison(local_has_raw_source() && ui.QNormCheckbox.Value);
            end
            local_update_all_views();
        catch ME
            local_show_error(ME, false);
        end
    end


    function local_on_display_change(~, ~)
        local_log_operation(local_describe_changes());
        local_sync_manual_y_state();
        local_update_all_views();
    end


    function local_on_norm_method_changed(~, ~)
        % Update lo/hi field values and labels when normalization mode changes
        method = char(ui.NormMethodDropdown.Value);
        if strcmpi(method, 'ZLP Peak')
            ui.AreaNormMinField.Value = -50;
            ui.AreaNormMaxField.Value = 50;
            ui.NormHiLabel.Text = "~";
        else
            ui.AreaNormMinField.Value = 200;
            ui.AreaNormMaxField.Value = 3800;
            ui.NormHiLabel.Text = "~";
        end
        local_on_display_change([], []);
    end


    function local_on_auto_y_changed(~, ~)
        local_log_operation(local_describe_changes());
        local_sync_manual_y_state();
        local_update_all_views();
    end



    function local_on_pick_peaks(~, ~)
        if isempty(state.physicalQE)
            return
        end
        state.manualPickArmed = true;
        ui.PickPtsButton.Text = "Picking...";
        ui.PickPtsButton.BackgroundColor = [0.85 1 0.85];
    end


    function local_on_undo_pt(~, ~)
        if ~isempty(state.manual_points)
            state.manual_points(end, :) = [];
            ui.PtsLabel.Text = sprintf("%d pts", size(state.manual_points, 1));
            local_update_all_views();
        end
    end


    function local_on_clear_pts(~, ~)
        state.manual_points = [];
        state.manual_branches = {};
        ui.PtsLabel.Text = "0 pts";
        local_update_all_views();
    end


    function local_on_save_pts(~, ~)
        if isempty(state.manual_points)
            return
        end
        [file, path] = uiputfile('*.mat', 'Save dispersion points', 'dispersion_points.mat');
        if isequal(file, 0)
            return
        end
        manual_dispersion = struct();
        manual_dispersion.abs_q_Ainv = state.manual_points(:, 1);
        manual_dispersion.energy_meV = state.manual_points(:, 2);
        manual_dispersion.timestamp = datestr(now);
        save(fullfile(path, file), 'manual_dispersion');
        ui.PtsLabel.Text = sprintf("%d pts (saved)", size(state.manual_points, 1));
    end


    function local_on_load_pts(~, ~)
        [file, path] = uigetfile('*.mat', 'Load dispersion points');
        if isequal(file, 0)
            return
        end
        loaded = load(fullfile(path, file));
        if isfield(loaded, 'manual_dispersion')
            md = loaded.manual_dispersion;
            state.manual_points = [md.abs_q_Ainv(:), md.energy_meV(:)];
        else
            warning('File does not contain manual_dispersion data.');
            return
        end
        state.manual_branches = {};
        ui.PtsLabel.Text = sprintf("%d pts (loaded)", size(state.manual_points, 1));
        local_update_all_views();
    end


    function local_on_split_branches(~, ~)
        if isempty(state.manual_points) || size(state.manual_points, 1) < 3
            return
        end
        pts = state.manual_points;
        energies = pts(:, 2);
        n_branches = 3;

        % k-means on energy to cluster branches
        % Manual k-means to avoid Statistics Toolbox dependency
        % Initialize centroids by evenly spacing in energy range
        e_min = min(energies);
        e_max = max(energies);
        centroids = linspace(e_min, e_max, n_branches)';

        for iter = 1:50
            % Assign each point to nearest centroid
            dists = abs(energies - centroids');
            [~, idx] = min(dists, [], 2);
            % Update centroids
            new_centroids = centroids;
            for k = 1:n_branches
                members = energies(idx == k);
                if ~isempty(members)
                    new_centroids(k) = mean(members);
                end
            end
            if max(abs(new_centroids - centroids)) < 1
                centroids = new_centroids;
                break
            end
            centroids = new_centroids;
        end

        % Final assignment
        dists = abs(energies - centroids');
        [~, idx] = min(dists, [], 2);

        % Sort centroids by energy (branch 1 = lowest energy)
        [~, sort_order] = sort(centroids);
        reorder = zeros(n_branches, 1);
        for k = 1:n_branches
            reorder(sort_order(k)) = k;
        end
        idx = reorder(idx);

        % Build branch arrays, sorted by |q|
        state.manual_branches = cell(1, n_branches);
        for k = 1:n_branches
            branch_pts = pts(idx == k, :);
            state.manual_branches{k} = sortrows(branch_pts, 1);
        end

        ui.PtsLabel.Text = sprintf("%d pts → %d branches", size(pts, 1), n_branches);
        local_log_operation(sprintf('Split → %d branches', n_branches));
        local_update_all_views();
    end


    function local_on_fit_model(~, ~)
        % Fit the da Jornada quasi-2D plasmon model to manual branches
        if isempty(state.manual_branches) && ...
                (isempty(state.manual_points) || size(state.manual_points, 1) < 3)
            ui.InfoLabel.Text = "Need ≥3 points or branches to fit";
            ui.InfoLabel.Visible = "on";
            return
        end

        % Collect branches to fit
        if ~isempty(state.manual_branches)
            branches = state.manual_branches;
        else
            branches = {state.manual_points};
        end

        ax = ui.DispersionAxes;
        hold(ax, 'on');

        % Color palette for branches
        colors = lines(numel(branches));
        fit_results = cell(1, numel(branches));
        info_parts = {};

        for b = 1:numel(branches)
            pts = branches{b};
            if isempty(pts) || size(pts, 1) < 2
                continue
            end

            q_vals = abs(pts(:, 1));
            e_vals = pts(:, 2);

            try
                result = fit_quasi2d_plasmon(q_vals, e_vals);
                fit_results{b} = result;

                % Plot fitted curve on dispersion axes
                plot(ax, result.q_fit, result.E_fit, '-', ...
                    'Color', colors(b, :), ...
                    'LineWidth', 2, ...
                    'DisplayName', sprintf( ...
                        'B%d: ρ₀=%.1fÅ  E_{flat}=%.0f meV', ...
                        b, result.rho0, result.E_flat_meV));

                info_parts{end+1} = sprintf( ...
                    'B%d: ρ₀=%.1fÅ  E_f=%.0fmeV  R²=%.3f', ...
                    b, result.rho0, result.E_flat_meV, result.R_squared); %#ok<AGROW>
            catch ME
                info_parts{end+1} = sprintf('B%d: fit failed (%s)', b, ME.message); %#ok<AGROW>
            end
        end

        hold(ax, 'off');

        % Update info label
        if ~isempty(info_parts)
            ui.InfoLabel.Text = strjoin(info_parts, '  |  ');
            ui.InfoLabel.Visible = "on";
        end

        % Store fit results in state for later use
        state.fitResults = fit_results;
        local_log_operation(sprintf('Fit model (%d branches)', numel(fit_results)));
    end


    function local_on_fit_spectrum(~, ~)
        % Fit the currently displayed single spectrum with Drude-Lorentz
        if isempty(state.physicalQE) && isempty(state.comparisonQE)
            ui.FitInfoLabel.Text = "Load data first";
            return
        end

        qe = local_get_active_qe();
        if isempty(qe); return; end

        % Apply same preprocessing as single spectrum display
        if ui.AreaNormCheckbox.Value
            qe = local_apply_area_normalization(qe);
        end
        if ui.DenoiseCheckbox.Value
            qe = local_apply_denoise(qe);
        end
        if ui.BgSubCheckbox.Value
            qe = local_apply_bg_subtraction(qe);
        end
        if ui.DeconvCheckbox.Value
            qe = local_apply_deconvolution(qe);
        end

        [mask, energy_axis] = local_energy_mask(qe);
        q_index = state.selectedQIndex;
        abs_q = abs(qe.q_Ainv(q_index));
        spectrum = double(qe.intensity(mask, q_index));

        if all(spectrum == 0) || all(isnan(spectrum))
            ui.FitInfoLabel.Text = "Empty spectrum";
            return
        end

        % Fit using current UI parameters
        E_min = max(ui.EnergyMinField.Value, 50);
        E_max = ui.EnergyMaxField.Value;

        % Parse manual peak guesses (comma-separated meV values)
        guess_str = strtrim(ui.GuessField.Value);
        if ~isempty(guess_str)
            guesses = str2double(strsplit(guess_str, {',', ' ', ';'}));
            guesses = guesses(~isnan(guesses) & guesses > 0);
        else
            guesses = [];
        end

        try
            result = fit_loss_function(energy_axis, spectrum, ...
                'E_min', E_min, 'E_max', E_max, ...
                'min_prominence', ui.PromField.Value, ...
                'smooth_width', ui.SmoothField.Value, ...
                'max_peaks', ui.MaxPeaksField.Value, ...
                'initial_guesses', guesses, ...
                'peak_model', ui.PeakModelDropdown.Value, ...
                'pre_subtracted', ui.BgSubCheckbox.Value);
        catch ME
            ui.FitInfoLabel.Text = sprintf("Fit failed: %s", ME.message);
            return
        end

        % Filter out peaks below E_min
        keep = result.omega_p >= E_min;
        if ~any(keep)
            ui.FitInfoLabel.Text = "No valid peaks found above E_min";
            return
        end

        % Store pending fit
        state.pendingFit = struct();
        state.pendingFit.abs_q = abs_q;
        state.pendingFit.q_index = q_index;
        state.pendingFit.result = result;
        state.pendingFit.keep = keep;

        % Overlay fit on spectrum axes
        ax = ui.SingleAxes;
        hold(ax, 'on');

        % Delete any previous fit overlay
        delete(findobj(ax, 'Tag', 'lorentz_fit'));

        % Total fit curve
        plot(ax, result.energy_fit, result.curve_fit, 'r-', ...
            'LineWidth', 2, 'Tag', 'lorentz_fit', ...
            'DisplayName', 'Lorentz fit');

        % Background curve (power-law)
        if isfield(result, 'bg_curve')
            plot(ax, result.energy_fit, result.bg_curve, ':', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
                'Tag', 'lorentz_fit', 'DisplayName', 'Background');
        end

        % Individual peak components (plotted ON TOP of background)
        peak_colors = [0.1 0.7 0.3; 0.9 0.5 0.1; 0.5 0.2 0.8];
        bg_vals = zeros(size(result.energy_fit));
        if isfield(result, 'bg_curve')
            bg_vals = result.bg_curve;
        end

        info_parts = {};
        for p = 1:result.n_peaks
            if ~keep(p); continue; end
            col = peak_colors(mod(p-1, size(peak_colors,1))+1, :);

            % Plot peak component + background so it's visible
            peak_on_bg = result.peak_curves{p} + bg_vals;
            plot(ax, result.energy_fit, peak_on_bg, '--', ...
                'Color', col, 'LineWidth', 1.5, 'Tag', 'lorentz_fit');

            % Peak marker at peak position on total fit
            [~, pk_idx] = min(abs(result.energy_fit - result.omega_p(p)));
            marker_y = result.curve_fit(pk_idx);
            plot(ax, result.omega_p(p), marker_y, 'v', ...
                'MarkerSize', 10, 'MarkerFaceColor', col, ...
                'MarkerEdgeColor', 'k', 'Tag', 'lorentz_fit');

            text(ax, result.omega_p(p), marker_y * 1.08, ...
                sprintf('%.0f meV', result.omega_p(p)), ...
                'FontSize', 9, 'FontWeight', 'bold', 'Color', col, ...
                'HorizontalAlignment', 'center', 'Tag', 'lorentz_fit');

            info_parts{end+1} = sprintf('ω_p=%.0f Γ=%.0f', ...
                result.omega_p(p), result.gamma(p)); %#ok<AGROW>
        end
        hold(ax, 'off');

        % Update fit info
        ui.FitInfoLabel.Text = sprintf('q=%.4f | %s | R²=%.3f', ...
            abs_q, strjoin(info_parts, ', '), result.R_squared);

        % Enable Accept button
        ui.AcceptFitButton.Enable = "on";
    end


    function local_on_accept_fit(~, ~)
        % Accept the pending fit and add peaks to dispersion
        if isempty(state.pendingFit)
            return
        end

        pf = state.pendingFit;
        result = pf.result;
        keep = pf.keep;
        abs_q = pf.abs_q;

        % Add accepted peaks to manual_points
        for p = 1:result.n_peaks
            if ~keep(p); continue; end
            new_pt = [abs_q, result.omega_p(p)];
            if isempty(state.manual_points)
                state.manual_points = new_pt;
            else
                state.manual_points(end+1, :) = new_pt;
            end
        end

        % Plot on dispersion axes
        ax = ui.DispersionAxes;
        hold(ax, 'on');
        for p = 1:result.n_peaks
            if ~keep(p); continue; end
            scatter(ax, abs_q, result.omega_p(p), 40, ...
                [0.1 0.5 0.9], 'filled', ...
                'HandleVisibility', 'off');
        end
        hold(ax, 'off');
        xlabel(ax, '|q| (Å^{-1})');
        ylabel(ax, 'Energy (meV)');
        title(ax, sprintf('Dispersion: %d pts', size(state.manual_points, 1)));
        grid(ax, 'on');

        % Update UI
        ui.PtsLabel.Text = sprintf('%d pts', size(state.manual_points, 1));
        ui.AcceptFitButton.Enable = "off";
        ui.FitInfoLabel.Text = sprintf('%s  ✓ Accepted', ui.FitInfoLabel.Text);

        local_log_operation(sprintf('Accept fit: q=%.4f, %d peaks', ...
            abs_q, sum(keep)));

        % Clear pending
        state.pendingFit = [];
    end


    function local_on_auto_fit_dispersion(~, ~)
        % Automatic multi-peak Drude-Lorentz loss function fitting
        if isempty(state.physicalQE) && isempty(state.comparisonQE)
            ui.InfoLabel.Text = "Load data first";
            ui.InfoLabel.Visible = "on";
            return
        end

        qe = local_get_active_qe();
        if isempty(qe)
            return
        end

        % Apply the same processing pipeline as single spectrum view
        if ui.AreaNormCheckbox.Value
            qe = local_apply_area_normalization(qe);
        end
        if ui.DenoiseCheckbox.Value
            qe = local_apply_denoise(qe);
        end
        if ui.BgSubCheckbox.Value
            qe = local_apply_bg_subtraction(qe);
        end
        if ui.DeconvCheckbox.Value
            qe = local_apply_deconvolution(qe);
        end

        % Also prepare raw (non-area-normalized) data for physical intensity
        qe_raw = local_get_active_qe();
        if ui.DenoiseCheckbox.Value
            qe_raw = local_apply_denoise(qe_raw);
        end
        if ui.BgSubCheckbox.Value
            qe_raw = local_apply_bg_subtraction(qe_raw);
        end
        if ui.DeconvCheckbox.Value
            qe_raw = local_apply_deconvolution(qe_raw);
        end
        % NOTE: qe_raw does NOT have area normalization applied

        [mask, energy_axis] = local_energy_mask(qe);
        q_axis = qe.q_Ainv;
        n_q = numel(q_axis);

        % Fitting window from current energy range
        E_min = max(ui.EnergyMinField.Value, 50);
        E_max = ui.EnergyMaxField.Value;

        % q range from display range controls (q start / q end)
        q_start = ui.QStartField.Value;
        q_end = ui.QEndField.Value;
        if q_start > q_end
            [q_start, q_end] = deal(q_end, q_start);
        end

        % Parse manual peak guesses (same as Fit Spectrum)
        guess_str = strtrim(ui.GuessField.Value);
        if ~isempty(guess_str)
            guesses = str2double(strsplit(guess_str, {',', ' ', ';'}));
            guesses = guesses(~isnan(guesses) & guesses > 0);
        else
            guesses = [];
        end

        % ═══════════ DUAL MODE: Seed vs Blind ═══════════
        if ~isempty(guesses)
            % ── SEED MODE: propagate from current q-channel ──
            ui.InfoLabel.Text = sprintf("Seed propagation (%d guesses, model=%s)...", ...
                numel(guesses), ui.PeakModelDropdown.Value);
            ui.InfoLabel.Visible = "on";
            drawnow;

            intensity = double(qe.intensity(mask, :));
            seed_idx = state.selectedQIndex;

            try
                prop = propagate_seed_peaks(intensity, energy_axis, q_axis(:), ...
                    'seed_guesses', guesses(:), ...
                    'seed_idx', seed_idx, ...
                    'direction', 'both', ...
                    'max_shift', ui.MaxShiftField.Value, ...
                    'peak_model', ui.PeakModelDropdown.Value, ...
                    'E_min', E_min, 'E_max', E_max, ...
                    'smooth_width', ui.SmoothField.Value, ...
                    'verbose', false);
            catch ME
                ui.InfoLabel.Text = sprintf("Seed propagation failed: %s", ME.message);
                return
            end

            if isempty(prop.peaks) || size(prop.peaks, 1) < 2
                ui.InfoLabel.Text = "Seed propagation: insufficient peaks found";
                return
            end

            all_peaks = prop.peaks;
            fit_details = prop.fit_details;
            n_success = prop.n_success;
            used_seed_mode = true;

        else
            % ── BLIND MODE: per-channel findpeaks (original) ──
            all_peaks = [];
            fit_details = cell(n_q, 1);
            used_seed_mode = false;

            n_in_range = sum(q_axis >= q_start & q_axis <= q_end);
            fprintf('Auto-fit: q range [%.4f, %.4f], %d/%d channels\n', q_start, q_end, n_in_range, n_q);
            ui.InfoLabel.Text = sprintf("Fitting %d channels (q: %.3f to %.3f)...", n_in_range, q_start, q_end);
            ui.InfoLabel.Visible = "on";
            drawnow;

            n_success = 0;
            for k = 1:n_q
                q_val = q_axis(k);
                if q_val < q_start || q_val > q_end
                    continue
                end

                spectrum = double(qe.intensity(mask, k));
                if all(spectrum == 0) || all(isnan(spectrum))
                    continue
                end

                try
                    result = fit_loss_function(energy_axis, spectrum, ...
                        'E_min', E_min, 'E_max', E_max, ...
                        'min_prominence', ui.PromField.Value, ...
                        'smooth_width', ui.SmoothField.Value, ...
                        'max_peaks', ui.MaxPeaksField.Value, ...
                        'initial_guesses', [], ...
                        'peak_model', ui.PeakModelDropdown.Value);
                    fit_details{k} = result;
                    n_success = n_success + 1;

                    % Also get raw (non-areanorm) peak heights at same positions
                    raw_spectrum = double(qe_raw.intensity(mask, k));
                    raw_peak_heights = NaN(result.n_peaks, 1);
                    for p = 1:result.n_peaks
                        % Evaluate the peak at E₀ using raw data scale
                        E0 = result.omega_p(p);
                        % Find the raw data value nearest the peak position
                        [~, e_idx] = min(abs(energy_axis - E0));
                        % Use a small window around E₀ to get robust peak height
                        hw = max(3, round(result.gamma(p) / mean(diff(energy_axis))));
                        win = max(1, e_idx-hw):min(numel(raw_spectrum), e_idx+hw);
                        raw_peak_heights(p) = max(raw_spectrum(win));
                    end

                    for p = 1:result.n_peaks
                        if result.omega_p(p) < E_min
                            continue
                        end
                        % Columns: [q, E, gamma, R², amp, E_ci_lo, E_ci_hi, G_ci_lo, G_ci_hi, A_ci_lo, A_ci_hi, raw_peak_height]
                        all_peaks(end+1, :) = [ ...
                            q_val, ...
                            result.omega_p(p), ...
                            result.gamma(p), ...
                            result.R_squared, ...
                            result.amplitude(p), ...
                            result.omega_p_ci(p, 1), result.omega_p_ci(p, 2), ...
                            result.gamma_ci(p, 1), result.gamma_ci(p, 2), ...
                            result.amplitude_ci(p, 1), result.amplitude_ci(p, 2), ...
                            raw_peak_heights(p)]; %#ok<AGROW>
                    end
                catch
                end

                if mod(k, 10) == 0
                    ui.InfoLabel.Text = sprintf("Auto-fitting... %d%%", round(k/n_q*100));
                    drawnow;
                end
            end
        end

        if isempty(all_peaks) || size(all_peaks, 1) < 2
            ui.InfoLabel.Text = "Auto-fit: insufficient peaks found";
            ui.InfoLabel.Visible = "on";
            return
        end

        % Filter by R² quality
        R2_threshold = 0.3;
        good = all_peaks(:, 4) > R2_threshold;
        peaks = all_peaks(good, :);

        if size(peaks, 1) < 2
            ui.InfoLabel.Text = sprintf('Auto-fit: only %d good peaks', size(peaks,1));
            ui.InfoLabel.Visible = "on";
            return
        end

        % ═══════════ Separate into branches ═══════════
        if used_seed_mode && size(peaks, 2) >= 6
            % --- SEED MODE: perfect assignment using branch_id ---
            unique_branches = unique(peaks(:, 6));
            n_branches = numel(unique_branches);
            branches = cell(n_branches, 1);
            for b = 1:n_branches
                b_id = unique_branches(b);
                branches{b} = peaks(peaks(:,6) == b_id, 1:5);  % drop branch_id for plot
                % Sort by q-value to fix display ordering
                [~, si] = sort(branches{b}(:,1));
                branches{b} = branches{b}(si, :);
            end
        else
            % --- BLIND MODE: use 1D gap-based clustering ---
            % Instead of fragile per-q sorting or kmeans (which needs a toolbox),
            % we find the (max_p - 1) largest gaps in the sorted energies.
            unique_q = unique(peaks(:,1));
            max_p = 1;
            for i = 1:numel(unique_q)
                max_p = max(max_p, sum(peaks(:,1) == unique_q(i)));
            end

            if max_p > 1
                energies = peaks(:,2);
                E_sorted = sort(energies);
                gaps = diff(E_sorted);
                [~, sort_gap_idx] = sort(gaps, 'descend');
                
                % The boundaries are the energies just before the largest gaps
                boundaries = sort(E_sorted(sort_gap_idx(1:max_p-1)));
                
                branches = cell(max_p, 1);
                for b = 1:max_p
                    if b == 1
                        mask = peaks(:,2) <= boundaries(1);
                    elseif b == max_p
                        mask = peaks(:,2) > boundaries(end);
                    else
                        mask = peaks(:,2) > boundaries(b-1) & peaks(:,2) <= boundaries(b);
                    end
                    
                    b_peaks = peaks(mask, :);
                    % Sort by q-value
                    [~, si] = sort(b_peaks(:,1));
                    branches{b} = b_peaks(si, :);
                end
            else
                % Only 1 branch detected
                branches = {peaks};
                [~, si] = sort(branches{1}(:,1));
                branches{1} = branches{1}(si, :);
            end
            
            % Remove empty branches just in case
            branches = branches(~cellfun('isempty', branches));
        end
        n_branches = numel(branches);

        autoResults = struct();
        autoResults.all_peaks = peaks;
        autoResults.branches = branches;
        autoResults.fit_details = fit_details;
        autoResults.n_success = n_success;
        state.autoFitResults = autoResults;

        % Use branch 1 (lowest energy) as manual_points
        state.manual_points = branches{1}(:, 1:2);

        % ═══════════ Plot ═══════════
        ax = ui.DispersionAxes;
        local_clear_axes(ax);
        hold(ax, 'on');

        branch_colors = [0.1 0.5 0.9; 0.9 0.3 0.1; 0.2 0.7 0.3; 0.6 0.2 0.8; 0.9 0.6 0.1];

        for b = 1:n_branches
            br = branches{b};
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
            branch_label = sprintf('Branch %d (%.0f-%.0f meV, %d pts)', ...
                b, min(br(:,2)), max(br(:,2)), size(br,1));

            % Plot with error bars if CI available (cols 6-7)
            if size(br, 2) >= 7 && ~all(isnan(br(:,6)))
                E_err = (br(:,7) - br(:,6)) / 2;  % half-width
                E_err(isnan(E_err)) = 0;
                errorbar(ax, br(:,1), br(:,2), E_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'LineWidth', 0.8, ...
                    'CapSize', 3, ...
                    'DisplayName', branch_label);
            else
                scatter(ax, br(:,1), br(:,2), 30, col, 'filled', ...
                    'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', branch_label);
            end

            % Try dispersion fit for branches with enough points
            disp_model_name = ui.DispModelDropdown.Value;
            if size(br, 1) >= 5
                try
                    disp_result = fit_dispersion_generic(br(:,1), br(:,2), ...
                        'model', disp_model_name);
                    plot(ax, disp_result.q_fit, disp_result.E_fit, '-', ...
                        'Color', col, 'LineWidth', 2, ...
                        'DisplayName', sprintf('Fit B%d: %s  R²=%.3f', ...
                            b, disp_result.model_label, disp_result.R_squared));
                    state.fitResults{b} = disp_result;
                catch
                end
            end
        end

        state.manual_branches = branches;

        hold(ax, 'off');
        legend(ax, 'Location', 'best', 'FontSize', 7);
        xlabel(ax, 'q (1/Å)');
        ylabel(ax, 'Energy (meV)');
        title(ax, sprintf('Auto-Fit: %d branches, %d peaks from %d ch', n_branches, size(peaks,1), n_success));
        grid(ax, 'on');

        ui.PtsLabel.Text = sprintf('%d pts', size(state.manual_points, 1));
        ui.InfoLabel.Text = sprintf('Auto-fit: %d peaks, %d branches', size(peaks,1), n_branches);
        ui.InfoLabel.Visible = "on";
        local_log_operation(sprintf('Auto fit: %d branches / %d peaks', n_branches, size(peaks,1)));
    end

    function local_on_reassign_points(~, ~)
        if isempty(state.manual_branches)
            ui.InfoLabel.Text = "No branches to reassign. Run Auto Fit first.";
            ui.InfoLabel.Visible = "on";
            return
        end
        
        ui.InfoLabel.Text = "Draw a region on the Dispersion map to select points to move...";
        ui.InfoLabel.Visible = "on";
        
        ax = ui.DispersionAxes;
        try
            roi = drawfreehand(ax, 'Color', 'r', 'LineWidth', 1.5);
        catch ME
            ui.InfoLabel.Text = "Selection cancelled or not supported.";
            return
        end
        
        if isempty(roi) || isempty(roi.Position)
            ui.InfoLabel.Text = "No region drawn.";
            return
        end
        
        % Check which points are inside
        in_b_idx = [];
        in_p_idx = [];
        poly_x = roi.Position(:,1);
        poly_y = roi.Position(:,2);
        
        for bi = 1:numel(state.manual_branches)
            bpts = state.manual_branches{bi};
            if isempty(bpts), continue; end
            in = inpolygon(bpts(:,1), bpts(:,2), poly_x, poly_y);
            
            p_idx = find(in);
            for k = 1:numel(p_idx)
                in_b_idx(end+1,1) = bi; %#ok<AGROW>
                in_p_idx(end+1,1) = p_idx(k); %#ok<AGROW>
            end
        end
        
        if isempty(in_b_idx)
            delete(roi);
            ui.InfoLabel.Text = "No points inside selected region.";
            return
        end
        
        n_sel = numel(in_b_idx);
        
        % Prompt for target branch
        cur_n_branches = numel(state.manual_branches);
        opts = cell(1, cur_n_branches + 2);
        for i = 1:cur_n_branches
            opts{i} = sprintf('Branch %d', i);
        end
        opts{cur_n_branches + 1} = sprintf('New Branch %d', cur_n_branches + 1);
        opts{end} = '❌ Delete Points';
        
        [sel_idx, ok] = listdlg('PromptString', sprintf('Move %d selected points to:', n_sel), ...
                                'SelectionMode', 'single', ...
                                'ListString', opts, ...
                                'Name', 'Reassign / Delete');
        
        if ~ok || isempty(sel_idx)
            delete(roi);
            ui.InfoLabel.Text = "Reassign cancelled.";
            return
        end
        
        % Delete points from original branches first
        pts_to_move = zeros(n_sel, size(state.manual_branches{1},2));
        for k = 1:n_sel
            pts_to_move(k,:) = state.manual_branches{in_b_idx(k)}(in_p_idx(k),:);
        end
        
        for bi = 1:cur_n_branches
            if ~any(in_b_idx == bi), continue; end
            keep_mask = true(size(state.manual_branches{bi}, 1), 1);
            keep_mask(in_p_idx(in_b_idx == bi)) = false;
            state.manual_branches{bi} = state.manual_branches{bi}(keep_mask, :);
        end
        
        % If not deleting, put them in target branch
        if sel_idx < numel(opts)
            if sel_idx > cur_n_branches
                state.manual_branches{sel_idx} = pts_to_move;
            else
                state.manual_branches{sel_idx} = [state.manual_branches{sel_idx}; pts_to_move];
            end
        end
        
        % Clean up empty branches and sort
        new_b = {};
        for bi = 1:numel(state.manual_branches)
            if ~isempty(state.manual_branches{bi})
                new_b{end+1} = sortrows(state.manual_branches{bi}, 1); %#ok<AGROW>
            end
        end
        state.manual_branches = new_b;
        
        delete(roi);
        
        if sel_idx == numel(opts)
            ui.InfoLabel.Text = sprintf("Successfully deleted %d points.", n_sel);
            local_log_operation(sprintf("Deleted %d points", n_sel));
        else
            ui.InfoLabel.Text = sprintf("Successfully moved %d points.", n_sel);
            local_log_operation(sprintf("Reassigned %d points to Branch %d", n_sel, sel_idx));
        end
        
        % Force update of overlay
        local_update_all_views();
    end

    function local_on_fit_dispersion(~, ~)
        % Re-fit dispersion to existing branches with the selected model
        if isempty(state.manual_branches)
            ui.DispInfoLabel.Text = "No branches available. Run Auto Fit first.";
            return
        end

        disp_model_name = ui.DispModelDropdown.Value;
        ax = ui.DispersionAxes;
        local_clear_axes(ax);
        hold(ax, 'on');

        branch_colors = [0.1 0.5 0.9; 0.9 0.3 0.1; 0.2 0.7 0.3; 0.6 0.2 0.8; 0.9 0.6 0.1];
        state.fitResults = {};
        n_branches = numel(state.manual_branches);

        for b = 1:n_branches
            br = state.manual_branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);

            % Plot points (with error bars if available)
            if size(br, 2) >= 7 && ~all(isnan(br(:,6)))
                E_err = (br(:,7) - br(:,6)) / 2;
                E_err(isnan(E_err)) = 0;
                errorbar(ax, br(:,1), br(:,2), E_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
                    'DisplayName', sprintf('Branch %d', b));
            else
                scatter(ax, br(:,1), br(:,2), 30, col, 'filled', ...
                    'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', sprintf('Branch %d', b));
            end

            % Fit dispersion
            if size(br, 1) >= 5
                try
                    disp_result = fit_dispersion_generic(br(:,1), br(:,2), ...
                        'model', disp_model_name);
                    plot(ax, disp_result.q_fit, disp_result.E_fit, '-', ...
                        'Color', col, 'LineWidth', 2, ...
                        'DisplayName', sprintf('Fit B%d: %s  R²=%.3f', ...
                            b, disp_result.model_label, disp_result.R_squared));
                    state.fitResults{b} = disp_result;
                catch ME
                    fprintf('Fit failed for branch %d: %s\n', b, ME.message);
                end
            end
        end

        hold(ax, 'off');
        grid(ax, 'on');
        xlabel(ax, 'q (1/Å)');
        ylabel(ax, 'Energy (meV)');
        title(ax, sprintf('Dispersion [%s] | %d branches', disp_model_name, n_branches));
        legend(ax, 'Location', 'best', 'FontSize', 7);

        % Summarize in info label
        n_fitted = sum(~cellfun('isempty', state.fitResults));
        ui.DispInfoLabel.Text = sprintf('Fitted %d/%d branches with %s', ...
            n_fitted, n_branches, disp_model_name);
        local_log_operation(sprintf('Fit dispersion: %s, %d branches', disp_model_name, n_fitted));
    end

    function local_on_export_dispersion(~, ~)
        % Export a standalone publication-quality dispersion figure
        if isempty(state.manual_branches)
            ui.DispInfoLabel.Text = "No branches to export.";
            return
        end

        % Create standalone figure
        fig_exp = figure('Name', 'Dispersion Export', ...
            'Position', [100 100 700 500], ...
            'Color', 'w');
        ax_exp = axes(fig_exp);
        hold(ax_exp, 'on');

        branch_colors = [0.1 0.5 0.9; 0.9 0.3 0.1; 0.2 0.7 0.3; 0.6 0.2 0.8; 0.9 0.6 0.1];

        for b = 1:numel(state.manual_branches)
            br = state.manual_branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);

            % Plot with error bars if available
            if size(br, 2) >= 7 && ~all(isnan(br(:,6)))
                E_err = (br(:,7) - br(:,6)) / 2;
                E_err(isnan(E_err)) = 0;
                errorbar(ax_exp, br(:,1), br(:,2), E_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 5, 'LineWidth', 1, 'CapSize', 4, ...
                    'DisplayName', sprintf('Branch %d', b));
            else
                scatter(ax_exp, br(:,1), br(:,2), 40, col, 'filled', ...
                    'MarkerFaceAlpha', 0.8, ...
                    'DisplayName', sprintf('Branch %d', b));
            end

            % Overlay fit curve with annotation
            if numel(state.fitResults) >= b && ~isempty(state.fitResults{b})
                fr = state.fitResults{b};
                plot(ax_exp, fr.q_fit, fr.E_fit, '-', ...
                    'Color', col, 'LineWidth', 2.5, ...
                    'DisplayName', sprintf('Fit: %s  R²=%.3f', ...
                        fr.model_label, fr.R_squared));
            end
        end

        hold(ax_exp, 'off');
        grid(ax_exp, 'on');
        box(ax_exp, 'on');
        xlabel(ax_exp, 'q (Å^{-1})', 'FontSize', 14);
        ylabel(ax_exp, 'Energy (meV)', 'FontSize', 14);
        title(ax_exp, 'Dispersion Relation', 'FontSize', 16);
        legend(ax_exp, 'Location', 'best', 'FontSize', 10);
        set(ax_exp, 'FontSize', 12, 'LineWidth', 1.2);

        % Save dialog
        [fname, fpath] = uiputfile({'*.png'; '*.pdf'; '*.svg'}, ...
            'Save Dispersion Figure', 'dispersion_export.png');
        if fname ~= 0
            full_path = fullfile(fpath, fname);
            exportgraphics(fig_exp, full_path, 'Resolution', 300);
            ui.DispInfoLabel.Text = sprintf('Exported to %s', fname);
            local_log_operation(sprintf('Exported dispersion figure: %s', full_path));
        end
    end

    function local_on_pick_guesses(~, ~)
        % Toggle guess-picking mode: click on spectrum to pick peak positions
        if state.guessPickArmed
            % Disarm
            state.guessPickArmed = false;
            ui.PickGuessesButton.Text = "Pick Guesses";
            ui.PickGuessesButton.BackgroundColor = [0.96 0.96 0.96];
            ui.SeedInfoLabel.Text = sprintf("Guesses: %s", ui.GuessField.Value);
        else
            % Arm — clear old guess markers and GuessField
            state.guessPickArmed = true;
            state.manualPickArmed = false;  % disarm the other picker
            ui.PickPtsButton.Text = "Pick Peaks";
            ui.PickPtsButton.BackgroundColor = [0.96 0.96 0.96];
            ui.GuessField.Value = "";
            delete(findobj(ui.SingleAxes, 'Tag', 'guess_marker'));
            ui.PickGuessesButton.Text = "Picking...";
            ui.PickGuessesButton.BackgroundColor = [1 0.9 0.85];
            ui.SeedInfoLabel.Text = "Click on spectrum peaks to add guesses";
        end
    end




    function local_on_show_gamma(~, ~)
        % Open a separate figure showing Γ(q), Γ(E), Amplitude(q), Amplitude(E)
        if ~isfield(state, 'autoFitResults') || isempty(state.autoFitResults)
            ui.InfoLabel.Text = "Run Auto Fit first";
            ui.InfoLabel.Visible = "on";
            return
        end

        peaks = state.autoFitResults.all_peaks;
        branches = state.autoFitResults.branches;
        n_branches = numel(branches);
        branch_colors = [0.1 0.5 0.9; 0.9 0.3 0.1; 0.2 0.7 0.3; 0.6 0.2 0.8; 0.9 0.6 0.1];

        % Check if CI columns are available
        % Cols: [q E gamma R² amp E_ci_lo E_ci_hi G_ci_lo G_ci_hi A_ci_lo A_ci_hi]
        %        1 2   3    4  5     6       7       8       9      10      11
        has_gamma_ci = size(peaks, 2) >= 9;
        has_amp_ci = size(peaks, 2) >= 11;

        fig = figure('Name', 'Linewidth Γ, Amplitude & EELS Prefactor', ...
            'NumberTitle', 'off', 'Color', 'w', ...
            'Position', [120 40 1100 1000]);

        % --- Subplot 1: Γ(q) ---
        ax1 = subplot(3, 2, 1, 'Parent', fig);
        hold(ax1, 'on');
        for b = 1:n_branches
            br = branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
            lbl = sprintf('B%d: Γ̄=%.0f meV', b, mean(br(:,3)));

            if has_gamma_ci && size(br, 2) >= 9 && ~all(isnan(br(:,8)))
                G_err = (br(:,9) - br(:,8)) / 2;
                G_err(isnan(G_err)) = 0;
                errorbar(ax1, br(:,1), br(:,3), G_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
                    'DisplayName', lbl);
            else
                scatter(ax1, br(:,1), br(:,3), 40, col, 'filled', ...
                    'MarkerFaceAlpha', 0.7, 'DisplayName', lbl);
            end
        end
        hold(ax1, 'off');
        xlabel(ax1, 'q (1/Å)'); ylabel(ax1, 'Γ (meV)');
        title(ax1, 'Linewidth Γ vs Momentum');
        legend(ax1, 'Location', 'best'); grid(ax1, 'on');

        % --- Subplot 2: Γ(E) ---
        ax2 = subplot(3, 2, 2, 'Parent', fig);
        hold(ax2, 'on');
        for b = 1:n_branches
            br = branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
            lbl = sprintf('B%d: Γ̄=%.0f meV', b, mean(br(:,3)));

            if has_gamma_ci && size(br, 2) >= 9 && ~all(isnan(br(:,8)))
                G_err = (br(:,9) - br(:,8)) / 2;
                G_err(isnan(G_err)) = 0;
                errorbar(ax2, br(:,2), br(:,3), G_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
                    'DisplayName', lbl);
            else
                scatter(ax2, br(:,2), br(:,3), 40, col, 'filled', ...
                    'MarkerFaceAlpha', 0.7, 'DisplayName', lbl);
            end
        end
        hold(ax2, 'off');
        xlabel(ax2, 'ω_p (meV)'); ylabel(ax2, 'Γ (meV)');
        title(ax2, 'Linewidth Γ vs Peak Energy');
        legend(ax2, 'Location', 'best'); grid(ax2, 'on');

        % --- Subplot 3: Amplitude(q) ---
        ax3 = subplot(3, 2, 3, 'Parent', fig);
        hold(ax3, 'on');
        for b = 1:n_branches
            br = branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
            lbl = sprintf('B%d', b);

            if has_amp_ci && size(br, 2) >= 11 && ~all(isnan(br(:,10)))
                A_err = (br(:,11) - br(:,10)) / 2;
                A_err(isnan(A_err)) = 0;
                errorbar(ax3, br(:,1), br(:,5), A_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
                    'DisplayName', lbl);
            else
                scatter(ax3, br(:,1), br(:,5), 40, col, 'filled', ...
                    'MarkerFaceAlpha', 0.7, 'DisplayName', lbl);
            end
        end
        hold(ax3, 'off');
        xlabel(ax3, 'q (1/Å)'); ylabel(ax3, 'Amplitude (a.u.)');
        title(ax3, 'Peak Amplitude vs Momentum');
        legend(ax3, 'Location', 'best'); grid(ax3, 'on');

        % --- Subplot 4: Amplitude(E) ---
        ax4 = subplot(3, 2, 4, 'Parent', fig);
        hold(ax4, 'on');
        for b = 1:n_branches
            br = branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
            lbl = sprintf('B%d', b);

            if has_amp_ci && size(br, 2) >= 11 && ~all(isnan(br(:,10)))
                A_err = (br(:,11) - br(:,10)) / 2;
                A_err(isnan(A_err)) = 0;
                errorbar(ax4, br(:,2), br(:,5), A_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
                    'DisplayName', lbl);
            else
                scatter(ax4, br(:,2), br(:,5), 40, col, 'filled', ...
                    'MarkerFaceAlpha', 0.7, 'DisplayName', lbl);
            end
        end
        hold(ax4, 'off');
        xlabel(ax4, 'ω_p (meV)'); ylabel(ax4, 'Amplitude (a.u.)');
        title(ax4, 'Peak Amplitude vs Energy');
        legend(ax4, 'Location', 'best'); grid(ax4, 'on');

        % --- Row 3: Raw Peak Height A(q) log-log with power law ---
        % Uses col 12 = raw (non-areanorm) peak height from original data
        ax5 = subplot(3, 2, 5, 'Parent', fig);
        hold(ax5, 'on');
        for b = 1:n_branches
            br = branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);

            q_abs = abs(br(:,1));
            if size(br, 2) >= 12 && ~all(isnan(br(:,12)))
                A_raw = br(:,12);
                src_label = 'raw';
            else
                A_raw = br(:,5);
                src_label = 'model';
            end
            valid = q_abs > 0 & A_raw > 0 & isfinite(A_raw);
            if sum(valid) < 2, continue; end

            loglog(ax5, q_abs(valid), A_raw(valid), 'o', ...
                'Color', col, 'MarkerFaceColor', col, ...
                'MarkerSize', 5, 'DisplayName', sprintf('B%d (%s)', b, src_label));

            log_q = log(q_abs(valid));
            log_A = log(A_raw(valid));
            P = polyfit(log_q, log_A, 1);
            q_fit_line = linspace(min(q_abs(valid)), max(q_abs(valid)), 100);
            loglog(ax5, q_fit_line, exp(P(2))*q_fit_line.^P(1), '--', ...
                'Color', col, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('B%d: q^{%.1f}', b, P(1)));

            q_median = median(q_abs(valid));
            low_mask = q_abs(valid) <= q_median;
            high_mask = q_abs(valid) > q_median;
            if sum(low_mask) >= 3 && sum(high_mask) >= 3
                P_lo = polyfit(log_q(low_mask), log_A(low_mask), 1);
                P_hi = polyfit(log_q(high_mask), log_A(high_mask), 1);
                q_lo = linspace(min(q_abs(valid)), q_median, 50);
                q_hi = linspace(q_median, max(q_abs(valid)), 50);
                loglog(ax5, q_lo, exp(P_lo(2))*q_lo.^P_lo(1), ':', ...
                    'Color', col*0.6, 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('lo: q^{%.1f}', P_lo(1)));
                loglog(ax5, q_hi, exp(P_hi(2))*q_hi.^P_hi(1), '-.', ...
                    'Color', col*0.6+0.4, 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('hi: q^{%.1f}', P_hi(1)));
            end
        end
        hold(ax5, 'off');
        xlabel(ax5, 'q (Å^{-1})'); ylabel(ax5, 'Peak Height (a.u.)');
        title(ax5, 'Raw Peak Height A(q) — log-log (no areanorm)');
        legend(ax5, 'Location', 'best', 'FontSize', 7); grid(ax5, 'on');

        % Subplot 6: Quality factor Q = E₀/Γ vs q
        ax6 = subplot(3, 2, 6, 'Parent', fig);
        hold(ax6, 'on');
        for b = 1:n_branches
            br = branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);
            Q_factor = br(:,2) ./ br(:,3);
            scatter(ax6, br(:,1), Q_factor, 30, col, 'filled', ...
                'MarkerFaceAlpha', 0.7, ...
                'DisplayName', sprintf('B%d: Q̄=%.1f', b, mean(Q_factor)));
        end
        hold(ax6, 'off');
        xlabel(ax6, 'q (Å^{-1})'); ylabel(ax6, 'Q = ω_p / Γ');
        title(ax6, 'Plasmon Quality Factor vs q');
        legend(ax6, 'Location', 'best'); grid(ax6, 'on');

        sgtitle(fig, sprintf('Γ, Amplitude & Quality Factor (%d branches, %d peaks)', ...
            n_branches, size(peaks, 1)), 'FontSize', 14);

        local_log_operation('Show Γ & Amplitude analysis');
    end




    function local_on_export(~, ~)
        % Export axes panels as publication-quality images.
        % Copies uiaxes content to a temporary invisible regular figure
        % for reliable export (avoids uiaxes rendering quirks).
        % Aspect ratio is applied via post-processing (white padding).
        panel_axes  = {ui.QEAxes, ui.ComparisonAxes, ui.SingleAxes, ui.DispersionAxes};
        panel_names = {'qE_map',  'comparison',      'spectrum',     'dispersion'};
        panel_labels = {'q-E Map', 'Comparison Map', 'Single Spectrum', 'Dispersion'};

        % Panel selection via listdlg (uiconfirm only supports up to 4 options)
        options = ['All panels', panel_labels];
        [idx, ok] = listdlg('ListString', options, ...
            'SelectionMode', 'single', ...
            'ListSize', [220 140], ...
            'Name', 'Export', ...
            'PromptString', 'Select panels to export:');
        if ~ok, return; end
        if idx == 1
            sel = 1:4;
        else
            sel = idx - 1;
        end

        % Use source data folder as default export location
        default_dir = pwd;
        if ~isempty(state.dataset) && isfield(state.dataset, 'source_path')
            [src_folder, ~, ~] = fileparts(char(state.dataset.source_path));
            if isfolder(src_folder)
                default_dir = src_folder;
            end
        end
        out_dir = uigetdir(default_dir, 'Select export folder');
        if isequal(out_dir, 0)
            return
        end

        % Parse aspect ratio from dropdown
        ratio_str = char(ui.ExportRatioDropdown.Value);
        % Sanitize for filenames: replace ':' with 'x' (e.g. '4:3' → '4x3')
        ratio_file = strrep(ratio_str, ':', 'x');
        ratio_parts = sscanf(ratio_str, '%d:%d');
        if numel(ratio_parts) == 2 && all(ratio_parts > 0)
            target_ratio = ratio_parts(1) / ratio_parts(2);  % W/H
        else
            target_ratio = 4 / 3;
        end

        ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

        exported = {};
        for k = sel
            src_ax = panel_axes{k};
            fname = panel_names{k};

            try
                % --- Create temp figure sized to target aspect ratio ---
                fig_h_px = 600;
                fig_w_px = round(fig_h_px * target_ratio);
                tmp_fig = figure('Visible', 'off', 'Color', 'w', ...
                    'Units', 'pixels', ...
                    'Position', [100 100 fig_w_px fig_h_px]);

                % Copy axes content to temp figure
                tmp_ax = copyobj(src_ax, tmp_fig);
                set(tmp_ax, 'Units', 'normalized', ...
                    'Position', [0.12 0.14 0.82 0.78], ...
                    'ActivePositionProperty', 'position');
                % Ensure colormap matches source
                try
                    colormap(tmp_ax, colormap(src_ax));
                catch
                end
                % Copy CLim if applicable
                try
                    tmp_ax.CLim = src_ax.CLim;
                catch
                end

                % --- Export PNG at 300 DPI ---
                png_path = fullfile(out_dir, sprintf('%s_%s_%s.png', ts, fname, ratio_file));
                exportgraphics(tmp_fig, png_path, 'Resolution', 300, 'BackgroundColor', 'white');

                exported{end+1} = fname; %#ok<AGROW>

                % --- Also export vector PDF ---
                try
                    pdf_path = fullfile(out_dir, sprintf('%s_%s_%s.pdf', ts, fname, ratio_file));
                    exportgraphics(tmp_fig, pdf_path, 'ContentType', 'vector', 'BackgroundColor', 'white');
                catch
                end

                % --- Clean up temp figure ---
                close(tmp_fig);
            catch ME
                fprintf('  Export failed for %s: %s\n', fname, ME.message);
                try close(tmp_fig); catch, end  %#ok<CTCH>
            end
        end

        if ~isempty(exported)
            msg = sprintf('Exported %d panels (%s) to %s', numel(exported), ratio_str, out_dir);
        else
            msg = 'Export failed — check MATLAB console for details';
        end
        ui.InfoLabel.Text = msg;
        ui.InfoLabel.Visible = "on";
        local_log_operation(sprintf('Export %d panels (%s): %s', numel(exported), ratio_str, strjoin(exported, ', ')));
        fprintf('  %s\n', msg);
    end


%<RENDER_CALLBACKS>
    function local_update_all_views()
        qe = local_get_active_qe();
        local_update_stage_label();
        local_update_info_label();

        if isempty(state.physicalQE)
            local_clear_axes(ui.QEAxes);
            local_clear_axes(ui.ComparisonAxes);
            local_clear_axes(ui.SingleAxes);
            local_clear_axes(ui.DispersionAxes);
            title(ui.QEAxes, "Physical q-E Map");
            title(ui.ComparisonAxes, "Normalized Off-axis Component");
            title(ui.SingleAxes, "Spectrum");
            title(ui.DispersionAxes, "Dispersion");
            return
        end

        state.selectedQIndex = min(max(1, state.selectedQIndex), size(state.physicalQE.intensity, 2));
        local_plot_qe_map(ui.QEAxes, state.physicalQE, "Physical");
        local_plot_normalized_map();
        local_plot_single_spectrum(qe);
        local_plot_dispersion_result();
    end


    function local_clear_axes(ax)
        if isempty(ax) || ~isgraphics(ax)
            return
        end

        hold(ax, "off");
        try
            legend(ax, "off");
        catch
        end
        delete(allchild(ax));
    end


    function local_plot_qe_map(ax, qe, map_label, display_style)
        if nargin < 4
            display_style = "physical";
        end
        local_clear_axes(ax);
        if ui.AreaNormCheckbox.Value
            qe = local_apply_area_normalization(qe);
        end
        if ui.DenoiseCheckbox.Value
            qe = local_apply_denoise(qe);
        end
        if ui.BgSubCheckbox.Value
            qe = local_apply_bg_subtraction(qe);
        end
        if ui.DeconvCheckbox.Value
            qe = local_apply_deconvolution(qe);
        end
        [energy_mask, energy_axis] = local_energy_mask(qe);
        [q_mask, q_axis] = local_q_mask(qe);
        map_values = local_build_map_values(qe, energy_mask, q_mask, "display");
        [display_map, display_clim] = local_prepare_map_display(map_values, "display", display_style);

        imagesc(ax, q_axis, energy_axis, display_map);
        axis(ax, "xy");
        colormap(ax, turbo);
        clim(ax, display_clim);
        xlabel(ax, "q (1/A)");
        ylabel(ax, "Energy relative to ZLP (meV)");
        title(ax, sprintf("%s q-E | %s", map_label, local_make_stage_title(qe)), ...
            "Interpreter", "none");
        ax.ButtonDownFcn = @(src, event)local_on_map_clicked(ax, event);
        xlim(ax, [q_axis(1), q_axis(end)]);
        ylim(ax, [energy_axis(1), energy_axis(end)]);

        image_handle = findobj(ax, "Type", "image");
        if ~isempty(image_handle)
            image_handle(1).PickableParts = "all";
            image_handle(1).ButtonDownFcn = @(src, event)local_on_map_clicked(ax, event);
        end

        hold(ax, "on");
        selected_q = qe.q_Ainv(state.selectedQIndex);
        xline(ax, selected_q, "-", ...
            "Color", [0.95 0.95 0.98], ...
            "LineWidth", 1.2, ...
            "Label", sprintf("q=%.4f", selected_q), ...
            "LabelHorizontalAlignment", "left");
        local_plot_dispersion_overlay(ax, qe);
        hold(ax, "off");
    end


    function local_plot_normalized_map()
        ax = ui.ComparisonAxes;
        local_clear_axes(ax);

        if isempty(state.comparisonQE)
            title(ax, "Normalized Off-axis Component | press Build Views");
            xlabel(ax, "q (1/A)");
            ylabel(ax, "Energy relative to ZLP (meV)");
            grid(ax, "on");
            return
        end

        local_plot_qe_map(ax, state.comparisonQE, "Normalized Off-axis Component", "normalized");
    end


    function local_plot_single_spectrum(qe)
        ax = ui.SingleAxes;
        local_clear_axes(ax);

        if ui.AreaNormCheckbox.Value
            qe = local_apply_area_normalization(qe);
        end
        if ui.DenoiseCheckbox.Value
            qe = local_apply_denoise(qe);
        end
        if ui.BgSubCheckbox.Value
            qe = local_apply_bg_subtraction(qe);
        end
        if ui.DeconvCheckbox.Value
            qe = local_apply_deconvolution(qe);
        end
        [mask, energy_axis] = local_energy_mask(qe);
        q_index = state.selectedQIndex;
        raw_spectrum = double(qe.intensity(mask, q_index));
        features = local_build_spectrum_features(energy_axis, raw_spectrum);
        mode_name = local_trace_mode();
        display_values = local_trace_values(features, mode_name);
        display_raw = raw_spectrum;
        display_smooth = features.smooth;
        y_scale = local_trace_y_scale();
        plot_handles = gobjects(0);

        if strcmpi(y_scale, "log")
            display_raw = max(display_raw, eps);
            display_smooth = max(display_smooth, eps);
            display_values = max(display_values, eps);
        end

        hold(ax, "on");
        switch mode_name
            case "display"
                if ~strcmpi(ui.SmoothingModeDropdown.Value, "off")
                    plot_handles(end + 1) = plot(ax, energy_axis, display_raw, "-", ... %#ok<AGROW>
                        "Color", [0.7 0.7 0.72], ...
                        "LineWidth", 0.9, ...
                        "DisplayName", "raw");
                end
                plot_handles(end + 1) = plot(ax, energy_axis, display_smooth, "-", ... %#ok<AGROW>
                    "Color", [0.05 0.32 0.75], ...
                    "LineWidth", 1.4, ...
                    "DisplayName", "display");

            case "background-subtracted"
                plot_handles(end + 1) = plot(ax, energy_axis, display_raw, "-", ... %#ok<AGROW>
                    "Color", [0.8 0.8 0.82], ...
                    "LineWidth", 0.8, ...
                    "DisplayName", "raw");
                plot_handles(end + 1) = plot(ax, energy_axis, features.baseline, "--", ... %#ok<AGROW>
                    "Color", [0.85 0.5 0.1], ...
                    "LineWidth", 1.1, ...
                    "DisplayName", "baseline");
                plot_handles(end + 1) = plot(ax, energy_axis, display_values, "-", ... %#ok<AGROW>
                    "Color", [0.05 0.32 0.75], ...
                    "LineWidth", 1.5, ...
                    "DisplayName", "residual");

            case "curvature"
                plot_handles(end + 1) = plot(ax, energy_axis, display_values, "-", ... %#ok<AGROW>
                    "Color", [0.4 0.1 0.75], ...
                    "LineWidth", 1.5, ...
                    "DisplayName", "curvature");
        end




        hold(ax, "off");

        grid(ax, "on");
        xlabel(ax, "Energy relative to ZLP (meV)");
        ylabel(ax, local_single_ylabel());
        title(ax, sprintf("q channel %.0f | q = %.4f 1/A | %s", ...
            qe.q_channel(q_index), qe.q_Ainv(q_index), mode_name), ...
            "Interpreter", "none");
        ax.YScale = y_scale;
        xlim(ax, [energy_axis(1), energy_axis(end)]);
        local_apply_y_limits(ax, [display_raw(:); display_smooth(:); display_values(:)]);
        ax.ButtonDownFcn = @local_on_single_axes_clicked;
        for handle_idx = 1:numel(plot_handles)
            if isgraphics(plot_handles(handle_idx))
                plot_handles(handle_idx).PickableParts = "all";
                plot_handles(handle_idx).ButtonDownFcn = @local_on_single_axes_clicked;
            end
        end
        if ~strcmpi(mode_name, "display") || ~strcmpi(ui.SmoothingModeDropdown.Value, "off")
            legend(ax, "Location", "best");
        end
    end


    function local_plot_waterfall(qe)
        ax = ui.WaterfallAxes;
        local_clear_axes(ax);

        [mask, energy_axis] = local_energy_mask(qe);
        q_indices = local_waterfall_indices(qe);
        if isempty(q_indices)
            title(ax, "Waterfall");
            return
        end

        offset_value = ui.OffsetField.Value;
        plotted_values = [];

        hold(ax, "on");
        for idx = 1:numel(q_indices)
            q_index = q_indices(idx);
            spectrum = double(qe.intensity(mask, q_index));
            features = local_build_spectrum_features(energy_axis, spectrum);
            spectrum = local_trace_values(features, local_trace_mode());
            if strcmpi(local_trace_y_scale(), "log")
                spectrum = max(spectrum, eps);
            end

            shifted = spectrum + (idx - 1) * offset_value;
            plotted_values = [plotted_values; shifted(:)]; %#ok<AGROW>

            line_width = 1.0;
            if q_index == state.selectedQIndex
                line_width = 1.8;
            end

            trace_color = local_signed_q_color(qe.q_Ainv(q_index), qe.q_Ainv);
            plot(ax, energy_axis, shifted, "-", ...
                "Color", trace_color, ...
                "LineWidth", line_width, ...
                "DisplayName", sprintf("q=%.4f", qe.q_Ainv(q_index)));
        end
        hold(ax, "off");

        grid(ax, "on");
        xlabel(ax, "Energy relative to ZLP (meV)");
        ylabel(ax, "Offset intensity");
        title(ax, sprintf("%s view waterfall | %s | q=0 auto-centered", ...
            local_current_view_name(), local_trace_mode()), ...
            "Interpreter", "none");
        ax.YScale = local_trace_y_scale();
        xlim(ax, [energy_axis(1), energy_axis(end)]);
        local_apply_y_limits(ax, plotted_values);
    end


    function local_plot_dispersion_overlay(ax, qe)
        % --- Manual points overlay ---
        if ~isempty(state.manual_branches)
            branch_colors = {[0.9 0.15 0.15], [0.15 0.45 0.9], [0.85 0.65 0.0]};
            for bi = 1:numel(state.manual_branches)
                bpts = state.manual_branches{bi};
                if isempty(bpts), continue; end
                bc = branch_colors{min(bi, numel(branch_colors))};
                scatter(ax, bpts(:,1), bpts(:,2), 36, 'filled', ...
                    'MarkerFaceColor', bc, 'MarkerEdgeColor', [1 1 1], ...
                    'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
        elseif ~isempty(state.manual_points)
            pts = sortrows(state.manual_points, 1);
            mc = [0.1 0.8 0.2];
            scatter(ax, pts(:,1), pts(:,2), 36, 'filled', ...
                'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', mc, ...
                'LineWidth', 1.2, 'HandleVisibility', 'off');
        end

    end


    function local_plot_dispersion_result()
        ax = ui.DispersionAxes;
        local_clear_axes(ax);

        has_manual = ~isempty(state.manual_points) || ~isempty(state.manual_branches);

        if ~has_manual
            title(ax, "Dispersion");
            xlabel(ax, "q (1/A)");
            ylabel(ax, "Peak energy (meV)");
            grid(ax, "on");
            return
        end

        hold(ax, "on");
        if ~isempty(state.manual_branches)
            branch_colors = {[0.9 0.15 0.15], [0.15 0.45 0.9], [0.85 0.65 0.0], [0.6 0.2 0.8], [0.9 0.6 0.1]};
            for bi = 1:numel(state.manual_branches)
                bpts = state.manual_branches{bi};
                if isempty(bpts), continue; end
                bc = branch_colors{min(bi, numel(branch_colors))};
                bn = sprintf('Branch %d', bi);

                if size(bpts, 2) >= 7 && ~all(isnan(bpts(:,6)))
                    E_err = (bpts(:,7) - bpts(:,6)) / 2;
                    E_err(isnan(E_err)) = 0;
                    errorbar(ax, bpts(:,1), bpts(:,2), E_err, ...
                        'o', 'Color', bc, 'MarkerFaceColor', bc, ...
                        'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
                        'DisplayName', bn);
                else
                    scatter(ax, bpts(:,1), bpts(:,2), 30, bc, 'filled', ...
                        'MarkerFaceAlpha', 0.7, ...
                        'DisplayName', bn);
                end
            end
        elseif ~isempty(state.manual_points)
            pts = sortrows(state.manual_points, 1);
            mc = [0.1 0.8 0.2];
            scatter(ax, pts(:,1), pts(:,2), 30, mc, 'filled', ...
                'MarkerFaceAlpha', 0.7, ...
                'DisplayName', sprintf('manual (%d pts)', size(pts,1)));
        end

        % Re-draw stored model fit curves
        if ~isempty(state.fitResults)
            fit_colors = {[0.8 0.0 0.4], [0.0 0.4 0.8], [0.7 0.5 0.0], [0.4 0.1 0.6], [0.8 0.5 0.0]};
            for fi = 1:numel(state.fitResults)
                fr = state.fitResults{fi};
                if isempty(fr), continue; end
                fc = fit_colors{min(fi, numel(fit_colors))};

                if isfield(fr, 'model_label')
                    lbl = sprintf('Fit B%d: %s  R²=%.3f', fi, fr.model_label, fr.R_squared);
                elseif isfield(fr, 'rho0')
                    lbl = sprintf('Fit B%d: ρ₀=%.1fÅ  E_f=%.0fmeV', fi, fr.rho0, fr.E_flat_meV);
                else
                    lbl = sprintf('Fit B%d: R²=%.3f', fi, fr.R_squared);
                end

                plot(ax, fr.q_fit, fr.E_fit, '-', ...
                    'Color', fc, 'LineWidth', 2.5, ...
                    'DisplayName', lbl);
            end
        end

        hold(ax, "off");

        grid(ax, "on");
        xlabel(ax, "q (1/A)");
        ylabel(ax, "Peak energy (meV)");
        if ~isempty(state.manual_branches)
            title(ax, sprintf("Manual dispersion | %d branches, %d pts", ...
                numel(state.manual_branches), size(state.manual_points, 1)));
        else
            title(ax, sprintf("Manual dispersion | %d points", size(state.manual_points, 1)));
        end
        xlim(ax, sort([ui.QStartField.Value, ui.QEndField.Value]));
        legend(ax, "Location", "best");
    end




    function local_plot_zlp_diagnostic()
        ax = ui.ZLPAxes;
        local_clear_axes(ax);

        [diag_q, raw_metric, scale_factor] = local_get_zlp_diagnostic_data();
        if isempty(diag_q)
            title(ax, "ZLP vs q");
            return
        end

        yyaxis(ax, "left");
        plot(ax, diag_q, raw_metric, "-", ...
            "Color", [0.8 0.25 0.1], ...
            "LineWidth", 1.3);
        ylabel(ax, "ZLP window metric");

        yyaxis(ax, "right");
        plot(ax, diag_q, scale_factor, "-", ...
            "Color", [0.12 0.45 0.25], ...
            "LineWidth", 1.2);
        ylabel(ax, "Normalization factor");

        selected_q = diag_q(min(max(state.selectedQIndex, 1), numel(diag_q)));
        state.zlpSelectionLine = xline(ax, selected_q, "--", ...
            "Color", [0.25 0.25 0.3], ...
            "LineWidth", 1.1);

        xlabel(ax, "q (1/A)");
        grid(ax, "on");
        title(ax, sprintf("ZLP diagnostic | window [%.1f, %.1f] meV", ...
            ui.ZlpMinField.Value, ui.ZlpMaxField.Value));
    end


    function local_on_map_clicked(ax, ~)
        qe = local_get_reference_qe();
        if isempty(qe) || isempty(state.physicalQE)
            return
        end

        current_point = ax.CurrentPoint;
        q_value = current_point(1, 1);
        energy_value = current_point(1, 2);
        [~, state.selectedQIndex] = min(abs(state.physicalQE.q_Ainv - q_value));

        local_update_all_views();
    end


    function local_on_single_axes_clicked(~, ~)
        if isempty(state.physicalQE)
            return
        end

        % --- Manual peak picking mode ---
        if state.manualPickArmed
            try
                current_point = ui.SingleAxes.CurrentPoint;
                energy_value = current_point(1, 1);
                [~, energy_axis] = local_energy_mask(state.physicalQE);
                if isempty(energy_axis)
                    return
                end
                energy_value = min(max(energy_value, energy_axis(1)), energy_axis(end));
                abs_q = abs(state.physicalQE.q_Ainv(state.selectedQIndex));
                state.manual_points(end+1, :) = [abs_q, energy_value];
                ui.PtsLabel.Text = sprintf("%d pts", size(state.manual_points, 1));
                local_update_all_views();
            catch ME
                local_show_error(ME, false);
            end
            return
        end

        % --- Guess picking mode: click on spectrum to add peak position ---
        if state.guessPickArmed
            try
                current_point = ui.SingleAxes.CurrentPoint;
                energy_value = round(current_point(1, 1));
                if energy_value > 0
                    % Append to GuessField
                    existing = strtrim(ui.GuessField.Value);
                    if isempty(existing)
                        ui.GuessField.Value = sprintf('%d', energy_value);
                    else
                        ui.GuessField.Value = sprintf('%s, %d', existing, energy_value);
                    end
                    % Mark on plot
                    ax = ui.SingleAxes;
                    hold(ax, 'on');
                    yl = ylim(ax);
                    plot(ax, [energy_value energy_value], yl, '--', ...
                        'Color', [0.9 0.2 0.2 0.6], 'LineWidth', 1.5, ...
                        'Tag', 'guess_marker');
                    hold(ax, 'off');
                    ui.SeedInfoLabel.Text = sprintf('Picked: %s', ui.GuessField.Value);
                end
            catch ME
                local_show_error(ME, false);
            end
            return
        end
    end


    function local_sync_action_buttons()
        has_dataset = ~isempty(state.dataset);
        has_raw = local_has_raw_source();

        ui.CropButton.Enable = local_on_off(has_raw);
        ui.AlignButton.Enable = local_on_off(has_raw);
        ui.BinButton.Enable = local_on_off(has_raw);
        ui.BuildViewsButton.Enable = local_on_off(has_dataset);
        ui.ResetButton.Enable = local_on_off(has_dataset);
        ui.UndoButton.Enable = local_on_off(numel(state.opHistory) >= 2);
        ui.QNormCheckbox.Enable = local_on_off(has_raw);

        if ~has_raw
            ui.QNormCheckbox.Value = false;
        end

        if ~has_dataset
            ui.SourceLabel.Text = "No data loaded";
            ui.DqDisplay.Value = "";
            ui.ResetStageDropdown.Items = {'raw'};
            ui.ResetStageDropdown.Value = 'raw';
            return
        end

        ui.SourceLabel.Text = state.dataset.label;
        ui.DqDisplay.Value = sprintf("%.6f", local_current_dq());
        local_sync_reset_stage_items();
    end


    function local_sync_reset_stage_items()
        items = strings(0, 1);
        if ~isempty(state.eq3dQE)
            items(end + 1, 1) = "eq3d"; %#ok<AGROW>
        else
            items(end + 1, 1) = "raw"; %#ok<AGROW>
            if ~isempty(state.cropped4d)
                items(end + 1, 1) = "crop"; %#ok<AGROW>
            end
            if ~isempty(state.aligned4d)
                items(end + 1, 1) = "align"; %#ok<AGROW>
            end
            if ~isempty(state.binned4d)
                items(end + 1, 1) = "bin"; %#ok<AGROW>
            end
        end

        item_cells = cellstr(items);
        ui.ResetStageDropdown.Items = item_cells;
        if ~any(strcmp(ui.ResetStageDropdown.Value, item_cells))
            ui.ResetStageDropdown.Value = item_cells{end};
        end
    end


    function local_sync_controls_from_qe(qe)
        if isempty(qe)
            return
        end

        [default_ref_q_min, default_ref_q_max] = local_default_reference_band(qe);
        ui.DqDisplay.Value = sprintf("%.6f", qe.dq_Ainv);
        ui.QStartField.Value = qe.q_Ainv(1);
        ui.QEndField.Value = qe.q_Ainv(end);
        ui.QStepField.Value = qe.dq_Ainv;
        ui.EnergyMinField.Value = qe.energy_meV(1);
        ui.EnergyMaxField.Value = qe.energy_meV(end);
        ui.OffsetField.Value = local_default_offset(qe);
        ui.RefAbsQMinField.Value = default_ref_q_min;
        ui.RefAbsQMaxField.Value = default_ref_q_max;
    end


    function [ref_q_min, ref_q_max] = local_default_reference_band(qe)
        q_abs_max = max(abs(double(qe.q_Ainv(:))));
        ref_q_min = min(max(2 * double(qe.dq_Ainv), double(qe.dq_Ainv)), q_abs_max);
        ref_q_max = min(max(2 * ref_q_min, ref_q_min + double(qe.dq_Ainv)), q_abs_max);

        if ref_q_max <= ref_q_min
            ref_q_min = min(max(double(qe.dq_Ainv), eps), q_abs_max);
            ref_q_max = min(max(ref_q_min + double(qe.dq_Ainv), 2 * ref_q_min), q_abs_max);
        end
    end


    function local_sync_manual_y_state()
        enabled = local_on_off(~ui.AutoYCheckbox.Value);
        ui.YMinField.Editable = enabled;
        ui.YMaxField.Editable = enabled;
    end


    function local_update_stage_label()
        if isempty(state.dataset)
            ui.StageLabel.Text = "Stage: none";
            return
        end
        ui.StageLabel.Text = sprintf("Stage: %s", char(state.activeStage));
    end


    function local_update_info_label()
        if isempty(state.dataset)
            ui.InfoLabel.Text = "Top-left shows the original Physical q-E map. Top-right shows the normalized off-axis component reconstructed from manual symmetric off-axis reference bands.";
            return
        end

        if strcmpi(ui.ViewModeDropdown.Value, "Normalized")
            if isempty(state.comparisonQE)
                ui.InfoLabel.Text = "Off-axis component not built yet. Press Build Views to reconstruct the center from manual symmetric off-axis reference bands.";
                return
            end

            if local_has_raw_source() && ui.QNormCheckbox.Value
                ui.InfoLabel.Text = sprintf("Off-axis component uses manual symmetric reference bands [%.4f, %.4f] 1/A, reconstructs |q| < %.4f 1/A, and optionally uses experimental q equalization as the base.", ...
                    ui.RefAbsQMinField.Value, ui.RefAbsQMaxField.Value, ui.RefAbsQMinField.Value);
            else
                ui.InfoLabel.Text = sprintf("Off-axis component uses manual symmetric reference bands [%.4f, %.4f] 1/A and reconstructs |q| < %.4f 1/A before q-wise normalization.", ...
                    ui.RefAbsQMinField.Value, ui.RefAbsQMaxField.Value, ui.RefAbsQMinField.Value);
            end
            return
        end

        if ~strcmp(local_trace_mode(), "display")
            ui.InfoLabel.Text = "Peak-hunting trace mode uses background removal or curvature. These views are for finding weak shoulders, not for reporting physical intensity.";
        else
            ui.InfoLabel.Text = "Physical view keeps q=0 on-axis signal stronger when present. That ZLP enhancement is expected and intentionally preserved.";
        end
    end


    function tf = local_has_raw_source()
        tf = ~isempty(state.dataset) ...
            && isfield(state.dataset, "source_kind") ...
            && strcmpi(char(state.dataset.source_kind), "raw") ...
            && ~isempty(state.raw4d);
    end





    function dq_Ainv = local_current_dq()
        dq_Ainv = ui.DqOverrideField.Value;
        if ~isfinite(dq_Ainv) || dq_Ainv <= 0
            if ~isempty(state.dataset) && isfield(state.dataset, "dq_Ainv")
                dq_Ainv = double(state.dataset.dq_Ainv);
            else
                dq_Ainv = 0.005;
            end
        end
    end





    function qe = local_get_reference_qe()
        if ~isempty(state.physicalQE)
            qe = state.physicalQE;
        elseif ~isempty(state.eq3dQE)
            qe = state.eq3dQE;
        elseif ~isempty(state.comparisonQE)
            qe = state.comparisonQE;
        else
            qe = [];
        end
    end


    function qe = local_get_active_qe()
        switch lower(char(ui.ViewModeDropdown.Value))
            case "normalized"
                if ~isempty(state.comparisonQE)
                    qe = state.comparisonQE;
                else
                    qe = state.physicalQE;
                end
            otherwise
                if ~isempty(state.physicalQE)
                    qe = state.physicalQE;
                else
                    qe = state.eq3dQE;
                end
        end
    end


    function qe = local_reaxis_qe(qe, dq_Ainv)
        if isempty(qe)
            return
        end
        qe.dq_Ainv = dq_Ainv;
        if isfield(qe, "q_zero_index") && isfinite(qe.q_zero_index)
            q_zero_index = double(qe.q_zero_index);
        else
            q_zero_index = round((size(qe.intensity, 2) + 1) / 2);
            qe.q_zero_index = q_zero_index;
        end
        qe.q_Ainv = ((1:size(qe.intensity, 2)) - q_zero_index) * dq_Ainv;
        qe.q_abs_Ainv = abs(qe.q_Ainv);
    end


    function working4d = local_get_working_4d()
        if ~isempty(state.binned4d)
            working4d = state.binned4d;
        elseif ~isempty(state.aligned4d)
            working4d = state.aligned4d;
        elseif ~isempty(state.cropped4d)
            working4d = state.cropped4d;
        else
            working4d = state.raw4d;
        end
    end


    function local_rebuild_comparison(show_progress)
        if nargin < 1
            show_progress = false;
        end

        if isempty(state.dataset)
            return
        end

        dialog = [];
        cleanup_handle = [];
        try
            if show_progress
                [dialog, progress_callback] = local_open_progress("Building physical and comparison views");
                cleanup_handle = onCleanup(@() local_close_dialog(dialog));
            else
                progress_callback = [];
            end

            dq_Ainv = local_current_dq();
            switch string(state.dataset.source_kind)
                case "eq3d"
                    state.eq3dQE = local_reaxis_qe(state.eq3dQE, dq_Ainv);
                    state.physicalQE = state.eq3dQE;
                    base_qe = state.eq3dQE;
                    state.activeStage = "eq3d";
                    local_report_progress(progress_callback, 0.35, "Preparing eq3D comparison");

                case "raw"
                    working4d = local_get_working_4d();
                    local_report_progress(progress_callback, 0.2, "Projecting physical q-E");
                    state.physicalQE = project_4d_to_qe_struct( ...
                        working4d, ...
                        dq_Ainv, ...
                        source_path=string(state.dataset.source_path), ...
                        source_kind="raw", ...
                        label=string(state.dataset.label), ...
                        view_kind="physical", ...
                        stage_name=string(state.activeStage));

                    base_qe = state.physicalQE;
                    if ui.QNormCheckbox.Value
                        local_report_progress(progress_callback, 0.45, "Applying experimental q equalization");
                        figures_before = findall(groot, "Type", "figure");
                        normalized4d = Normalize_qdim_4D(working4d);
                        local_close_new_figures(figures_before);
                        base_qe = project_4d_to_qe_struct( ...
                            normalized4d, ...
                            dq_Ainv, ...
                            source_path=string(state.dataset.source_path), ...
                            source_kind="raw", ...
                            label=string(state.dataset.label) + " | q-equalized base", ...
                            view_kind="comparison", ...
                            stage_name=string(state.activeStage), ...
                            note="comparison q-equalized base");
                    end

                otherwise
                    error("interactive_qe_browser:UnsupportedBuildKind", ...
                        "Unsupported source kind '%s'.", state.dataset.source_kind);
            end

            local_report_progress(progress_callback, 0.8, "Normalizing comparison view");
            state.comparisonQE = build_comparison_qe( ...
                base_qe, ...
                [ui.ZlpMinField.Value, ui.ZlpMaxField.Value], ...
                "zlp-area", ...
                string(ui.SmoothingModeDropdown.Value), ...
                ui.SmoothingWidthField.Value, ...
                ui.RefAbsQMinField.Value, ...
                ui.RefAbsQMaxField.Value);
            state.comparisonQE = local_reaxis_qe(state.comparisonQE, dq_Ainv);
            local_report_progress(progress_callback, 1.0, "Views ready");

            clear cleanup_handle
            local_close_dialog(dialog);
        catch ME
            clear cleanup_handle
            local_close_dialog(dialog);
            rethrow(ME);
        end
    end


    function [dialog, callback] = local_open_progress(title_text)
        dialog = uiprogressdlg(ui.Figure, ...
            "Title", char(title_text), ...
            "Message", "Working...", ...
            "Value", 0, ...
            "Cancelable", "off", ...
            "Indeterminate", "off");
        callback = @local_update_progress;

        function local_update_progress(fraction, message)
            if isempty(dialog) || ~isvalid(dialog)
                return
            end
            dialog.Value = min(max(fraction, 0), 1);
            dialog.Message = char(message);
            drawnow limitrate
        end
    end


    function local_close_dialog(dialog)
        if ~isempty(dialog) && isvalid(dialog)
            close(dialog);
        end
    end


    function local_close_new_figures(figures_before)
        figures_after = findall(groot, "Type", "figure");
        for idx = 1:numel(figures_after)
            fig_handle = figures_after(idx);
            if fig_handle == ui.Figure
                continue
            end
            if any(fig_handle == figures_before)
                continue
            end
            if isvalid(fig_handle)
                close(fig_handle);
            end
        end
    end


    function local_show_error(exception_object, throw_on_error)
        if nargin < 2
            throw_on_error = false;
        end

        if throw_on_error
            rethrow(exception_object);
        end

        if isvalid(ui.Figure)
            uialert(ui.Figure, exception_object.message, "q-E Browser Error");
        else
            warning("%s", exception_object.message);
        end
    end


    function view_name = local_current_view_name()
        qe = local_get_active_qe();
        if isempty(qe)
            view_name = char(ui.ViewModeDropdown.Value);
        elseif strcmpi(qe.view_kind, "comparison")
            view_name = "Normalized";
        else
            view_name = "Physical";
        end
    end


    function title_text = local_make_stage_title(qe)
        title_text = sprintf("%s | stage %s", qe.label, qe.stage_name);
    end


    function [mask, energy_axis] = local_energy_mask(qe)
        energy_min = min(ui.EnergyMinField.Value, ui.EnergyMaxField.Value);
        energy_max = max(ui.EnergyMinField.Value, ui.EnergyMaxField.Value);
        energy_axis_full = double(qe.energy_meV(:));
        mask = energy_axis_full >= energy_min & energy_axis_full <= energy_max;

        if ~any(mask)
            [~, nearest_idx] = min(abs(energy_axis_full - mean([energy_min, energy_max])));
            mask(nearest_idx) = true;
        end

        energy_axis = energy_axis_full(mask);
    end


    function [mask, q_axis] = local_q_mask(qe)
        q_min = min(ui.QStartField.Value, ui.QEndField.Value);
        q_max = max(ui.QStartField.Value, ui.QEndField.Value);
        q_axis_full = double(qe.q_Ainv(:)).';
        mask = q_axis_full >= q_min & q_axis_full <= q_max;

        if ~any(mask)
            [~, nearest_idx] = min(abs(q_axis_full - mean([q_min, q_max])));
            mask(nearest_idx) = true;
        end

        q_axis = q_axis_full(mask);
    end


    function map_values = local_build_map_values(qe, energy_mask, q_mask, mode_name)
        if nargin < 4
            mode_name = local_trace_mode();
        end
        spectra = double(qe.intensity(energy_mask, q_mask));
        energy_axis = double(qe.energy_meV(energy_mask));
        map_values = zeros(size(spectra));

        switch mode_name
            case "display"
                map_values = local_apply_smoothing(spectra);

            otherwise
                for col = 1:size(spectra, 2)
                    features = local_build_spectrum_features(energy_axis, spectra(:, col));
                    map_values(:, col) = local_trace_values(features, mode_name);
                end
        end
    end


    function [display_map, display_clim] = local_prepare_map_display(map_values, mode_name, display_style)
        if nargin < 2
            mode_name = local_trace_mode();
        end
        if nargin < 3
            display_style = "physical";
        end
        display_map = double(map_values);
        finite_values = display_map(isfinite(display_map));
        if isempty(finite_values)
            display_clim = [0 1];
            display_map(:) = 0;
            return
        end

        switch mode_name
            case "display"
                if strcmpi(display_style, "normalized")
                    lower_limit = max(min(finite_values), 0);
                    upper_limit = prctile(finite_values, 99.8);
                    if ~isfinite(upper_limit) || upper_limit <= lower_limit
                        upper_limit = max(finite_values);
                    end
                    if ~isfinite(upper_limit) || upper_limit <= lower_limit
                        upper_limit = lower_limit + 1;
                    end
                    display_map = max(display_map, lower_limit);
                    display_clim = [lower_limit, upper_limit];
                else
                    lower_limit = prctile(finite_values, 0.5);
                    upper_limit = prctile(finite_values, 99.8);
                    if ~isfinite(lower_limit) || ~isfinite(upper_limit) || upper_limit <= lower_limit
                        lower_limit = min(finite_values);
                        upper_limit = max(finite_values);
                    end
                    scale = max(upper_limit - lower_limit, eps);
                    clipped = max(display_map - lower_limit, 0);
                    stretch = 8;
                    display_map = asinh(stretch * clipped / scale) / asinh(stretch);
                    display_clim = [0 1];
                end

            otherwise
                center_value = median(finite_values, "omitnan");
                spread = prctile(abs(finite_values - center_value), 99);
                if ~isfinite(spread) || spread <= eps
                    spread = max(abs(finite_values - center_value));
                end
                if ~isfinite(spread) || spread <= eps
                    spread = 1;
                end
                stretch = 4;
                display_map = asinh(stretch * (display_map - center_value) / spread) / asinh(stretch);
                display_clim = [-1 1];
        end
    end


    function output = local_apply_smoothing(input_data)
        mode_name = lower(char(ui.SmoothingModeDropdown.Value));
        width = max(1, round(ui.SmoothingWidthField.Value));
        output = double(input_data);

        if strcmp(mode_name, "off")
            return
        end

        if isvector(output)
            output = local_smooth_vector(output, mode_name, width);
            return
        end

        for q_index = 1:size(output, 2)
            output(:, q_index) = local_smooth_vector(output(:, q_index), mode_name, width);
        end
    end


    function vector_out = local_smooth_vector(vector_in, mode_name, width)
        vector_out = double(vector_in(:));
        switch mode_name
            case "gaussian"
                window = min(numel(vector_out), max(1, width));
                vector_out = smoothdata(vector_out, "gaussian", window);

            case "sgolay"
                window = min(numel(vector_out), max(3, width));
                if mod(window, 2) == 0
                    window = max(3, window - 1);
                end
                if window < 3
                    window = min(numel(vector_out), max(1, width));
                    vector_out = smoothdata(vector_out, "movmean", window);
                else
                    vector_out = smoothdata(vector_out, "sgolay", window);
                end

            otherwise
                vector_out = vector_in(:);
        end
    end


    function y_label = local_single_ylabel()
        switch local_trace_mode()
            case "background-subtracted"
                y_label = "Residual intensity";
            case "curvature"
                y_label = "Curvature signal";
            otherwise
                if strcmpi(local_current_view_name(), "Normalized")
                    y_label = "Off-axis component intensity";
                else
                    y_label = "Physical intensity";
                end
        end
    end


    function mode_name = local_trace_mode()
        mode_name = lower(char(ui.TraceModeDropdown.Value));
    end


    function y_scale = local_trace_y_scale()
        if strcmp(local_trace_mode(), "display")
            y_scale = char(ui.YScaleDropdown.Value);
        else
            y_scale = 'linear';
        end
    end


    function features = local_build_spectrum_features(energy_axis, raw_spectrum)
        smooth_spectrum = local_apply_smoothing(raw_spectrum);
        baseline = local_asls_baseline(smooth_spectrum, ui.BaselineLambdaField.Value, 0.01, 12);
        residual = smooth_spectrum - baseline;
        curvature = -local_second_derivative(energy_axis, smooth_spectrum);

        features = struct();
        features.energy = energy_axis(:);
        features.raw = raw_spectrum(:);
        features.smooth = smooth_spectrum(:);
        features.baseline = baseline(:);
        features.residual = residual(:);
        features.curvature = curvature(:);
    end


    function values = local_trace_values(features, mode_name)
        switch mode_name
            case "background-subtracted"
                values = features.residual;
            case "curvature"
                values = features.curvature;
            otherwise
                values = features.smooth;
        end
    end





    function baseline = local_asls_baseline(y, lambda_value, asymmetry, iterations)
        y = double(y(:));
        n = numel(y);
        if n < 3
            baseline = y;
            return
        end

        lambda_value = max(double(lambda_value), 1);
        asymmetry = min(max(double(asymmetry), 1e-4), 0.49);
        iterations = max(1, round(iterations));

        d = diff(speye(n), 2);
        penalty = lambda_value * (d' * d);
        weights = ones(n, 1);

        for iter = 1:iterations %#ok<NASGU>
            w = spdiags(weights, 0, n, n);
            baseline = (w + penalty) \ (weights .* y);
            weights = asymmetry * (y > baseline) + (1 - asymmetry) * (y <= baseline);
        end
    end


    function second_derivative = local_second_derivative(x, y)
        x = double(x(:));
        y = double(y(:));
        if numel(y) < 3
            second_derivative = zeros(size(y));
            return
        end
        first_derivative = gradient(y, x);
        second_derivative = gradient(first_derivative, x);
    end


    function dE = local_energy_step(energy_axis)
        energy_axis = double(energy_axis(:));
        if numel(energy_axis) < 2
            dE = 1;
            return
        end
        dE = median(abs(diff(energy_axis)), "omitnan");
        if ~isfinite(dE) || dE <= 0
            dE = 1;
        end
    end





    function local_apply_y_limits(ax, values)
        finite_values = values(isfinite(values));
        if isempty(finite_values)
            return
        end

        if ui.AutoYCheckbox.Value
            if strcmpi(ax.YScale, "log")
                finite_values = finite_values(finite_values > 0);
                if isempty(finite_values)
                    finite_values = eps;
                end
            end

            lower_limit = min(finite_values);
            upper_limit = max(finite_values);
            if ~isfinite(lower_limit) || ~isfinite(upper_limit)
                return
            end
            if upper_limit <= lower_limit
                upper_limit = lower_limit + max(abs(lower_limit) * 0.05, 1);
            end
        else
            lower_limit = min(ui.YMinField.Value, ui.YMaxField.Value);
            upper_limit = max(ui.YMinField.Value, ui.YMaxField.Value);
            if strcmpi(ax.YScale, "log")
                lower_limit = max(lower_limit, eps);
                upper_limit = max(upper_limit, lower_limit * (1 + 1e-3));
            end
        end

        ylim(ax, [lower_limit, upper_limit]);
    end


    function q_indices = local_waterfall_indices(qe)
        q_axis = double(qe.q_Ainv(:)).';
        q_min = min(ui.QStartField.Value, ui.QEndField.Value);
        q_max = max(ui.QStartField.Value, ui.QEndField.Value);
        q_step = max(ui.QStepField.Value, eps);

        targets = q_min:q_step:q_max;
        if isempty(targets)
            targets = q_min;
        end

        q_indices = zeros(1, numel(targets));
        for idx = 1:numel(targets)
            [~, q_indices(idx)] = min(abs(q_axis - targets(idx)));
        end
        q_indices = unique(q_indices, "stable");

        if isempty(q_indices)
            [~, q_indices] = min(abs(q_axis - mean([q_min, q_max])));
        end
    end


    function [diag_q, raw_metric, scale_factor] = local_get_zlp_diagnostic_data()
        if ~isempty(state.comparisonQE) ...
                && isfield(state.comparisonQE, "raw_window_metric") ...
                && isfield(state.comparisonQE, "normalization_factor")
            diag_q = state.comparisonQE.q_Ainv;
            raw_metric = double(state.comparisonQE.raw_window_metric(:)).';
            scale_factor = double(state.comparisonQE.normalization_factor(:)).';
            return
        end

        qe = local_get_reference_qe();
        if isempty(qe)
            diag_q = [];
            raw_metric = [];
            scale_factor = [];
            return
        end

        diag_q = qe.q_Ainv;
        raw_metric = local_compute_window_metric(qe);
        scale_factor = ones(size(raw_metric));
    end


    function metric = local_compute_window_metric(qe)
        [mask, energy_axis] = local_energy_mask_for_window(qe, [ui.ZlpMinField.Value, ui.ZlpMaxField.Value]);
        metric = zeros(1, size(qe.intensity, 2));
        for q_index = 1:size(qe.intensity, 2)
            spectrum = max(double(qe.intensity(mask, q_index)), 0);
            if numel(energy_axis) > 1
                metric(q_index) = trapz(energy_axis, spectrum);
            else
                metric(q_index) = sum(spectrum);
            end
        end
    end


    function [mask, energy_axis] = local_energy_mask_for_window(qe, window)
        window = sort(double(window(:))).';
        energy_axis_full = double(qe.energy_meV(:));
        mask = energy_axis_full >= window(1) & energy_axis_full <= window(2);
        if ~any(mask)
            [~, nearest_idx] = min(abs(energy_axis_full - mean(window)));
            mask(nearest_idx) = true;
        end
        energy_axis = energy_axis_full(mask);
    end




    function value = local_default_offset(qe)
        sample = double(qe.intensity);
        finite_values = sample(isfinite(sample));
        if isempty(finite_values)
            value = 0.25;
            return
        end
        span = max(finite_values) - min(finite_values);
        if span <= 0 || ~isfinite(span)
            value = max(abs(max(finite_values)), 1) * 0.05;
        else
            value = span * 0.08;
        end
        value = max(value, eps);
    end


    function color_value = local_signed_q_color(q_value, q_axis)
        max_abs_q = max(abs(q_axis));
        if ~isfinite(max_abs_q) || max_abs_q <= eps
            color_value = [0 0 0];
            return
        end

        t = max(min(q_value / max_abs_q, 1), -1);
        if t >= 0
            color_value = (1 - t) * [0 0 0] + t * [0.05 0.9 0.05];
        else
            color_value = (1 + t) * [0 0 0] + (-t) * [0.95 0.05 0.05];
        end
    end


    function local_report_progress(callback, fraction, message)
        if isa(callback, "function_handle")
            callback(fraction, message);
        end
    end


    function value = local_on_off(tf)
        if tf
            value = "on";
        else
            value = "off";
        end
    end


    function qe_out = local_apply_area_normalization(qe_in)
        qe_out = qe_in;
        energy_axis = double(qe_in.energy_meV(:));
        norm_min = ui.AreaNormMinField.Value;
        norm_max = ui.AreaNormMaxField.Value;
        intensity = double(qe_in.intensity);
        method = char(ui.NormMethodDropdown.Value);

        if strcmpi(method, 'ZLP Peak')
            % ZLP Peak normalization: divide by ZLP peak height.
            % Preserves q-dependent spectral weight for A(q) studies.
            zlp_mask = energy_axis >= norm_min & energy_axis <= norm_max;
            if ~any(zlp_mask)
                [~, nearest] = min(abs(energy_axis));
                zlp_mask(nearest) = true;
            end
            for qi = 1:size(intensity, 2)
                zlp_peak = max(intensity(zlp_mask, qi));
                if isfinite(zlp_peak) && zlp_peak > eps
                    intensity(:, qi) = intensity(:, qi) ./ zlp_peak;
                end
            end
        else
            % Area normalization: divide by integrated spectral weight.
            % WARNING: destroys q-dependent A(q) scaling information.
            norm_mask = energy_axis >= norm_min & energy_axis <= norm_max;
            if ~any(norm_mask)
                norm_mask = true(size(energy_axis));
            end
            norm_energy = energy_axis(norm_mask);
            for qi = 1:size(intensity, 2)
                window_spec = max(intensity(norm_mask, qi), 0);
                if numel(norm_energy) > 1
                    area = trapz(norm_energy, window_spec);
                else
                    area = sum(window_spec);
                end
                if isfinite(area) && area > eps
                    intensity(:, qi) = intensity(:, qi) ./ area;
                end
            end
        end
        qe_out.intensity = intensity;
    end


    function qe_out = local_apply_bg_subtraction(qe_in)
        % EELS background subtraction with single low-energy window fitting.
        %
        % Uses ONLY the low-energy ZLP tail region [50, win1_hi] as anchor.
        % The high-energy window was removed because it can contain
        % interband transitions or plasmon tails, contaminating the
        % background estimate and causing over-subtraction.
        %
        % Safety cap: background is limited to at most 90% of the
        % smoothed local signal — prevents complete erasure at high |q|.
        qe_out = qe_in;
        energy_axis = double(qe_in.energy_meV(:));
        intensity = double(qe_in.intensity);
        method = char(ui.BgMethodDropdown.Value);

        % --- Define fit window ---
        % Low-energy window only: [fit_lo, win1_hi]
        % Between ZLP tail and the onset of loss features.
        % This region is physically clean (no plasmon/interband signal).
        fit_lo = 50;
        win1_hi = 300;  % conservative upper bound below typical Bi loss onset

        % Single window - low-energy ZLP tail only
        fit_mask = energy_axis >= fit_lo & energy_axis <= win1_hi;

        if nnz(fit_mask) < 5
            return
        end

        e_fit = energy_axis(fit_mask);

        % Only subtract for energies above fit_lo (avoid ZLP / negative E)
        sub_mask = energy_axis > fit_lo;
        e_sub = energy_axis(sub_mask);

        for qi = 1:size(intensity, 2)
            spec = intensity(:, qi);
            s_fit = spec(fit_mask);
            s_fit_pos = max(s_fit, eps);

            try
                switch method
                    case 'Power'
                        valid = e_fit > 0 & s_fit_pos > 0;
                        if nnz(valid) < 3, continue; end
                        p = polyfit(log(e_fit(valid)), log(s_fit_pos(valid)), 1);
                        bg = exp(polyval(p, log(e_sub)));

                    case 'ExpPoly3'
                        valid = s_fit_pos > 0;
                        if nnz(valid) < 5, continue; end
                        p = polyfit(e_fit(valid), log(s_fit_pos(valid)), 3);
                        bg = exp(polyval(p, e_sub));

                    case 'Pearson'
                        valid = e_fit > 0 & s_fit_pos > 0;
                        if nnz(valid) < 4, continue; end
                        p = polyfit(log(e_fit(valid)), log(s_fit_pos(valid)), 2);
                        bg = exp(polyval(p, log(e_sub)));

                    otherwise
                        continue
                end

                bg = real(bg(:));
                bg(~isfinite(bg)) = 0;
                bg = max(bg, 0);

                % --- Safety cap: prevent over-subtraction ---
                % Background should not exceed 90% of the local smoothed
                % signal. This prevents complete erasure of weak spectra
                % at high |q| where the signal is inherently faint.
                local_signal = max(spec(sub_mask), 0);
                local_smooth = smoothdata(local_signal, 'gaussian', ...
                    max(5, round(numel(local_signal) * 0.05)));
                bg = min(bg, local_smooth * 0.9);

                % Subtract — allow negative residuals for honest statistics
                intensity(sub_mask, qi) = spec(sub_mask) - bg;
            catch
            end
        end
        qe_out.intensity = intensity;
    end

    function qe_out = local_apply_deconvolution(qe_in)
        % Lucy-Richardson deconvolution using the ZLP (reference band) as PSF.
        % Only the narrow ZLP peak region is used, NOT the full spectrum,
        % to avoid including loss features in the instrument response.
        qe_out = qe_in;
        energy_axis = double(qe_in.energy_meV(:));
        intensity = double(qe_in.intensity);
        n_q = size(intensity, 2);
        n_iter = round(ui.DeconvIterField.Value);

        % --- Extract PSF from q≈0 channel only ---
        % Physics: the ZLP (instrument response) is independent of q.
        % Using only the q=0 channel (±1 pixel) avoids contaminating the
        % PSF with plasmon loss features present at finite q.
        if isfield(qe_in, 'q_zero_index') && isfinite(qe_in.q_zero_index)
            q0_idx = round(qe_in.q_zero_index);
        else
            % Fallback: find the channel with the highest ZLP peak
            [~, q0_idx] = max(max(intensity, [], 1));
        end
        q0_idx = max(1, min(q0_idx, n_q));
        % Average ±1 pixel around q=0 for noise reduction
        q0_range = max(1, q0_idx - 1):min(n_q, q0_idx + 1);
        ref_spectrum = mean(intensity(:, q0_range), 2, 'omitnan');
        ref_spectrum = max(ref_spectrum, 0);

        % --- Use ONLY the ZLP peak region as PSF ---
        % Find ZLP peak position (should be near 0 meV)
        [~, zlp_idx] = max(ref_spectrum);
        zlp_energy = energy_axis(zlp_idx);
        % Narrow window: ±100 meV around ZLP — pure instrument response
        zlp_half_width = 100;  % meV
        zlp_mask = energy_axis >= (zlp_energy - zlp_half_width) & ...
                   energy_axis <= (zlp_energy + zlp_half_width);
        psf = ref_spectrum(zlp_mask);
        psf = max(psf, 0);
        psf_sum = sum(psf);
        if psf_sum > 0
            psf = psf / psf_sum;
        else
            return  % cannot deconvolve with zero PSF
        end

        % --- Deconvolve each q-channel ---
        for qi = 1:n_q
            spec = intensity(:, qi);
            spec = max(spec, 0);  % LR requires non-negative
            if max(spec) <= 0
                continue
            end
            try
                intensity(:, qi) = deconvlucy(spec, psf, n_iter);
            catch
                % Fallback: manual Lucy-Richardson if toolbox unavailable
                intensity(:, qi) = local_lr_deconv(spec, psf, n_iter);
            end
        end
        qe_out.intensity = intensity;
    end


    function result = local_lr_deconv(signal, psf, n_iter)
        % Manual Lucy-Richardson deconvolution (no toolbox dependency)
        signal = double(signal(:));
        psf = double(psf(:));
        estimate = signal;
        psf_flip = flipud(psf);
        for iter = 1:n_iter
            blurred = conv(estimate, psf, 'same');
            blurred(blurred < eps) = eps;
            ratio = signal ./ blurred;
            correction = conv(ratio, psf_flip, 'same');
            estimate = estimate .* correction;
        end
        result = max(estimate, 0);
    end


    function qe_out = local_apply_denoise(qe_in)
        % Denoise the q-E intensity map.
        % Two methods:
        %   Wiener2D — 2D adaptive Wiener filter on the q-E map
        %   SavGol   — 1D Savitzky-Golay filter per q-channel along energy
        qe_out = qe_in;
        intensity = double(qe_in.intensity);
        method = char(ui.DenoiseMethodDropdown.Value);

        switch method
            case 'Wiener2D'
                % 2D adaptive Wiener filter (wiener2)
                sigma_input = ui.DenoiseSigmaField.Value;
                if sigma_input <= 0
                    % Auto-estimate noise via MAD (median absolute deviation)
                    noise_est = median(abs(intensity(:) - median(intensity(:)))) / 0.6745;
                else
                    noise_est = sigma_input;
                end
                try
                    % wiener2 from Image Processing Toolbox
                    denoised = wiener2(intensity, [3 5], noise_est^2);
                catch
                    % Fallback: simple 2D moving average
                    kernel = ones(3, 5) / 15;
                    denoised = conv2(intensity, kernel, 'same');
                end
                qe_out.intensity = denoised;

            case 'SavGol'
                % 1D Savitzky-Golay filter per q-channel
                sg_order = round(ui.SGOrderField.Value);
                sg_framelen = round(ui.SGFrameLenField.Value);
                % Frame length must be odd
                if mod(sg_framelen, 2) == 0
                    sg_framelen = sg_framelen + 1;
                end
                % Order must be < frame length
                sg_order = min(sg_order, sg_framelen - 1);
                n_q = size(intensity, 2);
                for qi = 1:n_q
                    spec = intensity(:, qi);
                    if numel(spec) >= sg_framelen
                        try
                            intensity(:, qi) = sgolayfilt(spec, sg_order, sg_framelen);
                        catch
                            % Fallback: moving average
                            intensity(:, qi) = movmean(spec, sg_framelen);
                        end
                    end
                end
                qe_out.intensity = intensity;
        end
    end

end
