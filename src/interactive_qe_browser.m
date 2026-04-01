function app = interactive_qe_browser(startPath, historyFile)
arguments
    startPath {mustBeTextScalar} = ""
    historyFile {mustBeTextScalar} = ""
end

state = local_initial_state();

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
        % Preprocessing cache (avoid redundant SVD/FFT on every click)
        state_out.ppCache = struct('hash', '', 'qe', [], 'qe_pre', [], 'bg_diag', []);
    end


    function ui_handles = local_build_ui()
        % Build callbacks struct — maps to local nested functions
        cb = struct();
        cb.on_load_data = @local_on_load_data;
        cb.on_build_views = @local_on_build_views;
        cb.on_undo = @local_on_undo;
        cb.on_reset_stage = @local_on_reset_stage;
        cb.on_display_change = @local_on_display_change;
        cb.on_normalization_change = @local_on_normalization_change;
        cb.on_dq_override_changed = @local_on_dq_override_changed;
        cb.on_qnorm_changed = @local_on_qnorm_changed;
        cb.on_norm_method_changed = @local_on_norm_method_changed;
        cb.on_auto_y_changed = @local_on_auto_y_changed;
        cb.on_pick_peaks = @local_on_pick_peaks;
        cb.on_undo_pt = @local_on_undo_pt;
        cb.on_clear_pts = @local_on_clear_pts;
        cb.on_save_pts = @local_on_save_pts;
        cb.on_load_pts = @local_on_load_pts;
        cb.on_split_branches = @local_on_split_branches;
        cb.on_fit_model = @local_on_fit_model;
        cb.on_export = @local_on_export;
        cb.on_auto_fit_dispersion = @local_on_auto_fit_dispersion;
        cb.on_fit_spectrum = @local_on_fit_spectrum;
        cb.on_accept_fit = @local_on_accept_fit;
        cb.on_show_gamma = @local_on_show_gamma;
        cb.on_pick_guesses = @local_on_pick_guesses;
        cb.on_reassign_points = @local_on_reassign_points;
        cb.on_fit_dispersion = @local_on_fit_dispersion;
        cb.on_export_dispersion = @local_on_export_dispersion;
        cb.on_history_select = @local_on_history_select;
        cb.on_save_history = @local_on_save_history;
        cb.on_load_history = @local_on_load_history;
        cb.on_clear_history = @local_on_clear_history;

        % Delegate all widget creation to external module
        ui_handles = qe_browser_ui(cb);
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
                dq_val = NaN;  % default: let load_raw_session prompt
                if contains(lower_path, "20w")
                    dq_val = 0.0025;
                    ui.DqOverrideField.Value = 0.0025;
                elseif contains(lower_path, "10w")
                    dq_val = 0.005;
                    ui.DqOverrideField.Value = 0.005;
                end
                dataset = load_qe_dataset(full_path, dq_val, ...
                    q_crop=[q_lo q_hi]);
                ui.DqOverrideField.Value = dataset.dq_Ainv;
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
        snap.bgWin1Lo = ui.BgWin1LoField.Value;
        snap.bgWin1Hi = ui.BgWin1HiField.Value;
        snap.bgDual = ui.BgDualCheckbox.Value;
        snap.bgWin2Lo = ui.BgWin2LoField.Value;
        snap.bgWin2Hi = ui.BgWin2HiField.Value;
        snap.bgIter = ui.BgIterCheckbox.Value;
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
        if isfield(snap, 'bgWin1Lo')
            ui.BgWin1LoField.Value = snap.bgWin1Lo;
            ui.BgWin1HiField.Value = snap.bgWin1Hi;
            ui.BgDualCheckbox.Value = snap.bgDual;
            ui.BgWin2LoField.Value = snap.bgWin2Lo;
            ui.BgWin2HiField.Value = snap.bgWin2Hi;
            ui.BgIterCheckbox.Value = snap.bgIter;
        end
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

        % Apply same preprocessing as single spectrum display (cached)
        [qe, ~, ~] = local_get_cached_preprocess(qe, local_preprocess_opts());

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
        % Delegates computation to qe_auto_fit.m
        if isempty(state.physicalQE) && isempty(state.comparisonQE)
            ui.InfoLabel.Text = "Load data first";
            ui.InfoLabel.Visible = "on";
            return
        end

        qe = local_get_active_qe();
        if isempty(qe), return; end

        % Preprocess
        qe = qe_preprocess(qe, local_preprocess_opts());

        % Raw (no area normalization) for physical intensity
        qe_raw = local_get_active_qe();
        raw_opts = local_preprocess_opts();
        raw_opts.do_normalize = false;
        qe_raw = qe_preprocess(qe_raw, raw_opts);

        [mask, energy_axis] = local_energy_mask(qe);

        % Parse manual peak guesses
        guess_str = strtrim(ui.GuessField.Value);
        if ~isempty(guess_str)
            guesses = str2double(strsplit(guess_str, {',', ' ', ';'}));
            guesses = guesses(~isnan(guesses) & guesses > 0);
        else
            guesses = [];
        end

        % Build auto-fit options from UI
        auto_opts = struct();
        auto_opts.E_min = max(ui.EnergyMinField.Value, 50);
        auto_opts.E_max = ui.EnergyMaxField.Value;
        auto_opts.q_start = ui.QStartField.Value;
        auto_opts.q_end = ui.QEndField.Value;
        if auto_opts.q_start > auto_opts.q_end
            [auto_opts.q_start, auto_opts.q_end] = deal(auto_opts.q_end, auto_opts.q_start);
        end
        auto_opts.prominence = ui.PromField.Value;
        auto_opts.smooth_width = ui.SmoothField.Value;
        auto_opts.max_peaks = ui.MaxPeaksField.Value;
        auto_opts.peak_model = ui.PeakModelDropdown.Value;
        auto_opts.guesses = guesses;
        auto_opts.seed_idx = state.selectedQIndex;
        auto_opts.max_shift = ui.MaxShiftField.Value;
        auto_opts.energy_mask = mask;
        auto_opts.energy_axis = energy_axis;
        auto_opts.progress_fn = @local_auto_fit_progress;

        ui.InfoLabel.Visible = "on";
        drawnow;

        try
            results = qe_auto_fit(qe, qe_raw, auto_opts);
        catch ME
            ui.InfoLabel.Text = ME.message;
            ui.InfoLabel.Visible = "on";
            return
        end

        % Store results
        peaks = results.all_peaks;
        branches = results.branches;
        n_branches = numel(branches);

        autoResults = struct();
        autoResults.all_peaks = peaks;
        autoResults.branches = branches;
        autoResults.fit_details = results.fit_details;
        autoResults.n_success = results.n_success;
        state.autoFitResults = autoResults;
        state.manual_points = branches{1}(:, 1:min(2, size(branches{1},2)));
        state.manual_branches = branches;

        % ═══════════ Plot ═══════════
        ax = ui.DispersionAxes;
        local_clear_axes(ax);
        hold(ax, 'on');

        disp_model_name = ui.DispModelDropdown.Value;
        for b = 1:n_branches
            br = branches{b};
            col = qe_plot_helpers.branch_color(b);
            branch_label = sprintf('Branch %d (%.0f-%.0f meV, %d pts)', ...
                b, min(br(:,2)), max(br(:,2)), size(br,1));

            qe_plot_helpers.plot_branch_scatter(ax, br, col, branch_label);

            if size(br, 1) >= 5
                try
                    disp_result = fit_dispersion_generic(br(:,1), br(:,2), ...
                        'model', disp_model_name);
                    qe_plot_helpers.plot_fit_curve(ax, disp_result, col, b);
                    state.fitResults{b} = disp_result;
                catch
                end
            end
        end

        hold(ax, 'off');
        legend(ax, 'Location', 'best', 'FontSize', 7);
        xlabel(ax, 'q (1/Å)');
        ylabel(ax, 'Energy (meV)');
        title(ax, sprintf('Auto-Fit: %d branches, %d peaks from %d ch', ...
            n_branches, size(peaks,1), results.n_success));
        grid(ax, 'on');

        ui.PtsLabel.Text = sprintf('%d pts', size(state.manual_points, 1));
        ui.InfoLabel.Text = sprintf('Auto-fit: %d peaks, %d branches', ...
            size(peaks,1), n_branches);
        ui.InfoLabel.Visible = "on";
        local_log_operation(sprintf('Auto fit: %d branches / %d peaks', ...
            n_branches, size(peaks,1)));

        function local_auto_fit_progress(~, message)
            ui.InfoLabel.Text = char(message);
            drawnow;
        end
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

        branch_colors = qe_plot_helpers.branch_colors();
        state.fitResults = {};
        n_branches = numel(state.manual_branches);

        for b = 1:n_branches
            br = state.manual_branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);

            % Plot points using shared helper
            qe_plot_helpers.plot_branch_scatter(ax, br, col, sprintf('Branch %d', b));

            % Fit dispersion
            if size(br, 1) >= 5
                try
                    disp_result = fit_dispersion_generic(br(:,1), br(:,2), ...
                        'model', disp_model_name);
                    qe_plot_helpers.plot_fit_curve(ax, disp_result, col, b);
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

        branch_colors = qe_plot_helpers.branch_colors();

        for b = 1:numel(state.manual_branches)
            br = state.manual_branches{b};
            if isempty(br), continue; end
            col = branch_colors(mod(b-1, size(branch_colors,1))+1, :);

            % Plot with error bars using shared helper
            qe_plot_helpers.plot_branch_scatter(ax_exp, br, col, sprintf('Branch %d', b));

            % Overlay fit curve
            if numel(state.fitResults) >= b && ~isempty(state.fitResults{b})
                qe_plot_helpers.plot_fit_curve(ax_exp, state.fitResults{b}, col, b);
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
        % Delegate to standalone qe_gamma_dashboard module
        if ~isfield(state, 'autoFitResults') || isempty(state.autoFitResults)
            ui.InfoLabel.Text = "Run Auto Fit first";
            ui.InfoLabel.Visible = "on";
            return
        end

        qe_gamma_dashboard(state.autoFitResults.all_peaks, ...
            state.autoFitResults.branches);

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
        qe = qe_preprocess(qe, local_preprocess_opts());
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
        pp_opts = local_preprocess_opts();

        % Use cached preprocessing result
        [qe, qe_pre, bg_diag_all] = local_get_cached_preprocess(qe, pp_opts);

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

        % --- BG overlay: use cached pre-BG spectrum and diagnostics ---
        bg_diag_qi = [];
        pre_bg_spectrum = [];
        if pp_opts.do_bg_sub && strcmpi(mode_name, "display")
            try
                if ~isempty(qe_pre)
                    pre_bg_spectrum = double(qe_pre.intensity(mask, q_index));
                end
                if ~isempty(bg_diag_all) && q_index <= numel(bg_diag_all)
                    bg_diag_qi = bg_diag_all(q_index);
                end
            catch
            end
        end

        hold(ax, "on");
        switch mode_name
            case "display"
                % --- Three-line BG overlay ---
                if ~isempty(bg_diag_qi) && ~isempty(pre_bg_spectrum)
                    % 1. Original pre-BG spectrum (blue)
                    pre_bg_disp = pre_bg_spectrum;
                    if strcmpi(y_scale, "log")
                        pre_bg_disp = max(pre_bg_disp, eps);
                    end
                    plot_handles(end + 1) = plot(ax, energy_axis, pre_bg_disp, "-", ... %#ok<AGROW>
                        "Color", [0.1 0.45 0.7], ...
                        "LineWidth", 1.3, ...
                        "DisplayName", "original");

                    % 2. Background fit curve (red dashed)
                    bg_curve = bg_diag_qi.bg_curve;
                    bg_in_mask = bg_curve(mask);
                    if strcmpi(y_scale, "log")
                        bg_in_mask = max(bg_in_mask, eps);
                    end
                    plot_handles(end + 1) = plot(ax, energy_axis, bg_in_mask, "--", ... %#ok<AGROW>
                        "Color", [0.85 0.24 0.24], ...
                        "LineWidth", 1.4, ...
                        "DisplayName", "BG fit");

                    % 3. Subtracted signal (green)
                    plot_handles(end + 1) = plot(ax, energy_axis, display_smooth, "-", ... %#ok<AGROW>
                        "Color", [0.22 0.65 0.22], ...
                        "LineWidth", 1.5, ...
                        "DisplayName", "signal");

                    % 4. Shade fit windows
                    yl = [min([pre_bg_disp(:); display_smooth(:)]), ...
                          max([pre_bg_disp(:); display_smooth(:)])];
                    if yl(2) <= yl(1), yl(2) = yl(1) + 1; end
                    win1 = pp_opts.bg_win_lo;
                    fill(ax, [win1(1) win1(2) win1(2) win1(1)], ...
                        [yl(1) yl(1) yl(2) yl(2)], ...
                        [0.996 0.914 0.914], ...
                        "EdgeColor", "none", "FaceAlpha", 0.5, ...
                        "HandleVisibility", "off");
                    if ~isempty(pp_opts.bg_win_hi)
                        win2 = pp_opts.bg_win_hi;
                        fill(ax, [win2(1) win2(2) win2(2) win2(1)], ...
                            [yl(1) yl(1) yl(2) yl(2)], ...
                            [0.996 0.914 0.914], ...
                            "EdgeColor", "none", "FaceAlpha", 0.5, ...
                            "HandleVisibility", "off");
                    end

                    % Move shading to back
                    ch = allchild(ax);
                    fills = findobj(ch, "Type", "patch");
                    for fi = 1:numel(fills)
                        uistack(fills(fi), "bottom");
                    end

                    % 5. Display diagnostics in FitInfoLabel
                    diag_parts = {};
                    if isfinite(bg_diag_qi.rsquare)
                        diag_parts{end+1} = sprintf('R²=%.4f', bg_diag_qi.rsquare);
                    end
                    if isfinite(bg_diag_qi.rmse)
                        diag_parts{end+1} = sprintf('RMSE=%.2g', bg_diag_qi.rmse);
                    end
                    if isfinite(bg_diag_qi.h_param)
                        diag_parts{end+1} = sprintf('h=%.1f', bg_diag_qi.h_param);
                    end
                    if isfinite(bg_diag_qi.snr)
                        diag_parts{end+1} = sprintf('SNR=%.1f', bg_diag_qi.snr);
                    end
                    if bg_diag_qi.iterations > 1
                        diag_parts{end+1} = sprintf('iter=%d', bg_diag_qi.iterations);
                    end
                    if ~isempty(diag_parts)
                        ui.FitInfoLabel.Text = sprintf('BG [%s]: %s', ...
                            pp_opts.bg_method, strjoin(diag_parts, ' | '));
                    end
                else
                    % No BG overlay — standard display
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
                end

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
        all_vals = [display_raw(:); display_smooth(:); display_values(:)];
        if ~isempty(pre_bg_spectrum)
            all_vals = [all_vals; pre_bg_spectrum(:)];
        end
        local_apply_y_limits(ax, all_vals);
        ax.ButtonDownFcn = @local_on_single_axes_clicked;
        for handle_idx = 1:numel(plot_handles)
            if isgraphics(plot_handles(handle_idx))
                plot_handles(handle_idx).PickableParts = "all";
                plot_handles(handle_idx).ButtonDownFcn = @local_on_single_axes_clicked;
            end
        end
        if ~strcmpi(mode_name, "display") || ~strcmpi(ui.SmoothingModeDropdown.Value, "off") ...
                || (~isempty(bg_diag_qi) && pp_opts.do_bg_sub)
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
            for bi = 1:numel(state.manual_branches)
                bpts = state.manual_branches{bi};
                if isempty(bpts), continue; end
                col = qe_plot_helpers.branch_color(bi);
                qe_plot_helpers.plot_branch_scatter(ax, bpts, col, sprintf('Branch %d', bi));
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
            for fi = 1:numel(state.fitResults)
                fr = state.fitResults{fi};
                if isempty(fr), continue; end
                col = qe_plot_helpers.branch_color(fi);
                qe_plot_helpers.plot_fit_curve(ax, fr, col, fi);
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
                case {"eq3d", "mat4d", "npy", "raw_import"}
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
        % Delegate to qe_plot_helpers (shared with other modules)
        if nargin < 2, mode_name = local_trace_mode(); end
        if nargin < 3, display_style = "physical"; end
        [display_map, display_clim] = qe_plot_helpers.prepare_map_display( ...
            map_values, mode_name, display_style);
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


    function opts = local_preprocess_opts()
        % Build preprocessing options struct from current UI state
        opts = struct();
        opts.do_normalize = ui.AreaNormCheckbox.Value;
        opts.norm_method = char(ui.NormMethodDropdown.Value);
        % Auto-sync: normalization window follows the energy display range
        opts.norm_min = ui.EnergyMinField.Value;
        opts.norm_max = ui.EnergyMaxField.Value;
        ui.AreaNormMinField.Value = opts.norm_min;
        ui.AreaNormMaxField.Value = opts.norm_max;
        opts.do_denoise = ui.DenoiseCheckbox.Value;
        opts.denoise_method = char(ui.DenoiseMethodDropdown.Value);
        opts.denoise_sigma = ui.DenoiseSigmaField.Value;
        opts.sg_order = ui.SGOrderField.Value;
        opts.sg_framelen = ui.SGFrameLenField.Value;
        opts.do_bg_sub = ui.BgSubCheckbox.Value;
        opts.bg_method = char(ui.BgMethodDropdown.Value);
        opts.do_deconv = ui.DeconvCheckbox.Value;
        opts.deconv_iter = ui.DeconvIterField.Value;

        % --- Enhanced background subtraction options ---
        opts.bg_win_lo = [ui.BgWin1LoField.Value, ui.BgWin1HiField.Value];
        if ui.BgDualCheckbox.Value
            opts.bg_win_hi = [ui.BgWin2LoField.Value, ui.BgWin2HiField.Value];
        else
            opts.bg_win_hi = [];
        end
        opts.bg_iterative = ui.BgIterCheckbox.Value;
    end


    function [qe_pp, qe_pre, bg_diag] = local_get_cached_preprocess(qe_raw, pp_opts)
        % Cached preprocessing to avoid redundant SVD/FFT on every click.
        %
        % Returns:
        %   qe_pp    — fully preprocessed data (despike+denoise+BG+deconv)
        %   qe_pre   — preprocessed WITHOUT BG and deconv (for overlay)
        %   bg_diag  — BG diagnostics struct array
        %
        % Cache is invalidated when opts change or the active dataset changes.

        % Build a hash from the opts and the dataset identity
        hash_str = local_opts_hash(pp_opts, qe_raw);

        if strcmp(state.ppCache.hash, hash_str) && ~isempty(state.ppCache.qe)
            % Cache hit — return stored results
            qe_pp   = state.ppCache.qe;
            qe_pre  = state.ppCache.qe_pre;
            bg_diag = state.ppCache.bg_diag;
            return
        end

        % Cache miss — compute everything once
        % 1. Full preprocessing
        [qe_pp, bg_diag] = qe_preprocess(qe_raw, pp_opts);

        % 2. Pre-BG version (for three-line overlay)
        if pp_opts.do_bg_sub
            pre_opts = pp_opts;
            pre_opts.do_bg_sub = false;
            pre_opts.do_deconv = false;
            qe_pre = qe_preprocess(qe_raw, pre_opts);

            % If BG diag wasn't returned from the full pass, get it
            if isempty(bg_diag)
                bg_opts = pre_opts;
                bg_opts.do_bg_sub = true;
                bg_opts.bg_method = pp_opts.bg_method;
                bg_opts.bg_win_lo = pp_opts.bg_win_lo;
                bg_opts.bg_win_hi = pp_opts.bg_win_hi;
                bg_opts.bg_iterative = pp_opts.bg_iterative;
                [~, bg_diag] = qe_preprocess(qe_pre, bg_opts);
            end
        else
            qe_pre  = [];
            bg_diag = [];
        end

        % Store in cache
        state.ppCache.hash    = hash_str;
        state.ppCache.qe      = qe_pp;
        state.ppCache.qe_pre  = qe_pre;
        state.ppCache.bg_diag = bg_diag;
    end


    function h = local_opts_hash(opts, qe)
        % Fast hash of preprocessing options + dataset identity.
        % Uses a simple string concatenation of all relevant fields.
        parts = {};
        parts{end+1} = sprintf('%d', opts.do_normalize);
        parts{end+1} = opts.norm_method;
        parts{end+1} = sprintf('%.1f', opts.norm_min);
        parts{end+1} = sprintf('%.1f', opts.norm_max);
        parts{end+1} = sprintf('%d', opts.do_denoise);
        parts{end+1} = opts.denoise_method;
        parts{end+1} = sprintf('%.1f', opts.denoise_sigma);
        parts{end+1} = sprintf('%d', opts.do_bg_sub);
        parts{end+1} = opts.bg_method;
        parts{end+1} = sprintf('%.1f_%.1f', opts.bg_win_lo(1), opts.bg_win_lo(2));
        if ~isempty(opts.bg_win_hi)
            parts{end+1} = sprintf('%.1f_%.1f', opts.bg_win_hi(1), opts.bg_win_hi(2));
        end
        parts{end+1} = sprintf('%d', opts.bg_iterative);
        parts{end+1} = sprintf('%d', opts.do_deconv);
        parts{end+1} = sprintf('%d', opts.deconv_iter);
        % Dataset identity: use size + first/last values as fingerprint
        parts{end+1} = sprintf('%dx%d', size(qe.intensity, 1), size(qe.intensity, 2));
        parts{end+1} = sprintf('%.6g', qe.intensity(1,1));
        h = strjoin(parts, '|');
    end


end
