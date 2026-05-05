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
        state_out.autoCorrectionArmed = false;
        state_out.manual_points = [];  % Nx2 [q, energy], kept in sync with manual_branches when available
        state_out.manual_branches = {};  % cell array of branch matrices; first two columns are [q, energy]
        state_out.branchCorrectionLog = [];  % struct array of auto-fit manual corrections
        state_out.fitResults = {};       % cell array of fit_quasi2d_plasmon results
        state_out.autoFitResults = [];   % struct array from auto Drude-Lorentz fit
        state_out.pendingFit = [];       % pending single-spectrum Lorentz fit result
        state_out.qeImage = gobjects(1);
        state_out.selectionMarker = gobjects(1);
        % Preprocessing cache (avoid redundant SVD/FFT on every click)
        state_out.ppCache = struct('hash', '', 'qe', [], 'qe_pre', [], 'bg_diag', []);
        state_out.singlePreviewCache = struct('hash', '', 'qe_pre', []);
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
        cb.on_fit_param_change = @local_on_fit_param_change;
        cb.on_correction_target_change = @local_on_correction_target_change;
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
        cb.on_correct_auto_peak = @local_on_correct_auto_peak;
        cb.on_undo_correction = @local_on_undo_correction;
        cb.on_fit_dispersion = @local_on_fit_dispersion;
        cb.on_export_dispersion = @local_on_export_dispersion;
        cb.on_show_loss_map = @local_on_show_loss_map;
        cb.on_export_data = @local_on_export_data;
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
                state.autoFitResults = [];   % clear stale fits (Issue #1)
                state.ppCache = struct('hash', '', 'qe', [], 'qe_pre', [], 'bg_diag', []);
                state.singlePreviewCache = struct('hash', '', 'qe_pre', []);
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
            state.autoFitResults = [];   % clear stale fit results (Issue #1)
            state.ppCache = struct('hash', '', 'qe', [], 'qe_pre', [], 'bg_diag', []);
            state.singlePreviewCache = struct('hash', '', 'qe_pre', []);

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


    function local_ensure_latest_history_snapshot(label)
        % Save History is the thesis-profile anchor.  Do not rely only on
        % per-widget callbacks: append the current live UI state if it differs
        % from the last logged snapshot.
        current = local_capture_snapshot();
        if isempty(state.opHistory) || ~isequaln(state.opHistory{end}.snapshot, current)
            local_log_operation(label);
        end
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
        snap.prominence = ui.PromField.Value;
        snap.autoFitSmoothWidth = ui.SmoothField.Value;
        snap.maxPeaks = ui.MaxPeaksField.Value;
        snap.peakModel = char(ui.PeakModelDropdown.Value);
        snap.maxShift = ui.MaxShiftField.Value;
        snap.guessText = char(string(ui.GuessField.Value));
        snap.branch1Min = ui.Branch1MinField.Value;
        snap.branch1Max = ui.Branch1MaxField.Value;
        snap.branch2Min = ui.Branch2MinField.Value;
        snap.branch2Max = ui.Branch2MaxField.Value;
        snap.branch3Min = ui.Branch3MinField.Value;
        snap.branch3Max = ui.Branch3MaxField.Value;
        snap.dispModel = char(ui.DispModelDropdown.Value);
        snap.fitBranchTarget = char(string(ui.FitBranchDropdown.Value));
        snap.correctionBranchTarget = char(string(ui.CorrectionBranchDropdown.Value));
        snap.exportRatio = char(ui.ExportRatioDropdown.Value);
        snap.selectedQIndex = state.selectedQIndex;
        snap.selectedQ_Ainv = NaN;
        qe_snap = local_get_reference_qe();
        if ~isempty(qe_snap) && isfield(qe_snap, 'q_Ainv') && ~isempty(qe_snap.q_Ainv)
            q_idx = min(max(1, state.selectedQIndex), numel(qe_snap.q_Ainv));
            snap.selectedQ_Ainv = qe_snap.q_Ainv(q_idx);
        end

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
        if isfield(snap, 'prominence')
            ui.PromField.Value = snap.prominence;
        end
        if isfield(snap, 'autoFitSmoothWidth')
            ui.SmoothField.Value = snap.autoFitSmoothWidth;
        end
        if isfield(snap, 'maxPeaks')
            ui.MaxPeaksField.Value = snap.maxPeaks;
        end
        if isfield(snap, 'peakModel') && any(strcmp(snap.peakModel, ui.PeakModelDropdown.ItemsData))
            ui.PeakModelDropdown.Value = snap.peakModel;
        end
        if isfield(snap, 'maxShift')
            ui.MaxShiftField.Value = snap.maxShift;
        end
        if isfield(snap, 'guessText')
            ui.GuessField.Value = string(snap.guessText);
        end
        if isfield(snap, 'branch1Min')
            ui.Branch1MinField.Value = snap.branch1Min;
            ui.Branch1MaxField.Value = snap.branch1Max;
            ui.Branch2MinField.Value = snap.branch2Min;
            ui.Branch2MaxField.Value = snap.branch2Max;
            ui.Branch3MinField.Value = snap.branch3Min;
            ui.Branch3MaxField.Value = snap.branch3Max;
        end
        if isfield(snap, 'dispModel') && any(strcmp(snap.dispModel, ui.DispModelDropdown.ItemsData))
            ui.DispModelDropdown.Value = snap.dispModel;
        end
        if isfield(snap, 'fitBranchTarget')
            local_sync_correction_controls();
            target = char(string(snap.fitBranchTarget));
            if any(strcmp(target, ui.FitBranchDropdown.Items))
                ui.FitBranchDropdown.Value = target;
            end
        end
        if isfield(snap, 'exportRatio') && any(strcmp(snap.exportRatio, ui.ExportRatioDropdown.Items))
            ui.ExportRatioDropdown.Value = snap.exportRatio;
        end
        if isfield(snap, 'selectedQ_Ainv') && isfinite(snap.selectedQ_Ainv)
            qe_snap = local_get_reference_qe();
            if ~isempty(qe_snap) && isfield(qe_snap, 'q_Ainv') && ~isempty(qe_snap.q_Ainv)
                [~, state.selectedQIndex] = min(abs(qe_snap.q_Ainv - snap.selectedQ_Ainv));
            end
        elseif isfield(snap, 'selectedQIndex')
            state.selectedQIndex = snap.selectedQIndex;
        end

        ui.ViewModeDropdown.Value = local_normalize_view_mode_value(snap.viewMode);
        if isfield(snap, 'correctionBranchTarget')
            local_sync_correction_controls();
            target = char(string(snap.correctionBranchTarget));
            if any(strcmp(target, ui.CorrectionBranchDropdown.Items))
                ui.CorrectionBranchDropdown.Value = target;
            end
        end
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
        % Save opHistory to a .mat file in the data folder.  Always append a
        % final live snapshot first so manually tuned controls are captured
        % even when an individual widget did not emit a history callback.
        local_ensure_latest_history_snapshot("Saved current GUI state");
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
            'prominence', 'Prominence', ...
            'autoFitSmoothWidth', 'Auto smooth W', ...
            'maxPeaks', 'Max peaks', ...
            'peakModel', 'Peak model', ...
            'maxShift', 'Max shift', ...
            'guessText', 'Guesses', ...
            'dispModel', 'Disp model', ...
            'correctionBranchTarget', 'Corr target', ...
            'exportRatio', 'Export ratio', ...
            'selectedQIndex', 'q index', ...
            'selectedQ_Ainv', 'selected q', ...
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
            if isequaln(old_val, new_val)
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

        local_log_operation(local_describe_changes());
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


    function local_on_fit_param_change(~, ~)
        % Auto-fit/dispersion/export parameters may not redraw the view, but
        % they are part of a reusable GUI profile and must enter history.
        local_log_operation(local_describe_changes());
    end


    function local_on_correction_target_change(~, ~)
        % Correction target affects where subsequent curated peak edits go.
        local_log_operation(local_describe_changes());
    end



    function local_on_pick_peaks(~, ~)
        if isempty(state.physicalQE)
            return
        end
        state.manualPickArmed = true;
        state.guessPickArmed = false;
        state.autoCorrectionArmed = false;
        ui.PickPtsButton.Text = "Picking...";
        ui.PickPtsButton.BackgroundColor = [0.85 1 0.85];
        ui.PickGuessesButton.Text = "Pick Guesses";
        ui.PickGuessesButton.BackgroundColor = [0.96 0.96 0.96];
        ui.CorrectAutoButton.Text = "Correct Peak";
        ui.CorrectAutoButton.BackgroundColor = [0.96 0.96 0.96];
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
        state.branchCorrectionLog = [];
        ui.PtsLabel.Text = "0 pts";
        local_sync_correction_controls();
        local_update_all_views();
    end


    function local_on_save_pts(~, ~)
        if isempty(state.manual_points) && isempty(state.manual_branches)
            return
        end
        [file, path] = uiputfile('*.mat', 'Save dispersion points', 'dispersion_points.mat');
        if isequal(file, 0)
            return
        end
        if ~isempty(state.manual_branches)
            pts = qe_flatten_branch_points(state.manual_branches);
        else
            pts = state.manual_points;
        end
        manual_dispersion = struct();
        manual_dispersion.q_Ainv = pts(:, 1);
        manual_dispersion.abs_q_Ainv = abs(pts(:, 1));  % backwards-compatible field name
        manual_dispersion.energy_meV = pts(:, 2);
        manual_dispersion.branches = state.manual_branches;
        manual_dispersion.corrections = state.branchCorrectionLog;
        manual_dispersion.snapshot = local_capture_snapshot();
        manual_dispersion.source = 'gui_auto_fit_with_manual_corrections';
        manual_dispersion.timestamp = datestr(now);
        save(fullfile(path, file), 'manual_dispersion');
        ui.PtsLabel.Text = sprintf("%d pts (saved)", size(pts, 1));
    end


    function local_on_load_pts(~, ~)
        [file, path] = uigetfile('*.mat', 'Load dispersion points');
        if isequal(file, 0)
            return
        end
        loaded = load(fullfile(path, file));
        if isfield(loaded, 'manual_dispersion')
            md = loaded.manual_dispersion;
            if isfield(md, 'branches') && ~isempty(md.branches)
                state.manual_branches = md.branches;
                state.manual_points = qe_flatten_branch_points(state.manual_branches);
            elseif isfield(md, 'q_Ainv')
                state.manual_points = [md.q_Ainv(:), md.energy_meV(:)];
                state.manual_branches = {};
            elseif isfield(md, 'abs_q_Ainv')
                state.manual_points = [md.abs_q_Ainv(:), md.energy_meV(:)];
                state.manual_branches = {};
            else
                warning('File does not contain q/energy dispersion data.');
                return
            end
            if isfield(md, 'corrections')
                state.branchCorrectionLog = md.corrections;
            else
                state.branchCorrectionLog = [];
            end
        else
            warning('File does not contain manual_dispersion data.');
            return
        end
        ui.PtsLabel.Text = sprintf("%d pts (loaded)", size(state.manual_points, 1));
        local_sync_correction_controls();
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
        q_value = qe.q_Ainv(q_index);
        abs_q = abs(q_value);
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

        % Filter out peaks below E_min for display, and keep only reliable
        % peaks for acceptance into curated branches.
        peak_energy = local_fit_peak_energy(result);
        display_keep = peak_energy >= E_min;
        if ~any(display_keep)
            ui.FitInfoLabel.Text = "No valid peaks found above E_min";
            return
        end
        keep = display_keep;
        if isfield(result, 'peak_valid')
            keep = keep & result.peak_valid(:);
        end

        % Store pending fit
        state.pendingFit = struct();
        state.pendingFit.abs_q = abs_q;
        state.pendingFit.q_Ainv = q_value;
        state.pendingFit.q_index = q_index;
        state.pendingFit.result = result;
        state.pendingFit.keep = keep;
        state.pendingFit.display_keep = display_keep;

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
            if ~display_keep(p); continue; end
            col = peak_colors(mod(p-1, size(peak_colors,1))+1, :);
            is_valid_peak = true;
            if isfield(result, 'peak_valid') && numel(result.peak_valid) >= p
                is_valid_peak = logical(result.peak_valid(p));
            end
            if ~is_valid_peak
                col = [0.95 0.55 0.05];
            end

            % Plot peak component + background so it's visible
            peak_on_bg = result.peak_curves{p} + bg_vals;
            plot(ax, result.energy_fit, peak_on_bg, '--', ...
                'Color', col, 'LineWidth', 1.5, 'Tag', 'lorentz_fit');

            % Peak marker at peak position on total fit
            [~, pk_idx] = min(abs(result.energy_fit - peak_energy(p)));
            marker_y = result.curve_fit(pk_idx);
            plot(ax, peak_energy(p), marker_y, 'v', ...
                'MarkerSize', 10, 'MarkerFaceColor', col, ...
                'MarkerEdgeColor', 'k', 'Tag', 'lorentz_fit');

            marker_label = sprintf('E0 %.0f', result.omega_p(p));
            if isfield(result, 'apex_energy_meV') && numel(result.apex_energy_meV) >= p
                marker_label = sprintf('E0 %.0f | apex %.0f', ...
                    result.omega_p(p), result.apex_energy_meV(p));
            end
            if ~is_valid_peak
                marker_label = sprintf('%s | low confidence', marker_label);
            end
            text(ax, peak_energy(p), marker_y * 1.08, ...
                marker_label, ...
                'FontSize', 9, 'FontWeight', 'bold', 'Color', col, ...
                'HorizontalAlignment', 'center', 'Tag', 'lorentz_fit');

            if isfield(result, 'gamma_ratio') && numel(result.gamma_ratio) >= p ...
                    && isfield(result, 'apex_energy_meV') && numel(result.apex_energy_meV) >= p
                quality_text = 'ok';
                if ~is_valid_peak
                    quality_text = 'low-confidence';
                end
                if isfield(result, 'fano_q') && numel(result.fano_q) >= p
                    info_parts{end+1} = sprintf('E0=%.0f apex=%.0f qF=%.2f Gamma/E0=%.2f %s', ...
                        result.omega_p(p), result.apex_energy_meV(p), ...
                        result.fano_q(p), result.gamma_ratio(p), quality_text); %#ok<AGROW>
                else
                    info_parts{end+1} = sprintf('E0=%.0f apex=%.0f Gamma/E0=%.2f %s', ...
                        result.omega_p(p), result.apex_energy_meV(p), ...
                        result.gamma_ratio(p), quality_text); %#ok<AGROW>
                end
            else
                info_parts{end+1} = sprintf('E0=%.0f Gamma=%.0f', ...
                    result.omega_p(p), result.gamma(p)); %#ok<AGROW>
            end
            continue
        end
        hold(ax, 'off');

        % Update fit info
        ui.FitInfoLabel.Text = sprintf('q=%.4f | %s | R²=%.3f', ...
            abs_q, strjoin(info_parts, ', '), result.R_squared);

        % Enable Accept button
        ui.AcceptFitButton.Enable = local_on_off(any(keep));
    end


    function peak_energy = local_fit_peak_energy(result)
        peak_energy = result.omega_p(:);
        if isfield(result, 'apex_energy_meV')
            apex = result.apex_energy_meV(:);
            use_apex = isfinite(apex);
            peak_energy(use_apex) = apex(use_apex);
        end
    end


    function local_on_accept_fit(~, ~)
        % Accept the pending fit and add peaks to dispersion
        if isempty(state.pendingFit)
            return
        end

        pf = state.pendingFit;
        result = pf.result;
        peak_energy = local_fit_peak_energy(result);
        keep = pf.keep;
        abs_q = pf.abs_q;
        if isfield(pf, 'q_Ainv')
            q_value = pf.q_Ainv;
        else
            q_value = abs_q;
        end

        if ~isempty(state.manual_branches)
            % Hybrid mode: use the accepted single-spectrum fit as a manual
            % correction of the existing auto-fit branch point at this q.
            n_corrected = 0;
            branch_spec = local_selected_correction_branch_spec();
            for p = 1:result.n_peaks
                if ~keep(p); continue; end
                [state.manual_branches, correction] = qe_apply_branch_correction( ...
                    state.manual_branches, q_value, peak_energy(p), ...
                    'Branch', branch_spec, ...
                    'QTolerance', local_branch_correction_q_tolerance());
                local_append_branch_correction(correction);
                n_corrected = n_corrected + 1;
            end
            local_sync_manual_points_from_branches();
            local_update_all_views();
            ui.AcceptFitButton.Enable = "off";
            ui.FitInfoLabel.Text = sprintf('%s  ✓ Corrected auto (%d)', ...
                ui.FitInfoLabel.Text, n_corrected);
            local_log_operation(sprintf('Correct auto from fit: q=%.4f, %d peaks', ...
                q_value, n_corrected));
            state.pendingFit = [];
            return
        end

        % Add accepted peaks to unbranched manual_points when no auto branches exist.
        for p = 1:result.n_peaks
            if ~keep(p); continue; end
            new_pt = [q_value, peak_energy(p)];
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
            scatter(ax, q_value, peak_energy(p), 40, ...
                [0.1 0.5 0.9], 'filled', ...
                'HandleVisibility', 'off');
        end
        hold(ax, 'off');
        xlabel(ax, 'q (Å^{-1})');
        ylabel(ax, 'Energy (meV)');
        title(ax, sprintf('Dispersion: %d pts', size(state.manual_points, 1)));
        grid(ax, 'on');

        % Update UI
        ui.PtsLabel.Text = sprintf('%d pts', size(state.manual_points, 1));
        ui.AcceptFitButton.Enable = "off";
        ui.FitInfoLabel.Text = sprintf('%s  ✓ Accepted', ui.FitInfoLabel.Text);

        local_log_operation(sprintf('Accept fit: q=%.4f, %d peaks', ...
            q_value, sum(keep)));

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

        try
            branch_specs = local_branch_specs_from_ui();
        catch ME
            ui.InfoLabel.Text = ME.message;
            ui.InfoLabel.Visible = "on";
            return
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
        % qe above is produced with the same background-subtraction checkbox,
        % so Auto Fit must not fit another power-law background when it is on.
        auto_opts.pre_subtracted = ui.BgSubCheckbox.Value;
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
        branch_filter_opts = struct();
        branch_filter_opts.q_skip_Ainv = 0.005;
        branch_filter_opts.min_R2 = 0.3;
        branch_filter_opts.max_gamma_ratio = 2.0;
        branch_filter_opts.score_column = 12;
        branch_assignment = qe_assign_peak_branches_by_windows(peaks, branch_specs, branch_filter_opts);
        branches = branch_assignment.branches;
        n_branches = numel(branches);

        autoResults = struct();
        autoResults.all_peaks = peaks;
        autoResults.raw_branches = results.branches;
        autoResults.branches = branches;
        autoResults.branch_assignment = branch_assignment;
        autoResults.branch_specs = branch_assignment.specs;
        autoResults.fit_details = results.fit_details;
        autoResults.n_success = results.n_success;
        state.autoFitResults = autoResults;
        state.autoFitResults.ppHash = local_opts_hash(local_preprocess_opts(), local_get_active_qe());
        state.manual_branches = branches;
        state.manual_points = qe_flatten_branch_points(branches);
        state.branchCorrectionLog = [];
        state.fitResults = {};
        local_sync_correction_controls();

        % ═══════════ Plot ═══════════
        ax = ui.DispersionAxes;
        local_clear_axes(ax);
        hold(ax, 'on');

        plotted = false;
        for b = 1:n_branches
            br = branches{b};
            col = qe_plot_helpers.branch_color(b);
            spec = branch_assignment.specs(b);
            win = spec.energy_window_meV;
            branch_label = sprintf('%s [%.0f-%.0f meV] (%d pts)', ...
                char(string(spec.name)), win(1), win(2), size(br,1));
            if isempty(br)
                continue
            end
            qe_plot_helpers.plot_branch_scatter(ax, br, col, branch_label);
            plotted = true;
        end

        hold(ax, 'off');
        if plotted
            legend(ax, 'Location', 'best', 'FontSize', 7);
        end
        xlabel(ax, 'q (1/A)');
        ylabel(ax, 'Energy (meV)');
        title(ax, sprintf('Auto candidates: %d branches, %d assigned / %d peaks', ...
            n_branches, size(state.manual_points, 1), size(peaks,1)));
        grid(ax, 'on');

        ui.PtsLabel.Text = sprintf('%d pts', size(state.manual_points, 1));
        ui.InfoLabel.Text = sprintf('Auto-fit candidates ready: %d assigned / %d peaks. Correct peaks, then Fit Curated.', ...
            size(state.manual_points, 1), size(peaks,1));
        ui.InfoLabel.Visible = "on";
        ui.DispInfoLabel.Text = "Candidates only. Use Correct Peak or Fit Curated.";
        local_log_operation(sprintf('Auto fit candidates: %d branches / %d assigned / %d peaks', ...
            n_branches, size(state.manual_points,1), size(peaks,1)));
        local_update_heatmap_branch_overlays();

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

        local_sync_manual_points_from_branches();

        % Force update of overlay
        local_update_all_views();
    end

    function local_on_correct_auto_peak(~, ~)
        % Toggle hybrid correction mode: click the corrected peak energy in
        % the single-spectrum panel to replace/add a point in auto branches.
        if isempty(state.manual_branches)
            ui.InfoLabel.Text = "Run Auto Fit first, then use Correct Peak.";
            ui.InfoLabel.Visible = "on";
            return
        end

        state.autoCorrectionArmed = ~state.autoCorrectionArmed;
        if state.autoCorrectionArmed
            state.manualPickArmed = false;
            state.guessPickArmed = false;
            ui.PickPtsButton.Text = "Pick Peaks";
            ui.PickPtsButton.BackgroundColor = [0.96 0.96 0.96];
            ui.PickGuessesButton.Text = "Pick Guesses";
            ui.PickGuessesButton.BackgroundColor = [0.96 0.96 0.96];
            ui.CorrectAutoButton.Text = "Click peak...";
            ui.CorrectAutoButton.BackgroundColor = [1.0 0.92 0.55];
            ui.InfoLabel.Text = sprintf("Correct Peak [%s]: click the corrected peak in the single-spectrum panel.", ...
                char(string(ui.CorrectionBranchDropdown.Value)));
            ui.InfoLabel.Visible = "on";
        else
            ui.CorrectAutoButton.Text = "Correct Peak";
            ui.CorrectAutoButton.BackgroundColor = [0.96 0.96 0.96];
            ui.InfoLabel.Text = "Correct Peak cancelled.";
            ui.InfoLabel.Visible = "on";
        end
    end


    function local_on_undo_correction(~, ~)
        if isempty(state.branchCorrectionLog)
            ui.InfoLabel.Text = "No manual correction to undo.";
            ui.InfoLabel.Visible = "on";
            local_sync_correction_controls();
            return
        end

        correction = state.branchCorrectionLog(end);
        try
            state.manual_branches = qe_revert_branch_correction(state.manual_branches, correction);
            state.branchCorrectionLog(end) = [];
            local_sync_manual_points_from_branches();
            local_update_all_views();
            ui.InfoLabel.Text = sprintf("Undid %s correction on Branch %d at q=%.4f.", ...
                char(correction.action), correction.branch_index, correction.new_q_Ainv);
            ui.InfoLabel.Visible = "on";
            local_log_operation(sprintf('Undo correction: B%d q=%.4f', ...
                correction.branch_index, correction.new_q_Ainv));
        catch ME
            local_show_error(ME, false);
        end
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
        fit_indices = qe_selected_branch_indices(n_branches, ui.FitBranchDropdown.Value);
        if isempty(fit_indices)
            hold(ax, 'off');
            ui.DispInfoLabel.Text = "No selected branch to fit.";
            return
        end

        for b = fit_indices
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
        title(ax, sprintf('Dispersion [%s] | %s', ...
            disp_model_name, char(string(ui.FitBranchDropdown.Value))));
        legend(ax, 'Location', 'best', 'FontSize', 7);

        % Summarize in info label
        n_fitted = sum(~cellfun('isempty', state.fitResults));
        ui.DispInfoLabel.Text = sprintf('Fitted %d selected branch(es) with %s', ...
            n_fitted, disp_model_name);
        local_log_operation(sprintf('Fit dispersion: %s, target=%s, %d branch(es)', ...
            disp_model_name, char(string(ui.FitBranchDropdown.Value)), n_fitted));
        local_update_heatmap_branch_overlays();
        local_update_selected_q_views();
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
            state.autoCorrectionArmed = false;
            ui.PickPtsButton.Text = "Pick Peaks";
            ui.PickPtsButton.BackgroundColor = [0.96 0.96 0.96];
            ui.CorrectAutoButton.Text = "Correct Peak";
            ui.CorrectAutoButton.BackgroundColor = [0.96 0.96 0.96];
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

        % Stale-fit check (same logic as I_kin-corrected map generation — Issue #3 round 3)
        qe = local_get_active_qe();
        if ~isempty(qe) && isfield(state.autoFitResults, 'ppHash')
            current_hash = local_opts_hash(local_preprocess_opts(), qe);
            if ~strcmp(state.autoFitResults.ppHash, current_hash)
                ui.InfoLabel.Text = "⚠ Fit results stale — re-run Auto Fit for accurate analysis.";
                ui.InfoLabel.Visible = "on";
            end
        end

        % Prefer manual_branches if user has reassigned (Issue #2 round 3)
        if ~isempty(state.manual_branches)
            branches = state.manual_branches;
        else
            branches = state.autoFitResults.branches;
        end

        % Reconstruct all_peaks from active branches so title count is accurate
        all_peaks = vertcat(branches{:});
        qe_gamma_dashboard(all_peaks, branches);

        local_log_operation('Show Γ & Amplitude analysis');
    end


    function local_on_show_loss_map(~, ~)
        % Generate an I_kin-corrected q-E map, or a loss-function-like
        % visualization when self-consistent physics extraction is available.
        qe = local_get_active_qe();
        if isempty(qe)
            ui.InfoLabel.Text = "Load data first";
            ui.InfoLabel.Visible = "on";
            return
        end

        % Preprocess the active qe data
        pp_opts = local_preprocess_opts();
        [qe_pp, ~, ~] = local_get_cached_preprocess(qe, pp_opts);

        % Check if we have fitting results for self-consistent mode
        phys = [];
        mode = 'simple';
        if isfield(state, 'autoFitResults') && ~isempty(state.autoFitResults)
            % Check if preprocessing has changed since the fit was run
            current_hash = local_opts_hash(local_preprocess_opts(), qe);
            if isfield(state.autoFitResults, 'ppHash') && ...
               ~strcmp(state.autoFitResults.ppHash, current_hash)
                fprintf('  I_kin Map: WARNING — preprocessing changed since Auto Fit.\n');
                fprintf('  Falling back to simple q⁻³ mode. Re-run Auto Fit for experimental mode.\n');
                ui.InfoLabel.Text = "⚠ Fit results stale — using simple mode. Re-run Auto Fit.";
                ui.InfoLabel.Visible = "on";
            else
                try
                    % Prefer manual_branches if user has reassigned (Issue #2 round 3)
                    if ~isempty(state.manual_branches)
                        branches_for_phys = state.manual_branches;
                    else
                        branches_for_phys = state.autoFitResults.branches;
                    end
                    phys = qe_physics_extract(branches_for_phys);
                    mode = 'experimental';
                    fprintf('  I_kin Map: using self-consistent I_kin (ρ₀=%.0fÅ)\n', phys.rho0);
                catch ME
                    fprintf('  I_kin Map: physics extraction failed (%s), using q⁻³ fallback\n', ME.message);
                end
            end
        else
            fprintf('  I_kin Map: no fit data, using simple q⁻³ mode\n');
            fprintf('  Tip: run Auto Fit first for self-consistent I_kin correction\n');
        end

        % Get display range from UI
        map_opts = struct();
        map_opts.mode = mode;
        map_opts.E_range = [ui.EnergyMinField.Value, ui.EnergyMaxField.Value];
        map_opts.q_range = [ui.QStartField.Value, ui.QEndField.Value];

        try
            qe_loss_map(qe_pp, phys, map_opts);
            local_log_operation(sprintf('I_kin Map (%s mode)', mode));
            ui.InfoLabel.Text = sprintf('I_kin-corrected map generated (%s mode)', mode);
            ui.InfoLabel.Visible = "on";
        catch ME
            local_show_error(ME, false);
        end
    end


    function local_on_export_data(~, ~)
        % Export preprocessed q-E data as a lossless .mat file.
        %
        % IMPORTANT: Always exports from physicalQE (the true q-E map),
        % regardless of the current view mode. This avoids exporting
        % the comparison view as if it were raw physical data.
        if isempty(state.physicalQE)
            ui.InfoLabel.Text = "Load data first";
            ui.InfoLabel.Visible = "on";
            return
        end

        % Always use physicalQE as the source — never the active view
        qe_physical = state.physicalQE;

        % Run full preprocessing (use cache if available)
        pp_opts = local_preprocess_opts();
        [qe_pp, qe_pre, bg_diag] = local_get_cached_preprocess(qe_physical, pp_opts);

        % Ask user for save location
        default_dir = pwd;
        src_name = 'qe_data';
        if ~isempty(state.dataset) && isfield(state.dataset, 'source_path')
            [src_folder, src_name, ~] = fileparts(char(state.dataset.source_path));
            if isfolder(src_folder)
                default_dir = src_folder;
            end
        end

        % Build default filename with preprocessing summary
        pp_summary = '';
        if pp_opts.do_normalize
            pp_summary = [pp_summary '_norm'];
        end
        if pp_opts.do_denoise
            pp_summary = [pp_summary '_' lower(pp_opts.denoise_method)];
        end
        if pp_opts.do_bg_sub
            pp_summary = [pp_summary '_bg' lower(pp_opts.bg_method)];
        end
        if pp_opts.do_deconv
            pp_summary = [pp_summary '_deconv'];
        end
        default_name = sprintf('%s_preprocessed%s.mat', src_name, pp_summary);

        [fname, fpath] = uiputfile({'*.mat', 'MATLAB data (*.mat)'}, ...
            'Export preprocessed data', fullfile(default_dir, default_name));
        if isequal(fname, 0)
            return
        end
        out_file = fullfile(fpath, fname);

        try
            % Pack export struct — all fields needed to reconstruct analysis
            export = struct();

            % Core data (double precision, lossless)
            export.intensity  = double(qe_pp.intensity);
            export.energy_meV = double(qe_pp.energy_meV(:));
            export.q_Ainv     = double(qe_pp.q_Ainv(:));

            % Pre-BG data (for reference / re-analysis)
            if ~isempty(qe_pre)
                export.intensity_pre_bg = double(qe_pre.intensity);
            end

            % Unprocessed physical data (always from physicalQE, never from
            % the active view — prevents Comparison view data from being
            % mislabeled as raw physical data)
            export.intensity_unprocessed = double(qe_physical.intensity);

            % Metadata — data provenance
            export.data_kind = 'physical';  % always physical, never 'normalized'
            export.view_mode_at_export = char(ui.ViewModeDropdown.Value);
            export.preprocessing = pp_opts;
            export.source_path = '';
            if ~isempty(state.dataset) && isfield(state.dataset, 'source_path')
                export.source_path = char(state.dataset.source_path);
            end
            export.dq_Ainv = qe_physical.dq_Ainv;
            export.export_timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
            export.toolbox_version = 'qe_browser v2.0';

            % BG fit diagnostics (R², model params per channel)
            if ~isempty(bg_diag)
                export.bg_diagnostics = bg_diag;
            end

            % Size info
            [nE, nQ] = size(export.intensity);
            export.dimensions = struct('n_energy', nE, 'n_q', nQ, ...
                'E_range_meV', [export.energy_meV(1), export.energy_meV(end)], ...
                'q_range_Ainv', [export.q_Ainv(1), export.q_Ainv(end)]);

            % Save with v7.3 for large data support
            save(out_file, '-struct', 'export', '-v7.3');

            file_info = dir(out_file);
            size_MB = file_info.bytes / 1e6;

            msg = sprintf('Exported %dx%d (%.1f MB): %s', nE, nQ, size_MB, fname);
            fprintf('  %s\n', msg);
            ui.InfoLabel.Text = msg;
            ui.InfoLabel.Visible = "on";
            local_log_operation(sprintf('Export data: %s', fname));
        catch ME
            local_show_error(ME, false);
        end
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
            local_clear_axes(ui.WaterfallAxes);
            title(ui.QEAxes, "Physical q-E Map");
            title(ui.ComparisonAxes, "Comparison Off-axis Component");
            title(ui.SingleAxes, "Spectrum");
            title(ui.DispersionAxes, "Dispersion");
            title(ui.WaterfallAxes, "Stacked Spectra");
            return
        end

        state.selectedQIndex = min(max(1, state.selectedQIndex), size(state.physicalQE.intensity, 2));

        % Render single spectrum FIRST — uses fast single-q preview if
        % the full batch cache is stale. This gives instant BG feedback.
        was_cache_miss = ~strcmp(state.ppCache.hash, local_opts_hash(local_preprocess_opts(), qe));
        local_plot_single_spectrum(qe);
        drawnow limitrate;  % flush the single spectrum to screen immediately

        % Then render the maps — this triggers full batch if cache misses,
        % but the user has already seen the single spectrum result.
        local_plot_qe_map(ui.QEAxes, state.physicalQE, "Physical");
        local_plot_normalized_map();
        local_plot_dispersion_result();
        local_plot_waterfall(qe);

        % If the cache was stale and deconv is on, the single spectrum
        % was drawn without deconv (fast preview skips it). Now that
        % the full batch is cached, redraw to ensure consistency. (Issue #2)
        if was_cache_miss && local_preprocess_opts().do_deconv
            local_plot_single_spectrum(qe);
        end


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
        [qe, ~, ~] = local_get_cached_preprocess(qe, local_preprocess_opts());
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
            "LabelHorizontalAlignment", "left", ...
            "Tag", "selected_q_marker");
        local_plot_dispersion_overlay(ax, qe);
        hold(ax, "off");
    end


    function local_plot_normalized_map()
        ax = ui.ComparisonAxes;
        local_clear_axes(ax);

        if isempty(state.comparisonQE)
            title(ax, "Comparison Off-axis Component | press Build Views");
            xlabel(ax, "q (1/A)");
            ylabel(ax, "Energy relative to ZLP (meV)");
            grid(ax, "on");
            return
        end

        local_plot_qe_map(ax, state.comparisonQE, "Comparison Off-axis Component", "normalized");
    end


    function local_plot_single_spectrum(qe)
        ax = ui.SingleAxes;
        local_clear_axes(ax);
        pp_opts = local_preprocess_opts();

        % Try cached full preprocessing first (O(1) if already computed)
        hash_str = local_opts_hash(pp_opts, qe);
        q_index = state.selectedQIndex;
        bg_diag_single = [];  % init before branch

        if strcmp(state.ppCache.hash, hash_str) && ~isempty(state.ppCache.qe)
            % Cache hit — use stored full batch result
            qe = state.ppCache.qe;
            qe_pre = state.ppCache.qe_pre;
            bg_diag_all = state.ppCache.bg_diag;
        else
            % Cache miss — use fast SINGLE-CHANNEL preview (instant)
            [qe, qe_pre, bg_diag_single] = local_get_single_q_bg_preview(qe, pp_opts, q_index);
            bg_diag_all = [];  % no batch diagnostics
        end

        [mask, energy_axis] = local_energy_mask(qe);
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

        % --- BG overlay: use cached or single-channel diagnostics ---
        bg_diag_qi = [];
        pre_bg_spectrum = [];
        if pp_opts.do_bg_sub && strcmpi(mode_name, "display")
            try
                if ~isempty(qe_pre)
                    pre_bg_spectrum = double(qe_pre.intensity(mask, q_index));
                end
                if ~isempty(bg_diag_all) && q_index <= numel(bg_diag_all)
                    bg_diag_qi = bg_diag_all(q_index);
                elseif ~isempty(bg_diag_single)
                    bg_diag_qi = bg_diag_single;
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
                    if isfield(bg_diag_qi, 'selected_method') && ~isempty(bg_diag_qi.selected_method)
                        selected_bg_method = char(bg_diag_qi.selected_method);
                    else
                        selected_bg_method = pp_opts.bg_method;
                    end
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
                    if isfield(bg_diag_qi, 'neg_fraction') && isfinite(bg_diag_qi.neg_fraction)
                        diag_parts{end+1} = sprintf('neg=%.1f%%', 100 * bg_diag_qi.neg_fraction);
                    end
                    if isfield(bg_diag_qi, 'neg_area_fraction') && isfinite(bg_diag_qi.neg_area_fraction)
                        diag_parts{end+1} = sprintf('negA=%.1f%%', 100 * bg_diag_qi.neg_area_fraction);
                    end
                    if isfield(bg_diag_qi, 'bg_fraction') && isfinite(bg_diag_qi.bg_fraction)
                        diag_parts{end+1} = sprintf('bg=%.1f%%', 100 * bg_diag_qi.bg_fraction);
                    end
                    if bg_diag_qi.iterations > 1
                        diag_parts{end+1} = sprintf('iter=%d', bg_diag_qi.iterations);
                    end
                    if ~isempty(diag_parts)
                        ui.FitInfoLabel.Text = sprintf('BG [%s]: %s', ...
                            selected_bg_method, strjoin(diag_parts, ' | '));
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




        local_plot_single_branch_fit_overlay(ax, qe, q_index, energy_axis, display_values);

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


    function local_plot_single_branch_fit_overlay(ax, qe, q_index, energy_axis, trace_values)
        if isempty(state.fitResults) || isempty(state.manual_branches)
            return
        end
        if isempty(qe) || q_index < 1 || q_index > numel(qe.q_Ainv)
            return
        end

        q_value = double(qe.q_Ainv(q_index));
        branch_colors = qe_plot_helpers.branch_colors();
        q_tol = local_branch_correction_q_tolerance();
        y_anchor = local_single_overlay_y_anchor(trace_values);

        for b = 1:numel(state.fitResults)
            if b > numel(state.manual_branches)
                continue
            end
            fit_result = state.fitResults{b};
            if isempty(fit_result)
                continue
            end

            col = branch_colors(mod(b-1, size(branch_colors, 1))+1, :);
            fit_energy = local_eval_branch_fit_at_q(fit_result, q_value);
            branch_energy = local_branch_energy_at_q(state.manual_branches{b}, q_value, q_tol);

            if isfinite(fit_energy) && fit_energy >= energy_axis(1) && fit_energy <= energy_axis(end)
                xline(ax, fit_energy, '--', ...
                    'Color', col, ...
                    'LineWidth', 1.6, ...
                    'Label', sprintf('B%d fit %.0f', b, fit_energy), ...
                    'LabelVerticalAlignment', 'bottom', ...
                    'LabelHorizontalAlignment', 'left', ...
                    'Tag', 'branch_fit_single_overlay');
            end

            if isfinite(branch_energy) && branch_energy >= energy_axis(1) && branch_energy <= energy_axis(end)
                plot(ax, branch_energy, y_anchor, 'o', ...
                    'MarkerSize', 7, ...
                    'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', col, ...
                    'LineWidth', 1.4, ...
                    'DisplayName', sprintf('B%d point %.0f', b, branch_energy), ...
                    'Tag', 'branch_fit_single_overlay');
            end
        end
    end


    function energy_meV = local_eval_branch_fit_at_q(fit_result, q_value)
        energy_meV = NaN;
        if isempty(fit_result) || ~isfield(fit_result, 'q_fit') || ~isfield(fit_result, 'E_fit')
            return
        end

        q_fit = double(fit_result.q_fit(:));
        e_fit = double(fit_result.E_fit(:));
        valid = isfinite(q_fit) & isfinite(e_fit);
        q_fit = q_fit(valid);
        e_fit = e_fit(valid);
        if numel(q_fit) < 2
            return
        end

        if all(q_fit >= 0) && q_value < 0
            query_q = abs(q_value);
        elseif all(q_fit <= 0) && q_value > 0
            query_q = -abs(q_value);
        else
            query_q = q_value;
        end

        [q_fit, order] = sort(q_fit);
        e_fit = e_fit(order);
        [q_unique, ia] = unique(q_fit, 'stable');
        e_unique = e_fit(ia);
        if query_q < min(q_unique) || query_q > max(q_unique)
            return
        end
        energy_meV = interp1(q_unique, e_unique, query_q, 'linear', NaN);
    end


    function energy_meV = local_branch_energy_at_q(branch_points, q_value, q_tol)
        energy_meV = NaN;
        if isempty(branch_points) || size(branch_points, 2) < 2
            return
        end
        [q_delta, idx] = min(abs(double(branch_points(:,1)) - q_value));
        if isempty(idx) || ~isfinite(q_delta) || q_delta > q_tol
            return
        end
        energy_meV = double(branch_points(idx, 2));
    end


    function y_anchor = local_single_overlay_y_anchor(values)
        finite_values = double(values(:));
        finite_values = finite_values(isfinite(finite_values));
        if isempty(finite_values)
            y_anchor = 1;
            return
        end

        if strcmpi(local_trace_y_scale(), "log")
            finite_values = finite_values(finite_values > 0);
            if isempty(finite_values)
                y_anchor = 1;
                return
            end
            y_min = min(finite_values);
            y_max = max(finite_values);
            y_anchor = exp(log(y_min) + 0.88 * (log(y_max) - log(y_min)));
            return
        end

        y_min = min(finite_values);
        y_max = max(finite_values);
        y_anchor = y_min + 0.88 * (y_max - y_min);
    end


    function local_plot_waterfall(qe)
        ax = ui.WaterfallAxes;
        local_clear_axes(ax);
        if isempty(ax) || ~isgraphics(ax)
            return
        end

        [qe, ~, ~] = local_get_cached_preprocess(qe, local_preprocess_opts());
        [mask, energy_axis] = local_energy_mask(qe);
        all_q_mask = true(1, size(qe.intensity, 2));
        spectra_matrix = local_build_map_values(qe, mask, all_q_mask, local_trace_mode());

        if strcmpi(local_trace_y_scale(), "log")
            spectra_matrix = max(spectra_matrix, eps);
        end

        stack_opts = struct( ...
            'q_start', ui.QStartField.Value, ...
            'q_end', ui.QEndField.Value, ...
            'q_step', ui.QStepField.Value, ...
            'offset', ui.OffsetField.Value);
        stack = qe_prepare_stacked_spectra(energy_axis, spectra_matrix, qe.q_Ainv, stack_opts);
        if isempty(stack.q_indices)
            title(ax, "Stacked Spectra");
            return
        end

        plotted_values = stack.shifted_traces(:);
        hold(ax, "on");
        for idx = 1:numel(stack.q_indices)
            q_index = stack.q_indices(idx);
            trace_color = local_signed_q_color(stack.q_values(idx), qe.q_Ainv);
            line_width = 1.0;
            if q_index == state.selectedQIndex
                trace_color = [0 0 0];
                line_width = 2.2;
            end

            plot(ax, stack.energy_meV, stack.shifted_traces(:, idx), "-", ...
                "Color", trace_color, ...
                "LineWidth", line_width, ...
                "HandleVisibility", "off");

            if numel(stack.q_indices) <= 12
                text(ax, stack.energy_meV(end), stack.shifted_traces(end, idx), ...
                    sprintf('  q=%.4f', stack.q_values(idx)), ...
                    'Color', trace_color, ...
                    'FontSize', 10, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle', ...
                    'Clipping', 'on');
            end
        end
        hold(ax, "off");

        grid(ax, "on");
        xlabel(ax, "Energy relative to ZLP (meV)");
        ylabel(ax, sprintf("%s + offset", local_single_ylabel()));
        title(ax, sprintf("Stacked spectra | %s | %s | %d traces", ...
            local_current_view_name(), local_trace_mode(), numel(stack.q_indices)), ...
            "Interpreter", "none");
        ax.YScale = local_trace_y_scale();
        xlim(ax, [stack.energy_meV(1), stack.energy_meV(end)]);
        local_apply_y_limits(ax, plotted_values);
    end


    function local_update_heatmap_branch_overlays()
        axes_to_update = [ui.QEAxes, ui.ComparisonAxes];
        for ax_idx = 1:numel(axes_to_update)
            ax = axes_to_update(ax_idx);
            if isempty(ax) || ~isgraphics(ax)
                continue
            end
            delete(findobj(ax, "Tag", "branch_fit_overlay"));
            hold(ax, "on");
            local_plot_dispersion_overlay(ax, []);
            hold(ax, "off");
        end
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
                    'LineWidth', 0.8, 'HandleVisibility', 'off', ...
                    'Tag', 'branch_fit_overlay');
            end
            local_plot_correction_markers(ax);
        elseif ~isempty(state.manual_points)
            pts = sortrows(state.manual_points, 1);
            mc = [0.1 0.8 0.2];
            scatter(ax, pts(:,1), pts(:,2), 36, 'filled', ...
                'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', mc, ...
                'LineWidth', 1.2, 'HandleVisibility', 'off', ...
                'Tag', 'branch_fit_overlay');
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
            local_plot_correction_markers(ax);
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
        [~, state.selectedQIndex] = min(abs(state.physicalQE.q_Ainv - q_value));

        local_update_selected_q_views();
    end


    function local_update_selected_q_views()
        qe = local_get_active_qe();
        if isempty(qe)
            return
        end
        local_plot_single_spectrum(qe);
        local_update_q_selection_markers();
        drawnow limitrate;
    end


    function local_update_q_selection_markers()
        if isempty(state.physicalQE)
            return
        end
        selected_q = state.physicalQE.q_Ainv(state.selectedQIndex);
        axes_to_update = [ui.QEAxes, ui.ComparisonAxes];
        for ax_idx = 1:numel(axes_to_update)
            ax = axes_to_update(ax_idx);
            if isempty(ax) || ~isgraphics(ax)
                continue
            end

            markers = findobj(ax, "Tag", "selected_q_marker");
            if isempty(markers)
                hold(ax, "on");
                xline(ax, selected_q, "-", ...
                    "Color", [0.95 0.95 0.98], ...
                    "LineWidth", 1.2, ...
                    "Label", sprintf("q=%.4f", selected_q), ...
                    "LabelHorizontalAlignment", "left", ...
                    "Tag", "selected_q_marker");
                hold(ax, "off");
                continue
            end

            for marker_idx = 1:numel(markers)
                try
                    markers(marker_idx).Value = selected_q;
                    markers(marker_idx).Label = sprintf("q=%.4f", selected_q);
                catch
                end
            end
        end
    end


    function local_on_single_axes_clicked(~, ~)
        if isempty(state.physicalQE)
            return
        end

        % --- Hybrid auto-fit correction mode ---
        if state.autoCorrectionArmed
            try
                if isempty(state.manual_branches)
                    ui.InfoLabel.Text = "Run Auto Fit first.";
                    ui.InfoLabel.Visible = "on";
                    return
                end
                current_point = ui.SingleAxes.CurrentPoint;
                energy_value = current_point(1, 1);
                [~, energy_axis] = local_energy_mask(state.physicalQE);
                if isempty(energy_axis)
                    return
                end
                energy_value = min(max(energy_value, energy_axis(1)), energy_axis(end));
                q_value = state.physicalQE.q_Ainv(state.selectedQIndex);
                branch_spec = local_selected_correction_branch_spec();

                [state.manual_branches, correction] = qe_apply_branch_correction( ...
                    state.manual_branches, q_value, energy_value, ...
                    'Branch', branch_spec, ...
                    'QTolerance', local_branch_correction_q_tolerance());
                local_append_branch_correction(correction);
                local_sync_manual_points_from_branches();

                ui.CorrectAutoButton.Text = "Correct Peak";
                ui.CorrectAutoButton.BackgroundColor = [0.96 0.96 0.96];
                state.autoCorrectionArmed = false;
                if isnan(correction.old_energy_meV)
                    ui.InfoLabel.Text = sprintf("Added Branch %d at q=%.4f: %.0f meV", ...
                        correction.branch_index, correction.new_q_Ainv, correction.new_energy_meV);
                else
                    ui.InfoLabel.Text = sprintf("Replaced Branch %d at q=%.4f: %.0f → %.0f meV", ...
                        correction.branch_index, correction.new_q_Ainv, ...
                        correction.old_energy_meV, correction.new_energy_meV);
                end
                ui.InfoLabel.Visible = "on";
                local_log_operation(sprintf('Correct auto peak: B%d q=%.4f %.0f meV', ...
                    correction.branch_index, correction.new_q_Ainv, correction.new_energy_meV));
                local_update_all_views();
            catch ME
                local_show_error(ME, false);
            end
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
                    local_log_operation(local_describe_changes());
                end
            catch ME
                local_show_error(ME, false);
            end
            return
        end
    end


    function local_sync_manual_points_from_branches()
        if ~isempty(state.manual_branches)
            state.manual_points = qe_flatten_branch_points(state.manual_branches);
        else
            state.manual_points = [];
        end
        ui.PtsLabel.Text = sprintf('%d pts', size(state.manual_points, 1));
        local_sync_correction_controls();
    end


    function local_append_branch_correction(correction)
        if isempty(state.branchCorrectionLog)
            state.branchCorrectionLog = correction;
        else
            state.branchCorrectionLog(end+1) = correction;
        end
        local_sync_correction_controls();
    end


    function branch_spec = local_selected_correction_branch_spec()
        branch_spec = 'auto';
        if ~isfield(ui, 'CorrectionBranchDropdown') || isempty(ui.CorrectionBranchDropdown)
            return
        end
        value = char(string(ui.CorrectionBranchDropdown.Value));
        value = strtrim(value);
        if isempty(value) || strcmpi(value, 'Auto')
            return
        end
        nums = regexp(value, '\d+', 'match');
        if ~isempty(nums)
            branch_spec = str2double(nums{end});
        end
    end

    function specs = local_branch_specs_from_ui()
        windows = [
            ui.Branch1MinField.Value, ui.Branch1MaxField.Value
            ui.Branch2MinField.Value, ui.Branch2MaxField.Value
            ui.Branch3MinField.Value, ui.Branch3MaxField.Value
            ];
        specs = repmat(struct('name', '', 'energy_window_meV', [0 0], 'enabled', true), 3, 1);
        for b = 1:3
            win = sort(double(windows(b, :)));
            if any(~isfinite(win)) || win(1) == win(2)
                error('interactive_qe_browser:invalidBranchWindow', ...
                    'Branch %d energy window must have two finite, different meV values.', b);
            end
            specs(b).name = sprintf('Branch %d', b);
            specs(b).energy_window_meV = win;
            specs(b).enabled = true;
        end
    end


    function local_sync_correction_controls()
        if ~isfield(ui, 'CorrectionBranchDropdown') || isempty(ui.CorrectionBranchDropdown)
            return
        end

        n_branches = numel(state.manual_branches);
        items = [{'Auto'}, arrayfun(@(b) sprintf('Branch %d', b), 1:n_branches, 'UniformOutput', false)];
        if n_branches > 0
            items{end+1} = sprintf('New Branch %d', n_branches + 1); %#ok<AGROW>
        end

        old_value = char(string(ui.CorrectionBranchDropdown.Value));
        ui.CorrectionBranchDropdown.Items = items;
        if any(strcmp(old_value, items))
            ui.CorrectionBranchDropdown.Value = old_value;
        else
            nums = regexp(old_value, 'New Branch\s+(\d+)', 'tokens', 'once');
            if ~isempty(nums)
                promoted = sprintf('Branch %s', nums{1});
                if any(strcmp(promoted, items))
                    ui.CorrectionBranchDropdown.Value = promoted;
                else
                    ui.CorrectionBranchDropdown.Value = 'Auto';
                end
            else
                ui.CorrectionBranchDropdown.Value = 'Auto';
            end
        end

        ui.CorrectionBranchDropdown.Enable = local_on_off(n_branches > 0);
        ui.CorrectAutoButton.Enable = local_on_off(n_branches > 0);
        ui.UndoCorrectionButton.Enable = local_on_off(~isempty(state.branchCorrectionLog));

        if isfield(ui, 'FitBranchDropdown') && ~isempty(ui.FitBranchDropdown)
            fit_items = [{'All'}, arrayfun(@(b) sprintf('Branch %d', b), ...
                1:n_branches, 'UniformOutput', false)];
            old_fit_value = char(string(ui.FitBranchDropdown.Value));
            ui.FitBranchDropdown.Items = fit_items;
            if any(strcmp(old_fit_value, fit_items))
                ui.FitBranchDropdown.Value = old_fit_value;
            else
                ui.FitBranchDropdown.Value = 'All';
            end
            ui.FitBranchDropdown.Enable = local_on_off(n_branches > 0);
        end
    end


    function q_tol = local_branch_correction_q_tolerance()
        q_tol = 1e-6;
        try
            dq = local_current_dq();
            if isfinite(dq) && dq > 0
                q_tol = max(0.55 * abs(dq), q_tol);
                return
            end
        catch
        end
        try
            qe = local_get_active_qe();
            if ~isempty(qe) && isfield(qe, 'q_Ainv') && numel(qe.q_Ainv) > 1
                dq_vals = diff(sort(unique(qe.q_Ainv(:))));
                dq_vals = dq_vals(isfinite(dq_vals) & dq_vals > 0);
                if ~isempty(dq_vals)
                    q_tol = max(0.55 * median(dq_vals), q_tol);
                end
            end
        catch
        end
    end


    function local_plot_correction_markers(ax)
        if isempty(state.branchCorrectionLog)
            return
        end
        try
            q_vals = [state.branchCorrectionLog.new_q_Ainv];
            e_vals = [state.branchCorrectionLog.new_energy_meV];
            if isempty(q_vals), return; end
            scatter(ax, q_vals, e_vals, 52, 'd', ...
                'MarkerFaceColor', [1.0 0.95 0.2], ...
                'MarkerEdgeColor', [0.05 0.05 0.05], ...
                'LineWidth', 0.9, ...
                'HandleVisibility', 'off', ...
                'Tag', 'branch_fit_overlay');
        catch
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
        local_sync_correction_controls();
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
            ui.InfoLabel.Text = "Top-left shows the original Physical q-E map. Top-right shows the comparison off-axis component reconstructed from manual symmetric off-axis reference bands.";
            return
        end

        if local_is_comparison_view_selected()
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
        switch lower(char(local_normalize_view_mode_value(ui.ViewModeDropdown.Value)))
            case {"comparison", "normalized"}
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


    function value = local_normalize_view_mode_value(value_in)
        value = char(string(value_in));
        if strcmpi(value, "Normalized")
            value = 'Comparison';
        elseif ~any(strcmp(value, {'Physical', 'Comparison'}))
            value = 'Physical';
        end
    end


    function tf = local_is_comparison_view_selected()
        tf = strcmpi(local_normalize_view_mode_value(ui.ViewModeDropdown.Value), 'Comparison');
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
            view_name = "Comparison";
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
                if strcmpi(local_current_view_name(), "Comparison")
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


    function [qe_qi, qe_pre_qi, bg_diag_qi] = local_get_single_q_bg_preview(qe_raw, pp_opts, q_index)
        % Fast single-channel BG preview — O(1) instead of O(n_q).
        %
        % Processes only the selected q-channel through the pipeline:
        %   normalize → denoise → BG subtract (single channel) + diagnostics
        %
        % Returns:
        %   qe_qi      — fully preprocessed data (only q_index modified by BG)
        %   qe_pre_qi  — preprocessed WITHOUT BG (for three-line overlay)
        %   bg_diag_qi — single BG diagnostic struct for q_index

        bg_diag_qi = [];
        qe_pre_qi = [];

        % 1. Pre-BG steps (norm + denoise) — these are fast, process all q
        pre_opts = pp_opts;
        pre_opts.do_bg_sub = false;
        pre_opts.do_deconv = false;
        pre_hash = local_opts_hash(pre_opts, qe_raw);
        if isfield(state, "singlePreviewCache") ...
                && strcmp(state.singlePreviewCache.hash, pre_hash) ...
                && ~isempty(state.singlePreviewCache.qe_pre)
            qe_pre_qi = state.singlePreviewCache.qe_pre;
        else
            qe_pre_qi = qe_preprocess(qe_raw, pre_opts);
            state.singlePreviewCache.hash = pre_hash;
            state.singlePreviewCache.qe_pre = qe_pre_qi;
        end

        % 2. BG subtraction on single q only
        if pp_opts.do_bg_sub
            bg_opts = pp_opts;
            bg_opts.do_normalize = false;  % already done
            bg_opts.do_denoise = false;    % already done
            bg_opts.do_deconv = false;
            bg_opts.q_indices = q_index;   % ← KEY: only process this channel
            [qe_qi, bg_diag_arr] = qe_preprocess(qe_pre_qi, bg_opts);
            if ~isempty(bg_diag_arr) && q_index <= numel(bg_diag_arr)
                bg_diag_qi = bg_diag_arr(q_index);
            end
        else
            qe_qi = qe_pre_qi;
        end
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
