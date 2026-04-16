function ULM_Master_GUI_v2()
% =========================================================================
% ULM_MASTER_GUI_V2 - Ultrasound Localization Microscopy Processing Suite
% =========================================================================
%
% DESCRIPTION:
%   Enhanced version with bug fixes, performance optimizations, and advanced features:
%   - Robust error handling and input validation
%   - Session save/load functionality
%   - Undo/Redo system for parameter changes
%   - Progress bars for long operations
%   - Memory profiling
%   - Parallel processing support
%   - Debounced display updates for smooth interaction
%   - Spatial Cropping to improve SVD performance
%
% FEATURES:
%   1. FILTER:   Spatial Crop, Clutter filtering (SVD, DCC, Butterworth)
%   2. DETECT:   Microbubble detection with ROI masking
%   3. LOCALIZE: Sub-pixel localization
%   4. TRACK:    Trajectory linking
%   5. PROCESS:  Track smoothing and interpolation
%   6. RENDER:   Super-resolution map generation
%
% AUTHOR: Grigori Shapiro 
% DATE:   February 2026
% =========================================================================

    % --- System Setup ---
    addpath(genpath(fileparts(mfilename('fullpath'))));
    disp('ULM System v2.0: Environment initialized.');
    
    % --- Initialize Application Structure FIRST ---
    app = struct();
    app.data = struct();
    app.ui = struct();
    app.state = struct();
    
    % --- Create Constants ---
    app.constants = ULM_Constants();
    
    % --- Initialize Main Figure ---
    fig = uifigure('Name', 'ULM Master GUI v2.0 - Professional Edition', ...
                   'Position', [50 50 1600 1000], ...
                   'CloseRequestFcn', @(src, evt) cleanupGUI(src));
    app.fig = fig;  % Store IMMEDIATELY
    
    % --- Initialize Session Manager ---
    app.sessionManager = SessionManager();
    
    % --- Initialize Undo/Redo System ---
    app.undoManager = UndoRedoManager();
    
    % --- Initialize Display Manager ---
    app.displayManager = DisplayManager();
    
    % --- Load Default Parameters with Safety Fallback ---
    isKidneyExperiment = true;
    try
        if exist('setDefaultParams', 'file')
            app.data.params = setDefaultParams(isKidneyExperiment);
            app.state.paramWarning = false;
        else
            error('setDefaultParams not found');
        end
    catch ME
        % If setDefaultParams is missing or crashes, use fallback GUI defaults
        app.state.paramWarning = true;
        app.data.params = createFallbackParams();
    end
    
    % --- Initialize Spatial Filter Params ---
    if ~isfield(app.data.params.filter, 'spatial_method')
        app.data.params.filter.spatial_method = 'Gaussian';
        app.data.params.filter.spatial_kernel = 3;
        app.data.params.filter.spatial_sigma1 = 1.0;
        app.data.params.filter.spatial_sigma2 = 2.0;
    end

    % Ensure crop_box exists in params
    if ~isfield(app.data.params, 'io')
        app.data.params.io = struct();
    end
    if ~isfield(app.data.params.io, 'crop_box')
        app.data.params.io.crop_box = [];
    end
    
    % --- Initialize Kalman Flags ---
    if ~isfield(app.data.params.track.kalman, 'use_direction')
        app.data.params.track.kalman.use_direction = true;
    end
    if ~isfield(app.data.params.track.kalman, 'use_angle')
        app.data.params.track.kalman.use_angle = true;
    end
    if ~isfield(app.data.params.track.kalman, 'use_brightness')
        app.data.params.track.kalman.use_brightness = true;
    end
    
    % --- Initialize Data Fields ---
    app.data.rawData = [];
    app.data.filteredData = [];
    app.data.candidateBubbles = [];
    app.data.localizations = [];
    app.data.tracks_raw = [];
    app.data.tracks_final = [];
    app.data.mask = [];
    app.data.vesselMap = [];
    app.data.baseVesselMap = [];
    
    % --- SVD Cache ---
    app.data.U = [];
    app.data.S_diag = [];
    app.data.V = [];
    app.data.svdDims = [];
    app.data.rawDataHash = '';  % Track data changes
    
    % --- DCC Components ---
    app.data.tissue_indices = [];
    app.data.blood_indices = [];
    app.data.noise_indices = [];
    
    % --- Visualization Properties ---
    app.data.rawClim = [0 1];
    app.data.filteredClim = [0 1];
    app.data.filteredMeanBG = [];
    app.data.MeanBG = [];
    app.data.filteredBGClim = [0 1];
    
    % --- State Machine ---
    app.state.currentState = -1;
    app.state.currentFrame = 1;
    app.state.maxFrame = 1;
    app.state.isPlaying = false;
    app.state.isROIPreview = false;
    app.state.isProcessing = false;
    
    % --- Timers ---
    app.state.playbackTimer = [];
    app.state.displayTimer = [];
    
    % --- Pre-computed Overlays ---
    app.data.redOverlayTemplate = [];
    
    % =========================================================================
    % --- GUI Layout Construction ---
    % =========================================================================
    
    % Save app before building (so functions can access fig)
    guidata(fig, app);
    
    % Build GUI (modifies app.ui)
    app = buildGUILayout(app);
    
    % --- Initialize State ---
    app = manageGUIState(app, -1);
    
    % --- Save app data BEFORE calling UI updaters ---
    guidata(fig, app);
    
    % --- Set up callbacks that need guidata ---
    app.ui.FilterMethodDropdown.ValueChangedFcn = @(s,e) updateFilterOptions(fig);
    app.ui.LocMethodDropdown.ValueChangedFcn = @(s,e) updateLocalizationOptions(fig);
    app.ui.TrackMethodDropdown.ValueChangedFcn = @(s,e) updateTrackingOptions(fig);
    app.ui.UseAdvancedCostCheckbox.ValueChangedFcn = @(s,e) updateTrackingOptions(fig);
    
    % Save again after callback setup
    guidata(fig, app);
    
    % --- Run UI updaters (AFTER guidata) ---
    updateFilterOptions(fig);
    updateLocalizationOptions(fig);
    updateTrackingOptions(fig);
    updateDetectionOptions(fig);

    % --- Set Ready Status ---
    setStatus(app, 'Ready', 'green');
    
    % --- Show Warning if Fallback was used ---
    if isfield(app.state, 'paramWarning') && app.state.paramWarning
        uialert(app.fig, ...
            'setDefaultParams.m is missing or failed to load. Using safe system GUI defaults (FPS: 200, Pixel Size: 0.05 mm).', ...
            'Warning: Missing Default Parameters', 'Icon', 'warning');
    end

    disp('GUI initialization complete.');
end

% =========================================================================
% --- GUI LAYOUT BUILDER ---
% =========================================================================

function app = buildGUILayout(app)
    fig = app.fig;
    
    % --- Main Grid: 3 rows (menu, main content, status) x 2 columns ---
    gl = uigridlayout(fig, [3, 2]);
    gl.RowHeight = {40, '1x', 40};  % Menu bar, main area, status bar
    gl.ColumnWidth = {'3x', '1x'};  % Visualization (wide), Controls (narrow)
    gl.Padding = [5 5 5 5];
    gl.RowSpacing = 5;
    gl.ColumnSpacing = 5;
    
    % --- Build components (all return modified app) ---
    app = buildMenuBar(app, gl);
    app = buildVisualizationPanel(app, gl);
    app = buildControlPanel(app, gl);
    app = buildStatusBar(app, gl);
end

function app = buildMenuBar(app, parentGrid)
    % Expanded grid layout to fit Fundamental Parameters on the right
    menuGrid = uigridlayout(parentGrid, [1, 11]);
    menuGrid.Layout.Row = 1;
    menuGrid.Layout.Column = [1 2];
    menuGrid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', '1x', 'fit', 60, 'fit', 60, 'fit', 60};
    menuGrid.Padding = [5 5 5 5];
    menuGrid.BackgroundColor = [0.94 0.94 0.94];  % Light gray background
    
    % Load/Save Session
    app.ui.btnLoadSession = uibutton(menuGrid, 'Text', '📁 Load Work Session', ...
        'ButtonPushedFcn', @(s,e) loadSession(app.fig), ...
        'Tooltip', 'Load previously saved GUI session');
    
    app.ui.btnSaveSession = uibutton(menuGrid, 'Text', '💾 Save Work Session', ...
        'ButtonPushedFcn', @(s,e) saveSession(app.fig), ...
        'Tooltip', 'Save current work for later');
    
    % Undo/Redo
    app.ui.btnUndo = uibutton(menuGrid, 'Text', '↶ Undo', ...
        'ButtonPushedFcn', @(s,e) performUndo(app.fig), ...
        'Enable', 'off');
    
    app.ui.btnRedo = uibutton(menuGrid, 'Text', '↷ Redo', ...
        'ButtonPushedFcn', @(s,e) performRedo(app.fig), ...
        'Enable', 'off');
    
    % Memory Monitor (Pushed to left by the '1x' space)
    app.ui.lblMemory = uilabel(menuGrid, 'Text', 'Memory: 0 MB', ...
        'HorizontalAlignment', 'center');
        
    % --- Fundamental Parameters (Right Side) ---
    uilabel(menuGrid, 'Text', 'FPS:', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    app.ui.TopFPSField = uieditfield(menuGrid, 'numeric', ...
        'Value', app.data.params.acq.framerate, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'framerate'));
        
    uilabel(menuGrid, 'Text', 'Px X (mm):', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    app.ui.TopPixelXField = uieditfield(menuGrid, 'numeric', ...
        'Value', app.data.params.track.pixel_X_size, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'pixel_x'));
        
    uilabel(menuGrid, 'Text', 'Px Z (mm):', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    app.ui.TopPixelZField = uieditfield(menuGrid, 'numeric', ...
        'Value', app.data.params.track.pixel_Z_size, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'pixel_z'));
    
    % Start memory monitor
    startMemoryMonitor(app.fig);
end

function app = buildVisualizationPanel(app, parentGrid)
    % Change grid to 2 columns: Left for visual controls, Right for the axes
    stageGrid = uigridlayout(parentGrid, [3, 2]);
    stageGrid.Layout.Row = 2;
    stageGrid.Layout.Column = 1;
    stageGrid.RowHeight = {'1x', 'fit', 30};
    stageGrid.ColumnWidth = {180, '1x'}; % 180px for the new visual controls
    stageGrid.Padding = [10 10 10 10];
    
    % --- Display Controls Panel (Left Column) ---
    app = buildDisplayControlsPanel(app, stageGrid);
    
    % --- Main Axes (Right Column) ---
    app.ui.ax = uiaxes(stageGrid);
    app.ui.ax.Layout.Row = 1;
    app.ui.ax.Layout.Column = 2;
    title(app.ui.ax, 'Please Load Raw Data to begin...');
    xlabel(app.ui.ax, 'X (pixels)');
    ylabel(app.ui.ax, 'Y (pixels)');
    axis(app.ui.ax, 'image');
    
    % --- Frame Slider (Right Column) ---
    app.ui.FrameSlider = uislider(stageGrid, 'Limits', [1 100], 'Value', 1, ...
        'ValueChangedFcn', @(s,e) onFrameSliderChanged(app.fig, e));
    app.ui.FrameSlider.Layout.Row = 2;
    app.ui.FrameSlider.Layout.Column = 2;
    
    % --- Bottom Controls (Spans both columns) ---
    bottomGrid = uigridlayout(stageGrid, [1, 3]);
    bottomGrid.Layout.Row = 3;
    bottomGrid.Layout.Column = [1 2];
    bottomGrid.ColumnWidth = {'2x', '1x', '1x'};
    bottomGrid.Padding = [0 0 0 0];
    
    app.ui.LoadButton = uibutton(bottomGrid, 'Text', '📂 Load Data (IQ / imageData)', ...
        'ButtonPushedFcn', @(s,e) loadRawData(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.3 0.7 1]);
    
    app.ui.PlayPauseButton = uibutton(bottomGrid, 'Text', 'Play', ...
        'ButtonPushedFcn', @(s,e) togglePlayback(app.fig));
    
    app.ui.FrameField = uieditfield(bottomGrid, 'numeric', 'Value', 1, ...
        'Limits', [1 Inf], 'RoundFractionalValues', 'on', ...
        'ValueChangedFcn', @(s,e) onFrameFieldChanged(app.fig, e));
end

function app = buildDisplayControlsPanel(app, parentGrid)
    % Create the side panel for visual adjustments
    p_disp = uipanel(parentGrid, 'Title', 'Visual Adjustments');
    p_disp.Layout.Row = 1;
    p_disp.Layout.Column = 1;

    g = uigridlayout(p_disp, [9, 1]);
    g.RowHeight = {'fit','fit','fit','fit','fit','fit','fit','fit','1x'};
    g.Padding = [5 5 5 5];

    % Callback to trigger real-time updates when sliders/fields are changed
    updateVis = @(s,e) displayCurrentFrame(app.fig);

    % Normalize toggle
    app.ui.disp_mat2gray = uicheckbox(g, 'Text', 'Normalize (mat2gray)', ...
        'Value', 1, 'ValueChangedFcn', updateVis);

    % Log Compression toggle
    app.ui.disp_log = uicheckbox(g, 'Text', 'Log Compression', ...
        'Value', 0, 'ValueChangedFcn', updateVis);

    % Gamma Slider
    uilabel(g, 'Text', 'Gamma (Stretch):', 'FontWeight', 'bold');
    app.ui.disp_gamma = uislider(g, 'Limits', [0.1 5.0], 'Value', 1.0, ...
        'ValueChangedFcn', updateVis);

    % Colormap Dropdown
    uilabel(g, 'Text', 'Colormap:', 'FontWeight', 'bold');
    app.ui.disp_cmap = uidropdown(g, 'Items', {'gray', 'hot', 'jet', 'parula'}, ...
        'Value', 'gray', 'ValueChangedFcn', updateVis);

    % CLim Settings
    app.ui.disp_clim_auto = uicheckbox(g, 'Text', 'Auto CLim', ...
        'Value', 1, 'ValueChangedFcn', updateVis);

    cg = uigridlayout(g, [1 2]);
    cg.Padding = [0 0 0 0];
    app.ui.disp_clim_min = uieditfield(cg, 'numeric', 'Value', 0, ...
        'Tooltip', 'Min CLim', 'ValueChangedFcn', updateVis);
    app.ui.disp_clim_max = uieditfield(cg, 'numeric', 'Value', 1, ...
        'Tooltip', 'Max CLim', 'ValueChangedFcn', updateVis);
end

function app = buildControlPanel(app, parentGrid)
    controlGrid = uigridlayout(parentGrid, [2, 1]);
    controlGrid.Layout.Row = 2;
    controlGrid.Layout.Column = 2;
    controlGrid.RowHeight = {'1x', 30};
    controlGrid.Padding = [10 10 10 10];
    
    % Tab Group
    app.ui.tabGroup = uitabgroup(controlGrid, ...
        'SelectionChangedFcn', @(s,e) onTabChanged(app.fig));
    app.ui.tabGroup.Layout.Row = 1;
    
    % Build all tabs
    app = buildFilterTab(app);
    app = buildDetectTab(app);
    app = buildLocalizeTab(app);
    app = buildTrackTab(app);
    app = buildPostProcessTab(app);
    app = buildRenderTab(app);
    
    % Reset Button
    app.ui.ResetParamsButton = uibutton(controlGrid, ...
        'Text', 'Reset All Params to Default', ...
        'ButtonPushedFcn', @(s,e) resetAllParams(app.fig));
    app.ui.ResetParamsButton.Layout.Row = 2;
end

function app = buildStatusBar(app, parentGrid)
    statusGrid = uigridlayout(parentGrid, [1, 3]);
    statusGrid.Layout.Row = 3;  % Bottom row
    statusGrid.Layout.Column = [1 2];  % Span both columns
    statusGrid.ColumnWidth = {'fit', '1x', 'fit'};
    statusGrid.Padding = [10 5 10 5];
    statusGrid.BackgroundColor = [0.94 0.94 0.94];  % Light gray background
    
    app.ui.lblStatusTitle = uilabel(statusGrid, 'Text', 'Status:', ...
        'FontWeight', 'bold');
    
    app.ui.StatusLabel = uilabel(statusGrid, 'Text', 'Ready');
    
    app.ui.StatusLamp = uilamp(statusGrid, 'Color', 'blue');
    
    % Progress dialog will be created on-demand in showProgress function
    app.ui.ProgressBar = [];
end

% =========================================================================
% --- TAB BUILDERS ---
% =========================================================================

function app = buildFilterTab(app)
    app.ui.tabFilter = uitab(app.ui.tabGroup, 'Title', '1. Filter');
    g_tab = uigridlayout(app.ui.tabFilter, [1 1]);
    app.ui.panel_filt = uipanel(g_tab, 'BorderType', 'none');
    
    % Increased to 11 rows to fit Spatial Filter
    g = uigridlayout(app.ui.panel_filt, [11, 1]);
    g.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
    
    % Method Selection
    uilabel(g, 'Text', 'Filter Method:', 'FontWeight', 'bold');
    app.ui.FilterMethodDropdown = uidropdown(g, ...
        'Items', {'svd_filter', 'svd_ssm', 'dcc_svd'}, ...
        'Value', app.data.params.filter.method);
        
    % Spatial Crop Panel
    app = buildCropPanel(app, g);
    app.ui.p_crop.Layout.Row = 4;
    
    % Masking Panel
    app.ui.p_mask = uipanel(g, 'Title', 'Masking (for Localization)');
    app.ui.p_mask.Layout.Row = 5;
    g_mask = uigridlayout(app.ui.p_mask, [2, 2]);
    app.ui.LoadMaskButton = uibutton(g_mask, 'Text', 'Load Mask', ...
        'ButtonPushedFcn', @(s,e) loadMask(app.fig));
    app.ui.CreateMaskButton = uibutton(g_mask, 'Text', 'Create New Mask', ...
        'ButtonPushedFcn', @(s,e) runCreateMask(app.fig));
    app.ui.ResetMaskButton = uibutton(g_mask, 'Text', 'Reset Mask', ...
        'ButtonPushedFcn', @(s,e) resetMask(app.fig));
    app.ui.maskStatusLabel = uilabel(g_mask, 'Text', 'Status: None', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % SVD Parameters
    app.ui.p_svd = uipanel(g, 'Title', 'Standard SVD Params');
    app.ui.p_svd.Layout.Row = 6;
    g_svd = uigridlayout(app.ui.p_svd, [2, 2]);
    uilabel(g_svd, 'Text', 'Cutoff Start:');
    app.ui.SVDCutoffStart = uieditfield(g_svd, 'numeric', ...
        'Value', app.data.params.filter.svd_cutoff(1), ...
        'ValueChangedFcn', @(s,e) runFilterWithValidation(app.fig));
    uilabel(g_svd, 'Text', 'Cutoff End:');
    app.ui.SVDCutoffEnd = uieditfield(g_svd, 'numeric', ...
        'Value', app.data.params.filter.svd_cutoff(2), ...
        'ValueChangedFcn', @(s,e) runFilterWithValidation(app.fig));
    
    % DCC Parameters
    app = buildDCCPanel(app, g);
    app.ui.p_dcc.Layout.Row = 7;
    
    % Butterworth Filter
    app = buildButterworthPanel(app, g);
    app.ui.p_butter.Layout.Row = 8;
    
    % Spatial Filter (NEW)
    app = buildSpatialFilterPanel(app, g);
    app.ui.p_spatial.Layout.Row = 9;
    
    % Run Button
    app.ui.RunFilterButton = uibutton(g, 'Text', 'Run Filter', ...
        'ButtonPushedFcn', @(s,e) runFilter(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunFilterButton.Layout.Row = 10;
end

function app = buildCropPanel(app, parentGrid)
    % Pre-Processing Crop Panel to reduce image dimensions and improve SVD
    app.ui.p_crop = uipanel(parentGrid, 'Title', 'Spatial Crop (Pre-Processing)');
    
    g_crop = uigridlayout(app.ui.p_crop, [2, 3]);
    g_crop.ColumnWidth = {'1x', '1x', '1x'};
    
    % Row 1
    uilabel(g_crop, 'Text', 'Crop Box [x y w h]:');
    app.ui.CropBoxField = uieditfield(g_crop, 'text', 'Value', '[]');
    app.ui.InteractiveCropBtn = uibutton(g_crop, 'Text', 'Interactive Crop', ...
        'ButtonPushedFcn', @(s,e) runInteractiveCrop(app.fig));
        
    % Row 2
    app.ui.LoadCropBtn = uibutton(g_crop, 'Text', 'Load Crop', ...
        'ButtonPushedFcn', @(s,e) loadCrop(app.fig));
    app.ui.SaveCropBtn = uibutton(g_crop, 'Text', 'Save Crop', ...
        'ButtonPushedFcn', @(s,e) saveCrop(app.fig));
    app.ui.ApplyCropBtn = uibutton(g_crop, 'Text', 'Apply Crop to Data', ...
        'ButtonPushedFcn', @(s,e) applyCrop(app.fig), ...
        'BackgroundColor', [1 0.8 0.6], 'FontWeight', 'bold');
end

function app = buildDCCPanel(app, parentGrid)
    app.ui.p_dcc = uipanel(parentGrid, 'Title', 'DCC-SVD Reconstruction');
    g = uigridlayout(app.ui.p_dcc, [3, 5]);
    g.ColumnWidth = {'fit', '1x', '0.5x', '1x', '0.5x'};
    
    % Tissue
    uilabel(g, 'Text', 'Tissue:', 'FontWeight', 'bold');
    uilabel(g, 'Text', 'Start (%):');
    app.ui.DCCTissueStart = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
    uilabel(g, 'Text', 'End (%):');
    app.ui.DCCTissueEnd = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
    
    % Blood
    uilabel(g, 'Text', 'Blood:', 'FontWeight', 'bold');
    uilabel(g, 'Text', 'Start (%):');
    app.ui.DCCBloodStart = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
    uilabel(g, 'Text', 'End (%):');
    app.ui.DCCBloodEnd = uieditfield(g, 'numeric', 'Value', 100, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
    
    % Noise
    uilabel(g, 'Text', 'Noise:', 'FontWeight', 'bold');
    uilabel(g, 'Text', 'Start (%):');
    app.ui.DCCNoiseStart = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
    uilabel(g, 'Text', 'End (%):');
    app.ui.DCCNoiseEnd = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
end

function app = buildButterworthPanel(app, parentGrid)
    app.ui.p_butter = uipanel(parentGrid, 'Title', 'Butterworth Filter (Optional)');
    g = uigridlayout(app.ui.p_butter, [3, 2]);
    
    app.ui.EnableButterworth = uicheckbox(g, 'Text', 'Enable', ...
        'Value', app.data.params.filter.enable_butterworth, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'butterworth'));
    app.ui.EnableButterworth.Layout.Column = [1 2];
    
    uilabel(g, 'Text', 'Cutoff [Low High]:');
    app.ui.ButterCutoff = uieditfield(g, 'text', ...
        'Value', mat2str(app.data.params.filter.butter_cutoff), ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'butter_cutoff'));
    
    uilabel(g, 'Text', 'Order:');
    app.ui.ButterOrder = uieditfield(g, 'numeric', ...
        'Value', app.data.params.filter.butter_order, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'butter_order'));
end

function app = buildSpatialFilterPanel(app, parentGrid)
    app.ui.p_spatial = uipanel(parentGrid, 'Title', 'Spatial Pre-filtering (Optional)');
    g = uigridlayout(app.ui.p_spatial, [2, 4]);
    g.ColumnWidth = {'fit', '1x', 'fit', '1x'};
    
    % Method Dropdown
    uilabel(g, 'Text', 'Method:');
    app.ui.SpatialMethodDrop = uidropdown(g, ...
        'Items', {'None', 'Gaussian', 'Median', 'DoG', 'Top-Hat'}, ...
        'Value', app.data.params.filter.spatial_method, ...
        'ValueChangedFcn', @(s,e) updateSpatialOptions(app.fig));
    
    % Kernel Size
    app.ui.lblSpatialKernel = uilabel(g, 'Text', 'Kernel (px):');
    app.ui.SpatialKernelField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.filter.spatial_kernel, ...
        'RoundFractionalValues', 'on', ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'spatial_kernel'));
    
    % Sigma 1
    app.ui.lblSpatialSigma1 = uilabel(g, 'Text', 'Sigma 1:');
    app.ui.SpatialSigma1Field = uieditfield(g, 'numeric', ...
        'Value', app.data.params.filter.spatial_sigma1, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'spatial_sigma1'));
        
    % Sigma 2
    app.ui.lblSpatialSigma2 = uilabel(g, 'Text', 'Sigma 2:');
    app.ui.SpatialSigma2Field = uieditfield(g, 'numeric', ...
        'Value', app.data.params.filter.spatial_sigma2, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'spatial_sigma2'));
        
    % Initial UI Sync
    updateSpatialOptions(app.fig, app);
end

function updateSpatialOptions(fig, app_in)
    if nargin < 2
        app = guidata(fig);
    else
        app = app_in;
    end
    
    method = app.ui.SpatialMethodDrop.Value;
    
    % Hide everything initially
    app.ui.lblSpatialKernel.Visible = 'off'; app.ui.SpatialKernelField.Visible = 'off';
    app.ui.lblSpatialSigma1.Visible = 'off'; app.ui.SpatialSigma1Field.Visible = 'off';
    app.ui.lblSpatialSigma2.Visible = 'off'; app.ui.SpatialSigma2Field.Visible = 'off';
    
    % Show based on method
    switch method
        case 'Gaussian'
            app.ui.lblSpatialKernel.Visible = 'on'; app.ui.SpatialKernelField.Visible = 'on';
            app.ui.lblSpatialSigma1.Text = 'Sigma:';
            app.ui.lblSpatialSigma1.Visible = 'on'; app.ui.SpatialSigma1Field.Visible = 'on';
        case 'Median'
            app.ui.lblSpatialKernel.Visible = 'on'; app.ui.SpatialKernelField.Visible = 'on';
        case 'Top-Hat'
            app.ui.lblSpatialKernel.Text = 'Radius:';
            app.ui.lblSpatialKernel.Visible = 'on'; app.ui.SpatialKernelField.Visible = 'on';
        case 'DoG'
            app.ui.lblSpatialSigma1.Text = 'Sigma 1:';
            app.ui.lblSpatialSigma1.Visible = 'on'; app.ui.SpatialSigma1Field.Visible = 'on';
            app.ui.lblSpatialSigma2.Visible = 'on'; app.ui.SpatialSigma2Field.Visible = 'on';
    end
    
    if nargin < 2
        guidata(fig, app);
        saveParamState(fig, 'spatial_method');
    end
    
end

function app = buildDetectTab(app)
    app.ui.tabDetect = uitab(app.ui.tabGroup, 'Title', '2. Detect');
    g_tab = uigridlayout(app.ui.tabDetect, [1 1]);
    app.ui.panel_detect = uipanel(g_tab, 'BorderType', 'none');
    
    g_main = uigridlayout(app.ui.panel_detect, [2, 1]);
    g_main.RowHeight = {'1.2x', '1x'};
    
    % ROI Panel
    app = buildROIPanel(app, g_main);
    
    % Detection Parameters
    app = buildDetectionParamsPanel(app, g_main);
end

function app = buildROIPanel(app, parentGrid)
    app.ui.p_roi = uipanel(parentGrid, 'Title', 'Step A: Define ROI Mask (Vessel Map)');
    app.ui.p_roi.Layout.Row = 1;
    
    g = uigridlayout(app.ui.p_roi, [2, 1]);
    g.RowHeight = {'1x', 'fit'};
    g.Padding = [5 5 5 5];
    
    % Histogram
    app.ui.axHist = uiaxes(g);
    app.ui.axHist.Layout.Row = 1;
    title(app.ui.axHist, 'Intensity Histogram (Log Scale)');
    app.ui.axHist.YScale = 'log';
    app.ui.axHist.XTickLabel = [];
    grid(app.ui.axHist, 'on');
    
    % Controls
    g_ctrl = uigridlayout(g, [4, 3]);
    g_ctrl.Layout.Row = 2;
    g_ctrl.RowHeight = {'fit', 'fit', 'fit', 30};
    g_ctrl.ColumnWidth = {110, 90, '1x'};
    
    % Enhancement
    uilabel(g_ctrl, 'Text', '1. Enhance:', 'FontWeight', 'bold');
    app.ui.EnhanceMethodDrop = uidropdown(g_ctrl, ...
        'Items', {'None', 'CLAHE (Local Contrast)', 'Top-Hat (Vesselness)', 'Sharpen'}, ...
        'Value', 'None', ...
        'ValueChangedFcn', @(s,e) applyVesselEnhancement(app.fig));
    app.ui.EnhanceAmountSlider = uislider(g_ctrl, 'Limits', [0 1], 'Value', 0.5, ...
        'ValueChangedFcn', @(s,e) applyVesselEnhancement(app.fig, e.Value));
    
    % Gamma
    uilabel(g_ctrl, 'Text', '2. Gamma:', 'FontWeight', 'bold');
    app.ui.ROIContrastField = uieditfield(g_ctrl, 'numeric', 'Value', 1, ...
        'ValueDisplayFormat', '%.3f', ...
        'ValueChangedFcn', @(s,e) onContrastChange(app.fig, e.Value));
    app.ui.ROIContrastSlider = uislider(g_ctrl, 'Limits', [0.1 3.0], 'Value', 1, ...
        'ValueChangedFcn', @(s,e) onContrastChange(app.fig, e.Value));
    
    % Threshold
    uilabel(g_ctrl, 'Text', '3. Threshold:', 'FontWeight', 'bold');
    app.ui.ROIThreshField = uieditfield(g_ctrl, 'numeric', 'Value', 0, ...
        'ValueDisplayFormat', '%.3f', ...
        'ValueChangedFcn', @(s,e) onROIChange(app.fig, e.Value));
    app.ui.ROIThreshSlider = uislider(g_ctrl, 'Limits', [0 1], 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onROIChange(app.fig, e.Value));
    
    % Save Button
    app.ui.SaveMaskButton = uibutton(g_ctrl, ...
        'Text', 'Save Current Mask (vesselMask)', ...
        'ButtonPushedFcn', @(s,e) saveCreatedMask(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.8 0.9 1]);
    app.ui.SaveMaskButton.Layout.Row = 4;
    app.ui.SaveMaskButton.Layout.Column = [1 3];
end

function app = buildDetectionParamsPanel(app, parentGrid)
    app.ui.p_detect_params = uipanel(parentGrid, 'Title', 'Step B: Detection');
    app.ui.p_detect_params.Layout.Row = 2;
    
    g = uigridlayout(app.ui.p_detect_params, [7, 3]);
    g.RowHeight = {'fit','fit','fit','fit','fit','fit'};
    g.ColumnWidth = {'fit', 60, '1x'};
    
    % --- NEW ROW 1: Detection Method Dropdown ---
    uilabel(g, 'Text', 'Detect Method:', 'FontWeight', 'bold');
    app.ui.DetectMethodDropdown = uidropdown(g, ...
        'Items', {'Intensity', 'NP', 'NCC'}, ...
        'Value', app.data.params.loc.DetectMethod, ...
        'ValueChangedFcn', @(s,e) updateDetectionOptions(app.fig));
    app.ui.DetectMethodDropdown.Layout.Column = [2 3];
    
    % --- NEW ROW 2: Method-Specific Parameters Panel ---
    app.ui.p_detect_method_params = uipanel(g, 'Title', 'Method Parameters');
    app.ui.p_detect_method_params.Layout.Row = 2;
    app.ui.p_detect_method_params.Layout.Column = [1 3];
    g_mp = uigridlayout(app.ui.p_detect_method_params, [2, 2]);
    g_mp.ColumnWidth = {'fit', '1x'};
    g_mp.RowHeight = {'fit', 'fit'};
    
    % NP params
    app.ui.lbl_NP_alpha = uilabel(g_mp, 'Text', 'NP alpha0:');
    app.ui.NP_AlphaField = uieditfield(g_mp, 'numeric', ...
        'Value', app.data.params.loc.NP_alpha0, ...
        'ValueDisplayFormat', '%.2e');
    
    % NCC params
    app.ui.lbl_NCC_thresh = uilabel(g_mp, 'Text', 'NCC tau:');
    app.ui.NCC_ThreshField = uieditfield(g_mp, 'numeric', ...
        'Value', app.data.params.loc.crosscor_threshold, ...
        'ValueDisplayFormat', '%.2f');
    
    % --- PSF FWHM Window Size ---
    uilabel(g, 'Text', 'PSF FWHM [x z] (px):');
    app.ui.DetectFWHMField = uieditfield(g, 'text', ...
        'Value', mat2str(app.data.params.loc.fwhm), ...
        'Tooltip', 'Bubble PSF window size, e.g. [3 3]. Used for localization ROI and NCC template.');
    app.ui.DetectFWHMField.Layout.Column = [2 3];

    % Intensity Threshold
    uilabel(g, 'Text', 'Intensity Thresh:');
    app.ui.LocThreshField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.loc.detection_threshold);
    app.ui.LocThreshSlider = uislider(g, 'Limits', [0.01 1], ...
        'Value', app.data.params.loc.detection_threshold);
    syncSliderField(app.ui.LocThreshSlider, app.ui.LocThreshField);
    
    % Max Bubbles
    uilabel(g, 'Text', 'Max Bubbles:');
    app.ui.LocMaxField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.loc.max_bubbles_per_frame);
    app.ui.LocMaxSlider = uislider(g, 'Limits', [0 5000], ...
        'Value', app.data.params.loc.max_bubbles_per_frame);
    syncSliderField(app.ui.LocMaxSlider, app.ui.LocMaxField);
    
    % Preview Checkbox
    app.ui.chkPreviewROI = uicheckbox(g, 'Text', 'Preview ROI Overlay', 'Value', 1, ...
        'ValueChangedFcn', @(s,e) onROIPreviewToggle(app.fig));
    app.ui.chkPreviewROI.Layout.Column = [1 3];
    
    % Run Button
    app.ui.RunDetectButton = uibutton(g, 'Text', 'Run Detection (Masked)', ...
        'ButtonPushedFcn', @(s,e) runDetection(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunDetectButton.Layout.Row = 7;
    app.ui.RunDetectButton.Layout.Column = [1 3];
end

function app = buildLocalizeTab(app)
    app.ui.tabLocalize = uitab(app.ui.tabGroup, 'Title', '3. Localize');
    g_tab = uigridlayout(app.ui.tabLocalize, [1 1]);
    app.ui.panel_loc = uipanel(g_tab, 'BorderType', 'none');
    g = uigridlayout(app.ui.panel_loc, [7, 1]);
    g.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
    
    % Method
    uilabel(g, 'Text', 'Localization Method:', 'FontWeight', 'bold');
    app.ui.LocMethodDropdown = uidropdown(g, ...
        'Items', {'radial', 'gaussian_fit'}, ...
        'Value', 'radial');  % Default to radial
    
    % QC Panels
    app = buildLocalizationQCPanels(app, g);
    
    % Run Button
    app.ui.RunLocalizationButton = uibutton(g, 'Text', 'Run Localization', ...
        'ButtonPushedFcn', @(s,e) runLocalization(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunLocalizationButton.Layout.Row = 6;
end

function app = buildLocalizationQCPanels(app, parentGrid)
    % Radial/Gaussian QC
    app.ui.p_loc_qc_radial = uipanel(parentGrid, 'Title', 'Radial Symmetry QC');
    app.ui.p_loc_qc_radial.Layout.Row = 3;
    g = uigridlayout(app.ui.p_loc_qc_radial, [3, 2]);
    g.ColumnWidth = {'fit', '1x'};
    
    app.ui.LocQCDivergence = uicheckbox(g, 'Text', 'Enable Divergence Check', ...
        'Value', app.data.params.loc.enable_divergence_check);
    app.ui.LocQCDivergence.Layout.Column = [1 2];
    
    uilabel(g, 'Text', 'Max Shift Factor:');
    app.ui.LocShiftFactor = uieditfield(g, 'numeric', ...
        'Value', app.data.params.loc.qc_max_shift_factor);
    
    app.ui.LocQCRoiMaxima = uicheckbox(g, 'Text', 'Enable ROI Maxima Check', ...
        'Value', app.data.params.loc.enable_roi_maxima_check);
    app.ui.LocQCRoiMaxima.Layout.Column = [1 2];
    
    % Gaussian-specific Params
    app.ui.p_loc_qc_gauss = uipanel(parentGrid, 'Title', 'Gaussian Fit QC');
    app.ui.p_loc_qc_gauss.Layout.Row = 4;
    g2 = uigridlayout(app.ui.p_loc_qc_gauss, [3, 2]);
    
    uilabel(g2, 'Text', 'FWHM [x z] (px):');
    app.ui.LocFWHM = uieditfield(g2, 'text', ...
        'Value', mat2str(app.data.params.loc.fwhm));
    
    uilabel(g2, 'Text', 'Box Radius (px):');
    app.ui.GaussBoxRadius = uieditfield(g2, 'numeric', ...
        'Value', app.data.params.loc.gauss_fit_box_radius, ...
        'RoundFractionalValues', 'on');
    
    uilabel(g2, 'Text', 'Min R² (fit quality):');
    app.ui.GaussMinRSquared = uieditfield(g2, 'numeric', ...
        'Value', app.data.params.loc.min_r_squared);
end

function app = buildTrackTab(app)
    app.ui.tabTrack = uitab(app.ui.tabGroup, 'Title', '4. Track');
    g_tab = uigridlayout(app.ui.tabTrack, [1 1]);
    app.ui.panel_track = uipanel(g_tab, 'BorderType', 'none');
    g = uigridlayout(app.ui.panel_track, [7, 1]);
    g.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
    
    % Method
    uilabel(g, 'Text', 'Tracking Method:', 'FontWeight', 'bold');
    app.ui.TrackMethodDropdown = uidropdown(g, ...
        'Items', {'Kalman', 'Hungarian', 'nn', 'Kalman_Advanced'}, ...
        'Value', app.data.params.track.method);
    
    % Core Parameters
    app = buildTrackCoreParams(app, g);
    
    % Kalman Settings
    app = buildKalmanPanel(app, g);
    
    % QC Filters
    app = buildTrackQCPanel(app, g);
    
    % Run Button
    app.ui.RunTrackingButton = uibutton(g, 'Text', 'Run Tracking', ...
        'ButtonPushedFcn', @(s,e) runTracking(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunTrackingButton.Layout.Row = 6;
end

function app = buildTrackCoreParams(app, parentGrid)
    p = uipanel(parentGrid, 'Title', 'Core Tracking Params');
    p.Layout.Row = 3;
    g = uigridlayout(p, [3, 3]);
    g.ColumnWidth = {'fit', '1x', '0.5x'};
    
    % Max Link Distance
    uilabel(g, 'Text', 'Max Link Dist (px):');
    app.ui.MaxDistSlider = uislider(g, 'Limits', [0.1 10], ...
        'Value', app.data.params.track.max_linking_distance);
    app.ui.MaxDistField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.max_linking_distance, ...
        'ValueDisplayFormat', '%.1f');
    syncSliderField(app.ui.MaxDistSlider, app.ui.MaxDistField);
    
    % Max Gap Frames
    uilabel(g, 'Text', 'Max Gap Frames:');
    app.ui.GapFramesSlider = uislider(g, 'Limits', [0 10], ...
        'Value', app.data.params.track.max_gap_closing_frames);
    app.ui.GapFramesField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.max_gap_closing_frames, ...
        'RoundFractionalValues', 'on');
    syncSliderField(app.ui.GapFramesSlider, app.ui.GapFramesField);
    
    % Min Track Length
    uilabel(g, 'Text', 'Min Track Length:');
    app.ui.MinLengthSlider = uislider(g, 'Limits', [2 20], ...
        'Value', app.data.params.track.min_track_length);
    app.ui.MinLengthField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.min_track_length, ...
        'RoundFractionalValues', 'on');
    syncSliderField(app.ui.MinLengthSlider, app.ui.MinLengthField);
end

function app = buildKalmanPanel(app, parentGrid)
    app.ui.p_kalman_adv = uipanel(parentGrid, 'Title', 'Kalman Settings');
    app.ui.p_kalman_adv.Layout.Row = 4;
    g = uigridlayout(app.ui.p_kalman_adv, [5, 3]);
    g.ColumnWidth = {'fit', '1x', '0.5x'};
    
    uilabel(g, 'Text', 'Model:');
    app.ui.KalmanModelDrop = uidropdown(g, ...
        'Items', {'ConstantVelocity', 'ConstantAcceleration'}, ...
        'Value', app.data.params.track.kalman.motion_model);
    app.ui.KalmanModelDrop.Layout.Column = [2 3];
    
    uilabel(g, 'Text', 'Process Noise:');
    app.ui.KalmanNoise = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.kalman.process_noise);
    app.ui.KalmanNoise.Layout.Column = [2 3];
    
    uilabel(g, 'Text', 'Assignment:');
    app.ui.AssignmentDrop = uidropdown(g, ...
        'Items', {'hungarian', 'nn'}, ...
        'Value', app.data.params.track.kalman.assignment_method);
    app.ui.AssignmentDrop.Layout.Column = [2 3];
    
    app.ui.UseAdvancedCostCheckbox = uicheckbox(g, ...
        'Text', 'Use Advanced Cost Matrix', ...
        'Value', app.data.params.track.use_advanced_cost_matrix);
    app.ui.UseAdvancedCostCheckbox.Layout.Row = 4;
    
    app.ui.BtnConfigCostMatrix = uibutton(g, ...
        'Text', 'Configure Advanced Cost Matrix', ...
        'ButtonPushedFcn', @(s,e) openCostMatrixGUI(app.fig));
    app.ui.BtnConfigCostMatrix.Layout.Row = 4;
    app.ui.BtnConfigCostMatrix.Layout.Column = [2 3];

    % Hierarchical Kalman Configuration Button ---
    app.ui.BtnConfigHK = uibutton(g, ...
        'Text', 'Configure Hierarchical Kalman (HK)', ...
        'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 1], ...
        'ButtonPushedFcn', @(s,e) openHKConfigGUI(app.fig));
    app.ui.BtnConfigHK.Layout.Row = 5;
    app.ui.BtnConfigHK.Layout.Column = [1 3];
end

function app = buildTrackQCPanel(app, parentGrid)
    p = uipanel(parentGrid, 'Title', 'Track QC Filters');
    p.Layout.Row = 5;
    g = uigridlayout(p, [3, 3]);
    g.ColumnWidth = {'fit', '1x', '0.5x'};
    
    app.ui.TrackQCDirection = uicheckbox(g, 'Text', 'Direction Constraint', ...
        'Value', app.data.params.track.qc.enable_direction_constraint);
    uilabel(g, 'Text', 'Max Angle (°):', 'HorizontalAlignment', 'right');
    app.ui.QCMaxAngle = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.qc.max_angle_change_deg);
    
    app.ui.TrackQCAcceleration = uicheckbox(g, 'Text', 'Accel Constraint', ...
        'Value', app.data.params.track.qc.enable_acceleration_constraint);
    uilabel(g, 'Text', 'Factor C:', 'HorizontalAlignment', 'right');
    app.ui.QCAccelFactor = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.qc.acceleration_C_factor);
    
    app.ui.TrackQCVD = uicheckbox(g, 'Text', 'VD (Jitter)', ...
        'Value', app.data.params.track.qc.enable_vd_constraint);
    uilabel(g, 'Text', 'Max Ratio:', 'HorizontalAlignment', 'right');
    app.ui.QCVDRatio = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.qc.max_vd_ratio);
end

function app = buildPostProcessTab(app)
    app.ui.tabPostProcess = uitab(app.ui.tabGroup, 'Title', '5. Post-Process');
    g_tab = uigridlayout(app.ui.tabPostProcess, [1 1]);
    app.ui.panel_post = uipanel(g_tab, 'BorderType', 'none');
    g = uigridlayout(app.ui.panel_post, [4, 1]);
    g.RowHeight = {'fit', 'fit', 'fit', '1x'};
    
    % Smoothing
    p_smooth = uipanel(g, 'Title', 'Track Smoothing');
    p_smooth.Layout.Row = 1;
    g_smooth = uigridlayout(p_smooth, [2, 3]);
    g_smooth.ColumnWidth = {'fit', '1x', '0.5x'};
    g_smooth.RowHeight = {22, 50};
    g_smooth.Padding = [5 5 5 12];
    
    app.ui.EnablePostProcessing = uicheckbox(g_smooth, 'Text', 'Enable Smoothing', ...
        'Value', app.data.params.track.enable_postprocessing);
    app.ui.EnablePostProcessing.Layout.Column = [1 3];
    
    uilabel(g_smooth, 'Text', 'Window Size:');
    app.ui.SmoothSlider = uislider(g_smooth, 'Limits', [3 21], ...
        'Value', app.data.params.track.smoothing_factor);
    app.ui.SmoothField = uieditfield(g_smooth, 'numeric', ...
        'Value', app.data.params.track.smoothing_factor, ...
        'RoundFractionalValues', 'on');
    syncSliderField(app.ui.SmoothSlider, app.ui.SmoothField);
    
    % Display Filter
    app.ui.p_display_filter = uipanel(g, 'Title', 'Final Display Filter (Live)');
    app.ui.p_display_filter.Layout.Row = 2;
    g_disp = uigridlayout(app.ui.p_display_filter, [1, 3]);
    g_disp.ColumnWidth = {'fit', '1x', '0.5x'};
    g_disp.RowHeight = {50};
    g_disp.Padding = [5 5 5 12]; 
    
    uilabel(g_disp, 'Text', 'Min Length:');
app.ui.DisplayMinLengthSlider = uislider(g_disp, 'Limits', [2 20], ...
    'Value', app.data.params.track.min_track_length, ...
    'ValueChangedFcn', @(s,e) onDisplayFilterChanged(app.fig, e));
app.ui.DisplayMinLengthField = uieditfield(g_disp, 'numeric', ...
    'Value', app.data.params.track.min_track_length, ...
    'RoundFractionalValues', 'on', ...
    'ValueChangedFcn', @(s,e) onDisplayFilterChanged(app.fig, e));
    
    % Run Button
    app.ui.RunPostProcessButton = uibutton(g, ...
        'Text', 'Run Post-Processing (Smoothing)', ...
        'ButtonPushedFcn', @(s,e) runPostProcessing(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunPostProcessButton.Layout.Row = 3;
end

function app = buildRenderTab(app)
    app.ui.tabRender = uitab(app.ui.tabGroup, 'Title', '6. Render');
    g_tab = uigridlayout(app.ui.tabRender, [1 1]);
    app.ui.panel_render = uipanel(g_tab, 'BorderType', 'none');
    g = uigridlayout(app.ui.panel_render, [4, 1]);
    g.RowHeight = {'fit', 'fit', 'fit', '1x'};
    
    % Settings
    p_settings = uipanel(g, 'Title', 'Rendering Settings');
    p_settings.Layout.Row = 1;
    g_set = uigridlayout(p_settings, [2, 2]);
    
    uilabel(g_set, 'Text', 'Upsampling Factor:');
    app.ui.UpsamplingField = uieditfield(g_set, 'numeric', ...
        'Value', app.data.params.render.upsampling_factor);
    
    uilabel(g_set, 'Text', 'Render Method:');
    app.ui.RenderMethodDrop = uidropdown(g_set, ...
        'Items', {'histogram', 'gaussian'}, ...
        'Value', app.data.params.render.method);
    
    % Generate Button
    app.ui.GenerateImagesButton = uibutton(g, ...
        'Text', 'Generate & Display Final Images (New Windows)', ...
        'ButtonPushedFcn', @(s,e) runRendering(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 0.8 1]);
    app.ui.GenerateImagesButton.Layout.Row = 2;
    
    uilabel(g, 'Text', 'Note: This will open 4 separate figure windows.', ...
        'FontAngle', 'italic');
end

% =========================================================================
% --- CROP CALLBACKS (Pre-SVD Data Reduction) ---
% =========================================================================

function runInteractiveCrop(fig)
    app = guidata(fig);
    
    if isempty(app.data.rawData)
        uialert(fig, 'Please load Raw Data first to define crop.', 'Error');
        return;
    end
    
    % Create temporary figure for crop selection
    f = figure('Name', 'Interactive Crop - Draw Rectangle');
    
    % Check if filtered data exists, otherwise fallback to raw data
    if ~isempty(app.data.filteredData)
        meanImg = mean(abs(app.data.filteredData), 3);
        dispTitle = 'Draw a rectangle to define the crop area (Filtered Data).';
    else
        meanImg = mean(abs(app.data.rawData), 3);
        dispTitle = 'Draw a rectangle to define the crop area (Raw Data).';
    end
    
    imagesc(meanImg);
    colormap gray;
    axis image;
    title([dispTitle, ' Close window to cancel.']);
    
    rect = drawrectangle('Label', 'Crop ROI', 'Color', 'r');
    wait(rect);
    
    if isvalid(rect)
        pos = round(rect.Position); % returns [x, y, w, h]
        app.ui.CropBoxField.Value = mat2str(pos);
        if ~isfield(app.data.params, 'io')
            app.data.params.io = struct();
        end
        app.data.params.io.crop_box = pos;
        guidata(fig, app);
    end
    
    % FIX: Always close the temporary figure to prevent zombie windows
    if isvalid(f)
        close(f);
    end
    
end

function applyCrop(fig)
    app = guidata(fig);
    if isempty(app.data.rawData)
        uialert(fig, 'No raw data loaded.', 'Error');
        return;
    end
    
    try
        pos = str2double(app.ui.CropBoxField.Value);
        if isempty(pos) || length(pos) ~= 4
            error('Invalid crop box format. Expected [x y width height].');
        end
        
        [H, W, ~] = size(app.data.rawData);
        
        x_start = max(1, round(pos(1)));
        y_start = max(1, round(pos(2)));
        w = round(pos(3));
        h = round(pos(4));
        x_end = min(W, x_start + w - 1);
        y_end = min(H, y_start + h - 1);
        
        if x_end <= x_start || y_end <= y_start
            error('Invalid dimensions for crop area. Width and height must be positive.');
        end
        
        if x_start == 1 && y_start == 1 && x_end == W && y_end == H
            uialert(fig, 'Crop area covers the full image. No changes made.', 'Info');
            return;
        end
        
        % Warn the user that this reduces the raw matrix permanently in this session
        selection = uiconfirm(fig, ...
            sprintf('Cropping will reduce data size from [%dx%d] to [%dx%d]. This action is permanent for the current session and will clear Undo history. Proceed?', W, H, x_end-x_start+1, y_end-y_start+1), ...
            'Confirm Spatial Crop', ...
            'Options', {'Proceed', 'Cancel'}, ...
            'DefaultOption', 1, 'CancelOption', 2);
            
        if strcmp(selection, 'Cancel')
            return;
        end
        
        showProgress(fig, 'Applying spatial crop to raw data...', true);
        
        % Crop operation
        app.data.rawData = app.data.rawData(y_start:y_end, x_start:x_end, :);
        
        % Invalidate caches & masks
        app.data.rawDataHash = DataHash(app.data.rawData);
        app.data.U = [];
        app.data.S_diag = [];
        app.data.V = [];
        app.data.svdDims = [];
        app.data.mask = [];
        app.ui.maskStatusLabel.Text = 'Status: None';
        app.data.vesselMap = [];
        app.data.baseVesselMap = [];
        
        % Clear downstream dependencies
        app.data = clearDownstreamData(app.data, 0);
        
        % FIX: Synchronize GUI state to disable downstream panels
        app = manageGUIState(app, 0);
        
        % FIX: Clear Undo history because data dimensions fundamentally changed
        app.undoManager.clear();
        app.ui.btnUndo.Enable = 'off';
        app.ui.btnRedo.Enable = 'off';
        
        if ~isfield(app.data.params, 'io')
            app.data.params.io = struct();
        end
        app.data.params.io.crop_box = [x_start, y_start, x_end-x_start+1, y_end-y_start+1];
        
        % Reset raw Clim based on new boundaries
        abs_data = abs(app.data.rawData(:));
        app.data.rawClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
        if app.data.rawClim(1) == app.data.rawClim(2)
            app.data.rawClim(2) = app.data.rawClim(1) + 1;
        end

        guidata(fig, app);
        
    catch ME
        hideProgress(fig);
        uialert(fig, ME.message, 'Crop Error');
    end
    
    uialert(fig, 'Crop applied successfully. Data dimensions reduced and Undo history cleared.', 'Success');

    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function saveCrop(fig)
    app = guidata(fig);
    pos = str2double(app.ui.CropBoxField.Value);
    if isempty(pos) || length(pos) ~= 4
        uialert(fig, 'Invalid crop box to save.', 'Error');
        return;
    end
    [file, path] = uiputfile('cropBox.mat', 'Save Crop Box As...');
    if file ~= 0
        crop_box = pos;
        save(fullfile(path, file), 'crop_box');
        uialert(fig, 'Crop box saved successfully.', 'Success');
    end
end

function loadCrop(fig)
    app = guidata(fig);
    [file, path] = uigetfile('*.mat', 'Load Crop Box');
    if file ~= 0
        data = load(fullfile(path, file));
        if isfield(data, 'crop_box')
            app.ui.CropBoxField.Value = mat2str(data.crop_box);
            if ~isfield(app.data.params, 'io')
                app.data.params.io = struct();
            end
            app.data.params.io.crop_box = data.crop_box;
            guidata(fig, app);
        else
            uialert(fig, 'No "crop_box" variable found in the selected file.', 'Error');
        end
    end
end

% =========================================================================
% --- MAIN PROCESSING CALLBACKS ---
% =========================================================================

function loadRawData(fig)
    app = guidata(fig);
    
    try
        defaultPath = app.data.params.io.data_folder;
        [file, path] = uigetfile(fullfile(defaultPath, '*.mat'), ...
            'Select Raw Data .mat file');
        
        if isequal(file, 0)
            return;
        end
        
        % Show progress
        showProgress(fig, 'Loading data...', true);
        
        fprintf('Loading %s...\n', fullfile(path, file));
        data = load(fullfile(path, file));
        
        % Smart data detection - find 3D matrix regardless of variable name
        fields = fieldnames(data);
        data_matrix = [];
        variable_name = '';
        
        fprintf('  Scanning file contents...\n');
        for i = 1:length(fields)
            field_content = data.(fields{i});
            
            % Check if it's a 3D numeric matrix
            if isnumeric(field_content) && ndims(field_content) == 3
                data_matrix = field_content;
                variable_name = fields{i};
                fprintf('  Found 3D data: "%s" [%d x %d x %d]\n', ...
                    variable_name, size(data_matrix, 1), size(data_matrix, 2), size(data_matrix, 3));
                break;
            end
        end
        
        if isempty(data_matrix)
            error('No 3D numeric data matrix found in file.\nFile contains: %s', strjoin(fields, ', '));
        end
        
        % Detect if data is IQ (complex) or already processed (real)
        if ~isreal(data_matrix)
            fprintf('  Data type: IQ (Complex) - will use abs() for display\n');
            app.data.isIQData = true;
        else
            fprintf('  Data type: Real (already processed)\n');
            app.data.isIQData = false;
        end
        
        % Store data
        app.data.rawData = data_matrix;
        app.data.rawDataVariableName = variable_name;
        [H, W, T] = size(app.data.rawData);
        
        fprintf('  Dimensions: %d x %d x %d frames\n', H, W, T);
        fprintf('  Data range: [%.2e, %.2e]\n', min(abs(data_matrix(:))), max(abs(data_matrix(:))));
        
        % Compute data hash for tracking changes
        app.data.rawDataHash = DataHash(app.data.rawData);
        
        % Clear SVD cache
        app.data.U = [];
        app.data.S_diag = [];
        app.data.V = [];
        app.data.svdDims = [];
        
        % Reset state
        app.state.maxFrame = T;
        app.state.currentFrame = 1;
        app.data.mask = [];
        app.ui.maskStatusLabel.Text = 'Status: None';
        
        % Calculate color limits (always use abs for display)
        abs_data = abs(app.data.rawData(:));
        app.data.rawClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
        if app.data.rawClim(1) == app.data.rawClim(2)
            app.data.rawClim(2) = app.data.rawClim(1) + 1;
        end
        
        % Update UI
        app.ui.FrameSlider.Limits = [1 T];
        app.ui.FrameSlider.Value = 1;
        app.ui.FrameField.Limits = [1 T];
        app.ui.FrameField.Value = 1;

        app.data.params.io.data_folder = path;
        
        % Clear downstream
        app.data = clearDownstreamData(app.data, 0);
        app = manageGUIState(app, 0);
        
        % Save state for undo
        saveParamState(fig, 'data_loaded');
        
        setStatus(app, sprintf('Loaded: %s (%dx%dx%d)', variable_name, H, W, T), 'green');
        guidata(fig, app);
        fprintf('✓ Data loaded successfully\n');
        
    catch ME
        errordlg(ME.message, 'Load Error');
        setStatus(app, 'Load failed', 'red');
    end
    
    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function runFilterWithValidation(fig)
    app = guidata(fig);
    
    % Validate SVD cutoffs
    startVal = app.ui.SVDCutoffStart.Value;
    endVal = app.ui.SVDCutoffEnd.Value;
    
    if startVal >= endVal
        uialert(fig, 'SVD Start must be less than End', 'Invalid Range');
        return;
    end
    
    runFilter(fig);
end

function runFilter(fig)
    app = guidata(fig);
    if isempty(app.data.rawData)
        errordlg('No raw data loaded.', 'Error');
        return;
    end
    app.state.isProcessing = true;
    guidata(fig, app);
    
    try
        showProgress(fig, 'Running filter...', true);
        
        % Update parameters
        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;
        
        % Save state
        saveParamState(fig, 'filter');
        
        rawData = app.data.rawData;
        currentHash = DataHash(rawData);
        filteredData = [];
        
        switch params.filter.method
            case 'svd_filter'
                % Check if SVD needs recalculation
                if isempty(app.data.U) || isempty(app.data.S_diag) || ...
                   ~strcmp(app.data.rawDataHash, currentHash)
                    
                    showProgress(fig, 'Computing SVD (first time only)...', true);
                    [app.data.U, app.data.S_diag, app.data.V, app.data.svdDims] = ...
                        run_SVD_Decomposition(rawData);
                    
                    % CRITICAL: Update hash after computing SVD
                    app.data.rawDataHash = currentHash;
                    
                    fprintf('  SVD computed and cached.\n');
                else
                    fprintf('  Using cached SVD (no recomputation needed).\n');
                end
                
                % Fast reconstruction
                showProgress(fig, 'Reconstructing filtered signal...', false);
                cutoff = params.filter.svd_cutoff;
                filteredData = reconstruct_SVD_Signal(app.data.U, ...
                    app.data.S_diag, app.data.V, app.data.svdDims, cutoff);
                
            case 'svd_ssm'
                filteredData = SVD_SSM(rawData, 'IndentPrefix', '  ');
                
            case 'dcc_svd'
                showProgress(app, 'Computing DCC-SVD...', true);
                
                [filteredData, dccInfo] = DCC_SVD(rawData, params.acq.framerate, ...
                    'ReconstructionMode', 'blood', ...
                    'DensityPercentile',  10, ...
                    'CanopySeparation',   2.0, ...
                    'PlotResults',        true, ...
                    'IndentPrefix',       '  ');
                
                % Cache the SVD decomposition for fast reconstruction when sliders change.
                app.data.U       = dccInfo.U;         
                app.data.S_diag  = dccInfo.singular_values;  
                app.data.V       = dccInfo.V;         
                app.data.svdDims = size(rawData);
                app.data.rawDataHash = currentHash;
                
                % Store cluster indices for interactive slider-based reconstruction.
                app.data.tissue_indices = dccInfo.tissue_indices;
                app.data.blood_indices  = dccInfo.blood_indices;
                app.data.noise_indices  = dccInfo.noise_indices;
                
                fprintf('  DCC clusters - Tissue: %d, Blood: %d, Noise: %d components.\n', ...
                    numel(dccInfo.tissue_indices), ...
                    numel(dccInfo.blood_indices), ...
                    numel(dccInfo.noise_indices));

                filteredData = reconstructDCCImage(app);
        end
        
        % Apply Butterworth if enabled
        if params.filter.enable_butterworth
            showProgress(fig, 'Applying Butterworth filter...', false);
            filteredData = Butterworth_bandpass_filter(filteredData, ...
                params.filter.butter_cutoff, params.acq.framerate, ...
                params.filter.butter_order);
        end
        
        % Convert to envelope (abs) BEFORE Spatial Filtering 
        % This is required because Median and Top-Hat are non-linear intensity operations
        filteredData = abs(filteredData);
        
        % Apply Spatial Conditioning
        if ~strcmp(params.filter.spatial_method, 'None')
            showProgress(fig, sprintf('Applying %s Spatial Filter...', params.filter.spatial_method), false);
            filteredData = applySpatialFilter(filteredData, params.filter);
        end
        
        % Store results
        app.data.filteredData = filteredData;
        
        % Calculate color limits
        abs_data = app.data.filteredData(:);
        app.data.filteredClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
        if app.data.filteredClim(1) == app.data.filteredClim(2)
            app.data.filteredClim(2) = app.data.filteredClim(1) + 1;
        end
        
        % Calculate background
        mean_bg = mean(app.data.filteredData, 3);
        app.data.filteredMeanBG = mean_bg .^ 0.5;
        app.data.MeanBG = mean(abs(rawData), 3) .^ 0.5;
        
        bg_clim_data = app.data.filteredMeanBG(:);
        app.data.filteredBGClim = [prctile(bg_clim_data, 1), prctile(bg_clim_data, 99)];
        
        % Clear ROI maps to force recalculation
        app.data.baseVesselMap = [];
        app.data.vesselMap = [];
        app.state.isROIPreview = false;
        
        % Clear downstream
        app.data = clearDownstreamData(app.data, 1);
        app = manageGUIState(app, 1);
        
        setStatus(app, 'Filter complete', 'green');
        guidata(fig, app);
        
    catch ME
        errordlg(sprintf('%s\n\n%s', ME.message, ME.getReport('basic')), 'Filter Error');
        setStatus(app, 'Filter failed', 'red');
    end
    
    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function runDetection(fig)
    app = guidata(fig);
    if isempty(app.data.filteredData)
        errordlg('No filtered data. Please run filtering first.', 'Error');
        return;
    end
    app.state.isProcessing = true;
    guidata(fig, app);
    
    % Ensure mask exists
    if isempty(app.data.mask)
        uialert(fig, 'No ROI Mask defined. Using full image.', 'Warning');
        [h, w, ~] = size(app.data.filteredData);
        app.data.mask = true(h, w);
    end
    
    try
        showProgress(fig, 'Detecting bubbles...', true);
        
        app.data.params = updateParamsFromGUI(app);
        saveParamState(fig, 'detection');
        
        % Read method-specific params from GUI into params before dispatching
        app.data.params.loc.DetectMethod      = app.ui.DetectMethodDropdown.Value;
        app.data.params.loc.NP_alpha0         = app.ui.NP_AlphaField.Value;
        app.data.params.loc.crosscor_threshold = app.ui.NCC_ThreshField.Value;
        
        % Dispatch to correct detector (mask is passed internally for NP/NCC,
        % and has already been applied to filteredData for Intensity — see
        % Step1_Filter note in ULM_Processor).
        roiMask = app.data.mask;  % may be empty — all detectors handle this
        
        switch upper(app.data.params.loc.DetectMethod)
            case 'INTENSITY'
                app.data.candidateBubbles = detectBubbles( ...
                    app.data.filteredData, app.data.params.loc, roiMask);
        
            case 'NP'
                app.data.candidateBubbles = detectBubbles_NP( ...
                    app.data.filteredData, app.data.params.loc, roiMask);
        
            case 'NCC'
                if ~isfield(app.data.params.loc, 'MB_image') || ...
                   isempty(app.data.params.loc.MB_image)
                   
                    answer = uiconfirm(fig, ...
                        ['No PSF template (MB_image) found. ' ...
                         'Auto-generate a Gaussian template from the current FWHM setting?'], ...
                        'NCC: Missing PSF Template', ...
                        'Options', {'Generate Gaussian PSF', 'Cancel'}, ...
                        'DefaultOption', 1, 'CancelOption', 2);
                    
                    if strcmp(answer, 'Cancel')
                        hideProgress(fig);
                        guidata(fig, app);
                        return;
                    end
                    
                    % Build Gaussian PSF from current FWHM params
                    fwhm   = app.data.params.loc.fwhm;
                    sz     = app.data.params.loc.psf_size;
                    sigma_x = fwhm(1) / 2.355;
                    sigma_z = fwhm(2) / 2.355;
                    % Enforce odd size
                    sz(mod(sz,2)==0) = sz(mod(sz,2)==0) + 1;
                    [Xg, Zg] = meshgrid(-(sz(1)-1)/2:(sz(1)-1)/2, -(sz(2)-1)/2:(sz(2)-1)/2);
                    app.data.params.loc.MB_image = exp(-(Xg.^2/(2*sigma_x^2) + Zg.^2/(2*sigma_z^2)));
                    fprintf('   -> Auto-generated Gaussian PSF [%dx%d], sigma=[%.2f, %.2f] px\n', ...
                        sz(1), sz(2), sigma_x, sigma_z);
                end
                
                app.data.candidateBubbles = detectBubbles_NCC( ...
                    app.data.filteredData, app.data.params.loc, roiMask);
        
            otherwise
                error('Unknown DetectMethod: %s', app.data.params.loc.DetectMethod);
        end
        
        if isempty(app.data.candidateBubbles)
            warning('No bubbles detected. Try adjusting parameters.');
        end
        
        app.data = clearDownstreamData(app.data, 2);
        app = manageGUIState(app, 2);
        
        % Turn off ROI preview
        app.state.isROIPreview = false;
        app.ui.chkPreviewROI.Value = 0;
        
        setStatus(app, sprintf('%d bubbles detected', height(app.data.candidateBubbles)), 'green');
        guidata(fig, app);
        
    catch ME
        errordlg(sprintf('%s\n\n%s', ME.message, ME.getReport('basic')), 'Detection Error');
        setStatus(app, 'Detection failed', 'red');
    end
    
    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function runLocalization(fig)
    app = guidata(fig);
    if isempty(app.data.candidateBubbles)
        errordlg('No candidates. Please run detection first.', 'Error');
        return;
    end
    app.state.isProcessing = true;
    guidata(fig, app);
    
    try
        showProgress(fig, 'Localizing particles...', true);
        
        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;
        saveParamState(fig, 'localization');
        
        dataToLocalize = app.data.filteredData;
        candidateBubbles = app.data.candidateBubbles;
        
        % Validate data
        if isempty(dataToLocalize)
            error('Filtered data is empty. Please run filtering first.');
        end
        
        % Capture QC output by redirecting stdout
        diary_file = [tempname '.txt'];
        diary(diary_file);
        try
            switch params.loc.method
                case 'radial'
                    locs = localizeRadialSymmetry(dataToLocalize, candidateBubbles, params.loc, '');
                case 'gaussian_fit'
                    locs = fit2DGaussian_Fast(dataToLocalize, candidateBubbles, params.loc, '');
                otherwise
                    error('Unknown localization method: %s', params.loc.method);
            end
            diary off;
        catch ME_loc
            diary off;  % FIX B-14: always close diary
            if exist(diary_file, 'file'), delete(diary_file); end
            rethrow(ME_loc);  % let outer catch handle display
        end
        
        % Read QC output (unchanged)
        fid = fopen(diary_file, 'r');
        if fid ~= -1
            qc_output = fread(fid, '*char')';
            fclose(fid);
            delete(diary_file);
        else
            qc_output = '';
        end
        
        % Validate output
        if isempty(locs)
            warning('No localizations found. Check parameters.');
        end
        
        app.data.localizations = locs;
        app.data = clearDownstreamData(app.data, 3);
        app = manageGUIState(app, 3);
        
        guidata(fig, app);
        % Hide progress before showing QC dialog
        hideProgress(fig);
        
        % Show QC summary dialog
        if ~isempty(qc_output)
            showQCDialog(fig, 'Localization QC Summary', qc_output);
        end
        
        setStatus(app, sprintf('%d particles localized', height(locs)), 'green');
        
        
    catch ME
        hideProgress(fig);
        errordlg(sprintf('%s\n\n%s', ME.message, ME.getReport('basic')), 'Localization Error');
        setStatus(app, 'Localization failed', 'red');
    end
    
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function locs = fit2DGaussianSafe(filteredData, candidateBubbles, locParams, indent_prefix)
    % Safe wrapper for fit2DGaussian that processes in batches to avoid memory issues
    
    num_candidates = height(candidateBubbles);
    batch_size = 500; % Process 500 at a time
    
    if num_candidates <= batch_size
        % Small enough, process normally
        locs = fit2DGaussian_Fast(filteredData, candidateBubbles, locParams, indent_prefix);
        return;
    end
    
    % Process in batches
    fprintf('Processing %d candidates in batches of %d...\n', num_candidates, batch_size);
    num_batches = ceil(num_candidates / batch_size);
    locs_cell = cell(num_batches, 1);
    
    for i = 1:num_batches
        start_idx = (i-1) * batch_size + 1;
        end_idx = min(i * batch_size, num_candidates);
        
        fprintf('  Batch %d/%d (candidates %d-%d)...\n', i, num_batches, start_idx, end_idx);
        
        batch_candidates = candidateBubbles(start_idx:end_idx, :);
        locs_cell{i} = fit2DGaussian_Fast(filteredData, batch_candidates, locParams, '    ');
        
        % Clear some memory
        if i < num_batches
            pause(0.1); % Brief pause to allow garbage collection
        end
    end
    
    % Combine all batches
    locs = vertcat(locs_cell{:});
    fprintf('Total localizations: %d\n', height(locs));
end

function showQCDialog(fig, titleStr, qc_text)
    d = uifigure('Name', titleStr, 'Position', [100 100 640 440]);
    
    gl = uigridlayout(d, [2, 1]);
    gl.RowHeight  = {'1x', 36};
    gl.Padding    = [10 10 10 10];
    gl.RowSpacing = 6;
    
    ta = uitextarea(gl, ...
        'Value',    strsplit(qc_text, newline), ...
        'Editable', 'off', ...
        'FontName', 'Courier New', ...
        'FontSize', 10);
    ta.Layout.Row = 1;
    
    btnGrid = uigridlayout(gl, [1, 3]);
    btnGrid.Layout.Row    = 2;
    btnGrid.ColumnWidth   = {'1x', 120, '1x'};
    btnGrid.Padding       = [0 0 0 0];
    uilabel(btnGrid, 'Text', '');   % spacer left
    uibutton(btnGrid, 'Text', 'OK', ...
        'ButtonPushedFcn', @(~,~) close(d));
    uilabel(btnGrid, 'Text', '');   % spacer right
end

function runTracking(fig)
    app = guidata(fig);
    if isempty(app.data.localizations)
        errordlg('No localizations. Please run localization first.', 'Error');
        return;
    end
    app.state.isProcessing = true;
    guidata(fig, app);
    
    try
        showProgress(fig, 'Tracking particles...', true);
        
        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;
        saveParamState(fig, 'tracking');
        
        % Enable parallel pool if available
        if ~isempty(gcp('nocreate'))
            fprintf('Parallel pool detected, using parallel processing...\n');
        end
        
        switch lower(params.track.method)
            case 'hungarian'
                tracks = trackHungarian(app.data.localizations, params, '  ');
            case 'nn'
                tracks = trackNearestNeighbor(app.data.localizations, params, '  ');
            case 'kalman'
                tracks = trackKalman(app.data.localizations, params, '  ');
            case 'kalman_advanced'
                tracks = trackKalman_Advanced(app.data.localizations, params, '  ');
            otherwise
                error('Unknown tracking method: %s', params.track.method);
        end
        
        if isempty(tracks)
            warning('No tracks generated. Check parameters.');
        end
        
        app.data.tracks_raw = applyQualityControl(tracks, params.track, '  ');
        app.data = clearDownstreamData(app.data, 4);
        app = manageGUIState(app, 4);
        
        setStatus(app, sprintf('%d tracks generated', length(app.data.tracks_raw)), 'green');
        guidata(fig, app);
        
    catch ME
        errordlg(sprintf('%s\n\n%s', ME.message, ME.getReport('basic')), 'Tracking Error');
        setStatus(app, 'Tracking failed', 'red');
    end
    
    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function runPostProcessing(fig)
    app = guidata(fig);
    if isempty(app.data.tracks_raw)
        errordlg('No tracks. Please run tracking first.', 'Error');
        return;
    end
    app.state.isProcessing = true;
    guidata(fig, app);
    
    try
        showProgress(fig, 'Post-processing tracks...', true);
        
        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;
        saveParamState(fig, 'postprocess');
        
        if ~params.track.enable_postprocessing
            % FIX B-04: Raw tracks need the same field set as processed tracks.
            % runRendering() and displayProcessedTracks() require original_length.
            raw = app.data.tracks_raw;
            numRaw = numel(raw);
            for k = 1:numRaw
                raw(k).original_length         = raw(k).length;
                raw(k).velocities_mm_s         = zeros(raw(k).length, 1);
                raw(k).average_velocity_mm_s   = 0;
            end
            app.data.tracks_final = raw;
        else
            raw_tracks = app.data.tracks_raw;
            smoothing_factor = params.track.smoothing_factor;
            interp_step = params.render.interpolation_step;
            
            % Parallel processing
            numTracks = length(raw_tracks);
            processed_tracks_cell = cell(1, numTracks);
            
            % Use parfor if parallel pool exists
            if ~isempty(gcp('nocreate'))
                parfor i = 1:numTracks
                    processed_tracks_cell{i} = processTrack(raw_tracks(i), ...
                        smoothing_factor, interp_step, params);
                end
            else
                for i = 1:numTracks
                    processed_tracks_cell{i} = processTrack(raw_tracks(i), ...
                        smoothing_factor, interp_step, params);
                    
                    % Update progress
                    if mod(i, 100) == 0
                        showProgress(fig, sprintf('Processing track %d/%d...', i, numTracks), false);
                    end
                end
            end
            
            app.data.tracks_final = [processed_tracks_cell{:}];
        end
        
        app = manageGUIState(app, 5);
        setStatus(app, 'Post-processing complete', 'green');
        guidata(fig, app);
        
    catch ME
        errordlg(sprintf('%s\n\n%s', ME.message, ME.getReport('basic')), 'Post-Process Error');
        setStatus(app, 'Post-processing failed', 'red');
    end
    
    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function processed = processTrack(track, smoothing_factor, interp_step, params)
    % Helper function for track processing
    path_to_process = double(track.path);
    act_win = min(smoothing_factor, track.length);
    if mod(act_win, 2) == 0
        act_win = act_win - 1;
    end
    
    % Safe Savitzky-Golay filtering
    poly_order = 3;
    min_window = poly_order + 2;
    
    if track.length >= min_window && act_win >= min_window
        path_to_process = [sgolayfilt(track.path(:,1), poly_order, act_win), ...
                          sgolayfilt(track.path(:,2), poly_order, act_win)];
    elseif track.length > 3
        % Fallback: moving average
        path_to_process = [movmean(track.path(:,1), min(3, track.length)), ...
                          movmean(track.path(:,2), min(3, track.length))];
    end
    
    % Interpolation
    orig_inds = 1:track.length;
    interp_inds = 1:interp_step:track.length;
    frames_interp = interp1(orig_inds, track.frames, interp_inds, 'linear');
    x_interp = fillmissing(interp1(orig_inds, path_to_process(:,1), interp_inds, 'pchip'), ...
        'linear', 'EndValues', 'nearest');
    y_interp = fillmissing(interp1(orig_inds, path_to_process(:,2), interp_inds, 'pchip'), ...
        'linear', 'EndValues', 'nearest');
    
    new_path = [x_interp', y_interp'];
    new_len = size(new_path, 1);
    
    % Velocity calculation
    if new_len > 1
        d_mm = sqrt((diff(new_path(:,1)) * params.track.pixel_X_size).^2 + ...
                    (diff(new_path(:,2)) * params.track.pixel_Z_size).^2);
        dt = diff(frames_interp') * params.track.dt;
        vels = zeros(size(dt));
        valid = dt > 1e-9;
        vels(valid) = d_mm(valid) ./ dt(valid);
        vels = [vels; vels(end)];
        avg_vel = repmat(mean(vels, 'omitnan'), new_len, 1);
    else
        vels = 0;
        avg_vel = 0;
    end
    
    processed = struct('id', track.id, 'path', new_path, 'frames', frames_interp', ...
        'length', new_len, 'original_length', track.length, ...
        'velocities_mm_s', vels, 'average_velocity_mm_s', avg_vel, ...
        'localizations', []);
end

function runRendering(fig)
    app = guidata(fig);
    app.state.isProcessing = true;
    guidata(fig, app);

    if isempty(app.data.tracks_final)
        errordlg('No processed tracks. Please run post-processing first.', 'Error');
        return;
    end
    
    try
        showProgress(fig, 'Generating maps...', true);
        
        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;
        
        tracks = app.data.tracks_final;
        min_len = app.ui.DisplayMinLengthField.Value;
        tracks = tracks([tracks.original_length] >= min_len);
        
        if isempty(tracks)
            error('No tracks meet minimum length criterion.');
        end
        
        upscale = params.render.upsampling_factor;
        [H, W, ~] = size(app.data.rawData);
        H_SR = round(H * upscale);
        W_SR = round(W * upscale);
        
        % Concatenate all track coordinate and velocity data at once.
        allX = round(vertcat_field(tracks, 'path', 1) * upscale);  % [N_total x 1]
        allY = round(vertcat_field(tracks, 'path', 2) * upscale);  % [N_total x 1]
        allV = vertcat_field(tracks, 'velocities_mm_s');            % [N_total x 1]
        
        valid = allX >= 1 & allX <= W_SR & allY >= 1 & allY <= H_SR;
        allX  = allX(valid);
        allY  = allY(valid);
        allV  = allV(valid);
        
        inds = sub2ind([H_SR, W_SR], allY, allX);   % linear indices
        N    = H_SR * W_SR;
        
        densityMap    = reshape(accumarray(inds, 1,    [N, 1], @sum, 0), H_SR, W_SR);
        velocityAccum = reshape(accumarray(inds, allV, [N, 1], @sum, 0), H_SR, W_SR);
        
        velocityMap = zeros(size(densityMap));
        mask = densityMap > 0;
        velocityMap(mask) = velocityAccum(mask) ./ densityMap(mask);
        
        % Rendering parameters
        res = params.render.upsampling_factor;
        res_pts = [1, 3, 5];
        font_pts = [6, 8, 10];
        line_pts = [2, 4, 5];
        res_c = max(min(res, res_pts(end)), res_pts(1));
        fs = round(interp1(res_pts, font_pts, res_c));
        lw = round(interp1(res_pts, line_pts, res_c), 1);
        
        % Close progress dialog BEFORE creating figures
        hideProgress(fig);
        drawnow;  % Ensure UI updates
        
        % Create figures
        createDensityFigure(densityMap, min_len, params, lw, fs);
        createVelocityFigures(velocityMap, params, lw, fs);
        createCombinedFigure(velocityMap, densityMap, params, lw, fs);
        
        setStatus(app, 'Rendering complete', 'green');
        guidata(fig, app);
        
    catch ME
        hideProgress(fig);  % Also close on error
        errordlg(sprintf('%s\n\n%s', ME.message, ME.getReport('basic')), 'Render Error');
        setStatus(app, 'Rendering failed', 'red');
    end
    
    % Save app state back to guidata
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function out = vertcat_field(tracks, fieldName, colIdx)
% VERTCAT_FIELD  Concatenate a field from a struct array vertically.
% colIdx is optional — if provided, selects that column from a matrix field.
    parts = cell(numel(tracks), 1);
    for k = 1:numel(tracks)
        v = tracks(k).(fieldName);
        if nargin >= 3
            v = v(:, colIdx);
        end
        parts{k} = v(:);   % ensure column vector
    end
    out = vertcat(parts{:});
end

% =========================================================================
% --- RENDERING HELPER FUNCTIONS ---
% =========================================================================

function createDensityFigure(densityMap, min_len, params, lw, fs)
    figure('Name', 'Density Map');
    d_proc = densityMap.^(1/3);
    imshow(d_proc, []);
    colormap(hot);
    if any(d_proc(:) > 0)
        clim([0 prctile(d_proc(d_proc>0), 99.5)]);
    end
    title(sprintf('Density (MinLen=%d)', min_len));
    colorbar;
    add_scale_bar(params, size(densityMap), lw, fs);
end

function createVelocityFigures(velocityMap, params, lw, fs)
    cm = [0 0 0; jet(256)];
    
    % Filtered
    figure('Name', 'Velocity (Filtered)');
    v_filt = imgaussfilt(velocityMap, 0.6);
    imshow(v_filt, []);
    colormap(cm);
    if any(v_filt(:) > 0)
        clim([0 prctile(v_filt(v_filt>0), 99.5)]);
    end
    title('Velocity Filtered');
    colorbar;
    add_scale_bar(params, size(velocityMap), lw, fs);
    
    % Raw
    figure('Name', 'Velocity (Raw)');
    imshow(velocityMap, []);
    colormap(cm);
    if any(velocityMap(:) > 0)
        clim([0 prctile(velocityMap(velocityMap>0), 99.5)]);
    end
    title('Velocity Raw');
    colorbar;
    add_scale_bar(params, size(velocityMap), lw, fs);
end

function createCombinedFigure(velocityMap, densityMap, params, lw, fs)
    figure('Name', 'Combined (Velocity × Density)');
    
    % Normalize velocity for color encoding
    if any(velocityMap(:) > 0)
        v_max = prctile(velocityMap(velocityMap > 0), 99.5);
    else
        v_max = 1;
    end
    v_norm = mat2gray(velocityMap, [0, v_max]);
    
    % Normalize density for brightness (alpha) encoding
    d_proc = densityMap .^ (1/3);
    if any(d_proc(:) > 0)
        d_max = prctile(d_proc(d_proc > 0), 99.0);
    else
        d_max = 1;
    end
    d_norm = mat2gray(d_proc, [0, d_max]) .^ 0.7;
    
    % Use imagesc (not imshow) so colormap + clim + colorbar all apply correctly.
    % Display velocity as the primary color-encoded channel.
    imagesc(v_norm, [0, 1]);
    colormap(jet(256));
    cb = colorbar;
    cb.Ticks = linspace(0, 1, 6);
    cb.TickLabels = arrayfun(@(v) sprintf('%.1f', v), linspace(0, v_max, 6), 'UniformOutput', false);
    ylabel(cb, 'Velocity (mm/s)', 'FontSize', max(fs, 7));
    clim([0, 1]);
    
    % Overlay density as alpha modulation (dark where no signal, bright where dense)
    hold on;
    hOv = image(repmat(d_norm, [1, 1, 3]));   % gray density image
    set(hOv, 'AlphaData', 1 - d_norm);       % transparent where dense → reveals velocity color
    
    % Mask background pixels completely
    bg_overlay = zeros([size(densityMap), 3]);   % black
    hBG = image(bg_overlay);
    set(hBG, 'AlphaData', double(densityMap == 0));
    hold off;
    
    axis image off;
    title('Combined: Velocity (color) × Density (brightness)');
    add_scale_bar(params, size(velocityMap), lw, fs);
end

% =========================================================================
% --- SESSION MANAGEMENT ---
% =========================================================================

function saveSession(fig)
    app = guidata(fig);
    
    [file, path] = uiputfile('*.mat', 'Save Session As');
    if isequal(file, 0)
        return;
    end
    
    try
        showProgress(fig, 'Saving session...', true);
        
        sessionData = app.sessionManager.createSessionData(app);
        save(fullfile(path, file), 'sessionData');
        
        setStatus(app, 'Session saved', 'green');
        guidata(fig, app);
        uialert(fig, 'Session saved successfully!', 'Success');
        
    catch ME
        errordlg(ME.message, 'Save Error');
        setStatus(app, 'Save failed', 'red');
    end
    
    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function loadSession(fig)
    app = guidata(fig);
    
    [file, path] = uigetfile('*.mat', 'Load Session');
    if isequal(file, 0)
        return;
    end
    
    try
        showProgress(fig, 'Loading session...', true);
        
        loaded = load(fullfile(path, file));
        if ~isfield(loaded, 'sessionData')
            error('Invalid session file.');
        end
        
        app = app.sessionManager.restoreSessionData(app, loaded.sessionData);
        
        % Update UI from loaded params
        populateGUIFromParams(app, app.data.params);
        
        % Update display
        app = manageGUIState(app, app.state.currentState);
        
        setStatus(app, 'Session loaded', 'green');
        guidata(fig, app);
        uialert(fig, 'Session loaded successfully!', 'Success');
        
    catch ME
        errordlg(ME.message, 'Load Error');
        setStatus(app, 'Load failed', 'red');
    end
    
    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

% =========================================================================
% --- UNDO/REDO SYSTEM ---
% =========================================================================

function saveParamState(fig, operation)
    app = guidata(fig);
    
    currentParams = updateParamsFromGUI(app);
    app.undoManager.saveState(currentParams, operation);
    
    % Update button states
    app.ui.btnUndo.Enable = app.undoManager.canUndo();
    app.ui.btnRedo.Enable = app.undoManager.canRedo();
    
    guidata(fig, app);
end

function performUndo(fig)
    app = guidata(fig);
    if app.undoManager.canUndo()
        prevParams = app.undoManager.undo();
        if isempty(prevParams), return; end  % reached beginning — nothing to apply
        app.data.params = prevParams;
        populateGUIFromParams(app, prevParams);
        app.ui.btnUndo.Enable = app.undoManager.canUndo();
        app.ui.btnRedo.Enable = app.undoManager.canRedo();
        setStatus(app, 'Undo complete', 'blue');
        guidata(fig, app);
    end
end

function performRedo(fig)
    app = guidata(fig);
    
    if app.undoManager.canRedo()
        nextParams = app.undoManager.redo();
        app.data.params = nextParams;
        populateGUIFromParams(app, nextParams);
        
        app.ui.btnUndo.Enable = app.undoManager.canUndo();
        app.ui.btnRedo.Enable = app.undoManager.canRedo();
        
        setStatus(app, 'Redo complete', 'blue');
        guidata(fig, app);
    end
end

% =========================================================================
% --- DISPLAY FUNCTIONS ---
% =========================================================================

function displayCurrentFrame(fig)
    app = guidata(fig);
    if ~isvalid(fig) || isempty(app)
        return;
    end
    
    % Use display manager
    app.displayManager.displayFrame(app);
end

% =========================================================================
% --- UTILITY FUNCTIONS ---
% =========================================================================

function setStatus(app, message, color)
    app.ui.StatusLabel.Text = message;
    app.ui.StatusLamp.Color = color;
    drawnow limitrate;
end

function showProgress(fig, message, indeterminate)
    % Accepts fig handle — reads/writes app via guidata to avoid stale struct bug
    if nargin < 3
        indeterminate = false;
    end
    app = guidata(fig);
    
    % Update status bar
    setStatus(app, message, 'red');
    
    % Skip progress dialog for rendering messages (uses status bar only)
    if contains(lower(message), 'generat') || contains(lower(message), 'map')
        drawnow;
        return;
    end
    
    % Create or update progress dialog
    if isempty(app.ui.ProgressBar) || ~isvalid(app.ui.ProgressBar)
        app.ui.ProgressBar = uiprogressdlg(app.fig, 'Title', 'Processing', ...
            'Message', message, 'Cancelable', 'off');
    end
    
    app.ui.ProgressBar.Message = message;
    app.ui.ProgressBar.Indeterminate = indeterminate;
    drawnow;
    
    % CRITICAL: save back so caller can get the handle
    guidata(fig, app);
end

function hideProgress(fig)
    % Accepts fig handle — reads fresh app from guidata
    if ~isvalid(fig), return; end
    app = guidata(fig);
    if ~isempty(app.ui.ProgressBar) && isvalid(app.ui.ProgressBar)
        try
            close(app.ui.ProgressBar);
            delete(app.ui.ProgressBar);
        catch
        end
        app.ui.ProgressBar = [];
        guidata(fig, app);
        drawnow;
    end
end

function startMemoryMonitor(fig)
    % Create timer for memory monitoring
    t = timer('ExecutionMode', 'fixedRate', 'Period', 2, ...
        'TimerFcn', @(~,~) updateMemoryDisplay(fig));
    start(t);
    
    % Store timer
    app = guidata(fig);
    app.state.memoryTimer = t;
    guidata(fig, app);
end

function updateMemoryDisplay(fig)
    if ~isvalid(fig), return; end
    app = guidata(fig);
    if isempty(app) || app.state.isProcessing, return; end
    
    if ispc()
        try
            memStats = memory();
            usedMB   = memStats.MemUsedMATLAB / 1024^2;
            app.ui.lblMemory.Text = sprintf('Memory: %.0f MB', usedMB);
        catch
            app.ui.lblMemory.Text = 'Memory: N/A';
        end
    else
        % memory() is only available on Windows — use process RSS on Linux/macOS
        try
            if ismac()
                [~, txt] = system('ps -o rss= -p ' + string(feature('getpid')));
                usedMB   = str2double(strtrim(txt)) / 1024;
            else  % Linux
                fid = fopen(sprintf('/proc/%d/status', feature('getpid')), 'r');
                txt = fread(fid, '*char')';
                fclose(fid);
                tok = regexp(txt, 'VmRSS:\s*(\d+)', 'tokens', 'once');
                usedMB = str2double(tok{1}) / 1024;
            end
            app.ui.lblMemory.Text = sprintf('Memory: %.0f MB', usedMB);
        catch
            app.ui.lblMemory.Text = 'Memory: N/A';
        end
    end
end

function cleanupGUI(fig)
    if ~isvalid(fig), return; end
    app = guidata(fig);
    if isempty(app), delete(fig); return; end
    
    % Stop and delete all timers
    timerFields = {'playbackTimer', 'displayTimer', 'memoryTimer'};
    for k = 1:numel(timerFields)
        fn = timerFields{k};
        if isfield(app.state, fn) && ~isempty(app.state.(fn)) && isvalid(app.state.(fn))
            stop(app.state.(fn));
            delete(app.state.(fn));
        end
    end
    
    if isfield(app, 'displayManager') && ~isempty(app.displayManager)
        app.displayManager.cleanup();
    end
    
    % Close progress dialog
    if isfield(app.ui, 'ProgressBar') && ~isempty(app.ui.ProgressBar) && isvalid(app.ui.ProgressBar)
        close(app.ui.ProgressBar);
    end
    
    fprintf('ULM GUI closed. Resources cleaned up.\n');
    delete(fig);
end

function app = manageGUIState(app, ns)
    app.state.currentState = ns;
    
    % Disable all panels
    app.ui.panel_filt.Enable = 'off';
    app.ui.panel_detect.Enable = 'off';
    app.ui.panel_loc.Enable = 'off';
    app.ui.panel_track.Enable = 'off';
    app.ui.panel_post.Enable = 'off';
    app.ui.panel_render.Enable = 'off';
    app.ui.PlayPauseButton.Enable = 'off';
    app.ui.FrameSlider.Enable = 'off';
    app.ui.FrameField.Enable = 'off';
    
    % Enable based on state
    if ns >= 0
        app.ui.panel_filt.Enable = 'on';
        app.ui.PlayPauseButton.Enable = 'on';
        app.ui.FrameSlider.Enable = 'on';
        app.ui.FrameField.Enable = 'on';
    end
    if ns >= 1, app.ui.panel_detect.Enable = 'on'; end
    if ns >= 2, app.ui.panel_loc.Enable = 'on'; end
    if ns >= 3, app.ui.panel_track.Enable = 'on'; end
    if ns >= 4, app.ui.panel_post.Enable = 'on'; end
    if ns >= 5
        app.ui.panel_render.Enable = 'on';
        app.ui.PlayPauseButton.Enable = 'off';
        app.ui.FrameSlider.Enable = 'off';
        app.ui.FrameField.Enable = 'off';
    end
end

function data = clearDownstreamData(data, lvl)
    if lvl < 5, data.tracks_final = []; end
    if lvl < 4, data.tracks_raw = []; end
    if lvl < 3, data.localizations = []; end
    if lvl < 2, data.candidateBubbles = []; end
    if lvl < 1
        data.filteredData = [];
        data.filteredMeanBG = [];
        data.U = [];
        data.S_diag = [];
        data.V = [];
    end
end

% --- Parameter Management ---
function params = updateParamsFromGUI(app)
    params = app.data.params;
    
    % === Filter / Crop ===
    params.filter.method             = app.ui.FilterMethodDropdown.Value;
    params.filter.svd_cutoff         = [app.ui.SVDCutoffStart.Value, app.ui.SVDCutoffEnd.Value];
    params.filter.enable_butterworth = app.ui.EnableButterworth.Value;
    params.filter.spatial_method     = app.ui.SpatialMethodDrop.Value;
    params.filter.spatial_kernel     = app.ui.SpatialKernelField.Value;
    params.filter.spatial_sigma1     = app.ui.SpatialSigma1Field.Value;
    params.filter.spatial_sigma2     = app.ui.SpatialSigma2Field.Value;
    
    try
        params.io.crop_box = str2double(app.ui.CropBoxField.Value);
    catch
        params.io.crop_box = [];
    end
    
    % Validate butter_cutoff with feedback (not silent)
    vals = str2double(strsplit(strtrim(app.ui.ButterCutoff.Value)));
    if numel(vals) == 2 && all(~isnan(vals)) && vals(1) < vals(2) && all(vals > 0)
        params.filter.butter_cutoff = vals;
    else
        app.ui.ButterCutoff.Value = mat2str(params.filter.butter_cutoff); % revert bad input
    end
    params.filter.butter_order = app.ui.ButterOrder.Value;
    
    % === Top Menu Bar (fundamental params) ===
    params.acq.framerate        = app.ui.TopFPSField.Value;
    params.track.pixel_X_size   = app.ui.TopPixelXField.Value;
    params.track.pixel_Z_size   = app.ui.TopPixelZField.Value;
    params.track.dt             = 1 / params.acq.framerate;
    
    % === Detection / Localization ===
    params.loc.DetectMethod           = app.ui.DetectMethodDropdown.Value;
    params.loc.NP_alpha0              = app.ui.NP_AlphaField.Value;
    params.loc.crosscor_threshold     = app.ui.NCC_ThreshField.Value;
    params.loc.detection_threshold    = app.ui.LocThreshField.Value;
    params.loc.max_bubbles_per_frame  = app.ui.LocMaxField.Value;
    params.loc.method                 = app.ui.LocMethodDropdown.Value;
    params.loc.enable_divergence_check  = app.ui.LocQCDivergence.Value;
    params.loc.enable_roi_maxima_check  = app.ui.LocQCRoiMaxima.Value;
    params.loc.qc_max_shift_factor    = app.ui.LocShiftFactor.Value;
    try
        params.loc.fwhm = str2num(app.ui.DetectFWHMField.Value);
    catch
    end
    params.loc.gauss_fit_box_radius = app.ui.GaussBoxRadius.Value;
    params.loc.min_r_squared        = app.ui.GaussMinRSquared.Value;
    
    % === Tracking ===
    params.track.method                    = app.ui.TrackMethodDropdown.Value;
    params.track.max_linking_distance      = app.ui.MaxDistField.Value;
    params.track.max_gap_closing_frames    = app.ui.GapFramesField.Value;
    params.track.min_track_length          = app.ui.MinLengthField.Value;
    params.track.kalman.motion_model       = app.ui.KalmanModelDrop.Value;
    params.track.kalman.process_noise      = app.ui.KalmanNoise.Value;
    params.track.kalman.assignment_method  = app.ui.AssignmentDrop.Value;
    params.track.use_advanced_cost_matrix  = app.ui.UseAdvancedCostCheckbox.Value;
    
    % === QC ===
    params.track.qc.enable_direction_constraint    = app.ui.TrackQCDirection.Value;
    params.track.qc.max_angle_change_deg           = app.ui.QCMaxAngle.Value;
    params.track.qc.enable_acceleration_constraint = app.ui.TrackQCAcceleration.Value;
    params.track.qc.acceleration_C_factor          = app.ui.QCAccelFactor.Value;
    params.track.qc.enable_vd_constraint           = app.ui.TrackQCVD.Value;
    params.track.qc.max_vd_ratio                   = app.ui.QCVDRatio.Value;
    
    % === Post-processing ===
    params.track.enable_postprocessing = app.ui.EnablePostProcessing.Value;
    params.track.smoothing_factor      = app.ui.SmoothField.Value;
    
    % === Rendering ===
    params.render.upsampling_factor  = app.ui.UpsamplingField.Value;
    params.render.method             = app.ui.RenderMethodDrop.Value;
    params.render.interpolation_step = 0.5;  % FIX B-03: no UI control; use constant
end

function populateGUIFromParams(app, p)
    % Filter / Crop
    if isfield(p, 'io') && isfield(p.io, 'crop_box') && ~isempty(p.io.crop_box)
        app.ui.CropBoxField.Value = mat2str(p.io.crop_box);
    else
        app.ui.CropBoxField.Value = '[]';
    end

    % Top Menu Bar parameters
    if isfield(p, 'acq') && isfield(p.acq, 'framerate')
        app.ui.TopFPSField.Value = p.acq.framerate;
    end
    if isfield(p, 'track')
        if isfield(p.track, 'pixel_X_size'), app.ui.TopPixelXField.Value = p.track.pixel_X_size; end
        if isfield(p.track, 'pixel_Z_size'), app.ui.TopPixelZField.Value = p.track.pixel_Z_size; end
    end
    
    app.ui.FilterMethodDropdown.Value = p.filter.method;
    app.ui.SVDCutoffStart.Value = p.filter.svd_cutoff(1);
    app.ui.SVDCutoffEnd.Value = p.filter.svd_cutoff(2);
    app.ui.EnableButterworth.Value = p.filter.enable_butterworth;
    app.ui.ButterCutoff.Value = mat2str(p.filter.butter_cutoff);
    app.ui.ButterOrder.Value = p.filter.butter_order;

    if isfield(p.filter, 'spatial_method')
        app.ui.SpatialMethodDrop.Value = p.filter.spatial_method;
        app.ui.SpatialKernelField.Value = p.filter.spatial_kernel;
        app.ui.SpatialSigma1Field.Value = p.filter.spatial_sigma1;
        app.ui.SpatialSigma2Field.Value = p.filter.spatial_sigma2;
        updateSpatialOptions(app.fig, app);
    end
    
    % Detection
    app.ui.LocThreshField.Value = p.loc.detection_threshold;
    app.ui.LocThreshSlider.Value = p.loc.detection_threshold;
    app.ui.LocMaxField.Value = p.loc.max_bubbles_per_frame;
    app.ui.LocMaxSlider.Value = p.loc.max_bubbles_per_frame;
    if isfield(p.loc, 'DetectMethod')
        app.ui.DetectMethodDropdown.Value = p.loc.DetectMethod;
        updateDetectionOptions(app.fig);
    end
    if isfield(p.loc, 'NP_alpha0')
        app.ui.NP_AlphaField.Value = p.loc.NP_alpha0;
    end
    if isfield(p.loc, 'crosscor_threshold')
        app.ui.NCC_ThreshField.Value = p.loc.crosscor_threshold;
    end
    
    % Localization
    app.ui.LocMethodDropdown.Value = p.loc.method;
    app.ui.LocQCDivergence.Value = p.loc.enable_divergence_check;
    app.ui.LocQCRoiMaxima.Value = p.loc.enable_roi_maxima_check;
    app.ui.LocShiftFactor.Value = p.loc.qc_max_shift_factor;
    app.ui.LocFWHM.Value = mat2str(p.loc.fwhm);
    app.ui.DetectFWHMField.Value = mat2str(p.loc.fwhm);
    app.ui.GaussBoxRadius.Value = p.loc.gauss_fit_box_radius;
    app.ui.GaussMinRSquared.Value = p.loc.min_r_squared;
    
    % Tracking
    app.ui.TrackMethodDropdown.Value = p.track.method;
    app.ui.MaxDistField.Value = p.track.max_linking_distance;
    app.ui.MaxDistSlider.Value = p.track.max_linking_distance;
    app.ui.GapFramesField.Value = p.track.max_gap_closing_frames;
    app.ui.GapFramesSlider.Value = p.track.max_gap_closing_frames;
    app.ui.MinLengthField.Value = p.track.min_track_length;
    app.ui.MinLengthSlider.Value = p.track.min_track_length;
    app.ui.KalmanModelDrop.Value = p.track.kalman.motion_model;
    app.ui.KalmanNoise.Value = p.track.kalman.process_noise;
    app.ui.AssignmentDrop.Value = p.track.kalman.assignment_method;
    app.ui.UseAdvancedCostCheckbox.Value = p.track.use_advanced_cost_matrix;
    
    % QC
    app.ui.TrackQCDirection.Value = p.track.qc.enable_direction_constraint;
    app.ui.QCMaxAngle.Value = p.track.qc.max_angle_change_deg;
    app.ui.TrackQCAcceleration.Value = p.track.qc.enable_acceleration_constraint;
    app.ui.QCAccelFactor.Value = p.track.qc.acceleration_C_factor;
    app.ui.TrackQCVD.Value = p.track.qc.enable_vd_constraint;
    app.ui.QCVDRatio.Value = p.track.qc.max_vd_ratio;
    
    % Post-processing
    app.ui.EnablePostProcessing.Value = p.track.enable_postprocessing;
    app.ui.SmoothField.Value = p.track.smoothing_factor;
    app.ui.SmoothSlider.Value = p.track.smoothing_factor;
    app.ui.DisplayMinLengthField.Value = p.track.min_track_length;
    app.ui.DisplayMinLengthSlider.Value = p.track.min_track_length;
    
    % Rendering
    app.ui.UpsamplingField.Value = p.render.upsampling_factor;
    app.ui.RenderMethodDrop.Value = p.render.method;

    updateFilterOptions(app.fig);
    updateLocalizationOptions(app.fig);
    updateDetectionOptions(app.fig);
    updateTrackingOptions(app.fig);
end

function resetAllParams(fig)
    app = guidata(fig);
    try
        app.data.params = setDefaultParams(true, struct('apply_roi_mask', false));
    catch
        app.data.params = createFallbackParams();
        uialert(fig, 'setDefaultParams not found — reset to built-in defaults.', 'Warning', 'Icon', 'warning');
    end
    populateGUIFromParams(app, app.data.params);
    guidata(fig, app);
end

% --- UI Callbacks ---
function syncSliderField(slider, field)
    slider.ValueChangedFcn = @(o,e) set(field, 'Value', e.Value);
    field.ValueChangedFcn = @(o,e) set(slider, 'Value', max(min(e.Value, slider.Limits(2)), slider.Limits(1)));
end

function onFrameSliderChanged(fig, e)
    app = guidata(fig);
    app.state.currentFrame = round(e.Value);
    app.ui.FrameField.Value = app.state.currentFrame;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function onFrameFieldChanged(fig, e)
    app = guidata(fig);
    app.state.currentFrame = round(e.Value);
    app.ui.FrameSlider.Value = app.state.currentFrame;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function onDisplayFilterChanged(fig, e)
    app = guidata(fig);
    newVal = round(e.Value);
    newVal = max(2, min(20, newVal));  % clamp to slider limits
    
    % Cross-sync: update whichever control was NOT changed
    if isa(e.Source, 'matlab.ui.control.Slider')
        app.ui.DisplayMinLengthField.Value  = newVal;
        app.ui.DisplayMinLengthSlider.Value = newVal;
    else
        app.ui.DisplayMinLengthSlider.Value = newVal;
        app.ui.DisplayMinLengthField.Value  = newVal;
    end
    
    guidata(fig, app);
    if app.state.currentState >= 5
        displayCurrentFrame(fig);
    end
end

function togglePlayback(fig)
    app = guidata(fig);
    if app.state.currentState == 5, return; end
    
    % isempty() misses deleted-but-non-empty handles. Use isvalid().
    if isempty(app.state.playbackTimer) || ~isvalid(app.state.playbackTimer)
        app.state.playbackTimer = timer('ExecutionMode', 'fixedRate', ...
            'Period', 1 / ULM_Constants.PLAYBACK_FPS, ...
            'TimerFcn', @(~,~) timerCallback(fig));
    end
    
    if app.state.isPlaying
        stop(app.state.playbackTimer);
        app.state.isPlaying          = false;
        app.ui.PlayPauseButton.Text  = 'Play';
    else
        start(app.state.playbackTimer);
        app.state.isPlaying          = true;
        app.ui.PlayPauseButton.Text  = 'Pause';
    end
    guidata(fig, app);
end

function timerCallback(fig)
    if ~isvalid(fig), return; end
    app = guidata(fig);
    if isempty(app) || app.state.isProcessing, return; end  % FIX B-09
    if ~app.state.isPlaying, return; end
    
    app.state.currentFrame = mod(app.state.currentFrame, app.state.maxFrame) + 1;
    app.ui.FrameSlider.Value = app.state.currentFrame;
    app.ui.FrameField.Value  = app.state.currentFrame;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function updateFilterOptions(fig)
    app = guidata(fig);
    s = app.ui.FilterMethodDropdown.Value;
    app.ui.p_svd.Visible = strcmp(s, 'svd_filter');
    app.ui.p_dcc.Visible = strcmp(s, 'dcc_svd');
end

function updateDetectionOptions(fig)
    % Shows/hides method-specific parameter fields based on selected
    % detection method. Called when DetectMethodDropdown changes.
    app = guidata(fig);
    method = app.ui.DetectMethodDropdown.Value;
    
    % Hide all method-specific controls first
    app.ui.lbl_NP_alpha.Visible  = 'off';
    app.ui.NP_AlphaField.Visible = 'off';
    app.ui.lbl_NCC_thresh.Visible  = 'off';
    app.ui.NCC_ThreshField.Visible = 'off';
    
    switch upper(method)
        case 'NP'
            app.ui.lbl_NP_alpha.Visible  = 'on';
            app.ui.NP_AlphaField.Visible = 'on';
        case 'NCC'
            app.ui.lbl_NCC_thresh.Visible  = 'on';
            app.ui.NCC_ThreshField.Visible = 'on';
    end
    
    guidata(fig, app);
end

function updateLocalizationOptions(fig)
    app = guidata(fig);
    s = app.ui.LocMethodDropdown.Value;
    
    % Show appropriate QC panel based on method
    app.ui.p_loc_qc_radial.Visible = 'on';  % Always visible (common params)
    
    if strcmp(s, 'gaussian_fit')
        app.ui.p_loc_qc_gauss.Visible = 'on';
    else
        app.ui.p_loc_qc_gauss.Visible = 'off';
    end
end

function updateTrackingOptions(fig)
    app = guidata(fig);
    trackMethod = app.ui.TrackMethodDropdown.Value;
    
    isKalmanStandard = strcmpi(trackMethod, 'Kalman');
    isKalmanAdv = strcmpi(trackMethod, 'Kalman_Advanced');
    isAnyKalman = isKalmanStandard || isKalmanAdv;
    
    app.ui.KalmanModelDrop.Enable = isAnyKalman;
    
    % Standard process noise only matters for standard Kalman
    app.ui.KalmanNoise.Enable = isKalmanStandard; 
    
    % Enable HK button only for Advanced Kalman
    if isfield(app.ui, 'BtnConfigHK')
        app.ui.BtnConfigHK.Enable = isKalmanAdv;
    end
    
    isAdvCost = app.ui.UseAdvancedCostCheckbox.Value;
    app.ui.BtnConfigCostMatrix.Enable = isAdvCost;
end

function onDCCSliderChanged(fig)
    app = guidata(fig);
    
    % Only act if SVD has already been computed (U, S, V cached)
    % and DCC clustering has been run (indices exist).
    if app.state.currentState < 1 || isempty(app.data.U) || ...
       isempty(app.data.tissue_indices)
        return;
    end
    
    try
        % Reconstruct using only the slider-selected component sub-ranges.
        % No SVD recomputation needed — uses the cached U, S, V directly.
        filteredMatrix = reconstructDCCImage(app);
        app.data.filteredData = abs(filteredMatrix);
        
        % Recompute color limits for display
        abs_data = app.data.filteredData(:);
        app.data.filteredClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
        if app.data.filteredClim(1) == app.data.filteredClim(2)
            app.data.filteredClim(2) = app.data.filteredClim(1) + 1;
        end
        
        % Recompute background for the Detect tab
        app.data.filteredMeanBG = mean(app.data.filteredData, 3) .^ 0.5;
        app.data.baseVesselMap  = [];   % force recalculation on next tab switch
        
        guidata(fig, app);
        displayCurrentFrame(fig);
        
    catch ME
        setStatus(app, sprintf('DCC slider error: %s', ME.message), 'red');
    end
end

function onTabChanged(fig)
    app = guidata(fig);
    currentTab = app.ui.tabGroup.SelectedTab.Title;
    
    if contains(currentTab, 'Detect')
        prepareROITab(fig);
    else
        app.state.isROIPreview = false;
        guidata(fig, app);
        displayCurrentFrame(fig);
    end
end

% --- ROI Functions (from original code, adapted) ---
function prepareROITab(fig)
    app = guidata(fig);
    if isempty(app.data.filteredData), return; end
    
    if isempty(app.data.baseVesselMap)
        showProgress(fig, 'Calculating vessel map...', true);  % FIX B-01 sig
        
        if ~isreal(app.data.filteredData)
            absData = abs(app.data.filteredData);
        else
            absData = app.data.filteredData;
        end
        
        vMap = mean(absData, 3);
        vMap = vMap .^ 0.5;
        mx = max(vMap(:));
        if mx > 0
            vMap = vMap / mx;
        end
        
        app.data.baseVesselMap = vMap;
        app.data.vesselMap     = vMap;
        
        % threshold=1 means "exclude everything".
        % Use 0 so the full map is shown on first open.
        app.ui.ROIContrastSlider.Value = 1;   % gamma=1: no stretch (correct)
        app.ui.ROIThreshSlider.Value   = 0;   % threshold=0: show all (fixed)
        app.ui.ROIThreshField.Value    = 0;   % threshold=0: show all (fixed)
        
        hideProgress(fig);  % FIX B-01 sig
    end
    
    app.ui.EnhanceMethodDrop.Value    = 'None';
    app.ui.EnhanceAmountSlider.Value  = 0.5;
    app.ui.axHist.YScale              = 'log';
    app.state.isROIPreview            = true;
    app.ui.chkPreviewROI.Value        = 1;
    
    guidata(fig, app);
    updateHistogram(app);
    displayCurrentFrame(fig);
end

function applyVesselEnhancement(fig, sliderVal)
    app = guidata(fig);
    if isempty(app.data.baseVesselMap)
        return;
    end
    
    if nargin < 2
        val = app.ui.EnhanceAmountSlider.Value;
    else
        val = sliderVal;
    end
    
    method = app.ui.EnhanceMethodDrop.Value;
    baseImg = app.data.baseVesselMap;
    
    switch method
        case 'None'
            procImg = baseImg;
        case 'CLAHE (Local Contrast)'
            clipLim = 0.001 + (val * 0.04);
            procImg = adapthisteq(baseImg, 'ClipLimit', clipLim, 'Distribution', 'rayleigh');
        case 'Top-Hat (Vesselness)'
            radius = 1 + round(val * 10);
            se = strel('disk', radius);
            procImg = imtophat(baseImg, se);
            mx = max(procImg(:));
            if mx > 0
                procImg = procImg / mx;
            end
        case 'Sharpen'
            amount = val * 2;
            procImg = imsharpen(baseImg, 'Radius', 1, 'Amount', amount);
    end
    
    currentGamma = app.ui.ROIContrastSlider.Value;
    procImg = procImg .^ currentGamma;
    mn = min(procImg(:));
    mx = max(procImg(:));
    if mx > mn
        procImg = (procImg - mn) / (mx - mn);
    else
        procImg = zeros(size(procImg), 'like', procImg);
    end
    
    app.data.vesselMap = procImg;
    currThresh = app.ui.ROIThreshSlider.Value;
    app.data.mask = procImg >= currThresh;
    
    app.state.isROIPreview = true;
    app.ui.chkPreviewROI.Value = 1;
    guidata(fig, app);
    
    updateHistogram(app);
    displayCurrentFrame(fig);
end

function onContrastChange(fig, gammaVal)
    app = guidata(fig);
    gammaVal = max(0.1, min(3, gammaVal));
    app.ui.ROIContrastSlider.Value = gammaVal;
    app.ui.ROIContrastField.Value = gammaVal;
    guidata(fig, app);
    applyVesselEnhancement(fig);
end

function updateHistogram(app)
    vMap = app.data.vesselMap;
    validPixels = vMap(vMap > 0);
    
    cla(app.ui.axHist);
    histogram(app.ui.axHist, validPixels, 100, 'EdgeColor', 'none', 'FaceColor', [0.2 0.4 0.8]);
    app.ui.axHist.YScale = 'log';
    app.ui.axHist.XLim = [0 1];
    
    xline(app.ui.axHist, app.ui.ROIThreshSlider.Value, 'r-', 'LineWidth', 2, 'Tag', 'threshLine');
end

function onROIChange(fig, val)
    app = guidata(fig);
    val = max(0, min(1, val));
    app.ui.ROIThreshSlider.Value = val;
    app.ui.ROIThreshField.Value = val;
    
    l = findobj(app.ui.axHist, 'Tag', 'threshLine');
    if ~isempty(l)
        l.Value = val;
    end
    
    if ~isempty(app.data.vesselMap)
        app.data.mask = app.data.vesselMap >= val;
    end
    
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function onROIPreviewToggle(fig)
    app = guidata(fig);
    app.state.isROIPreview = app.ui.chkPreviewROI.Value;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function saveCreatedMask(fig)
    app = guidata(fig);
    
    if isempty(app.data.mask)
        uialert(fig, 'No mask created yet.', 'Save Error');
        return;
    end
    
    [file, path] = uiputfile('vesselMask.mat', 'Save Mask As...');
    if isequal(file, 0)
        return;
    end
    
    try
        vesselMask = app.data.mask;
        save(fullfile(path, file), 'vesselMask');
        uialert(fig, sprintf('Mask saved as ''vesselMask'' in:\n%s', file), 'Success');
    catch ME
        uialert(fig, ME.message, 'Error Saving File');
    end
end

function loadMask(fig)
    app = guidata(fig);
    
    if isempty(app.data.rawData)
        uialert(fig, 'Please load Raw Data first.', 'Sequence Error');
        return;
    end
    
    defaultPath = app.data.params.io.data_folder;
    [f, p] = uigetfile(fullfile(defaultPath, '*.mat'), 'Select Mask File');
    if isequal(f, 0)
        return;
    end
    
    try
        loadedStruct = load(fullfile(p, f));
        fieldNames = fieldnames(loadedStruct);
        maskData = [];
        
        if isfield(loadedStruct, 'mask')
            maskData = loadedStruct.mask;
        elseif length(fieldNames) == 1
            maskData = loadedStruct.(fieldNames{1});
        else
            for i = 1:length(fieldNames)
                tempData = loadedStruct.(fieldNames{i});
                if ismatrix(tempData) && (islogical(tempData) || isnumeric(tempData))
                    maskData = tempData;
                    break;
                end
            end
        end
        
        if isempty(maskData)
            uialert(fig, 'No valid mask variable found.', 'Load Error');
            return;
        end
        
        [h, w, ~] = size(app.data.rawData);
        [mh, mw] = size(maskData);
        
        if mh ~= h || mw ~= w
            uialert(fig, sprintf('Mask size (%d rows × %d cols) does not match data (%d rows × %d cols).\nMaybe this mask is from a different crop or experiment?', ...
    mh, mw, h, w), 'Dimension Mismatch');
            return;
        end
        
        app.data.mask = logical(maskData);
        app.ui.maskStatusLabel.Text = 'Status: Loaded';
        app.data = clearDownstreamData(app.data, 1);
        
        guidata(fig, app);
        displayCurrentFrame(fig);
        
    catch ME
        uialert(fig, ME.message, 'File Load Error');
    end
end

function runCreateMask(fig)
    app = guidata(fig);
    if isempty(app.data.filteredData)
        return;
    end
    
    showProgress(fig, 'Drawing mask...', true);
    f = figure;
    imagesc(mean(abs(app.data.filteredData), 3));
    colormap gray;
    axis image;
    roi = drawfreehand;
    wait(roi);
    
    app.data.mask = createMask(roi);
    delete(f);
    app.ui.maskStatusLabel.Text = 'Created';
    app.data = clearDownstreamData(app.data, 1);

    guidata(fig, app);

    hideProgress(fig);
    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function resetMask(fig)
    app = guidata(fig);
    app.data.mask = [];
    app.ui.maskStatusLabel.Text = 'None';
    app.data = clearDownstreamData(app.data, 1);
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function openCostMatrixGUI(mainFig)
    app = guidata(mainFig);
    p = app.data.params.track.kalman;
    
    f = uifigure('Name', 'Advanced Cost Matrix Configuration', ...
        'Position', [300 300 400 300], 'WindowStyle', 'modal');
    gl = uigridlayout(f, [5, 3]);
    gl.RowHeight = {'fit', 'fit', 'fit', '1x', 'fit'};
    gl.ColumnWidth = {'fit', '1x', '1x'};
    
    chkDir = uicheckbox(gl, 'Text', 'Enable Direction Penalty', 'Value', p.use_direction);
    uilabel(gl, 'Text', 'Weight:');
    efDirW = uieditfield(gl, 'numeric', 'Value', p.direction_penalty_weight);
    
    chkAng = uicheckbox(gl, 'Text', 'Enable Angle Penalty', 'Value', p.use_angle);
    uilabel(gl, 'Text', 'Slope:');
    efAngS = uieditfield(gl, 'numeric', 'Value', p.angle_penalty_slope);
    
    chkBri = uicheckbox(gl, 'Text', 'Enable Brightness Penalty', 'Value', p.use_brightness);
    uilabel(gl, 'Text', 'Weight:');
    efBriW = uieditfield(gl, 'numeric', 'Value', p.brightness_penalty_weight);
    
    uilabel(gl, 'Text', 'Dir History Pts:');
    efHist = uieditfield(gl, 'numeric', 'Value', p.direction_history_points);
    
    btnSave = uibutton(gl, 'Text', 'Save & Close', ...
        'ButtonPushedFcn', @(s,e) saveAndClose());
    btnSave.Layout.Row = 5;
    btnSave.Layout.Column = [2 3];
    
    function saveAndClose()
        app.data.params.track.kalman.use_direction = chkDir.Value;
        app.data.params.track.kalman.direction_penalty_weight = efDirW.Value;
        app.data.params.track.kalman.use_angle = chkAng.Value;
        app.data.params.track.kalman.angle_penalty_slope = efAngS.Value;
        app.data.params.track.kalman.use_brightness = chkBri.Value;
        app.data.params.track.kalman.brightness_penalty_weight = efBriW.Value;
        app.data.params.track.kalman.direction_history_points = efHist.Value;
        guidata(mainFig, app);
        close(f);
    end
end

function openHKConfigGUI(mainFig)
    app = guidata(mainFig);
    p = app.data.params.track.kalman;
    
    f = uifigure('Name', 'Hierarchical Kalman (HK) Settings', ...
        'Position', [300 200 450 540], 'WindowStyle', 'modal');
    gl = uigridlayout(f, [11, 2]);
    gl.RowHeight = {'fit','fit','fit','fit','fit','fit','fit','fit','fit','1x','fit'};
    gl.ColumnWidth = {'1x', '1x'};
    
    uilabel(gl, 'Text', 'Alpha (Process Noise multiplier):', 'FontWeight', 'bold');
    efAlpha = uieditfield(gl, 'numeric', 'Value', p.hk_alpha, 'ValueDisplayFormat', '%.4f');
    
    uilabel(gl, 'Text', 'Beta (Measurement Noise multiplier):', 'FontWeight', 'bold');
    efBeta = uieditfield(gl, 'numeric', 'Value', p.hk_beta, 'ValueDisplayFormat', '%.4f');
    
    chkFB = uicheckbox(gl, 'Text', 'Enable Forward-Backward Pass', 'Value', p.hk_forward_backward);
    chkFB.Layout.Column = [1 2];
    
    uilabel(gl, 'Text', 'Max Velocity (v_max) [mm/s]:');
    efVmax = uieditfield(gl, 'numeric', 'Value', p.hk_v_max);
    
    uilabel(gl, 'Text', 'Number of Levels:');
    efLevels = uieditfield(gl, 'numeric', 'Value', p.hk_num_levels, 'RoundFractionalValues', 'on');
    
    uilabel(gl, 'Text', 'Spacing Power (1.0 = linear):');
    efPower = uieditfield(gl, 'numeric', 'Value', p.hk_spacing_power);
    
    chkOverlap = uicheckbox(gl, 'Text', 'Enable Level Overlap', 'Value', p.hk_enable_overlap);
    
    uilabel(gl, 'Text', 'Overlap Band Width [mm/s]:');
    efOverlap = uieditfield(gl, 'numeric', 'Value', p.hk_overlap_mm_s);
    
    % Toggle overlap field based on checkbox
    efOverlap.Enable = chkOverlap.Value;
    chkOverlap.ValueChangedFcn = @(s,e) onOverlapToggle(e);
    
    % --- Velocity Ladder Preview ---
    lblLadder = uilabel(gl, 'Text', 'Velocity Ladder Preview:', 'FontWeight', 'bold');
    lblLadder.Layout.Column = [1 2];
    
    lbLadder = uilistbox(gl, 'Items', {}, 'Enable', 'off');
    lbLadder.Layout.Column = [1 2];
    
    btnSave = uibutton(gl, 'Text', 'Save & Close', ...
        'ButtonPushedFcn', @(s,e) saveHKAndClose(), ...
        'BackgroundColor', [0.8 1 0.8], 'FontWeight', 'bold');
    btnSave.Layout.Row = 11;
    btnSave.Layout.Column = [1 2];
    
    % Initial ladder render
    refreshLadder();
    
    % Hook live-update callbacks
    efVmax.ValueChangedFcn   = @(s,e) refreshLadder();
    efLevels.ValueChangedFcn = @(s,e) refreshLadder();
    efPower.ValueChangedFcn  = @(s,e) refreshLadder();
    efOverlap.ValueChangedFcn = @(s,e) refreshLadder();
    
    % ---- Nested helpers ----
    function onOverlapToggle(e)
        efOverlap.Enable = e.Value;
        refreshLadder();
    end
    
    function refreshLadder()
        v_max        = efVmax.Value;
        n            = max(1, round(efLevels.Value));
        power        = max(0.01, efPower.Value);
        overlap_mm_s = 0;
        if chkOverlap.Value
            overlap_mm_s = efOverlap.Value;
        end
        
        norm_boundaries = linspace(0, 1, n + 1) .^ power;
        boundaries      = v_max * norm_boundaries;
        
        lower_bounds = boundaries(1:end-1);
        upper_bounds = boundaries(2:end);
        
        if overlap_mm_s > 0
            half_overlap = overlap_mm_s / 2;
            lower_bounds(2:end)   = lower_bounds(2:end)   - half_overlap;
            upper_bounds(1:end-1) = upper_bounds(1:end-1) + half_overlap;
            lower_bounds = max(lower_bounds, 0);
            upper_bounds = min(upper_bounds, v_max);
        end
        
        levels_matrix = round([lower_bounds', upper_bounds'], 1);
        
        items = cell(1, n);
        for k = 1:n
            items{k} = sprintf('Level %d:  %.1f  –  %.1f  mm/s', k, levels_matrix(k,1), levels_matrix(k,2));
        end
        lbLadder.Items = items;
    end

    function saveHKAndClose()
        % 1. Save the scalar UI values
        app.data.params.track.kalman.hk_alpha            = efAlpha.Value;
        app.data.params.track.kalman.hk_beta             = efBeta.Value;
        app.data.params.track.kalman.hk_forward_backward = chkFB.Value;
        app.data.params.track.kalman.hk_v_max            = efVmax.Value;
        app.data.params.track.kalman.hk_num_levels       = efLevels.Value;
        app.data.params.track.kalman.hk_spacing_power    = efPower.Value;
        app.data.params.track.kalman.hk_enable_overlap   = chkOverlap.Value;
        app.data.params.track.kalman.hk_overlap_mm_s     = efOverlap.Value;
        
        v_max        = efVmax.Value;
        n            = max(1, round(efLevels.Value));
        power        = max(0.01, efPower.Value);
        overlap_mm_s = 0;
        if chkOverlap.Value
            overlap_mm_s = efOverlap.Value;
        end
        
        norm_boundaries = linspace(0, 1, n + 1) .^ power;
        boundaries      = v_max * norm_boundaries;
        
        lower_bounds = boundaries(1:end-1);
        upper_bounds = boundaries(2:end);
        
        if overlap_mm_s > 0
            half_overlap = overlap_mm_s / 2;
            lower_bounds(2:end)   = lower_bounds(2:end)   - half_overlap;
            upper_bounds(1:end-1) = upper_bounds(1:end-1) + half_overlap;
            lower_bounds = max(lower_bounds, 0);
            upper_bounds = min(upper_bounds, v_max);
        end
        
        levels_matrix = round([lower_bounds', upper_bounds'], 1);
        
        % Convert matrix to cell array of [min max] vectors just like setDefaultParams does
        app.data.params.track.kalman.velocity_levels = mat2cell(levels_matrix, ones(n, 1), 2);
        % ---------------------------------------------------------

        % 3. Save undo state
        app.undoManager.saveState(app.data.params, 'Update HK Params');
        app.ui.btnUndo.Enable = app.undoManager.canUndo();
        app.ui.btnRedo.Enable = app.undoManager.canRedo();
        
        % 4. Push back to main app and close
        guidata(mainFig, app);
        close(f);
    end
end

function f = reconstructDCCImage(app)
    % Selects a percentage-based sub-range from each cluster's indices,
    % then reconstructs the signal using an efficient low-rank product -
    % avoiding the full T x T diag(S) matrix.
    
    if isempty(app.data.tissue_indices) && ...
       isempty(app.data.blood_indices)  && ...
       isempty(app.data.noise_indices)
        warning('reconstructDCCImage: No DCC indices found. Returning raw data.');
        f = app.data.rawData;
        return;
    end
    
    % Helper: select a percentage sub-range from a sorted index list.
    % s and e are percentages [0-100] from the GUI sliders.
    calc_idx = @(list, s, e) list( ...
        max(1,          floor(numel(list) * (s / 100)) + 1) : ...
        min(numel(list), round(numel(list) * (e / 100))));
    
    t_idx = calc_idx(app.data.tissue_indices, ...
        app.ui.DCCTissueStart.Value, app.ui.DCCTissueEnd.Value);
    b_idx = calc_idx(app.data.blood_indices,  ...
        app.ui.DCCBloodStart.Value,  app.ui.DCCBloodEnd.Value);
    n_idx = calc_idx(app.data.noise_indices,  ...
        app.ui.DCCNoiseStart.Value,  app.ui.DCCNoiseEnd.Value);
    
    active_components = unique([t_idx(:); b_idx(:); n_idx(:)]);
    
    fprintf('[DCC Reconstruction] Using %d active components: %s\n', ...
        numel(active_components), formatIndicesToRanges(active_components));
    
    if isempty(active_components)
        warning('reconstructDCCImage: No active components selected. Returning zeros.');
        f = zeros(size(app.data.rawData), 'like', app.data.rawData);
        return;
    end
    
    % Efficient low-rank reconstruction - avoids building a full T x T
    % diagonal matrix. Equivalent to sum_k( S_k * U_k * V_k' ).
    U_k = app.data.U(:, active_components);          % [H*W x K]
    S_k = app.data.S_diag(active_components);        % [K x 1]
    V_k = app.data.V(:, active_components);          % [T x K]
    
    filteredMatrix = (U_k .* S_k') * V_k';           % [H*W x T]
    f = reshape(filteredMatrix, size(app.data.rawData));
end

function str = formatIndicesToRanges(inds)
    if isempty(inds)
        str = 'None';
        return;
    end
    inds = sort(inds(:)');
    d = diff(inds);
    break_pts = [0, find(d > 1), length(inds)];
    parts = {};
    for i = 1:length(break_pts)-1
        startVal = inds(break_pts(i)+1);
        endVal = inds(break_pts(i+1));
        if startVal == endVal
            parts{end+1} = sprintf('%d', startVal);
        else
            parts{end+1} = sprintf('%d-%d', startVal, endVal);
        end
    end
    str = strjoin(parts, ', ');
end

function add_scale_bar(params, mapSize, linewidth, fontsize)
    if ~isfield(params, 'expParams') || ~isfield(params.expParams, 'fovX')
        return;
    end
    
    hold on;
    try
        px_mm = params.expParams.fovX / mapSize(2);
        if px_mm <= 0
            return;
        end
        len_px = 1 / px_mm;
        x = 0.05 * mapSize(2);
        y = 0.95 * mapSize(1);
        plot([x, x+len_px], [y, y], 'w-', 'LineWidth', linewidth);
        text(x+len_px/2, y-4*linewidth, '1 mm', 'Color', 'w', ...
            'FontSize', fontsize, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom');
    catch
    end
    hold off;
end

function outData = applySpatialFilter(inData, filterParams)
    % Applies frame-by-frame 2D spatial filtering to 3D matrix
    [H, W, T] = size(inData);
    outData = zeros(size(inData), 'like', inData);
    
    method = filterParams.spatial_method;
    kSize = max(1, round(filterParams.spatial_kernel));
    sig1 = filterParams.spatial_sigma1;
    sig2 = filterParams.spatial_sigma2;
    
    % Ensure kernel is odd for Gaussian and Median
    if (strcmp(method, 'Gaussian') || strcmp(method, 'Median')) && mod(kSize, 2) == 0
        kSize = kSize + 1; 
    end
    
    % Pre-compute structuring element for Top-Hat to save time
    if strcmp(method, 'Top-Hat')
        se = strel('disk', kSize);
    end
    
    for t = 1:T
        frame = inData(:,:,t);
        
        switch method
            case 'Gaussian'
                outData(:,:,t) = imgaussfilt(frame, sig1, 'FilterSize', kSize);
            case 'Median'
                outData(:,:,t) = medfilt2(frame, [kSize kSize]);
            case 'DoG'
                g1 = imgaussfilt(frame, sig1);
                g2 = imgaussfilt(frame, sig2);
                outData(:,:,t) = max(0, g1 - g2); % Clamp negative values
            case 'Top-Hat'
                outData(:,:,t) = imtophat(frame, se);
        end
    end
end

function params = createFallbackParams()
    params = struct();
    
    params.acq.framerate        = 200;
    params.track.pixel_X_size   = 0.05;
    params.track.pixel_Z_size   = 0.05;
    params.track.dt             = 1 / params.acq.framerate;  % FIX B-02
    
    params.io.crop_box          = [];
    params.io.data_folder       = pwd;
    
    params.filter.method              = 'svd_filter';
    params.filter.svd_cutoff          = [1, 100];
    params.filter.enable_butterworth  = false;
    params.filter.butter_cutoff       = [10, 100];
    params.filter.butter_order        = 4;
    params.filter.spatial_method      = 'Gaussian';
    params.filter.spatial_kernel      = 3;
    params.filter.spatial_sigma1      = 1.0;
    params.filter.spatial_sigma2      = 2.0;
    
    params.loc.method                 = 'radial';
    params.loc.detection_threshold    = 0.3;
    params.loc.max_bubbles_per_frame  = 2000;
    params.loc.enable_divergence_check  = true;
    params.loc.enable_roi_maxima_check  = false;
    params.loc.qc_max_shift_factor    = 2.0;
    params.loc.fwhm                   = [1.5, 1.5];
    params.loc.gauss_fit_box_radius   = 3;
    params.loc.min_r_squared          = 0.5;
    params.loc.DetectMethod           = 'Intensity';
    params.loc.NP_alpha0              = 1e-4;
    params.loc.crosscor_threshold     = 0.6;
    params.loc.MB_image               = [];
    params.loc.psf_size               = [7, 7];  % FIX B-05: NCC PSF template size
    
    params.track.method                    = 'Kalman';
    params.track.max_linking_distance      = 2.0;
    params.track.max_gap_closing_frames    = 2;
    params.track.min_track_length          = 3;
    params.track.use_advanced_cost_matrix  = false;
    
    params.track.kalman.motion_model              = 'ConstantVelocity';
    params.track.kalman.process_noise             = 10;
    params.track.kalman.assignment_method         = 'hungarian';
    params.track.kalman.use_direction             = true;
    params.track.kalman.direction_penalty_weight  = 1;
    params.track.kalman.use_angle                 = true;
    params.track.kalman.angle_penalty_slope       = 1;
    params.track.kalman.use_brightness            = true;
    params.track.kalman.brightness_penalty_weight = 1;
    params.track.kalman.direction_history_points  = 3;
    params.track.kalman.hk_alpha                  = 0.01;
    params.track.kalman.hk_beta                   = 0.25;
    params.track.kalman.hk_forward_backward       = true;
    params.track.kalman.hk_v_max                  = 30;
    params.track.kalman.hk_num_levels             = 5;
    params.track.kalman.hk_spacing_power          = 1.0;
    params.track.kalman.hk_enable_overlap         = true;
    params.track.kalman.hk_overlap_mm_s           = 1.0;
    
    params.track.qc.enable_direction_constraint    = false;
    params.track.qc.max_angle_change_deg           = 90;
    params.track.qc.enable_acceleration_constraint = false;
    params.track.qc.acceleration_C_factor          = 3.0;
    params.track.qc.enable_vd_constraint           = false;
    params.track.qc.max_vd_ratio                   = 2.0;
    
    params.track.enable_postprocessing = true;
    params.track.smoothing_factor      = 5;
    
    params.render.upsampling_factor    = 3;
    params.render.method               = 'histogram';
    params.render.interpolation_step   = 0.5;  % FIX B-03: used in runPostProcessing
end