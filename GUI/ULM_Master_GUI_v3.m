function ULM_Master_GUI_v3()
% =========================================================================
% ULM_MASTER_GUI_V3 - Ultrasound Localization Microscopy Processing Suite
% =========================================================================
%
% DESCRIPTION:
%   Enhanced version with bug fixes, performance optimizations, and advanced
%   features driven from a single central algorithm registry.
%
% CORE INFRASTRUCTURE:
%   - Robust error handling and input validation
%   - Session save/load functionality
%   - Undo/Redo system for parameter changes
%   - Progress bars for long operations
%   - Memory profiling
%   - Parallel processing support
%   - Debounced display updates for smooth interaction
%   - Spatial cropping to improve SVD performance
%
% PIPELINE STAGES:
%   1. FILTER:   Spatial Crop, Clutter filtering (SVD, SSM, DCC, Block-wise SVD,
%                optional Butterworth, optional spatial conditioning)
%   2. DETECT:   Microbubble detection (Intensity / Neyman-Pearson / NCC)
%                with ROI masking
%   3. LOCALIZE: Sub-pixel localization (Radial Symmetry, Gaussian Fit,
%                Gaussian Fit Fast)
%   4. TRACK:    Trajectory linking (Kalman, Kalman_Advanced (HK),
%                Hungarian, Nearest Neighbor). Includes Kalman Gain
%                diagnostics panel.
%   5. PROCESS:  Track smoothing (Savitzky-Golay / LOESS / Gaussian / movmean)
%                and interpolation (spline / pchip / linear / makima)
%   6. RENDER:   Super-resolution density & velocity maps (histogram /
%                Gaussian splat)
%
% NEW IN THIS VERSION:
%   - Central algorithm registry (getAlgorithmRegistry) drives every dropdown
%     and dispatch; adding a new tracker is a one-line change.
%   - Block-wise SVD (Song 2017) integrated as filter option with its own
%     parameter panel.
%   - "Advanced..." modal dialogs on Localize / Track / Render tabs expose
%     every parameter from setDefaultParams.m — each field is grouped by
%     pipeline stage and carries a plain-English descriptor.
%   - Global hover tooltip system (applyTooltips + getTooltipDictionary)
%     gives new users a built-in mini user-guide.
%   - Kalman Gain diagnostics panel: effective gain summary, histogram,
%     spatial-influence heatmap, and per-track gain trace. Works for any
%     tracker; automatically uses tracker-reported gains if available.
%
% AUTHOR: Grigori Shapiro
% DATE:   February 2026 (registry + advanced modals + Kalman gain: April 2026)
% =========================================================================

    % --- System Setup ---
    addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
    disp('ULM System v3.0: Environment initialized.');

    % --- Initialize Application Structure FIRST ---
    app = struct();
    app.data = struct();
    app.ui = struct();
    app.state = struct();

    % --- Create Constants ---
    app.constants = ULM_Constants();

    % --- Initialize Main Figure ---
    fig = uifigure('Name', 'ULM Master GUI v3.0 - Professional Edition', ...
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

    % --- Ensure every field the GUI expects actually exists ---
    % (defensive: older sessions or minimal fallbacks may omit some fields)
    app.data.params = ensureAllParamFields(app.data.params);

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

    % --- Apply hover tooltips to every control (global helper) ---
    applyTooltips(guidata(fig));

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
% --- CENTRAL ALGORITHM REGISTRY ------------------------------------------
% Single source of truth for every algorithm choice in the GUI. To add a
% new method (e.g. a new Kalman variant), add ONE entry here; dropdowns,
% option panels, tooltips and dispatch switches pick it up automatically.
% =========================================================================

function reg = getAlgorithmRegistry()
    % ---------------- CLUTTER FILTERS ----------------
    reg.filter = struct( ...
        'id',          {'svd_filter',    'svd_ssm',        'dcc_svd',       'svd_blockwise'}, ...
        'display',     {'Global SVD',    'SVD (SSM-Auto)', 'DCC-SVD',       'Block-wise SVD'}, ...
        'func',        {'run_SVD_Decomposition', 'SVD_SSM', 'DCC_SVD',      'SVD_blockwise'}, ...
        'panel_field', {'p_svd',         '',               'p_dcc',         'p_blockwise'}, ...
        'tooltip',     { ...
            'Single global SVD over the whole FOV. Fast, but assumes tissue motion is uniform.', ...
            'Spatial-similarity auto-thresholding (Baranger 2018). Picks cutoffs automatically.', ...
            'Density-clustering of SVD components into tissue/blood/noise (interactive sliders).', ...
            'Block-wise adaptive SVD (Song 2017). Local cutoffs per block - best for in vivo.' ...
        });

    % ---------------- TRACKERS ----------------
    reg.track = struct( ...
        'id',         {'Kalman',         'Kalman_Advanced',     'Hungarian',            'nn'}, ...
        'display',    {'Kalman',         'Hierarchical Kalman', 'Hungarian',            'Nearest Neighbor'}, ...
        'func',       {'trackKalman',    'trackKalman_Advanced','trackHungarian',       'trackNearestNeighbor'}, ...
        'isKalman',   {true,                true,                  false,                  false}, ...
        'showsGain',  {true,                true,                  false,                  false}, ...
        'usesHK',     {false,               true,                  false,                  false}, ...
        'tooltip',    { ...
            'Classic Kalman filter with fixed process noise and refined gating. Solid baseline.', ...
            'Hierarchical Kalman: multiple velocity bands tracked in parallel. Best for mixed flow speeds.', ...
            'Linear-assignment solver, no motion model. Robust but ignores direction.', ...
            'Greedy nearest-neighbor linking. Fastest, most fragile in dense regions.' ...
        });

    % ---------------- DETECTORS ----------------
    reg.detect = struct( ...
        'id',            {'Intensity',           'NP',             'NCC'}, ...
        'display',       {'Intensity Threshold', 'Neyman-Pearson', 'Normalized Cross-Correlation'}, ...
        'needsTemplate', {false,                 false,            true}, ...
        'tooltip',       { ...
            'Local maxima above a fixed intensity fraction. Simple and fast.', ...
            'Statistical detector with controlled false-alarm rate (alpha0). Robust to variable noise.', ...
            'Template matching against a PSF kernel. Best shape selectivity - slower.' ...
        });

    % ---------------- LOCALIZERS ----------------
    reg.loc = struct( ...
        'id',         {'radial',                         'gaussian_fit',                      'gaussian_fit_fast'}, ...
        'display',    {'Radial Symmetry',                'Gaussian Fit',                      'Gaussian Fit (Fast)'}, ...
        'isGaussian', {false,                            true,                                true}, ...
        'tooltip',    { ...
            'Parsakhoo-Parthasarathy radial symmetry. Sub-pixel accurate, no parametric fit.', ...
            'Full 2D Gaussian nonlinear fit per ROI. Most accurate, reports R^2.', ...
            'Linearized Gaussian fit in closed form. ~10x faster, slightly less accurate.' ...
        });

    % ---------------- RENDERERS ----------------
    reg.render = struct( ...
        'id',      {'histogram',                                'gaussian'}, ...
        'display', {'Histogram (hard bins)',                    'Gaussian Splat'}, ...
        'tooltip', { ...
            'Accumulates one count per track point into the nearest super-resolution pixel.', ...
            'Splats each track point as a small Gaussian - smoother, slightly blurrier.' ...
        });

    % ---------------- SMOOTHERS ----------------
    reg.smoothing = struct( ...
        'id',      {'sgolay',                                    'rloess',                                 'gaussian',                   'movmean'}, ...
        'display', {'Savitzky-Golay',                            'Robust LOESS',                           'Gaussian',                   'Moving Mean'}, ...
        'tooltip', { ...
            'Polynomial smoothing that preserves peaks. Best default for tracks.', ...
            'Outlier-resistant local regression. Slower, better when localizations are noisy.', ...
            'Simple Gaussian kernel smoothing.', ...
            'Plain moving average - fastest, most aggressive.' ...
        });
end

% ---- Registry convenience accessors ----
function ids = registryIds(entries)
    ids = {entries.id};
end

function idx = registryIndexById(entries, id)
    idx = find(strcmpi({entries.id}, id), 1);
    if isempty(idx), idx = 1; end
end

function e = registryEntryById(entries, id)
    e = entries(registryIndexById(entries, id));
end

% ---- Generic safe nested-field accessor ----
function v = getDefault(s, path, fallback)
% Traverse a dot-separated path into struct `s`.  Returns `fallback` only
% when an intermediate field is missing — never for intentionally empty
% leaf values (e.g. cropBox = []).
    parts = strsplit(path, '.');
    cur = s;
    for i = 1:numel(parts)
        if isstruct(cur) && isfield(cur, parts{i})
            cur = cur.(parts{i});
        else
            v = fallback; return;
        end
    end
    v = cur;
end

% ---- Fill in any missing fields expected by the GUI ----
function params = ensureAllParamFields(params)
    % Filter defaults
    if ~isfield(params, 'filter'), params.filter = struct(); end
    if ~isfield(params.filter, 'spatial_method'), params.filter.spatial_method = 'Gaussian'; end
    if ~isfield(params.filter, 'spatial_kernel'), params.filter.spatial_kernel = 3; end
    if ~isfield(params.filter, 'spatial_sigma1'), params.filter.spatial_sigma1 = 1.0; end
    if ~isfield(params.filter, 'spatial_sigma2'), params.filter.spatial_sigma2 = 2.0; end

    % Block-wise SVD defaults
    if ~isfield(params.filter, 'blockwise'), params.filter.blockwise = struct(); end
    bwd = { ...
        'threshold_method',    'DopplerGradient'; ...
        'block_size_mm',        4.0; ...
        'overlap_pct',          75.0; ...
        'manual_cutoff',        [10 200]; ...
        'tissue_freq_hz',       -1; ...
        'mp_deviation_sigma',   2.0; ...
        'gradient_pct',         0.10; ...
        'min_blood_comps',      3; ...
        'max_tissue_frac',      0.60; ...
        'plot_maps',            false};
    for i = 1:size(bwd,1)
        if ~isfield(params.filter.blockwise, bwd{i,1})
            params.filter.blockwise.(bwd{i,1}) = bwd{i,2};
        end
    end

    % IO
    if ~isfield(params, 'io'), params.io = struct(); end
    if ~isfield(params.io, 'cropBox'), params.io.cropBox = []; end

    % Kalman flags
    if ~isfield(params, 'track'), params.track = struct(); end
    if ~isfield(params.track, 'kalman'), params.track.kalman = struct(); end
    kdef = { ...
        'use_direction',                  true; ...
        'use_angle',                      true; ...
        'use_brightness',                 true; ...
        'direction_penalty_weight',       2; ...
        'angle_penalty_slope',            0.3; ...
        'brightness_penalty_weight',      2; ...
        'direction_history_points',       4; ...
        'max_angle_change_deg',           70; ...
        'gating_max_angle_change_deg',    90; ...
        'hk_alpha',                       0.01; ...
        'hk_beta',                        0.025; ...
        'hk_forward_backward',            true; ...
        'hk_v_max',                       20; ...
        'hk_num_levels',                  5; ...
        'hk_spacing_power',               1.0; ...
        'hk_enable_overlap',              true; ...
        'hk_overlap_mm_s',                2.0};
    for i = 1:size(kdef,1)
        if ~isfield(params.track.kalman, kdef{i,1})
            params.track.kalman.(kdef{i,1}) = kdef{i,2};
        end
    end

    % Localization extras
    if ~isfield(params, 'loc'), params.loc = struct(); end
    ldef = { ...
        'psf_type',             'Gaussian'; ...
        'psf_size',             [5 5]; ...
        'psf_file_path',        ''; ...
        'qc_max_roi_maxima',    3; ...
        'min_gradient_squared', 1e-6; ...
        'min_determinant',      1e-6};
    for i = 1:size(ldef,1)
        if ~isfield(params.loc, ldef{i,1})
            params.loc.(ldef{i,1}) = ldef{i,2};
        end
    end

    if ~isfield(params.loc, 'h_contrast')
        params.loc.h_contrast = 0;  % 0 = classic mode
    end

    % Smoothing / interpolation / render extras
    if ~isfield(params.track, 'smoothing_method'), params.track.smoothing_method = 'sgolay'; end
    if ~isfield(params, 'render'), params.render = struct(); end
    if ~isfield(params.render, 'interpolation_method'), params.render.interpolation_method = 'spline'; end
    if ~isfield(params.render, 'gaussian_sigma'), params.render.gaussian_sigma = 0.3; end

    % Analysis
    if ~isfield(params, 'analysis'), params.analysis = struct(); end
    if ~isfield(params.analysis, 'tortuosity_bins'), params.analysis.tortuosity_bins = 0:0.05:8; end
    if ~isfield(params.analysis, 'velocity_hist_num_bins'), params.analysis.velocity_hist_num_bins = 60; end
    if ~isfield(params.analysis, 'density_grid_size_mm'), params.analysis.density_grid_size_mm = 0.5; end

    % Initialize transient data fields on params struct consumers ---------
    % (the GUI also initializes these on app.data, below)
end

% =========================================================================
% --- GUI LAYOUT BUILDER --------------------------------------------------
% =========================================================================

function app = buildGUILayout(app)
    fig = app.fig;

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
    app.data.vesselMapSource = 'none';

    % --- SVD Cache ---
    app.data.U = [];
    app.data.S_diag = [];
    app.data.V = [];
    app.data.svdDims = [];
    app.data.rawDataHash = '';

    % --- DCC Components ---
    app.data.tissue_indices = [];
    app.data.blood_indices = [];
    app.data.noise_indices = [];

    % --- Block-wise diagnostics slot ---
    app.data.blockwiseDiag = [];

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

    % --- Main Grid: 3 rows (menu, main content, status) x 2 columns ---
    gl = uigridlayout(fig, [3, 2]);
    gl.RowHeight = {40, '1x', 40};
    gl.ColumnWidth = {'1x', 340};
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
    menuGrid.BackgroundColor = [0.94 0.94 0.94];

    % Load/Save Session
    app.ui.btnLoadSession = uibutton(menuGrid, 'Text', 'Load Work Session', ...
        'ButtonPushedFcn', @(s,e) loadSession(app.fig));

    app.ui.btnSaveSession = uibutton(menuGrid, 'Text', 'Save Work Session', ...
        'ButtonPushedFcn', @(s,e) saveSession(app.fig));

    % Undo/Redo
    app.ui.btnUndo = uibutton(menuGrid, 'Text', 'Undo', ...
        'ButtonPushedFcn', @(s,e) performUndo(app.fig), ...
        'Enable', 'off');

    app.ui.btnRedo = uibutton(menuGrid, 'Text', 'Redo', ...
        'ButtonPushedFcn', @(s,e) performRedo(app.fig), ...
        'Enable', 'off');

    % Memory Monitor
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
    stageGrid = uigridlayout(parentGrid, [3, 2]);
    stageGrid.Layout.Row = 2;
    stageGrid.Layout.Column = 1;
    stageGrid.RowHeight = {'1x', 'fit', 30};
    stageGrid.ColumnWidth = {180, '1x'};
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

    app.ui.LoadButton = uibutton(bottomGrid, 'Text', 'Load Data (IQ / imageData)', ...
        'ButtonPushedFcn', @(s,e) loadRawData(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.3 0.7 1]);

    app.ui.PlayPauseButton = uibutton(bottomGrid, 'Text', 'Play', ...
        'ButtonPushedFcn', @(s,e) togglePlayback(app.fig));

    app.ui.FrameField = uieditfield(bottomGrid, 'numeric', 'Value', 1, ...
        'Limits', [1 Inf], 'RoundFractionalValues', 'on', ...
        'ValueChangedFcn', @(s,e) onFrameFieldChanged(app.fig, e));
end

function app = buildDisplayControlsPanel(app, parentGrid)
    p_disp = uipanel(parentGrid, 'Title', 'Visual Adjustments', 'Scrollable', 'on');
    p_disp.Layout.Row = 1;
    p_disp.Layout.Column = 1;

    g = uigridlayout(p_disp, [9, 1]);
    g.RowHeight = {'fit','fit','fit','fit','fit','fit','fit','fit','1x'};
    g.Padding = [5 5 5 5];

    updateVis = @(s,e) displayCurrentFrame(app.fig);

    app.ui.disp_mat2gray = uicheckbox(g, 'Text', 'Normalize (mat2gray)', ...
        'Value', 1, 'ValueChangedFcn', updateVis);

    app.ui.disp_log = uicheckbox(g, 'Text', 'Log Compression', ...
        'Value', 0, 'ValueChangedFcn', updateVis);

    uilabel(g, 'Text', 'Gamma (Stretch):', 'FontWeight', 'bold');
    app.ui.disp_gamma = uislider(g, 'Limits', [0.1 5.0], 'Value', 1.0, ...
        'ValueChangedFcn', updateVis);

    uilabel(g, 'Text', 'Colormap:', 'FontWeight', 'bold');
    app.ui.disp_cmap = uidropdown(g, 'Items', {'gray', 'hot', 'jet', 'parula'}, ...
        'Value', 'gray', 'ValueChangedFcn', updateVis);

    app.ui.disp_clim_auto = uicheckbox(g, 'Text', 'Auto CLim', ...
        'Value', 1, 'ValueChangedFcn', updateVis);

    cg = uigridlayout(g, [1 2]);
    cg.Padding = [0 0 0 0];
    app.ui.disp_clim_min = uieditfield(cg, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', updateVis);
    app.ui.disp_clim_max = uieditfield(cg, 'numeric', 'Value', 1, ...
        'ValueChangedFcn', updateVis);
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
    statusGrid.Layout.Row = 3;
    statusGrid.Layout.Column = [1 2];
    statusGrid.ColumnWidth = {'fit', '1x', 'fit'};
    statusGrid.Padding = [10 5 10 5];
    statusGrid.BackgroundColor = [0.94 0.94 0.94];

    app.ui.lblStatusTitle = uilabel(statusGrid, 'Text', 'Status:', 'FontWeight', 'bold');
    app.ui.StatusLabel = uilabel(statusGrid, 'Text', 'Ready');
    app.ui.lblStatus = app.ui.StatusLabel;   % alias for setStatus()
    app.ui.StatusLamp = uilamp(statusGrid, 'Color', 'blue');

    app.ui.ProgressBar = [];
end

% =========================================================================
% --- TAB BUILDERS --------------------------------------------------------
% =========================================================================

function scrollPanel = createScrollableTabPanel(tab)
% CREATESCROLLABLETABPANEL  Fills a uitab with a scrollable panel.
%   Uses a 1x1 grid wrapper (required because uitab's AutoResizeChildren
%   overrides Units/Position), then places a Scrollable panel inside it.
    g = uigridlayout(tab, [1 1]);
    g.Padding = [0 0 0 0];
    scrollPanel = uipanel(g, ...
        'BorderType', 'none', ...
        'Scrollable', 'on');
end

function app = buildFilterTab(app)
    reg = getAlgorithmRegistry();
    app.ui.tabFilter = uitab(app.ui.tabGroup, 'Title', '1. Filter');
    
    % 1. MAIN TAB LAYOUT: 2 Rows
    % Row 1: Scrollable area ('1x' - takes all available space)
    % Row 2: Pinned Run Button (45px - fixed at the bottom)
    g_tab = uigridlayout(app.ui.tabFilter, [2, 1]);
    g_tab.RowHeight = {'1x', 45};
    g_tab.ColumnWidth = {'1x'};
    g_tab.Padding = [5 5 5 5];
    g_tab.RowSpacing = 10;

    % 2. INNER GRID (SCROLLABLE)
    % We apply 'Scrollable' DIRECTLY to the grid layout. No uipanel needed!
    g = uigridlayout(g_tab, [10, 1], 'Scrollable', 'on');
    g.Layout.Row = 1;
    g.Layout.Column = 1;
    app.ui.g_filter = g; % Save reference for dynamic height adjustments
    
    g.ColumnWidth = {'1x'};
    g.RowSpacing = 8;
    g.Padding = [5 5 20 5]; % Extra right padding to avoid hiding text behind scrollbar
    
    % STRICT ABSOLUTE HEIGHTS IN PIXELS. NO '1x' ALLOWED HERE!
    % Order: 1.Label, 2.Drop, 3.Crop, 4.Mask, 5.SVD, 6.DCC, 7.Butter, 8.Spatial, 9.Blockwise, 10.Spacer
    g.RowHeight = {25, 30, 140, 140, 0, 0, 130, 110, 0, 20};

    % Row 1 & 2: Method Selection
    lbl = uilabel(g, 'Text', 'Filter Method:', 'FontWeight', 'bold');
    lbl.Layout.Row = 1;

    filterItems = registryIds(reg.filter);
    matchIdx = find(strcmpi(filterItems, app.data.params.filter.method), 1);
    if isempty(matchIdx), matchIdx = 1; end
    
    app.ui.FilterMethodDropdown = uidropdown(g, ...
        'Items', filterItems, ...
        'Value', filterItems{matchIdx});
    app.ui.FilterMethodDropdown.Layout.Row = 2;

    % Row 3: Spatial Crop Panel
    app = buildCropPanel(app, g);
    app.ui.p_crop.Layout.Row = 3;

    % Row 4: Masking Panel
    app.ui.p_mask = uipanel(g, 'Title', 'Masking (for Localization)');
    app.ui.p_mask.Layout.Row = 4;
    g_mask = uigridlayout(app.ui.p_mask, [2, 2]);
    app.ui.LoadMaskButton = uibutton(g_mask, 'Text', 'Load Mask', ...
        'ButtonPushedFcn', @(s,e) loadMask(app.fig));
    app.ui.CreateMaskButton = uibutton(g_mask, 'Text', 'Create New Mask', ...
        'ButtonPushedFcn', @(s,e) runCreateMask(app.fig));
    app.ui.ResetMaskButton = uibutton(g_mask, 'Text', 'Reset Mask', ...
        'ButtonPushedFcn', @(s,e) resetMask(app.fig));
    app.ui.maskStatusLabel = uilabel(g_mask, 'Text', 'Status: None', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Row 5: SVD Parameters
    app.ui.p_svd = uipanel(g, 'Title', 'Standard SVD Params');
    app.ui.p_svd.Layout.Row = 5;
    g_svd = uigridlayout(app.ui.p_svd, [2, 2]);
    uilabel(g_svd, 'Text', 'Cutoff Start:');
    app.ui.SVDCutoffStart = uieditfield(g_svd, 'numeric', ...
        'Value', app.data.params.filter.svd_cutoff(1), ...
        'ValueChangedFcn', @(s,e) runFilterWithValidation(app.fig));
    uilabel(g_svd, 'Text', 'Cutoff End:');
    app.ui.SVDCutoffEnd = uieditfield(g_svd, 'numeric', ...
        'Value', app.data.params.filter.svd_cutoff(2), ...
        'ValueChangedFcn', @(s,e) runFilterWithValidation(app.fig));

    % Row 6: DCC Parameters
    app = buildDCCPanel(app, g);
    app.ui.p_dcc.Layout.Row = 6;

    % Row 7: Butterworth Filter
    app = buildButterworthPanel(app, g);
    app.ui.p_butter.Layout.Row = 7;

    % Row 8: Spatial Filter
    app = buildSpatialFilterPanel(app, g);
    app.ui.p_spatial.Layout.Row = 8;

    % Row 9: Block-wise SVD Panel
    app = buildBlockwisePanel(app, g);
    app.ui.p_blockwise.Layout.Row = 9;

    % Row 10: Empty label to provide bottom padding inside the scroll area
    lblPad = uilabel(g, 'Text', '');
    lblPad.Layout.Row = 10;

    % 3. RUN BUTTON: Pinned to the outer tab grid (Row 2), ALWAYS visible
    app.ui.RunFilterButton = uibutton(g_tab, 'Text', 'Run Filter', ...
        'ButtonPushedFcn', @(s,e) runFilter(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunFilterButton.Layout.Row = 2;
    app.ui.RunFilterButton.Layout.Column = 1;
end

function app = buildCropPanel(app, parentGrid)
    app.ui.p_crop = uipanel(parentGrid, 'Title', 'Spatial Crop (Pre-Processing)');

    g_crop = uigridlayout(app.ui.p_crop, [2, 3]);
    g_crop.ColumnWidth = {'1x', '1x', '1x'};
    g_crop.RowHeight = {'fit', 'fit'};

    uilabel(g_crop, 'Text', 'Crop Box [x y w h]:');
    app.ui.CropBoxField = uieditfield(g_crop, 'text', 'Value', '[]');
    app.ui.InteractiveCropBtn = uibutton(g_crop, 'Text', 'Interactive Crop', ...
        'ButtonPushedFcn', @(s,e) runInteractiveCrop(app.fig));

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
    g.RowHeight = {'fit', 'fit', 'fit'};

    uilabel(g, 'Text', 'Tissue:', 'FontWeight', 'bold');
    uilabel(g, 'Text', 'Start (%):');
    app.ui.DCCTissueStart = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
    uilabel(g, 'Text', 'End (%):');
    app.ui.DCCTissueEnd = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));

    uilabel(g, 'Text', 'Blood:', 'FontWeight', 'bold');
    uilabel(g, 'Text', 'Start (%):');
    app.ui.DCCBloodStart = uieditfield(g, 'numeric', 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));
    uilabel(g, 'Text', 'End (%):');
    app.ui.DCCBloodEnd = uieditfield(g, 'numeric', 'Value', 100, ...
        'ValueChangedFcn', @(s,e) onDCCSliderChanged(app.fig));

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
    g.RowHeight = {'fit', 'fit', 'fit'};

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
    g.RowHeight = {'fit', 'fit'};

    uilabel(g, 'Text', 'Method:');
    app.ui.SpatialMethodDrop = uidropdown(g, ...
        'Items', {'None', 'Gaussian', 'Median', 'DoG', 'Top-Hat'}, ...
        'Value', app.data.params.filter.spatial_method, ...
        'ValueChangedFcn', @(s,e) updateSpatialOptions(app.fig));

    app.ui.lblSpatialKernel = uilabel(g, 'Text', 'Kernel (px):');
    app.ui.SpatialKernelField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.filter.spatial_kernel, ...
        'RoundFractionalValues', 'on', ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'spatial_kernel'));

    app.ui.lblSpatialSigma1 = uilabel(g, 'Text', 'Sigma 1:');
    app.ui.SpatialSigma1Field = uieditfield(g, 'numeric', ...
        'Value', app.data.params.filter.spatial_sigma1, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'spatial_sigma1'));

    app.ui.lblSpatialSigma2 = uilabel(g, 'Text', 'Sigma 2:');
    app.ui.SpatialSigma2Field = uieditfield(g, 'numeric', ...
        'Value', app.data.params.filter.spatial_sigma2, ...
        'ValueChangedFcn', @(s,e) saveParamState(app.fig, 'spatial_sigma2'));

    updateSpatialOptions(app.fig, app);
end

function updateSpatialOptions(fig, app_in)
    if nargin < 2
        app = guidata(fig);
    else
        app = app_in;
    end

    method = app.ui.SpatialMethodDrop.Value;

    app.ui.lblSpatialKernel.Visible = 'off'; app.ui.SpatialKernelField.Visible = 'off';
    app.ui.lblSpatialSigma1.Visible = 'off'; app.ui.SpatialSigma1Field.Visible = 'off';
    app.ui.lblSpatialSigma2.Visible = 'off'; app.ui.SpatialSigma2Field.Visible = 'off';

    switch method
        case 'Gaussian'
            app.ui.lblSpatialKernel.Text = 'Kernel (px):';
            app.ui.lblSpatialKernel.Visible = 'on'; app.ui.SpatialKernelField.Visible = 'on';
            app.ui.lblSpatialSigma1.Text = 'Sigma:';
            app.ui.lblSpatialSigma1.Visible = 'on'; app.ui.SpatialSigma1Field.Visible = 'on';
        case 'Median'
            app.ui.lblSpatialKernel.Text = 'Kernel (px):';
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

% -------------------------------------------------------------------------
% Block-wise SVD panel (Song et al. 2017). Appears when filter method is
% 'svd_blockwise'. All parameter semantics match SVD_blockwise.m exactly.
% -------------------------------------------------------------------------
function app = buildBlockwisePanel(app, parent)
    % 1. Create the main wrapper panel
    app.ui.p_blockwise = uipanel(parent, 'Title', 'Block-wise SVD Parameters');
    
    % 2. Create a strict internal grid (3 Rows, 2 Columns) matching your original UI
    g_bw = uigridlayout(app.ui.p_blockwise, [3, 2]);
    g_bw.RowHeight = {25, 25, 25};
    g_bw.ColumnWidth = {'fit', '1x'};
    g_bw.RowSpacing = 8;
    g_bw.Padding = [10 10 10 10];

    % 3. Ensure fallback parameters exist to prevent crashes
    if ~isfield(app.data.params.filter, 'block_size')
        app.data.params.filter.block_size = [32 32];
    end
    if ~isfield(app.data.params.filter, 'overlap_pct')
        app.data.params.filter.overlap_pct = 0.5;
    end
    if ~isfield(app.data.params.filter, 'svd_cutoff')
        app.data.params.filter.svd_cutoff = [5 100];
    end

    % --- ROW 1: Block Size [Z, X] (Original Text Array) ---
    lblZ = uilabel(g_bw, 'Text', 'Block Size [Z, X]:');
    lblZ.Layout.Row = 1; lblZ.Layout.Column = 1;
    app.ui.BlockSize = uieditfield(g_bw, 'text', ...
        'Value', num2str(app.data.params.filter.block_size), ...
        'ValueChangedFcn', @(s,e) runFilterWithValidation(app.fig));
    app.ui.BlockSize.Layout.Row = 1; app.ui.BlockSize.Layout.Column = 2;

    % --- ROW 2: Overlap [%] (Original Numeric) ---
    lblOv = uilabel(g_bw, 'Text', 'Overlap [%]:');
    lblOv.Layout.Row = 2; lblOv.Layout.Column = 1;
    app.ui.BlockOverlap = uieditfield(g_bw, 'numeric', ...
        'Value', app.data.params.filter.overlap_pct * 100, ...
        'ValueChangedFcn', @(s,e) runFilterWithValidation(app.fig));
    app.ui.BlockOverlap.Layout.Row = 2; app.ui.BlockOverlap.Layout.Column = 2;

    % --- ROW 3: SVD Cutoffs [Start, End] (Original Text Array) ---
    lblCut = uilabel(g_bw, 'Text', 'SVD Cutoffs [Start, End]:');
    lblCut.Layout.Row = 3; lblCut.Layout.Column = 1;
    app.ui.BlockCutoff = uieditfield(g_bw, 'text', ...
        'Value', num2str(app.data.params.filter.svd_cutoff), ...
        'ValueChangedFcn', @(s,e) runFilterWithValidation(app.fig));
    app.ui.BlockCutoff.Layout.Row = 3; app.ui.BlockCutoff.Layout.Column = 2;
end

function updateBlockwiseOptions(fig, app_in)
    if nargin < 2, app = guidata(fig); else, app = app_in; end
    if ~isfield(app.ui, 'BWThresholdMethod'), return; end
    isManual = strcmpi(app.ui.BWThresholdMethod.Value, 'Manual');
    app.ui.lblBWManualCutoff.Visible = isManual;
    app.ui.BWManualCutoff.Visible    = isManual;
    if nargin < 2, guidata(fig, app); end
end

function app = buildDetectTab(app)
    app.ui.tabDetect = uitab(app.ui.tabGroup, 'Title', '2. Detect');
    
    % 1. MAIN TAB LAYOUT: 2 Rows (Scrollable content area + Pinned Run Button)
    g_tab = uigridlayout(app.ui.tabDetect, [2, 1]);
    g_tab.RowHeight = {'1x', 45};
    g_tab.ColumnWidth = {'1x'};
    g_tab.Padding = [5 5 5 5];
    g_tab.RowSpacing = 10;

    % 2. INNER GRID (SCROLLABLE) - Holds the sequential step panels
    g_scroll = uigridlayout(g_tab, [2, 1], 'Scrollable', 'on');
    g_scroll.Layout.Row = 1;
    g_scroll.Layout.Column = 1;
    g_scroll.RowHeight = {380, 'fit'}; % ROI Panel takes 380px, Detect parameters panel auto-fits
    g_scroll.Padding = [5 5 20 5]; % Extra right padding to prevent clipping behind the scrollbar

    % Step A: Define ROI Mask
    app = buildROIPanel(app, g_scroll);

    % Step B: Detection Parameters
    app = buildDetectionParamsPanel(app, g_scroll);

    % 3. RUN DETECTION BUTTON: Pinned permanently to the bottom of the tab (Row 2)
    app.ui.RunDetectButton = uibutton(g_tab, 'Text', 'Run Detection (Masked)', ...
        'ButtonPushedFcn', @(s,e) runDetection(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunDetectButton.Layout.Row = 2;
    app.ui.RunDetectButton.Layout.Column = 1;
end

function app = buildROIPanel(app, parentGrid)
    app.ui.p_roi = uipanel(parentGrid, 'Title', 'Step A: Define ROI Mask (Vessel Map)');
    app.ui.p_roi.Layout.Row = 1;

    g = uigridlayout(app.ui.p_roi, [2, 1]);
    g.RowHeight = {'1x', 'fit'};
    g.Padding = [5 5 5 5];

    % Histogram only — the image is shown on the main central axes
    app.ui.axHist = uiaxes(g);
    app.ui.axHist.Layout.Row = 1;
    title(app.ui.axHist, 'Intensity Histogram (Log Scale)');
    app.ui.axHist.YScale = 'log';
    app.ui.axHist.XTickLabel = [];
    grid(app.ui.axHist, 'on');

    g_ctrl = uigridlayout(g, [4, 3]);
    g_ctrl.Layout.Row = 2;
    % Strict pixel heights prevent MATLAB from aggressively crushing sliders/buttons
    g_ctrl.RowHeight = {40, 40, 40, 35};
    g_ctrl.RowSpacing = 5;
    g_ctrl.ColumnWidth = {80, 70, '1x'};

    uilabel(g_ctrl, 'Text', '1. Enhance:', 'FontWeight', 'bold');
    app.ui.EnhanceMethodDrop = uidropdown(g_ctrl, ...
        'Items', {'None', 'CLAHE (Local Contrast)', 'Top-Hat (Vesselness)', 'Sharpen'}, ...
        'Value', 'None', ...
        'ValueChangedFcn', @(s,e) applyVesselEnhancement(app.fig));
    app.ui.EnhanceAmountSlider = uislider(g_ctrl, 'Limits', [0 1], 'Value', 0.5, ...
        'ValueChangedFcn', @(s,e) applyVesselEnhancement(app.fig, e.Value));

    uilabel(g_ctrl, 'Text', '2. Gamma:', 'FontWeight', 'bold');
    app.ui.ROIContrastField = uieditfield(g_ctrl, 'numeric', 'Value', 1, ...
        'ValueDisplayFormat', '%.3f', ...
        'ValueChangedFcn', @(s,e) onContrastChange(app.fig, e.Value));
    app.ui.ROIContrastSlider = uislider(g_ctrl, 'Limits', [0.1 3.0], 'Value', 1, ...
        'ValueChangedFcn', @(s,e) onContrastChange(app.fig, e.Value));

    uilabel(g_ctrl, 'Text', '3. Threshold:', 'FontWeight', 'bold');
    app.ui.ROIThreshField = uieditfield(g_ctrl, 'numeric', 'Value', 0, ...
        'ValueDisplayFormat', '%.3f', ...
        'ValueChangedFcn', @(s,e) onROIChange(app.fig, e.Value));
    app.ui.ROIThreshSlider = uislider(g_ctrl, 'Limits', [0 1.1], 'Value', 0, ...
        'ValueChangedFcn', @(s,e) onROIChange(app.fig, e.Value));

    app.ui.SaveMaskButton = uibutton(g_ctrl, ...
        'Text', 'Save Mask', ...
        'ButtonPushedFcn', @(s,e) saveCreatedMask(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.8 0.9 1]);
    app.ui.SaveMaskButton.Layout.Row = 4;
    app.ui.SaveMaskButton.Layout.Column = [1 2];

    app.ui.ResetROIButton = uibutton(g_ctrl, ...
        'Text', 'Reset All', ...
        'ButtonPushedFcn', @(s,e) resetROIPanel(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [1.0 0.85 0.85]);
    app.ui.ResetROIButton.Layout.Row = 4;
    app.ui.ResetROIButton.Layout.Column = 3;
end

function app = buildDetectionParamsPanel(app, parentGrid)
    reg = getAlgorithmRegistry();
    
    app.ui.p_detect_params = uipanel(parentGrid, 'Title', 'Step B: Detection');
    app.ui.p_detect_params.Layout.Row = 2;
   
    g = uigridlayout(app.ui.p_detect_params, [9, 3]);
    g.RowHeight = {'fit','fit','fit','fit','fit','fit','fit','fit','fit'};
    g.ColumnWidth = {'fit', 60, '1x'};
    
    % Row 1: Method
    lbl_method = uilabel(g, 'Text', 'Detect Method:', 'FontWeight', 'bold');
    lbl_method.Layout.Row = 1;
    lbl_method.Layout.Column = 1;
    
    app.ui.DetectMethodDropdown = uidropdown(g, ...
        'Items', registryIds(reg.detect), ...
        'Value', app.data.params.loc.DetectMethod, ...
        'ValueChangedFcn', @(s,e) updateDetectionOptions(app.fig));
    app.ui.DetectMethodDropdown.Layout.Row = 1;
    app.ui.DetectMethodDropdown.Layout.Column = [2 3];
    
    % Row 2: Method-specific sub-panel
    app.ui.p_detect_method_params = uipanel(g, 'Title', 'Method Parameters');
    app.ui.p_detect_method_params.Layout.Row = 2;
    app.ui.p_detect_method_params.Layout.Column = [1 3];
    
    g_mp = uigridlayout(app.ui.p_detect_method_params, [2, 2]);
    g_mp.ColumnWidth = {'fit', '1x'};
    g_mp.RowHeight = {'fit', 'fit'};
    
    app.ui.lbl_NP_alpha = uilabel(g_mp, 'Text', 'NP alpha0:');
    app.ui.NP_AlphaField = uieditfield(g_mp, 'numeric', ...
        'Value', app.data.params.loc.NP_alpha0, ...
        'ValueDisplayFormat', '%.2e');
        
    app.ui.lbl_NCC_thresh = uilabel(g_mp, 'Text', 'NCC tau:');
    app.ui.NCC_ThreshField = uieditfield(g_mp, 'numeric', ...
        'Value', app.data.params.loc.crosscor_threshold, ...
        'ValueDisplayFormat', '%.2f');
        
    % Row 3: Peak Contrast
    lbl_ncc_peak = uilabel(g, 'Text', 'Peak Contrast (h):');
    lbl_ncc_peak.Layout.Row = 3;
    lbl_ncc_peak.Layout.Column = 1;
    
    app.ui.PeakContrastField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.loc.h_contrast, ...
        'ValueDisplayFormat', '%.3f', ...
        'Limits', [0 0.5]);
    app.ui.PeakContrastField.Layout.Row = 3;
    app.ui.PeakContrastField.Layout.Column = [2 3];
    
    % Row 4: PSF FWHM
    lbl_fwhm = uilabel(g, 'Text', 'PSF FWHM [x z] (px):');
    lbl_fwhm.Layout.Row = 4;
    lbl_fwhm.Layout.Column = 1;
    
    app.ui.DetectFWHMField = uieditfield(g, 'text', ...
        'Value', mat2str(app.data.params.loc.fwhm));
    app.ui.DetectFWHMField.Layout.Row = 4;
    app.ui.DetectFWHMField.Layout.Column = [2 3];
    
    % Row 5: Intensity threshold
    lbl_thresh = uilabel(g, 'Text', 'Intensity Thresh:');
    lbl_thresh.Layout.Row = 5;
    lbl_thresh.Layout.Column = 1;
    
    app.ui.LocThreshField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.loc.detection_threshold);
    app.ui.LocThreshField.Layout.Row = 5;
    app.ui.LocThreshField.Layout.Column = 2;
    
    app.ui.LocThreshSlider = uislider(g, 'Limits', [0.001 1], ...
        'Value', app.data.params.loc.detection_threshold);
    app.ui.LocThreshSlider.Layout.Row = 5;
    app.ui.LocThreshSlider.Layout.Column = 3;
    
    syncSliderField(app.ui.LocThreshSlider, app.ui.LocThreshField);
    
    % Row 6: Max bubbles
    lbl_max = uilabel(g, 'Text', 'Max Bubbles:');
    lbl_max.Layout.Row = 6;
    lbl_max.Layout.Column = 1;
    
    app.ui.LocMaxField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.loc.max_bubbles_per_frame);
    app.ui.LocMaxField.Layout.Row = 6;
    app.ui.LocMaxField.Layout.Column = 2;
    
    app.ui.LocMaxSlider = uislider(g, 'Limits', [0 5000], ...
        'Value', app.data.params.loc.max_bubbles_per_frame);
    app.ui.LocMaxSlider.Layout.Row = 6;
    app.ui.LocMaxSlider.Layout.Column = 3;
    
    syncSliderField(app.ui.LocMaxSlider, app.ui.LocMaxField);
    
    % Row 7: Preview ROI checkbox
    app.ui.chkPreviewROI = uicheckbox(g, 'Text', 'Preview ROI Overlay', 'Value', 1, ...
        'ValueChangedFcn', @(s,e) onROIPreviewToggle(app.fig));
    app.ui.chkPreviewROI.Layout.Row = 7;
    app.ui.chkPreviewROI.Layout.Column = [1 3];
    
    % Row 8: Advanced detection modal
    app.ui.BtnAdvDetect = uibutton(g, ...
        'Text', 'Advanced Detection / PSF Parameters...', ...
        'ButtonPushedFcn', @(s,e) openAdvancedLocGUI(app.fig));
    app.ui.BtnAdvDetect.Layout.Row = 8;
    app.ui.BtnAdvDetect.Layout.Column = [1 3];
    
    % Row 9: Detection Hint Label (For PSF template warnings)
    app.ui.LblDetectHint = uilabel(g, 'Text', '', 'WordWrap', 'on');
    app.ui.LblDetectHint.Layout.Row = 9;
    app.ui.LblDetectHint.Layout.Column = [1 3];
end

function app = buildLocalizeTab(app)
    reg = getAlgorithmRegistry();
    app.ui.tabLocalize = uitab(app.ui.tabGroup, 'Title', '3. Localize');

    % 1. MAIN TAB LAYOUT: 2 Rows (Scrollable content area + Pinned Run Button)
    g_tab = uigridlayout(app.ui.tabLocalize, [2, 1]);
    g_tab.RowHeight = {'1x', 45};
    g_tab.ColumnWidth = {'1x'};
    g_tab.Padding = [5 5 5 5];
    g_tab.RowSpacing = 10;

    % 2. INNER GRID (SCROLLABLE) - Holds the localization panels
    g = uigridlayout(g_tab, [7, 1], 'Scrollable', 'on');
    g.Layout.Row = 1;
    g.Layout.Column = 1;
    g.RowSpacing = 8;
    g.Padding = [5 5 20 5]; % Padding to avoid overlap with scrollbar
    g.RowHeight = {25, 30, 110, 110, 40, 0, 40}; % Strict height allocation

    % Method Selection
    lblMethod = uilabel(g, 'Text', 'Localization Method:', 'FontWeight', 'bold');
    lblMethod.Layout.Row = 1;
    app.ui.LocMethodDropdown = uidropdown(g, ...
        'Items', registryIds(reg.loc), ...
        'Value', app.data.params.loc.method);
    app.ui.LocMethodDropdown.Layout.Row = 2;

    % QC Panels (Row 3 & 4)
    % Note: These build functions create their own internal uipanels
    app = buildLocalizationQCPanels(app, g);
    app.ui.p_loc_qc_radial.Layout.Row = 3;
    app.ui.p_loc_qc_gauss.Layout.Row = 4;

    % Advanced Options
    app.ui.BtnAdvLoc = uibutton(g, ...
        'Text', 'Advanced Localization / Detection Parameters...', ...
        'ButtonPushedFcn', @(s,e) openAdvancedLocGUI(app.fig));
    app.ui.BtnAdvLoc.Layout.Row = 5;

    % Density Preview
    app.ui.ShowLocDensityButton = uibutton(g, ...
        'Text', 'Show Localization Density Map', ...
        'ButtonPushedFcn', @(s,e) showLocalizationDensity(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 0.8 1]);
    app.ui.ShowLocDensityButton.Layout.Row = 7;

    % 3. RUN LOCALIZATION BUTTON: Pinned permanently to the bottom (Row 2)
    app.ui.RunLocalizationButton = uibutton(g_tab, 'Text', 'Run Localization', ...
        'ButtonPushedFcn', @(s,e) runLocalization(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunLocalizationButton.Layout.Row = 2;
    app.ui.RunLocalizationButton.Layout.Column = 1;
end

function app = buildLocalizationQCPanels(app, parentGrid)
    % ---- Common QC panel (used by ALL localization methods) ----
    app.ui.p_loc_qc_radial = uipanel(parentGrid, 'Title', 'Localization QC');
    app.ui.p_loc_qc_radial.Layout.Row = 3;
    g = uigridlayout(app.ui.p_loc_qc_radial, [3, 2]);
    g.ColumnWidth = {'fit', '1x'};
    g.RowHeight = {'fit', 'fit', 'fit'};

    app.ui.LocQCDivergence = uicheckbox(g, 'Text', 'Enable Divergence Check', ...
        'Value', app.data.params.loc.enable_divergence_check);
    app.ui.LocQCDivergence.Layout.Column = [1 2];

    uilabel(g, 'Text', 'Max Shift Factor:');
    app.ui.LocShiftFactor = uieditfield(g, 'numeric', ...
        'Value', app.data.params.loc.qc_max_shift_factor);

    app.ui.LocQCRoiMaxima = uicheckbox(g, 'Text', 'Enable ROI Maxima Check', ...
        'Value', app.data.params.loc.enable_roi_maxima_check);
    app.ui.LocQCRoiMaxima.Layout.Column = [1 2];

    % ---- Gaussian Fit QC panel (only for Gaussian methods) ----
    app.ui.p_loc_qc_gauss = uipanel(parentGrid, 'Title', 'Gaussian Fit QC');
    app.ui.p_loc_qc_gauss.Layout.Row = 4;
    g2 = uigridlayout(app.ui.p_loc_qc_gauss, [3, 2]);
    g2.ColumnWidth = {'fit', '1x'};
    g2.RowHeight = {'fit', 'fit', 'fit'};

    uilabel(g2, 'Text', 'FWHM [x z] (px):');
    app.ui.LocFWHM = uieditfield(g2, 'text', ...
        'Value', mat2str(app.data.params.loc.fwhm));

    uilabel(g2, 'Text', 'Min R-squared:');
    app.ui.GaussMinRSquared = uieditfield(g2, 'numeric', ...
        'Value', app.data.params.loc.min_r_squared);
end

function app = buildTrackTab(app)
    reg = getAlgorithmRegistry();
    app.ui.tabTrack = uitab(app.ui.tabGroup, 'Title', '4. Track');
    
    % 1. MAIN TAB LAYOUT: 2 Rows (Scrollable content area + Pinned Run Button)
    g_tab = uigridlayout(app.ui.tabTrack, [2, 1]);
    g_tab.RowHeight = {'1x', 45};
    g_tab.ColumnWidth = {'1x'};
    g_tab.Padding = [5 5 5 5];
    g_tab.RowSpacing = 10;

    % 2. INNER GRID (SCROLLABLE) - 6 rows only, no run button here
    g = uigridlayout(g_tab, [6, 1], 'Scrollable', 'on');
    g.Layout.Row = 1;
    g.Layout.Column = 1;
    g.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    g.Padding = [5 5 20 5]; % Extra right padding for scrollbar visibility

    uilabel(g, 'Text', 'Tracking Method:', 'FontWeight', 'bold');
    
    % Safely find the matching item ignoring case
    trackItems = registryIds(reg.track);
    matchIdx = find(strcmpi(trackItems, app.data.params.track.method), 1);
    if isempty(matchIdx)
        matchIdx = 1; % Fallback to the first item if not found at all
    end
    
    app.ui.TrackMethodDropdown = uidropdown(g, ...
        'Items', trackItems, ...
        'Value', trackItems{matchIdx});

    % Row 3: Core params
    app = buildTrackCoreParams(app, g);

    % Row 4: Kalman settings
    app = buildKalmanPanel(app, g);

    % Row 5: QC filters
    app = buildTrackQCPanel(app, g);

    % Row 6: Kalman gain diagnostics (NEW)
    app = buildKalmanGainPanel(app, g);

    % 3. RUN TRACKING BUTTON: Pinned permanently to the bottom of the tab (Row 2)
    app.ui.RunTrackingButton = uibutton(g_tab, 'Text', 'Run Tracking', ...
        'ButtonPushedFcn', @(s,e) runTracking(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunTrackingButton.Layout.Row = 2;
    app.ui.RunTrackingButton.Layout.Column = 1;
end

function app = buildTrackCoreParams(app, parentGrid)
    p = uipanel(parentGrid, 'Title', 'Core Tracking Params');
    p.Layout.Row = 3;
    g = uigridlayout(p, [3, 3]);
    g.ColumnWidth = {'fit', '1x', '0.5x'};
    g.RowHeight = {50, 50, 50};

    uilabel(g, 'Text', 'Max Link Dist (px):');
    app.ui.MaxDistSlider = uislider(g, 'Limits', [0.1 10], ...
        'Value', app.data.params.track.max_linking_distance);
    app.ui.MaxDistField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.max_linking_distance, ...
        'ValueDisplayFormat', '%.1f');
    syncSliderField(app.ui.MaxDistSlider, app.ui.MaxDistField);

    uilabel(g, 'Text', 'Max Gap Frames:');
    app.ui.GapFramesSlider = uislider(g, 'Limits', [0 10], ...
        'Value', app.data.params.track.max_gap_closing_frames);
    app.ui.GapFramesField = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.max_gap_closing_frames, ...
        'RoundFractionalValues', 'on');
    syncSliderField(app.ui.GapFramesSlider, app.ui.GapFramesField);

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
    g = uigridlayout(app.ui.p_kalman_adv, [6, 3]);
    g.ColumnWidth = {'fit', '1x', '0.5x'};
    g.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};

    uilabel(g, 'Text', 'Model:');
    app.ui.KalmanModelDrop = uidropdown(g, ...
        'Items', {'ConstantVelocity', 'ConstantAcceleration'}, ...
        'Value', app.data.params.track.kalman.motion_model);
    app.ui.KalmanModelDrop.Layout.Column = [2 3];

    uilabel(g, 'Text', 'Process Noise:');
    app.ui.KalmanNoise = uieditfield(g, 'numeric', ...
        'Value', app.data.params.track.kalman.process_noise, ...
        'ValueChangedFcn', @(s,e) updateKalmanGainSummary(app.fig));
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

    % Hierarchical Kalman modal (for Kalman_Advanced)
    app.ui.BtnConfigHK = uibutton(g, ...
        'Text', 'Configure Hierarchical Kalman (HK)', ...
        'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 1], ...
        'ButtonPushedFcn', @(s,e) openHKConfigGUI(app.fig));
    app.ui.BtnConfigHK.Layout.Row = 5;
    app.ui.BtnConfigHK.Layout.Column = [1 3];

    % Advanced Kalman modal (cost-matrix weights, gating, HK noise)
    app.ui.BtnAdvKalman = uibutton(g, ...
        'Text', 'Advanced Kalman Parameters...', ...
        'ButtonPushedFcn', @(s,e) openAdvancedKalmanGUI(app.fig));
    app.ui.BtnAdvKalman.Layout.Row = 6;
    app.ui.BtnAdvKalman.Layout.Column = [1 3];
end

function app = buildTrackQCPanel(app, parentGrid)
    p = uipanel(parentGrid, 'Title', 'Track QC Filters');
    p.Layout.Row = 5;
    g = uigridlayout(p, [3, 3]);
    g.ColumnWidth = {'fit', '1x', '0.5x'};
    g.RowHeight = {'fit', 'fit', 'fit'};

    app.ui.TrackQCDirection = uicheckbox(g, 'Text', 'Direction Constraint', ...
        'Value', app.data.params.track.qc.enable_direction_constraint);
    uilabel(g, 'Text', 'Max Angle (deg):', 'HorizontalAlignment', 'right');
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

% -------------------------------------------------------------------------
% Kalman Trust Balance panel.  Shows the theoretical split between
% trusting the motion model (algorithm) vs. trusting the localizations,
% computed from the user's current noise parameters — no tracking needed.
%
%   Standard Kalman:      K = Q / (Q + R),   Q = process_noise,
%                                             R = (mean(FWHM)/2.355)^2
%   Hierarchical Kalman:  K = (alpha*v_max) / (alpha*v_max + beta)
%
%   K -> 0  =>  trusts the motion model  (smooth, predictive)
%   K -> 1  =>  trusts the localizations (follows raw data)
% -------------------------------------------------------------------------
function app = buildKalmanGainPanel(app, parentGrid)
    app.ui.p_kgain = uipanel(parentGrid, 'Title', 'Kalman Trust Balance');
    app.ui.p_kgain.Layout.Row = 6;

    g = uigridlayout(app.ui.p_kgain, [3, 1]);
    g.RowHeight = {'fit', 28, 'fit'};
    g.Padding = [8 4 8 4];
    g.RowSpacing = 4;

    % Row 1: K value + formula
    app.ui.KGainSummary = uilabel(g, ...
        'Text', 'K = —', ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');

    % Row 2: Colored split bar  (Model | Localizations)
    gBar = uigridlayout(g, [1, 2]);
    gBar.ColumnWidth = {'5x', '5x'};   % 50/50 default; updated dynamically
    gBar.Padding = [0 0 0 0];
    gBar.ColumnSpacing = 2;
    gBar.RowHeight = {26};

    app.ui.KGainModelBar = uilabel(gBar, ...
        'Text', 'Model 50%', ...
        'BackgroundColor', [0.25 0.65 0.35], ...
        'FontColor', 'w', 'FontWeight', 'bold', 'FontSize', 11, ...
        'HorizontalAlignment', 'center');

    app.ui.KGainMeasBar = uilabel(gBar, ...
        'Text', 'Localizations 50%', ...
        'BackgroundColor', [0.25 0.48 0.85], ...
        'FontColor', 'w', 'FontWeight', 'bold', 'FontSize', 11, ...
        'HorizontalAlignment', 'center');

    app.ui.KGainBarGrid = gBar;  % store handle for dynamic column update

    % Row 3: Formula description
    app.ui.KGainFormula = uilabel(g, ...
        'Text', 'K = Q / (Q + R)', ...
        'FontAngle', 'italic', 'FontSize', 10, ...
        'FontColor', [0.35 0.35 0.35], ...
        'HorizontalAlignment', 'center', ...
        'WordWrap', 'on');
end

function app = buildPostProcessTab(app)
    app.ui.tabPostProcess = uitab(app.ui.tabGroup, 'Title', '5. Post-Process');
    app.ui.panel_post = createScrollableTabPanel(app.ui.tabPostProcess);

    g = uigridlayout(app.ui.panel_post, [3, 1]);
    g.RowHeight = {'fit', 'fit', 'fit'};

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

    app.ui.RunPostProcessButton = uibutton(g, ...
        'Text', 'Run Post-Processing (Smoothing)', ...
        'ButtonPushedFcn', @(s,e) runPostProcessing(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 1 0.6]);
    app.ui.RunPostProcessButton.Layout.Row = 3;
end

function app = buildRenderTab(app)
    reg = getAlgorithmRegistry();
    app.ui.tabRender = uitab(app.ui.tabGroup, 'Title', '6. Render');
    app.ui.panel_render = createScrollableTabPanel(app.ui.tabRender);

    g = uigridlayout(app.ui.panel_render, [4, 1]);
    g.RowHeight = {'fit', 'fit', 'fit', 'fit'};

    p_settings = uipanel(g, 'Title', 'Rendering Settings');
    p_settings.Layout.Row = 1;
    g_set = uigridlayout(p_settings, [2, 2]);
    g_set.RowHeight = {'fit', 'fit'};

    uilabel(g_set, 'Text', 'Upsampling Factor:');
    app.ui.UpsamplingField = uieditfield(g_set, 'numeric', ...
        'Value', app.data.params.render.upsampling_factor);

    uilabel(g_set, 'Text', 'Render Method:');
    app.ui.RenderMethodDrop = uidropdown(g_set, ...
        'Items', registryIds(reg.render), ...
        'Value', app.data.params.render.method);

    % Row 2: Advanced
    app.ui.BtnAdvRender = uibutton(g, ...
        'Text', 'Advanced Rendering / Smoothing / Analysis Parameters...', ...
        'ButtonPushedFcn', @(s,e) openAdvancedRenderGUI(app.fig));
    app.ui.BtnAdvRender.Layout.Row = 2;

    % Row 3: Generate
    app.ui.GenerateImagesButton = uibutton(g, ...
        'Text', 'Generate & Display Final Images (New Windows)', ...
        'ButtonPushedFcn', @(s,e) runRendering(app.fig), ...
        'FontWeight', 'bold', 'BackgroundColor', [0.6 0.8 1]);
    app.ui.GenerateImagesButton.Layout.Row = 3;

    uilabel(g, 'Text', 'Note: This will open 4 separate figure windows.', ...
        'FontAngle', 'italic');
end

% =========================================================================
% --- CROP CALLBACKS (Pre-SVD Data Reduction) -----------------------------
% =========================================================================

function runInteractiveCrop(fig)
    app = guidata(fig);

    if isempty(app.data.rawData)
        uialert(fig, 'Please load Raw Data first to define crop.', 'Error');
        return;
    end

    f = figure('Name', 'Interactive Crop - Draw Rectangle', 'Tag', 'ULM_InteractiveCrop');

    if ~isempty(app.data.filteredData)
        meanImg = mean(abs(app.data.filteredData), 3);
        dispTitle = 'Draw a rectangle to define the crop area (Filtered Data).';
    else
        meanImg = mean(abs(app.data.rawData), 3);
        dispTitle = 'Draw a rectangle to define the crop area (Raw Data).';
    end

    imagesc(meanImg); colormap gray; axis image;
    title([dispTitle, ' Close window to cancel.']);

    rect = drawrectangle('Label', 'Crop ROI', 'Color', 'r');
    wait(rect);

    if isvalid(rect)
        pos = round(rect.Position);
        app.ui.CropBoxField.Value = mat2str(pos);
        if ~isfield(app.data.params, 'io')
            app.data.params.io = struct();
        end
        app.data.params.io.cropBox = pos;
        guidata(fig, app);
    end

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
        pos = str2num(app.ui.CropBoxField.Value); %#ok<ST2NM>
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

        selection = uiconfirm(fig, ...
            sprintf('Cropping will reduce data size from [%dx%d] to [%dx%d]. This action is permanent for the current session and will clear Undo history. Proceed?', ...
            W, H, x_end-x_start+1, y_end-y_start+1), ...
            'Confirm Spatial Crop', ...
            'Options', {'Proceed', 'Cancel'}, ...
            'DefaultOption', 1, 'CancelOption', 2);

        if strcmp(selection, 'Cancel')
            return;
        end

        showProgress(fig, 'Applying spatial crop to raw data...', true);

        app.data.rawData = app.data.rawData(y_start:y_end, x_start:x_end, :);

        app.data.rawDataHash = DataHash(app.data.rawData);
        app.data.U = []; app.data.S_diag = []; app.data.V = [];
        app.data.svdDims = [];
        app.data.mask = [];
        app.ui.maskStatusLabel.Text = 'Status: None';
        app.data.vesselMap = [];
        app.data.baseVesselMap = [];

        app.data = clearDownstreamData(app.data, 0);
        app = manageGUIState(app, 0);

        app.undoManager.clear();
        app.ui.btnUndo.Enable = 'off';
        app.ui.btnRedo.Enable = 'off';

        if ~isfield(app.data.params, 'io')
            app.data.params.io = struct();
        end
        app.data.params.io.cropBox = [x_start, y_start, x_end-x_start+1, y_end-y_start+1];

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
    pos = str2num(app.ui.CropBoxField.Value); %#ok<ST2NM>
    if isempty(pos) || length(pos) ~= 4
        uialert(fig, 'Invalid crop box to save.', 'Error');
        return;
    end
    [file, path] = uiputfile('cropBox.mat', 'Save Crop Box As...');
    if file ~= 0
        cropBox = pos; %#ok<NASGU>
        save(fullfile(path, file), 'cropBox');
        uialert(fig, 'Crop box saved successfully.', 'Success');
    end
end

function loadCrop(fig)
    app = guidata(fig);
    [file, path] = uigetfile('*.mat', 'Load Crop Box');
    if file ~= 0
        data = load(fullfile(path, file));
        if isfield(data, 'cropBox')
            app.ui.CropBoxField.Value = mat2str(data.cropBox);
            if ~isfield(app.data.params, 'io')
                app.data.params.io = struct();
            end
            app.data.params.io.cropBox = data.cropBox;
            guidata(fig, app);
        else
            uialert(fig, 'No "cropBox" variable found in the selected file.', 'Error');
        end
    end
end

% =========================================================================
% --- MAIN PROCESSING CALLBACKS -------------------------------------------
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

        showProgress(fig, 'Loading data...', true);

        fprintf('Loading %s...\n', fullfile(path, file));
        data = load(fullfile(path, file));

        fields = fieldnames(data);
        data_matrix = [];
        variable_name = '';

        fprintf('  Scanning file contents...\n');
        for i = 1:length(fields)
            field_content = data.(fields{i});
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

        if ~isreal(data_matrix)
            fprintf('  Data type: IQ (Complex) - will use abs() for display\n');
            app.data.isIQData = true;
        else
            fprintf('  Data type: Real (already processed)\n');
            app.data.isIQData = false;
        end

        app.data.rawData = data_matrix;
        app.data.rawDataVariableName = variable_name;
        [H, W, T] = size(app.data.rawData);

        fprintf('  Dimensions: %d x %d x %d frames\n', H, W, T);
        fprintf('  Data range: [%.2e, %.2e]\n', min(abs(data_matrix(:))), max(abs(data_matrix(:))));

        app.data.rawDataHash = DataHash(app.data.rawData);

        app.data.U = []; app.data.S_diag = []; app.data.V = [];
        app.data.svdDims = [];

        app.state.maxFrame = T;
        app.state.currentFrame = 1;
        app.data.mask = [];
        app.ui.maskStatusLabel.Text = 'Status: None';

        abs_data = abs(app.data.rawData(:));
        app.data.rawClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
        if app.data.rawClim(1) == app.data.rawClim(2)
            app.data.rawClim(2) = app.data.rawClim(1) + 1;
        end

        app.ui.FrameSlider.Limits = [1 T];
        app.ui.FrameSlider.Value = 1;
        app.ui.FrameField.Limits = [1 T];
        app.ui.FrameField.Value = 1;

        app.data.params.io.data_folder = path;

        app.data = clearDownstreamData(app.data, 0);
        app = manageGUIState(app, 0);

        saveParamState(fig, 'data_loaded');

        setStatus(app, sprintf('Loaded: %s (%dx%dx%d)', variable_name, H, W, T), 'green');
        guidata(fig, app);
        fprintf('  Data loaded successfully\n');

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
    startVal = app.ui.SVDCutoffStart.Value;
    endVal   = app.ui.SVDCutoffEnd.Value;
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

        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;

        saveParamState(fig, 'filter');

        rawData = app.data.rawData;
        currentHash = DataHash(rawData);
        filteredData = [];

        switch params.filter.method
            case 'svd_filter'
                cacheKey = [currentHash '_svd_filter'];
                if isempty(app.data.U) || isempty(app.data.S_diag) || ...
                   ~strcmp(app.data.rawDataHash, cacheKey)

                    showProgress(fig, 'Computing SVD (first time only)...', true);
                    [app.data.U, app.data.S_diag, app.data.V, app.data.svdDims] = ...
                        run_SVD_Decomposition(rawData);
                    app.data.rawDataHash = cacheKey;
                    fprintf('  SVD computed and cached.\n');
                else
                    fprintf('  Using cached SVD (no recomputation needed).\n');
                end

                showProgress(fig, 'Reconstructing filtered signal...', false);
                cutoff = params.filter.svd_cutoff;
                filteredData = reconstruct_SVD_Signal(app.data.U, ...
                    app.data.S_diag, app.data.V, app.data.svdDims, cutoff);

            case 'svd_ssm'
                filteredData = SVD_SSM(rawData, 'IndentPrefix', '  ');

            case 'dcc_svd'
                showProgress(app.fig, 'Computing DCC-SVD...', true);

                [filteredData, dccInfo] = DCC_SVD(rawData, params.acq.framerate, ...
                    'ReconstructionMode', 'blood', ...
                    'DensityPercentile',  10, ...
                    'CanopySeparation',   2.0, ...
                    'PlotResults',        true, ...
                    'IndentPrefix',       '  ');

                app.data.U       = dccInfo.U;
                app.data.S_diag  = dccInfo.singular_values;
                app.data.V       = dccInfo.V;
                app.data.svdDims = size(rawData);
                app.data.rawDataHash = [currentHash '_dcc_svd'];

                app.data.tissue_indices = dccInfo.tissue_indices;
                app.data.blood_indices  = dccInfo.blood_indices;
                app.data.noise_indices  = dccInfo.noise_indices;

                fprintf('  DCC clusters - Tissue: %d, Blood: %d, Noise: %d components.\n', ...
                    numel(dccInfo.tissue_indices), ...
                    numel(dccInfo.blood_indices), ...
                    numel(dccInfo.noise_indices));

            case 'svd_blockwise'
                showProgress(fig, 'Running block-wise SVD (Song 2017)...', true);
                bw = params.filter.blockwise;
                opts = { ...
                    'ThresholdMethod',       bw.threshold_method, ...
                    'OverlapPct',            bw.overlap_pct, ...
                    'ManualCutoff',          bw.manual_cutoff, ...
                    'TissueFreqThreshHz',    bw.tissue_freq_hz, ...
                    'MPDeviationSigma',      bw.mp_deviation_sigma, ...
                    'GradientInflectionPct', bw.gradient_pct, ...
                    'MinBloodComponents',    bw.min_blood_comps, ...
                    'MaxTissueFraction',     bw.max_tissue_frac, ...
                    'PlotThresholdMaps',     bw.plot_maps, ...
                    'IndentPrefix',          '  '};
                if isnumeric(bw.block_size_mm) && isscalar(bw.block_size_mm) && bw.block_size_mm > 0
                    opts = [opts, {'BlockSizeMm', [1 1] * bw.block_size_mm}];
                end

                [filteredData, bwDiag] = SVD_blockwise(rawData, params, opts{:});
                app.data.blockwiseDiag = bwDiag;
                fprintf('  Block-wise SVD: %d blocks, %.2f s elapsed.\n', ...
                    bwDiag.n_blocks, bwDiag.elapsed_sec);

            otherwise
                error('Unknown filter method: %s', params.filter.method);
        end

        % Butterworth (optional) - works on IQ or real
        if params.filter.enable_butterworth
            showProgress(fig, 'Applying Butterworth filter...', false);
            filteredData = Butterworth_bandpass_filter(filteredData, ...
                params.filter.butter_cutoff, params.acq.framerate, ...
                params.filter.butter_order);
        end

        % Envelope before spatial filters (needed for non-linear ops)
        filteredData = abs(filteredData);

        if ~strcmp(params.filter.spatial_method, 'None')
            showProgress(fig, sprintf('Applying %s Spatial Filter...', params.filter.spatial_method), false);
            filteredData = applySpatialFilter(filteredData, params.filter);
        end

        app.data.filteredData = filteredData;

        abs_data = app.data.filteredData(:);
        app.data.filteredClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
        if app.data.filteredClim(1) == app.data.filteredClim(2)
            app.data.filteredClim(2) = app.data.filteredClim(1) + 1;
        end

        mean_bg = mean(app.data.filteredData, 3);
        app.data.filteredMeanBG = mean_bg .^ 0.5;
        app.data.MeanBG = mean(abs(rawData), 3) .^ 0.5;

        bg_clim_data = app.data.filteredMeanBG(:);
        app.data.filteredBGClim = [prctile(bg_clim_data, 1), prctile(bg_clim_data, 99)];

        app.data.baseVesselMap = [];
        app.data.vesselMap = [];
        app.state.isROIPreview = false;

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

    % Use committed mask, fall back to preview mask, then full-image default
    if isfield(app.data.params, 'proc') && ...
            isfield(app.data.params.proc, 'ROIMask') && ...
            ~isempty(app.data.params.proc.ROIMask)
        app.data.mask = app.data.params.proc.ROIMask;
    elseif isempty(app.data.mask)
        uialert(fig, 'No ROI Mask defined. Using full image.', 'Warning');
        [h, w, ~] = size(app.data.filteredData);
        app.data.mask = true(h, w);
    end

    try
        showProgress(fig, 'Detecting bubbles...', true);

        app.data.params = updateParamsFromGUI(app);
        saveParamState(fig, 'detection');

        app.data.params.loc.DetectMethod       = app.ui.DetectMethodDropdown.Value;
        app.data.params.loc.NP_alpha0          = app.ui.NP_AlphaField.Value;
        app.data.params.loc.crosscor_threshold = app.ui.NCC_ThreshField.Value;

        roiMask = app.data.mask;

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
                        app.state.isProcessing = false;
                        guidata(fig, app);
                        return;
                    end

                    fwhm    = app.data.params.loc.fwhm;
                    sz      = app.data.params.loc.psf_size;
                    sigma_x = fwhm(1) / 2.355;
                    sigma_z = fwhm(2) / 2.355;
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

        dataToLocalize   = app.data.filteredData;
        candidateBubbles = app.data.candidateBubbles;

        if isempty(dataToLocalize)
            error('Filtered data is empty. Please run filtering first.');
        end

        diary_file = [tempname '.txt'];
        diary(diary_file);
        try
            switch params.loc.method
                case 'radial'
                    locs = localizeRadialSymmetry(dataToLocalize, candidateBubbles, params.loc, '');
                case 'gaussian_fit'
                    locs = fit2DGaussian(dataToLocalize, candidateBubbles, params.loc, '');
                case 'gaussian_fit_fast'
                    locs = fit2DGaussianSafe(dataToLocalize, candidateBubbles, params.loc, '');
                otherwise
                    error('Unknown localization method: %s', params.loc.method);
            end
            diary off;
        catch ME_loc
            diary off;
            if exist(diary_file, 'file'), delete(diary_file); end
            rethrow(ME_loc);
        end

        fid = fopen(diary_file, 'r');
        if fid ~= -1
            qc_output = fread(fid, '*char')';
            fclose(fid);
            delete(diary_file);
        else
            qc_output = '';
        end

        if isempty(locs)
            warning('No localizations found. Check parameters.');
        end

        app.data.localizations = locs;
        app.data = clearDownstreamData(app.data, 3);
        app = manageGUIState(app, 3);

        guidata(fig, app);
        hideProgress(fig);

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
    num_candidates = height(candidateBubbles);
    batch_size = 3000;

    if num_candidates <= batch_size
        % --- Slim the data: send only needed frames to parfor workers ---
        batch_frames = unique(candidateBubbles.Frame);
        batch_data   = filteredData(:,:,batch_frames);
        remapped     = candidateBubbles;
        [~, remapped.Frame] = ismember(candidateBubbles.Frame, batch_frames);

        locs = fit2DGaussian_Fast(batch_data, remapped, locParams, indent_prefix);

        % Restore original frame numbers
        if ~isempty(locs) && height(locs) > 0
            locs.Frame = batch_frames(locs.Frame);
        end
        return;
    end

    fprintf('Processing %d candidates in batches of %d...\n', num_candidates, batch_size);
    num_batches = ceil(num_candidates / batch_size);
    locs_cell = cell(num_batches, 1);

    for i = 1:num_batches
        start_idx = (i-1) * batch_size + 1;
        end_idx   = min(i * batch_size, num_candidates);

        fprintf('  Batch %d/%d (candidates %d-%d)...\n', i, num_batches, start_idx, end_idx);

        batch_candidates = candidateBubbles(start_idx:end_idx, :);

        % --- Slim the data: extract only frames this batch needs ---
        batch_frames = unique(batch_candidates.Frame);
        batch_data   = filteredData(:,:,batch_frames);
        [~, batch_candidates.Frame] = ismember(batch_candidates.Frame, batch_frames);

        locs_cell{i} = fit2DGaussian_Fast(batch_data, batch_candidates, locParams, '    ');

        % Restore original frame numbers in the output
        if ~isempty(locs_cell{i}) && height(locs_cell{i}) > 0
            locs_cell{i}.Frame = batch_frames(locs_cell{i}.Frame);
        end

        if i < num_batches
            pause(0.1);
        end
    end

    locs = vertcat(locs_cell{:});
    fprintf('Total localizations: %d\n', height(locs));
end

function showQCDialog(fig, titleStr, qc_text)
    d = uifigure('Name', titleStr, 'Position', [100 100 640 440], 'Tag', 'ULM_QCDialog');

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
    uilabel(btnGrid, 'Text', '');
    uibutton(btnGrid, 'Text', 'OK', ...
        'ButtonPushedFcn', @(~,~) close(d));
    uilabel(btnGrid, 'Text', '');
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
            case 'kalman_v2'
                tracks = trackKalman_v2(app.data.localizations, params, '  ');
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
        updateKalmanGainSummary(fig);

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
            raw = app.data.tracks_raw;
            numRaw = numel(raw);
            for k = 1:numRaw
                raw(k).original_length         = raw(k).length;
                raw(k).velocities_mm_s         = zeros(raw(k).length, 1);
                raw(k).average_velocity_mm_s   = 0;
            end
            app.data.tracks_final = raw;
        else
            raw_tracks       = app.data.tracks_raw;
            smoothing_factor = params.track.smoothing_factor;
            interp_step      = params.render.interpolation_step;

            numTracks = length(raw_tracks);
            processed_tracks_cell = cell(1, numTracks);

            if ~isempty(gcp('nocreate'))
                parfor i = 1:numTracks
                    processed_tracks_cell{i} = processTrack(raw_tracks(i), ...
                        smoothing_factor, interp_step, params);
                end
            else
                for i = 1:numTracks
                    processed_tracks_cell{i} = processTrack(raw_tracks(i), ...
                        smoothing_factor, interp_step, params);

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
    path_to_process = double(track.path);
    act_win = min(smoothing_factor, track.length);
    if mod(act_win, 2) == 0
        act_win = act_win - 1;
    end

    poly_order = ULM_Constants.SGOLAY_POLY_ORDER;
    min_window = poly_order + 2;

    if track.length >= min_window && act_win >= min_window
        path_to_process = [sgolayfilt(track.path(:,1), poly_order, act_win), ...
                          sgolayfilt(track.path(:,2), poly_order, act_win)];
    elseif track.length > 3
        path_to_process = [movmean(track.path(:,1), min(3, track.length)), ...
                          movmean(track.path(:,2), min(3, track.length))];
    end

    orig_inds   = 1:track.length;
    interp_inds = 1:interp_step:track.length;
    frames_interp = interp1(orig_inds, track.frames, interp_inds, 'linear');
    x_interp = fillmissing(interp1(orig_inds, path_to_process(:,1), interp_inds, 'pchip'), ...
        'linear', 'EndValues', 'nearest');
    y_interp = fillmissing(interp1(orig_inds, path_to_process(:,2), interp_inds, 'pchip'), ...
        'linear', 'EndValues', 'nearest');

    new_path = [x_interp', y_interp'];
    new_len  = size(new_path, 1);

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

    if isempty(app.data.tracks_final)
        errordlg('No processed tracks. Please run post-processing first.', 'Error');
        return;
    end

    app.state.isProcessing = true;
    guidata(fig, app);

    try
        showProgress(fig, 'Generating maps...', true);

        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;

        tracks  = app.data.tracks_final;
        min_len = app.ui.DisplayMinLengthField.Value;
        tracks  = tracks([tracks.original_length] >= min_len);

        if isempty(tracks)
            error('No tracks meet minimum length criterion.');
        end

        upscale = params.render.upsampling_factor;
        [H, W, ~] = size(app.data.rawData);
        H_SR = round(H * upscale);
        W_SR = round(W * upscale);

        allX = round(vertcat_field(tracks, 'path', 1) * upscale);
        allY = round(vertcat_field(tracks, 'path', 2) * upscale);
        allV = vertcat_field(tracks, 'velocities_mm_s');

        valid = allX >= 1 & allX <= W_SR & allY >= 1 & allY <= H_SR;
        allX  = allX(valid);
        allY  = allY(valid);
        allV  = allV(valid);

        inds = sub2ind([H_SR, W_SR], allY, allX);
        N    = H_SR * W_SR;

        densityMap    = reshape(accumarray(inds, 1,    [N, 1], @sum, 0), H_SR, W_SR);
        velocityAccum = reshape(accumarray(inds, allV, [N, 1], @sum, 0), H_SR, W_SR);

        velocityMap = zeros(size(densityMap));
        mask = densityMap > 0;
        velocityMap(mask) = velocityAccum(mask) ./ densityMap(mask);

        res = params.render.upsampling_factor;
        res_pts  = [1, 3, 5];
        font_pts = [6, 8, 10];
        line_pts = [2, 4, 5];
        res_c = max(min(res, res_pts(end)), res_pts(1));
        fs = round(interp1(res_pts, font_pts, res_c));
        lw = round(interp1(res_pts, line_pts, res_c), 1);

        hideProgress(fig);
        drawnow;

        createDensityFigure(densityMap, min_len, params, lw, fs);
        createVelocityFigures(velocityMap, params, lw, fs);
        createCombinedFigure(velocityMap, densityMap, params, lw, fs);

        setStatus(app, 'Rendering complete', 'green');
        guidata(fig, app);

    catch ME
        hideProgress(fig);
        errordlg(sprintf('%s\n\n%s', ME.message, ME.getReport('basic')), 'Render Error');
        setStatus(app, 'Rendering failed', 'red');
    end

    app = guidata(fig);
    app.state.isProcessing = false;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function out = vertcat_field(tracks, fieldName, colIdx)
    parts = cell(numel(tracks), 1);
    for k = 1:numel(tracks)
        v = tracks(k).(fieldName);
        if nargin >= 3
            v = v(:, colIdx);
        end
        parts{k} = v(:);
    end
    out = vertcat(parts{:});
end

% =========================================================================
% --- RENDERING HELPER FUNCTIONS ------------------------------------------
% =========================================================================

function enforceMinimumSize(fig)
% ENFORCEMINIMUMSIZE  Prevents the figure from being resized below usable limits.
    MIN_W = 900;
    MIN_H = 600;
    pos = fig.Position;
    if pos(3) < MIN_W || pos(4) < MIN_H
        fig.Position = [pos(1), pos(2), max(pos(3), MIN_W), max(pos(4), MIN_H)];
    end
end

function showLocalizationDensity(fig)
    app = guidata(fig);

    if isempty(app.data.localizations) || height(app.data.localizations) == 0
        uialert(fig, 'No localizations available. Please run Localization first.', 'No Data');
        return;
    end

    if isempty(app.data.rawData)
        uialert(fig, 'Raw data not available - cannot determine image dimensions.', 'Error');
        return;
    end

    try
        app.data.params = updateParamsFromGUI(app);
        params = app.data.params;

        locs = app.data.localizations;

        upscale = params.render.upsampling_factor;
        [H, W, ~] = size(app.data.rawData);
        H_SR = round(H * upscale);
        W_SR = round(W * upscale);

        allX = round(locs.X * upscale);
        allY = round(locs.Y * upscale);

        valid = allX >= 1 & allX <= W_SR & allY >= 1 & allY <= H_SR;
        allX = allX(valid);
        allY = allY(valid);

        n_total   = height(locs);
        n_in_grid = numel(allX);
        if n_in_grid == 0
            uialert(fig, 'No localizations fall within the image grid after upsampling.', 'Empty Map');
            return;
        end

        inds = sub2ind([H_SR, W_SR], allY, allX);
        N    = H_SR * W_SR;
        densityMap = reshape(accumarray(inds, 1, [N, 1], @sum, 0), H_SR, W_SR);

        res      = params.render.upsampling_factor;
        res_pts  = [1, 3, 5];
        font_pts = [6, 8, 10];
        line_pts = [2, 4, 5];
        res_c    = max(min(res, res_pts(end)), res_pts(1));
        fs       = round(interp1(res_pts, font_pts, res_c));
        lw       = round(interp1(res_pts, line_pts, res_c), 1);

        figure('Name', 'Localization Density Map (pre-tracking)', 'Tag', 'ULM_LocDensity');
        d_proc = densityMap .^ (1/3);
        imshow(d_proc, []);
        colormap(hot);
        if any(d_proc(:) > 0)
            clim([0 prctile(d_proc(d_proc > 0), 99.5)]);
        end
        title(sprintf('Localization Density (N = %d points, %.1f%% in grid)', ...
            n_in_grid, 100 * n_in_grid / n_total));
        colorbar;
        add_scale_bar(params, size(densityMap), lw, fs);

        setStatus(app, sprintf('Localization density rendered (%d points)', n_in_grid), 'green');

    catch ME
        uialert(fig, sprintf('Failed to render localization density:\n\n%s', ME.message), 'Render Error');
        setStatus(app, 'Localization density render failed', 'red');
    end
end

function createDensityFigure(densityMap, min_len, params, lw, fs)
    figure('Name', 'Density Map', 'Tag', 'ULM_DensityMap');
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

    figure('Name', 'Velocity (Filtered)', 'Tag', 'ULM_VelocityFiltered');
    v_filt = imgaussfilt(velocityMap, 0.6);
    imshow(v_filt, []);
    colormap(cm);
    if any(v_filt(:) > 0)
        clim([0 prctile(v_filt(v_filt>0), 99.5)]);
    end
    title('Velocity Filtered');
    colorbar;
    add_scale_bar(params, size(velocityMap), lw, fs);

    figure('Name', 'Velocity (Raw)', 'Tag', 'ULM_VelocityRaw');
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
    figure('Name', 'Combined (Velocity x Density)', 'Tag', 'ULM_Combined');

    if any(velocityMap(:) > 0)
        v_max = prctile(velocityMap(velocityMap > 0), 99.5);
    else
        v_max = 1;
    end
    v_norm = mat2gray(velocityMap, [0, v_max]);

    d_proc = densityMap .^ (1/3);
    if any(d_proc(:) > 0)
        d_max = prctile(d_proc(d_proc > 0), 99.0);
    else
        d_max = 1;
    end
    d_norm = mat2gray(d_proc, [0, d_max]) .^ 0.7;

    imagesc(v_norm, [0, 1]);
    colormap(jet(256));
    cb = colorbar;
    cb.Ticks = linspace(0, 1, 6);
    cb.TickLabels = arrayfun(@(v) sprintf('%.1f', v), linspace(0, v_max, 6), 'UniformOutput', false);
    ylabel(cb, 'Velocity (mm/s)', 'FontSize', max(fs, 7));
    clim([0, 1]);

    hold on;
    hOv = image(repmat(d_norm, [1, 1, 3]));
    set(hOv, 'AlphaData', 1 - d_norm);

    bg_overlay = zeros([size(densityMap), 3]);
    hBG = image(bg_overlay);
    set(hBG, 'AlphaData', double(densityMap == 0));
    hold off;

    axis image off;
    title('Combined: Velocity (color) x Density (brightness)');
    add_scale_bar(params, size(velocityMap), lw, fs);
end

% =========================================================================
% SECTION: Session Management
% =========================================================================

function saveSession(fig)
    app = guidata(fig);
    [file, path] = uiputfile('*.mat', 'Save Session');
    if isequal(file, 0), return; end
    try
        session = app.sessionManager.createSession(app.data, app.data.params);
        save(fullfile(path, file), 'session', '-v7.3');
        setStatus(app, sprintf('Session saved: %s', file), 'green');
    catch ME
        uialert(app.fig, sprintf('Save failed: %s', ME.message), 'Error');
    end
end

function loadSession(fig)
    app = guidata(fig);
    [file, path] = uigetfile('*.mat', 'Load Session');
    if isequal(file, 0), return; end
    try
        loaded = load(fullfile(path, file));
        if ~isfield(loaded, 'session')
            error('Invalid session file.');
        end

        [restoredData, restoredParams] = app.sessionManager.restoreSession(loaded.session);
        restoredParams = ensureAllParamFields(restoredParams);

        % Separate state fields from data fields
        stateFields = {'currentState', 'currentFrame', 'maxFrame'};
        for i = 1:numel(stateFields)
            if isfield(restoredData, stateFields{i})
                app.state.(stateFields{i}) = restoredData.(stateFields{i});
                restoredData = rmfield(restoredData, stateFields{i});
            end
        end

        if isfield(restoredData, 'params')
            restoredData = rmfield(restoredData, 'params');
        end

        app.data = restoredData;
        app.data.params = restoredParams;

        % Initialise transient fields that buildGUILayout normally sets
        if ~isfield(app.data, 'rawDataHash'),       app.data.rawDataHash = ''; end
        if ~isfield(app.data, 'U'),                  app.data.U = []; end
        if ~isfield(app.data, 'S_diag'),             app.data.S_diag = []; end
        if ~isfield(app.data, 'V'),                  app.data.V = []; end
        if ~isfield(app.data, 'svdDims'),            app.data.svdDims = []; end
        if ~isfield(app.data, 'redOverlayTemplate'), app.data.redOverlayTemplate = []; end
        if ~isfield(app.data, 'vesselMap'),           app.data.vesselMap = []; end
        if ~isfield(app.data, 'baseVesselMap'),       app.data.baseVesselMap = []; end
        if ~isfield(app.data, 'rawClim'),             app.data.rawClim = [0 1]; end
        if ~isfield(app.data, 'filteredClim'),        app.data.filteredClim = [0 1]; end
        if ~isfield(app.data, 'filteredMeanBG'),      app.data.filteredMeanBG = []; end
        if ~isfield(app.data, 'MeanBG'),              app.data.MeanBG = []; end
        if ~isfield(app.data, 'filteredBGClim'),      app.data.filteredBGClim = [0 1]; end
        if ~isfield(app.data, 'blockwiseDiag'),       app.data.blockwiseDiag = []; end
        if ~isfield(app.data, 'tissue_indices'),      app.data.tissue_indices = []; end
        if ~isfield(app.data, 'blood_indices'),       app.data.blood_indices = []; end
        if ~isfield(app.data, 'noise_indices'),       app.data.noise_indices = []; end

        % Recalculate colour limits from restored data
        if ~isempty(app.data.rawData)
            abs_data = abs(app.data.rawData(:));
            app.data.rawClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
            if app.data.rawClim(1) == app.data.rawClim(2)
                app.data.rawClim(2) = app.data.rawClim(1) + 1;
            end
        end
        if isfield(app.data, 'filteredData') && ~isempty(app.data.filteredData)
            abs_filt = abs(app.data.filteredData(:));
            app.data.filteredClim = [prctile(abs_filt, 1), prctile(abs_filt, 99.9)];
            if app.data.filteredClim(1) == app.data.filteredClim(2)
                app.data.filteredClim(2) = app.data.filteredClim(1) + 1;
            end
            mean_bg = mean(abs(app.data.filteredData), 3);
            app.data.filteredMeanBG = mean_bg .^ 0.5;
            bg_clim_data = app.data.filteredMeanBG(:);
            app.data.filteredBGClim = [prctile(bg_clim_data, 1), prctile(bg_clim_data, 99)];
        end
        if ~isempty(app.data.rawData)
            app.data.MeanBG = mean(abs(app.data.rawData), 3) .^ 0.5;
        end

        % Determine correct pipeline stage from data presence
        if ~isempty(app.data.tracks_final),          stage = 5;
        elseif ~isempty(app.data.tracks_raw),        stage = 4;
        elseif ~isempty(app.data.localizations),     stage = 3;
        elseif ~isempty(app.data.candidateBubbles),  stage = 2;
        elseif ~isempty(app.data.filteredData),      stage = 1;
        elseif ~isempty(app.data.rawData),           stage = 0;
        else,                                         stage = -1;
        end

        guidata(fig, app);
        populateGUIFromParams(fig);
        app = guidata(fig);
        app = manageGUIState(app, stage);
        guidata(fig, app);
        displayCurrentFrame(fig);
        setStatus(app, sprintf('Session loaded: %s (stage %d)', file, stage), 'green');
    catch ME
        uialert(app.fig, sprintf('Load failed: %s', ME.message), 'Error');
    end
end

% =========================================================================
% SECTION: Undo / Redo
% =========================================================================

function saveParamState(fig, ~)
    % Second argument is an optional textual tag describing which parameter
    % changed. It is accepted for call-site compatibility but not required.
    app = guidata(fig);
    if ~isfield(app, 'undoManager') || isempty(app.undoManager), return; end
    % Sync the changed control into params BEFORE pushing the snapshot so the
    % undo stack actually reflects the new value rather than the previous one.
    try
        updateParamsFromGUI(fig);
        app = guidata(fig);
    catch
    end
    try
        app.undoManager.push(app.data.params);
    catch
    end
    guidata(fig, app);
    updateUndoRedoButtons(fig);
end

function performUndo(fig)
    app = guidata(fig);
    if ~isfield(app, 'undoManager'), return; end
    if ~app.undoManager.canUndo(), return; end
    try
        app.data.params = app.undoManager.undo(app.data.params);
        app.data.params = ensureAllParamFields(app.data.params);
    catch ME
        warning(ME.identifier, 'Undo failed: %s', ME.message);
        return;
    end
    guidata(fig, app);
    populateGUIFromParams(fig);
    updateUndoRedoButtons(fig);
    setStatus(guidata(fig), 'Undo.', 'blue');
end

function performRedo(fig)
    app = guidata(fig);
    if ~isfield(app, 'undoManager'), return; end
    if ~app.undoManager.canRedo(), return; end
    try
        app.data.params = app.undoManager.redo(app.data.params);
        app.data.params = ensureAllParamFields(app.data.params);
    catch ME
        warning(ME.identifier, 'Redo failed: %s', ME.message);
        return;
    end
    guidata(fig, app);
    populateGUIFromParams(fig);
    updateUndoRedoButtons(fig);
    setStatus(guidata(fig), 'Redo.', 'blue');
end

function updateUndoRedoButtons(fig)
    app = guidata(fig);
    if ~isfield(app, 'undoManager'), return; end
    if isfield(app.ui, 'btnUndo') && isvalid(app.ui.btnUndo)
        try
            app.ui.btnUndo.Enable = matlab.lang.OnOffSwitchState(app.undoManager.canUndo());
        catch
        end
    end
    if isfield(app.ui, 'btnRedo') && isvalid(app.ui.btnRedo)
        try
            app.ui.btnRedo.Enable = matlab.lang.OnOffSwitchState(app.undoManager.canRedo());
        catch
        end
    end
end

% =========================================================================
% SECTION: Status bar, progress bar, memory display
% =========================================================================
% setStatus accepts either an app struct OR a figure handle as the first
% argument so existing 3-arg call sites (setStatus(app,msg,color)) and new
% callbacks that only have `fig` can both use it.

function setStatus(ref, msg, color)
    if nargin < 3 || isempty(color), color = 'black'; end
    if isgraphics(ref)
        app = guidata(ref);
    else
        app = ref;
    end
    if isfield(app.ui, 'lblStatus') && isvalid(app.ui.lblStatus)
        app.ui.lblStatus.Text = msg;
        try
            if ischar(color) || isstring(color)
                switch lower(char(color))
                    case 'green',  app.ui.lblStatus.FontColor = [0.0 0.55 0.0];
                    case 'red',    app.ui.lblStatus.FontColor = [0.75 0.0 0.0];
                    case 'blue',   app.ui.lblStatus.FontColor = [0.0 0.3 0.7];
                    case 'orange', app.ui.lblStatus.FontColor = [0.8 0.45 0.0];
                    otherwise,     app.ui.lblStatus.FontColor = [0.1 0.1 0.1];
                end
            elseif isnumeric(color) && numel(color) == 3
                app.ui.lblStatus.FontColor = color;
            end
        catch
        end
    end
    % Update the StatusLamp color
    if isfield(app.ui, 'StatusLamp') && isvalid(app.ui.StatusLamp)
        switch lower(char(color))
            case 'green',  app.ui.StatusLamp.Color = [0 0.8 0];
            case 'red',    app.ui.StatusLamp.Color = [0.9 0 0];
            case 'blue',   app.ui.StatusLamp.Color = [0 0.4 0.9];
            case 'orange', app.ui.StatusLamp.Color = [1 0.6 0];
            otherwise,     app.ui.StatusLamp.Color = [0.5 0.5 0.5];
        end
    end
    drawnow limitrate;
end

function showProgress(ref, msg, ~)
% SHOWPROGRESS  Create or update a modal indeterminate progress dialog.
%
%   On the FIRST call inside a pipeline run, this creates a uiprogressdlg
%   that is modal over the parent uifigure.  Subsequent calls merely update
%   the Message text, so that each sub-step (e.g. "Computing SVD…",
%   "Applying Butterworth…") appears in the same dialog without flicker.
%
%   The dialog handle is stored via setappdata (key: 'ULM_ModalDlg')
%   rather than in guidata, which avoids any data-integrity conflicts with
%   the app struct that pipeline functions continuously read and write.
%
%   The optional 3rd argument is accepted for backward compatibility and
%   is silently ignored.

    % --- Resolve figure handle ---
    if isgraphics(ref)
        fig = ref;
        app = guidata(ref);
    else
        app = ref;
        fig = app.fig;
    end

    % --- Modal dialog: create or update ---
    dlg = [];
    if isappdata(fig, 'ULM_ModalDlg')
        dlg = getappdata(fig, 'ULM_ModalDlg');
    end

    if ~isempty(dlg) && isobject(dlg) && isvalid(dlg)
        % Dialog already visible — just update the message text
        dlg.Message = msg;
    else
        % First call in this pipeline run — create the dialog
        try
            dlgTitle = deriveProgressTitle(msg);
            dlg = uiprogressdlg(fig, ...
                'Title',         dlgTitle, ...
                'Message',       msg, ...
                'Indeterminate', 'on');      % animated bar, no % needed
            setappdata(fig, 'ULM_ModalDlg', dlg);
        catch
            % Fail silently — the pipeline must not abort because of UI
        end
    end

    % --- Existing status-bar update (non-modal, always runs) ---
    if isfield(app.ui, 'ProgressBar') && ~isempty(app.ui.ProgressBar)
        try
            if isobject(app.ui.ProgressBar) && isvalid(app.ui.ProgressBar)
                app.ui.ProgressBar.Visible = 'on';
                if isprop(app.ui.ProgressBar, 'Text')
                    app.ui.ProgressBar.Text = msg;
                end
            end
        catch
        end
    end
    setStatus(app, msg, 'blue');
end

function hideProgress(ref)
% HIDEPROGRESS  Close the modal progress dialog and re-enable the parent figure.
%
%   This function is safe to call multiple times.  If the dialog has
%   already been closed (or was never created), it silently does nothing.
%   The fail-safe try-catch ensures the parent window is never permanently
%   locked, even if the dialog object is in a bad state.

    % --- Resolve figure handle ---
    if isgraphics(ref)
        fig = ref;
        app = guidata(ref);
    else
        app = ref;
        fig = app.fig;
    end

    % --- Close modal dialog ---
    if isappdata(fig, 'ULM_ModalDlg')
        dlg = getappdata(fig, 'ULM_ModalDlg');
        try
            if ~isempty(dlg) && isobject(dlg) && isvalid(dlg)
                close(dlg);
            end
        catch
            % Force-delete if close() fails
            try delete(dlg); catch, end
        end
        rmappdata(fig, 'ULM_ModalDlg');
    end

    % --- Existing ProgressBar logic (backward compat) ---
    if isfield(app.ui, 'ProgressBar') && ~isempty(app.ui.ProgressBar)
        try
            if isobject(app.ui.ProgressBar) && isvalid(app.ui.ProgressBar)
                app.ui.ProgressBar.Visible = 'off';
            end
        catch
        end
    end
end

function t = deriveProgressTitle(msg)
% DERIVEPROGRESSTITLE  Map a progress message to a human-readable stage title.
%
%   This keeps the dialog window title informative without requiring any
%   change to the calling pipeline functions.  The mapping is based on
%   keywords that already appear in the existing showProgress calls.

    msg_lower = lower(msg);
    if     contains(msg_lower, 'filter') || contains(msg_lower, 'svd') || contains(msg_lower, 'butterworth') || contains(msg_lower, 'dcc')
        t = 'Step 1 — Clutter Filtering';
    elseif contains(msg_lower, 'detect') || contains(msg_lower, 'bubble')
        t = 'Step 2 — Bubble Detection';
    elseif contains(msg_lower, 'localiz')
        t = 'Step 3 — Sub-pixel Localization';
    elseif contains(msg_lower, 'track')
        t = 'Step 4 — Particle Tracking';
    elseif contains(msg_lower, 'post-process') || contains(msg_lower, 'processing track')
        t = 'Step 5 — Track Post-Processing';
    elseif contains(msg_lower, 'render') || contains(msg_lower, 'generat')
        t = 'Step 6 — Rendering Maps';
    elseif contains(msg_lower, 'saving') || contains(msg_lower, 'loading') || contains(msg_lower, 'session')
        t = 'Session I/O';
    else
        t = 'ULM Processing';
    end
end

function startMemoryMonitor(fig)
    app = guidata(fig);
    try
        app.memTimer = timer( ...
            'ExecutionMode', 'fixedSpacing', ...
            'Period', ULM_Constants.MEMORY_UPDATE_INTERVAL, ...
            'TimerFcn', @(~,~) updateMemoryDisplay(fig));
        guidata(fig, app);
        start(app.memTimer);
    catch
    end
end

function updateMemoryDisplay(fig)
    if ~isvalid(fig), return; end
    app = guidata(fig);
    if ~isfield(app.ui, 'lblMemory') || ~isvalid(app.ui.lblMemory), return; end
    try
        if ispc
            [~, sys] = memory;
            used_gb = (sys.PhysicalMemory.Total - sys.PhysicalMemory.Available) / 1e9;
            app.ui.lblMemory.Text = sprintf('Memory: %.1f GB', used_gb);
        else
            % Estimate from the known large data fields
            totalMB = 0;
            largeFields = {'rawData','filteredData','U','S_diag','V', ...
                           'candidateBubbles','localizations', ...
                           'tracks_raw','tracks_final','mask', ...
                           'vesselMap','baseVesselMap','redOverlayTemplate'};
            for k = 1:numel(largeFields)
                if isfield(app.data, largeFields{k}) && ~isempty(app.data.(largeFields{k}))
                    d = app.data.(largeFields{k});
                    if isnumeric(d) || islogical(d)
                        info = whos('d');
                        totalMB = totalMB + info.bytes / 1e6;
                    end
                end
            end
            app.ui.lblMemory.Text = sprintf('Data: %.0f MB', totalMB);
        end
    catch
        app.ui.lblMemory.Text = 'Memory: N/A';
    end
end

function cleanupGUI(fig)
    if ~isvalid(fig), return; end
    app = guidata(fig);
    if isempty(app), delete(fig); return; end

    % ---- 1. Stop & delete ALL timers ----
    % memTimer lives on app directly
    if isfield(app, 'memTimer') && ~isempty(app.memTimer) && isvalid(app.memTimer)
        stop(app.memTimer);
        delete(app.memTimer);
    end
    % displayTimer lives on app directly (created in playback)
    if isfield(app, 'displayTimer') && ~isempty(app.displayTimer) && isvalid(app.displayTimer)
        stop(app.displayTimer);
        delete(app.displayTimer);
    end
    % playbackTimer in app.state
    if isfield(app.state, 'playbackTimer') && ~isempty(app.state.playbackTimer) ...
            && isvalid(app.state.playbackTimer)
        stop(app.state.playbackTimer);
        delete(app.state.playbackTimer);
    end
    % Catch any orphaned MATLAB timers created by this GUI
    allTimers = timerfindall;
    for k = 1:numel(allTimers)
        try
            cb = func2str(allTimers(k).TimerFcn);
            if contains(cb, 'updateMemoryDisplay') || contains(cb, 'timerCallback')
                stop(allTimers(k));
                delete(allTimers(k));
            end
        catch
        end
    end

    % ---- 2. DisplayManager debounce timer ----
    if isfield(app, 'displayManager') && ~isempty(app.displayManager)
        try app.displayManager.cleanup(); catch, end
    end

    % ---- 3. UndoRedoManager history ----
    if isfield(app, 'undoManager') && ~isempty(app.undoManager)
        try app.undoManager.clear(); catch, end
    end

    % ---- 3b. Close modal progress dialog (if still open) ----
    if isappdata(fig, 'ULM_ModalDlg')
        try
            dlg = getappdata(fig, 'ULM_ModalDlg');
            if ~isempty(dlg) && isvalid(dlg), close(dlg); end
        catch
        end
        try rmappdata(fig, 'ULM_ModalDlg'); catch, end
    end

    % ---- 4. Close progress dialog ----
    if isfield(app.ui, 'ProgressBar') && ~isempty(app.ui.ProgressBar)
        try
            if isobject(app.ui.ProgressBar) && isvalid(app.ui.ProgressBar)
                close(app.ui.ProgressBar);
                delete(app.ui.ProgressBar);
            end
        catch
        end
    end

    % ---- 5. Release large data arrays ----
    largeFields = {'rawData', 'filteredData', 'candidateBubbles', ...
                   'localizations', 'tracks_raw', 'tracks_final', ...
                   'U', 'S_diag', 'V', 'mask', 'vesselMap', ...
                   'baseVesselMap', 'redOverlayTemplate', ...
                   'filteredMeanBG', 'MeanBG', 'blockwiseDiag'};
    for k = 1:numel(largeFields)
        if isfield(app.data, largeFields{k})
            app.data.(largeFields{k}) = [];
        end
    end

    % ---- 6. Close any child figures opened by the GUI ----
    allFigs = findall(0, 'Type', 'figure');
    for k = 1:numel(allFigs)
        if allFigs(k) ~= fig && isvalid(allFigs(k))
            tag = get(allFigs(k), 'Tag');
            if startsWith(tag, 'ULM_')
                delete(allFigs(k));
            end
        end
    end

    % ---- 7. Final cleanup ----
    fprintf('ULM GUI closed. All resources released.\n');
    delete(fig);
end

% =========================================================================
% SECTION: GUI State Machine
% =========================================================================
% Stage codes:  -1 initial | 0 raw loaded | 1 filter | 2 detect | 3 loc
%                4 track    | 5 post      | 6 render

function app = manageGUIState(app, stage)
    if nargin < 2, stage = -1; end

    hasRaw  = isfield(app.data, 'rawData')          && ~isempty(app.data.rawData);
    hasFilt = isfield(app.data, 'filteredData')     && ~isempty(app.data.filteredData);
    hasDet  = isfield(app.data, 'candidateBubbles') && ~isempty(app.data.candidateBubbles);
    hasLoc  = isfield(app.data, 'localizations')    && ~isempty(app.data.localizations);
    hasTrk  = isfield(app.data, 'tracks_raw')       && ~isempty(app.data.tracks_raw);

    ST = @(b) matlab.lang.OnOffSwitchState(logical(b));

    if isfield(app.ui, 'RunFilterButton'),  app.ui.RunFilterButton.Enable  = ST(hasRaw);  end
    if isfield(app.ui, 'RunDetectButton'),  app.ui.RunDetectButton.Enable  = ST(hasFilt); end
    if isfield(app.ui, 'btnRunLoc'),        app.ui.btnRunLoc.Enable        = ST(hasDet);  end
    if isfield(app.ui, 'btnRunTrack'),      app.ui.btnRunTrack.Enable      = ST(hasLoc);  end
    if isfield(app.ui, 'btnRunPost'),       app.ui.btnRunPost.Enable       = ST(hasTrk);  end
    if isfield(app.ui, 'btnRunRender'),     app.ui.btnRunRender.Enable     = ST(hasTrk);  end
    if isfield(app.ui, 'btnRunAll'),        app.ui.btnRunAll.Enable        = ST(hasRaw);  end
    if isfield(app.ui, 'GenerateImagesButton'), app.ui.GenerateImagesButton.Enable = ST(hasTrk); end

    if isfield(app.ui, 'InteractiveCropBtn'), app.ui.InteractiveCropBtn.Enable = ST(hasRaw); end
    if isfield(app.ui, 'ApplyCropBtn'),       app.ui.ApplyCropBtn.Enable       = ST(hasRaw); end

    if hasRaw && isfield(app.ui, 'FrameSlider') && isvalid(app.ui.FrameSlider)
        n = size(app.data.rawData, 3);
        app.ui.FrameSlider.Limits = [1, max(2, n)];
        app.ui.FrameSlider.Enable = 'on';
        if isfield(app.ui, 'FrameField') && isvalid(app.ui.FrameField)
            app.ui.FrameField.Limits = [1, n];
        end
    end

    if isfield(app, 'undoManager') && isfield(app.ui, 'btnUndo')
        try
            app.ui.btnUndo.Enable = ST(app.undoManager.canUndo());
            app.ui.btnRedo.Enable = ST(app.undoManager.canRedo());
        catch
        end
    end

    app.state.currentState = stage;
    app.state.currentStage = stage;  % alias for any code reading the old name
end

function data = clearDownstreamData(data, stage)
    % Clears cached outputs downstream of a given pipeline stage.
    % Stage codes: 0 = raw loaded, 1 = filter, 2 = detect, 3 = localize,
    %              4 = track,    5 = post,   6 = render.
    % Works on the `app.data` struct directly so legacy call sites like
    %   app.data = clearDownstreamData(app.data, 0)
    % keep working.
    if isnumeric(stage)
        if stage < 1
            data.filteredData     = [];
            data.candidateBubbles = [];
            data.localizations    = [];
            data.tracks_raw       = [];
            data.tracks_final     = [];
        elseif stage < 2
            data.candidateBubbles = [];
            data.localizations    = [];
            data.tracks_raw       = [];
            data.tracks_final     = [];
        elseif stage < 3
            data.localizations    = [];
            data.tracks_raw       = [];
            data.tracks_final     = [];
        elseif stage < 4
            data.tracks_raw       = [];
            data.tracks_final     = [];
        elseif stage < 5
            data.tracks_final     = [];
        end
    else
        % String interface used by my callbacks
        switch lower(string(stage))
            case "filter"
                data.filteredData     = [];
                data.candidateBubbles = [];
                data.localizations    = [];
                data.tracks_raw       = [];
                data.tracks_final     = [];
            case {"detect","detection"}
                data.candidateBubbles = [];
                data.localizations    = [];
                data.tracks_raw       = [];
                data.tracks_final     = [];
            case {"localize","localization"}
                data.localizations    = [];
                data.tracks_raw       = [];
                data.tracks_final     = [];
            case {"track","tracking"}
                data.tracks_raw       = [];
                data.tracks_final     = [];
            case {"post","postprocess"}
                data.tracks_final     = [];
        end
    end
end
% =========================================================================
% SECTION: Parameter <-> GUI Synchronization  (actual setDefaultParams schema)
% =========================================================================

function p = updateParamsFromGUI(ref)
    % Polymorphic parameter sync:
    %   p = updateParamsFromGUI(fig)  - read from guidata, write back, return params
    %   p = updateParamsFromGUI(app)  - read from app struct, return params only
    % Both call-site styles used across the codebase are supported.
    if isgraphics(ref)
        fig = ref;
        app = guidata(fig);
        writeBack = true;
    else
        app = ref;
        fig = [];
        writeBack = false;
    end
    p = app.data.params;

    % -------- Acquisition / pixel size (Top Row) --------
    if isfield(app.ui, 'TopFPSField'),     p.acq.framerate         = app.ui.TopFPSField.Value;    end
    if isfield(app.ui, 'TopPixelXField'),  p.track.pixel_X_size    = app.ui.TopPixelXField.Value; end
    if isfield(app.ui, 'TopPixelZField'),  p.track.pixel_Z_size    = app.ui.TopPixelZField.Value; end

    % -------- Stage 1: Filter --------
    if isfield(app.ui, 'FilterMethodDropdown')
        p.filter.method = app.ui.FilterMethodDropdown.Value;
    end
    if isfield(app.ui, 'SVDCutoffStart') && isfield(app.ui, 'SVDCutoffEnd')
        p.filter.svd_cutoff = [app.ui.SVDCutoffStart.Value, app.ui.SVDCutoffEnd.Value];
    end

    % Block-wise SVD
    if isfield(p.filter, 'blockwise')
        if isfield(app.ui, 'BWThresholdMethod'), p.filter.blockwise.threshold_method = app.ui.BWThresholdMethod.Value; end
        if isfield(app.ui, 'BWBlockSize'),       p.filter.blockwise.block_size_mm   = app.ui.BWBlockSize.Value; end
        if isfield(app.ui, 'BWOverlapPct'),      p.filter.blockwise.overlap_pct     = app.ui.BWOverlapPct.Value; end
        if isfield(app.ui, 'BWTissueFreq'),      p.filter.blockwise.tissue_freq_hz  = app.ui.BWTissueFreq.Value; end
        if isfield(app.ui, 'BWMPSigma'),         p.filter.blockwise.mp_deviation_sigma = app.ui.BWMPSigma.Value; end
        if isfield(app.ui, 'BWGradientPct'),     p.filter.blockwise.gradient_pct    = app.ui.BWGradientPct.Value; end
        if isfield(app.ui, 'BWMinBlood'),        p.filter.blockwise.min_blood_comps = app.ui.BWMinBlood.Value; end
        if isfield(app.ui, 'BWMaxTissueFrac'),   p.filter.blockwise.max_tissue_frac = app.ui.BWMaxTissueFrac.Value; end
        if isfield(app.ui, 'BWPlotMaps'),        p.filter.blockwise.plot_maps       = app.ui.BWPlotMaps.Value; end
        if isfield(app.ui, 'BWManualCutoff')
            raw = app.ui.BWManualCutoff.Value;
            try
                v = str2num(raw); %#ok<ST2NM>
                if numel(v) == 2
                    p.filter.blockwise.manual_cutoff = v;
                end
            catch
            end
        end
    end

    % Butterworth
    if isfield(app.ui, 'EnableButterworth'), p.filter.enable_butterworth = app.ui.EnableButterworth.Value; end
    if isfield(app.ui, 'ButterCutoff')
        try
            v = str2num(app.ui.ButterCutoff.Value); %#ok<ST2NM>
            if numel(v) == 1 || numel(v) == 2
                p.filter.butter_cutoff = v;
            end
        catch
        end
    end
    if isfield(app.ui, 'ButterOrder'), p.filter.butter_order = app.ui.ButterOrder.Value; end

    % Spatial conditioning
    if isfield(app.ui, 'SpatialMethodDrop'),  p.filter.spatial_method = app.ui.SpatialMethodDrop.Value; end

    % -------- Stage 2: Detection --------
    if isfield(app.ui, 'DetectMethodDropdown'), p.loc.DetectMethod = app.ui.DetectMethodDropdown.Value; end
    if isfield(app.ui, 'DetectFWHMField')
        try
            v = str2num(app.ui.DetectFWHMField.Value); %#ok<ST2NM>
            if ~isempty(v)
                if isscalar(v), v = [v v]; end
                p.loc.fwhm = v;
            end
        catch
        end
    end
    if isfield(app.ui, 'LocThreshField'),   p.loc.detection_threshold    = app.ui.LocThreshField.Value; end
    if isfield(app.ui, 'LocMaxField'),      p.loc.max_bubbles_per_frame  = app.ui.LocMaxField.Value; end
    if isfield(app.ui, 'LocShiftFactor'),   p.loc.qc_max_shift_factor    = app.ui.LocShiftFactor.Value; end
    if isfield(app.ui, 'LocQCDivergence'),  p.loc.enable_divergence_check= app.ui.LocQCDivergence.Value; end
    if isfield(app.ui, 'LocQCRoiMaxima'),   p.loc.enable_roi_maxima_check= app.ui.LocQCRoiMaxima.Value; end
    if isfield(app.ui, 'NP_AlphaField'),    p.loc.NP_alpha0              = app.ui.NP_AlphaField.Value; end
    if isfield(app.ui, 'NCC_ThreshField'),  p.loc.crosscor_threshold     = app.ui.NCC_ThreshField.Value; end
    if isfield(app.ui, 'PeakContrastField'), p.loc.h_contrast            = app.ui.PeakContrastField.Value; end

    % -------- Stage 3: Localization --------
    if isfield(app.ui, 'LocMethodDropdown'), p.loc.method = app.ui.LocMethodDropdown.Value; end
    if isfield(app.ui, 'GaussMinRSquared'),  p.loc.min_r_squared        = app.ui.GaussMinRSquared.Value; end

    % -------- Stage 4: Tracking --------
    if isfield(app.ui, 'TrackMethodDropdown'),      p.track.method                 = app.ui.TrackMethodDropdown.Value; end
    if isfield(app.ui, 'MaxDistField'),             p.track.max_linking_distance   = app.ui.MaxDistField.Value; end
    if isfield(app.ui, 'GapFramesField'),           p.track.max_gap_closing_frames = app.ui.GapFramesField.Value; end
    if isfield(app.ui, 'MinLengthField'),           p.track.min_track_length       = app.ui.MinLengthField.Value; end
    if isfield(app.ui, 'UseAdvancedCostCheckbox'),  p.track.use_advanced_cost_matrix = app.ui.UseAdvancedCostCheckbox.Value; end

    if ~isfield(p.track, 'kalman'), p.track.kalman = struct(); end
    if isfield(app.ui, 'KalmanModelDrop'), p.track.kalman.motion_model      = app.ui.KalmanModelDrop.Value; end
    if isfield(app.ui, 'KalmanNoise'),     p.track.kalman.process_noise     = app.ui.KalmanNoise.Value; end
    if isfield(app.ui, 'AssignmentDrop'),  p.track.kalman.assignment_method = app.ui.AssignmentDrop.Value; end

    % Track QC
    if ~isfield(p.track, 'qc'), p.track.qc = struct(); end
    if isfield(app.ui, 'TrackQCDirection'),    p.track.qc.enable_direction_constraint    = app.ui.TrackQCDirection.Value; end
    if isfield(app.ui, 'QCMaxAngle'),          p.track.qc.max_angle_change_deg           = app.ui.QCMaxAngle.Value; end
    if isfield(app.ui, 'TrackQCAcceleration'), p.track.qc.enable_acceleration_constraint = app.ui.TrackQCAcceleration.Value; end
    if isfield(app.ui, 'QCAccelFactor'),       p.track.qc.acceleration_C_factor          = app.ui.QCAccelFactor.Value; end
    if isfield(app.ui, 'TrackQCVD'),           p.track.qc.enable_vd_constraint           = app.ui.TrackQCVD.Value; end
    if isfield(app.ui, 'QCVDRatio'),           p.track.qc.max_vd_ratio                   = app.ui.QCVDRatio.Value; end

    % -------- Stage 5: Post-process --------
    if isfield(app.ui, 'EnablePostProcessing'),  p.track.enable_postprocessing = app.ui.EnablePostProcessing.Value; end
    if isfield(app.ui, 'SmoothField'),           p.track.smoothing_factor      = app.ui.SmoothField.Value; end
    if isfield(app.ui, 'DisplayMinLengthField'), p.track.display_min_length    = app.ui.DisplayMinLengthField.Value; end

    % -------- Stage 6: Rendering --------
    if isfield(app.ui, 'UpsamplingField'),  p.render.upsampling_factor = app.ui.UpsamplingField.Value; end
    if isfield(app.ui, 'RenderMethodDrop'), p.render.method            = app.ui.RenderMethodDrop.Value; end

    % Crop box
    if isfield(app.ui, 'CropBoxField')
        try
            v = str2num(app.ui.CropBoxField.Value); %#ok<ST2NM>
            if numel(v) == 4 || isempty(v)
                p.io.cropBox = v;
            end
        catch
        end
    end

    if writeBack
        app.data.params = p;
        guidata(fig, app);
    end
end

function populateGUIFromParams(fig)
    app = guidata(fig);
    p   = app.data.params;

    sv = @(name, val) setSafe(app.ui, name, val);

    % Acquisition / pixel size
    sv('TopFPSField',    getDefault(p, 'acq.framerate', 200));
    sv('TopPixelXField', getDefault(p, 'track.pixel_X_size', 0.05));
    sv('TopPixelZField', getDefault(p, 'track.pixel_Z_size', 0.05));

    % Filter
    sv('FilterMethodDropdown', getDefault(p, 'filter.method', 'svd'));
    svd_cut = getDefault(p, 'filter.svd_cutoff', [5 100]);
    if numel(svd_cut) >= 2
        sv('SVDCutoffStart', svd_cut(1));
        sv('SVDCutoffEnd',   svd_cut(2));
    end

    % Block-wise
    if isfield(p.filter, 'blockwise')
        bw = p.filter.blockwise;
        sv('BWThresholdMethod', getDefault(bw, 'threshold_method', 'DopplerGradient'));
        sv('BWBlockSize',       getDefault(bw, 'block_size_mm', 4.0));
        sv('BWOverlapPct',      getDefault(bw, 'overlap_pct', 75));
        sv('BWTissueFreq',      getDefault(bw, 'tissue_freq_hz', -1));
        sv('BWMPSigma',         getDefault(bw, 'mp_deviation_sigma', 2.0));
        sv('BWGradientPct',     getDefault(bw, 'gradient_pct', 0.10));
        sv('BWMinBlood',        getDefault(bw, 'min_blood_comps', 3));
        sv('BWMaxTissueFrac',   getDefault(bw, 'max_tissue_frac', 0.6));
        sv('BWPlotMaps',        getDefault(bw, 'plot_maps', false));
        mc = getDefault(bw, 'manual_cutoff', [10 200]);
        if numel(mc) == 2
            sv('BWManualCutoff', sprintf('[%g %g]', mc(1), mc(2)));
        end
        updateBlockwiseOptions(fig, app);
    end

    % Butterworth
    sv('EnableButterworth', getDefault(p, 'filter.enable_butterworth', false));
    bc = getDefault(p, 'filter.butter_cutoff', [50 250]);
    if numel(bc) == 2
        sv('ButterCutoff', sprintf('[%g %g]', bc(1), bc(2)));
    else
        sv('ButterCutoff', sprintf('%g', bc));
    end
    sv('ButterOrder', getDefault(p, 'filter.butter_order', 2));

    % Spatial
    sv('SpatialMethodDrop', getDefault(p, 'filter.spatial_method', 'Gaussian'));

    % Detection
    sv('DetectMethodDropdown', getDefault(p, 'loc.DetectMethod', 'Intensity'));
    fwhm_val = getDefault(p, 'loc.fwhm', [3 3]);
    sv('DetectFWHMField', mat2str(fwhm_val));
    sv('LocThreshField',  getDefault(p, 'loc.detection_threshold', 0.5));
    sv('LocMaxField',     getDefault(p, 'loc.max_bubbles_per_frame', 200));
    sv('LocShiftFactor',  getDefault(p, 'loc.qc_max_shift_factor', 1));
    sv('LocQCDivergence', getDefault(p, 'loc.enable_divergence_check', true));
    sv('LocQCRoiMaxima',  getDefault(p, 'loc.enable_roi_maxima_check', true));
    sv('NP_AlphaField',   getDefault(p, 'loc.NP_alpha0', 0.01));
    sv('NCC_ThreshField', getDefault(p, 'loc.crosscor_threshold', 0.7));

    % Localization
    sv('LocMethodDropdown', getDefault(p, 'loc.method', 'radial'));
    sv('GaussMinRSquared',  getDefault(p, 'loc.min_r_squared', 0.3));

    % Tracking
    sv('TrackMethodDropdown',       getDefault(p, 'track.method', 'kalman'));
    sv('MaxDistField',              getDefault(p, 'track.max_linking_distance', 5));
    sv('GapFramesField',            getDefault(p, 'track.max_gap_closing_frames', 2));
    sv('MinLengthField',            getDefault(p, 'track.min_track_length', 8));
    sv('UseAdvancedCostCheckbox',   getDefault(p, 'track.use_advanced_cost_matrix', false));
    sv('KalmanModelDrop',           getDefault(p, 'track.kalman.motion_model', 'ConstantVelocity'));
    sv('KalmanNoise',               getDefault(p, 'track.kalman.process_noise', 0.1));
    sv('AssignmentDrop',            getDefault(p, 'track.kalman.assignment_method', 'hungarian'));

    % Track QC
    sv('TrackQCDirection',    getDefault(p, 'track.qc.enable_direction_constraint', true));
    sv('QCMaxAngle',          getDefault(p, 'track.qc.max_angle_change_deg', 60));
    sv('TrackQCAcceleration', getDefault(p, 'track.qc.enable_acceleration_constraint', true));
    sv('QCAccelFactor',       getDefault(p, 'track.qc.acceleration_C_factor', 3));
    sv('TrackQCVD',           getDefault(p, 'track.qc.enable_vd_constraint', false));
    sv('QCVDRatio',           getDefault(p, 'track.qc.max_vd_ratio', 0.5));

    % Post
    sv('EnablePostProcessing',  getDefault(p, 'track.enable_postprocessing', true));
    sv('SmoothField',           getDefault(p, 'track.smoothing_factor', 5));
    sv('DisplayMinLengthField', getDefault(p, 'track.display_min_length', 10));

    % Render
    sv('UpsamplingField',  getDefault(p, 'render.upsampling_factor', 10));
    sv('RenderMethodDrop', getDefault(p, 'render.method', 'histogram'));

    % Crop
    cb = getDefault(p, 'io.cropBox', []);
    if isempty(cb)
        sv('CropBoxField', '[]');
    else
        sv('CropBoxField', mat2str(cb));
    end

    % Mirror paired sliders where they exist
    mirrorPair(app.ui, 'LocThreshField',       'LocThreshSlider');
    mirrorPair(app.ui, 'LocMaxField',          'LocMaxSlider');
    mirrorPair(app.ui, 'MaxDistField',         'MaxDistSlider');
    mirrorPair(app.ui, 'GapFramesField',       'GapFramesSlider');
    mirrorPair(app.ui, 'MinLengthField',       'MinLengthSlider');
    mirrorPair(app.ui, 'DisplayMinLengthField','DisplayMinLengthSlider');

    updateFilterOptions(fig);
    updateDetectionOptions(fig);
    updateLocalizationOptions(fig);
    updateTrackingOptions(fig);
end

function setSafe(ui, name, val)
    if ~isfield(ui, name) || isempty(ui.(name)) || ~isvalid(ui.(name)), return; end
    h = ui.(name);
    try
        cls = class(h);
        if contains(cls, 'CheckBox')
            h.Value = logical(val);
        elseif contains(cls, 'DropDown')
            if isnumeric(val), val = num2str(val); end
            v = char(val);
            if any(strcmp(h.Items, v))
                h.Value = v;
            elseif ~isempty(h.Items)
                % try case-insensitive match
                idx = find(strcmpi(h.Items, v), 1);
                if ~isempty(idx), h.Value = h.Items{idx}; end
            end
        elseif contains(cls, 'Slider')
            if isnumeric(val)
                v = max(h.Limits(1), min(h.Limits(2), double(val)));
                h.Value = v;
            end
        elseif contains(cls, 'EditField')
            if contains(cls, 'Numeric')
                if isnumeric(val) && isscalar(val)
                    h.Value = double(val);
                elseif ischar(val) || isstring(val)
                    try, h.Value = str2double(val); catch, end
                end
            else
                if isnumeric(val), val = num2str(val); end
                h.Value = char(val);
            end
        else
            h.Value = val;
        end
    catch
    end
end

function mirrorPair(ui, fieldName, sliderName)
    if isfield(ui, fieldName) && isfield(ui, sliderName) ...
            && isvalid(ui.(fieldName)) && isvalid(ui.(sliderName))
        try
            v = ui.(fieldName).Value;
            if isnumeric(v)
                lim = ui.(sliderName).Limits;
                ui.(sliderName).Value = max(lim(1), min(lim(2), v));
            end
        catch
        end
    end
end

function resetAllParams(fig)
    app = guidata(fig);
    saveParamState(fig);
    try
        if exist('setDefaultParams', 'file') == 2
            app.data.params = setDefaultParams(true);
        else
            app.data.params = createFallbackParams();
        end
    catch
        app.data.params = createFallbackParams();
    end
    app.data.params = ensureAllParamFields(app.data.params);
    guidata(fig, app);
    populateGUIFromParams(fig);
    app = guidata(fig);
    app = manageGUIState(app, -1);
    guidata(fig, app);
    setStatus(app, 'All parameters reset to defaults.', 'green');
end

% =========================================================================
% SECTION: Option visibility (driven by the algorithm registry)
% =========================================================================

function updateFilterOptions(fig)
    app = guidata(fig);
    if ~isfield(app.ui, 'FilterMethodDropdown'), return; end
    method = app.ui.FilterMethodDropdown.Value;

    off = matlab.lang.OnOffSwitchState.off;
    on  = matlab.lang.OnOffSwitchState.on;

    % 1. Hide all conditional panels by default
    panels = {'p_svd', 'p_dcc', 'p_blockwise'};
    for i = 1:numel(panels)
        if isfield(app.ui, panels{i}) && isvalid(app.ui.(panels{i}))
            app.ui.(panels{i}).Visible = off;
        end
    end

    % 2. Define the baseline absolute row heights in pixels
    % Index mapping: 5=SVD, 6=DCC, 9=Blockwise (start hidden at 0)
    rowHeights = {25, 30, 140, 140, 0, 0, 130, 110, 0, 20};

    % 3. Reveal the active panel and assign its strict pixel height
    switch lower(char(method))
        case {'svd_filter', 'svd_ssm'}
            if isfield(app.ui, 'p_svd'), app.ui.p_svd.Visible = on; end
            rowHeights{5} = 110; % Expanded SVD Panel (was 90)
            
        case 'dcc_svd'
            if isfield(app.ui, 'p_dcc'), app.ui.p_dcc.Visible = on; end
            rowHeights{6} = 130; % Expanded DCC Panel (was 110)
            
        case 'svd_blockwise'
            if isfield(app.ui, 'p_blockwise')
                app.ui.p_blockwise.Visible = on;
            end
            updateBlockwiseOptions(fig, app);
            rowHeights{9} = 220; % Generously expanded Blockwise Panel (was 190)
    end

    % 4. Apply the exact pixel heights to the inner grid
    % Because NO row is '1x' or 'fit', MATLAB is forced to overflow and show the scrollbar!
    if isfield(app.ui, 'g_filter') && isvalid(app.ui.g_filter)
        app.ui.g_filter.RowHeight = rowHeights;
    end

    updateSpatialOptions(fig, app);
end

function updateDetectionOptions(fig)
    app = guidata(fig);
    if ~isfield(app.ui, 'DetectMethodDropdown'), return; end
    method = app.ui.DetectMethodDropdown.Value;
    reg = getAlgorithmRegistry();
    entry = registryEntryById(reg.detect, method);

    % Show/hide NP and NCC parameter sub-panels
    needsTemplate = ~isempty(entry) && isfield(entry, 'needsTemplate') && entry.needsTemplate;
    isNP  = any(strcmpi(method, {'np','neyman-pearson','neyman_pearson','neymanpearson'}));
    isNCC = any(strcmpi(method, {'ncc','cross-correlation','crosscorrelation','crosscor'}));

    % NP alpha field + label
    for name = {'NP_AlphaField','lbl_NP_alpha'}
        if isfield(app.ui, name{1}) && isvalid(app.ui.(name{1}))
            app.ui.(name{1}).Visible = matlab.lang.OnOffSwitchState(isNP);
        end
    end

    % NCC threshold field + label
    for name = {'NCC_ThreshField','lbl_NCC_thresh'}
        if isfield(app.ui, name{1}) && isvalid(app.ui.(name{1}))
            app.ui.(name{1}).Visible = matlab.lang.OnOffSwitchState(isNCC);
        end
    end

    % Hide the entire parent panel if neither NP nor NCC is selected
    if isfield(app.ui, 'p_detect_method_params') && isvalid(app.ui.p_detect_method_params)
        app.ui.p_detect_method_params.Visible = matlab.lang.OnOffSwitchState(isNP || isNCC);
    end

    % Hint label if any template-requiring method is active
    if isfield(app.ui, 'LblDetectHint') && isvalid(app.ui.LblDetectHint)
        if needsTemplate
            app.ui.LblDetectHint.Text = 'This method uses a PSF template. Configure via "Advanced..." button.';
            app.ui.LblDetectHint.FontColor = [0.7 0.35 0];
        else
            app.ui.LblDetectHint.Text = '';
        end
    end
end

function updateLocalizationOptions(fig)
    app = guidata(fig);
    if ~isfield(app.ui, 'LocMethodDropdown'), return; end
    method = app.ui.LocMethodDropdown.Value;
    reg = getAlgorithmRegistry();
    entry = registryEntryById(reg.loc, method);
    isGauss = ~isempty(entry) && isfield(entry,'isGaussian') && entry.isGaussian;

    st = matlab.lang.OnOffSwitchState(isGauss);
    if isfield(app.ui, 'GaussMinRSquared'), app.ui.GaussMinRSquared.Enable = st; end

    % Gaussian Fit QC panel: show only for Gaussian localizers
    if isfield(app.ui, 'p_loc_qc_gauss') && isvalid(app.ui.p_loc_qc_gauss)
        app.ui.p_loc_qc_gauss.Visible = st;
    end

    % Common Localization QC panel: ALWAYS visible (used by all methods)
    if isfield(app.ui, 'p_loc_qc_radial') && isvalid(app.ui.p_loc_qc_radial)
        app.ui.p_loc_qc_radial.Visible = 'on';
    end
end

function updateTrackingOptions(fig)
    app = guidata(fig);
    if ~isfield(app.ui, 'TrackMethodDropdown'), return; end
    method = app.ui.TrackMethodDropdown.Value;
    reg = getAlgorithmRegistry();
    entry = registryEntryById(reg.track, method);
    if isempty(entry)
        entry = struct('isKalman', false, 'usesHK', false, 'showsGain', false);
    end

    stK    = matlab.lang.OnOffSwitchState(getDefault(entry, 'isKalman',  false));
    stHK   = matlab.lang.OnOffSwitchState(getDefault(entry, 'usesHK',    false));
    stGain = matlab.lang.OnOffSwitchState(getDefault(entry, 'showsGain', false));

    % Allow advanced cost matrix for Kalman AND Hungarian variants
    isHungarian = strcmpi(method, 'Hungarian');
    stAdvancedCostAllowed = matlab.lang.OnOffSwitchState(getDefault(entry, 'isKalman', false) || isHungarian);
    
    stCost = matlab.lang.OnOffSwitchState( ...
        stAdvancedCostAllowed && ...
        isfield(app.ui,'UseAdvancedCostCheckbox') && ...
        app.ui.UseAdvancedCostCheckbox.Value);

    if isfield(app.ui, 'KalmanModelDrop'),     app.ui.KalmanModelDrop.Enable     = stK;  end
    if isfield(app.ui, 'KalmanNoise'),         app.ui.KalmanNoise.Enable         = stK;  end
    if isfield(app.ui, 'AssignmentDrop'),      app.ui.AssignmentDrop.Enable      = stK;  end
    
    % Use the new combined flag for the checkbox
    if isfield(app.ui, 'UseAdvancedCostCheckbox'), app.ui.UseAdvancedCostCheckbox.Enable = stAdvancedCostAllowed; end
    if isfield(app.ui, 'BtnConfigCostMatrix'), app.ui.BtnConfigCostMatrix.Enable = stCost; end
    if isfield(app.ui, 'BtnConfigHK'),         app.ui.BtnConfigHK.Enable         = stHK; end
    if isfield(app.ui, 'BtnAdvKalman'),        app.ui.BtnAdvKalman.Enable        = stK;  end

    if isfield(app.ui, 'p_kgain') && isvalid(app.ui.p_kgain)
        app.ui.p_kgain.Visible = stGain;
    end

    updateKalmanGainSummary(fig);
end
% =========================================================================
% SECTION: Frame Display
% =========================================================================

function displayCurrentFrame(fig)
    app = guidata(fig);
    if ~isfield(app.data, 'rawData') || isempty(app.data.rawData)
        if isfield(app.ui, 'ax') && isvalid(app.ui.ax)
            cla(app.ui.ax);
            title(app.ui.ax, 'Please Load Raw Data to begin...');
        end
        return;
    end
    if ~isfield(app, 'displayManager') || isempty(app.displayManager), return; end

    % Sync the current frame index from the slider/field so DisplayManager
    % picks up the right frame.
    if isfield(app.ui, 'FrameSlider') && isvalid(app.ui.FrameSlider)
        app.state.currentFrame = round(app.ui.FrameSlider.Value);
    end
    guidata(fig, app);

    try
        app.displayManager.displayFrame(app);
    catch ME
        warning(ME.identifier, 'Display error: %s', ME.message);
    end
end

function onFrameSliderChanged(fig, ~)
    app = guidata(fig);
    if ~isfield(app.ui, 'FrameSlider'), return; end
    f = round(app.ui.FrameSlider.Value);
    app.ui.FrameSlider.Value = f;
    if isfield(app.ui, 'FrameField') && isvalid(app.ui.FrameField)
        app.ui.FrameField.Value = f;
    end
    app.state.currentFrame = f;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function onFrameFieldChanged(fig, ~)
    app = guidata(fig);
    if ~isfield(app.ui, 'FrameField'), return; end
    f = app.ui.FrameField.Value;
    if isfield(app.ui, 'FrameSlider')
        f = max(app.ui.FrameSlider.Limits(1), min(app.ui.FrameSlider.Limits(2), f));
        app.ui.FrameSlider.Value = f;
    end
    app.ui.FrameField.Value  = f;
    app.state.currentFrame = f;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function onDisplayFilterChanged(fig, ~)
    displayCurrentFrame(fig);
end

function togglePlayback(fig)
    app = guidata(fig);
    if isfield(app, 'displayTimer') && ~isempty(app.displayTimer) ...
            && isvalid(app.displayTimer) && strcmp(app.displayTimer.Running, 'on')
        stop(app.displayTimer);
        if isfield(app.ui, 'PlayPauseButton') && isvalid(app.ui.PlayPauseButton)
            app.ui.PlayPauseButton.Text = 'Play';
        end
    else
        if ~isfield(app, 'displayTimer') || isempty(app.displayTimer) || ~isvalid(app.displayTimer)
            app.displayTimer = timer( ...
                'ExecutionMode', 'fixedRate', ...
                'Period', ULM_Constants.PLAYBACK_TIMER_PERIOD, ...
                'TimerFcn', @(~,~) timerCallback(fig));
        end
        guidata(fig, app);
        start(app.displayTimer);
        if isfield(app.ui, 'PlayPauseButton') && isvalid(app.ui.PlayPauseButton)
            app.ui.PlayPauseButton.Text = 'Pause';
        end
    end
    guidata(fig, app);
end

function timerCallback(fig)
    if ~isvalid(fig), return; end
    app = guidata(fig);
    if ~isfield(app.ui, 'FrameSlider') || ~isvalid(app.ui.FrameSlider), return; end
    f = round(app.ui.FrameSlider.Value);
    lim = app.ui.FrameSlider.Limits;
    f_next = f + 1;
    if f_next > lim(2), f_next = lim(1); end
    app.ui.FrameSlider.Value = f_next;
    if isfield(app.ui, 'FrameField') && isvalid(app.ui.FrameField)
        app.ui.FrameField.Value = f_next;
    end
    app.state.currentFrame = f_next;
    guidata(fig, app);
    displayCurrentFrame(fig);
end

function syncSliderField(slider, field)
% Bidirectional link between a uislider and a numeric uieditfield.
    if ~isvalid(slider) || ~isvalid(field), return; end
    field.Value = slider.Value;
    slider.ValueChangedFcn = @(~,~) updateFieldFromSlider(slider, field);
    field.ValueChangedFcn  = @(~,~) updateSliderFromField(slider, field);
end

function updateFieldFromSlider(slider, field)
    if isvalid(slider) && isvalid(field)
        field.Value = slider.Value;
    end
end

function updateSliderFromField(slider, field)
    if isvalid(slider) && isvalid(field)
        v = max(slider.Limits(1), min(slider.Limits(2), field.Value));
        slider.Value = v;
        field.Value  = v;
    end
end

function onDCCSliderChanged(fig)
    % DCC band sliders are paired with editfields; re-build the band
    % visualization to give immediate feedback on cluster boundaries.
    app = guidata(fig);
    saveParamState(fig);

    % Copy DCC slider values into params
    p = app.data.params;
    if ~isfield(p, 'filter'), p.filter = struct(); end
    if ~isfield(p.filter, 'dcc'), p.filter.dcc = struct(); end
    if isfield(app.ui, 'DCCTissueStart'), p.filter.dcc.tissue_start = app.ui.DCCTissueStart.Value; end
    if isfield(app.ui, 'DCCTissueEnd'),   p.filter.dcc.tissue_end   = app.ui.DCCTissueEnd.Value; end
    if isfield(app.ui, 'DCCBloodStart'),  p.filter.dcc.blood_start  = app.ui.DCCBloodStart.Value; end
    if isfield(app.ui, 'DCCBloodEnd'),    p.filter.dcc.blood_end    = app.ui.DCCBloodEnd.Value; end
    if isfield(app.ui, 'DCCNoiseStart'),  p.filter.dcc.noise_start  = app.ui.DCCNoiseStart.Value; end
    if isfield(app.ui, 'DCCNoiseEnd'),    p.filter.dcc.noise_end    = app.ui.DCCNoiseEnd.Value; end
    app.data.params = p;
    guidata(fig, app);

    app = guidata(fig); app.data = clearDownstreamData(app.data, 1); app = manageGUIState(app, 0); guidata(fig, app);
end

function onTabChanged(fig)
    app = guidata(fig);
    if ~isfield(app.ui, 'tabGroup'), return; end
    tab = app.ui.tabGroup.SelectedTab;
    if isempty(tab), return; end
    if contains(tab.Title, 'ROI', 'IgnoreCase', true) || contains(tab.Title, 'Mask', 'IgnoreCase', true) || contains(tab.Title, 'Detect', 'IgnoreCase', true)
        prepareROITab(fig);
    end
end

% =========================================================================
% SECTION: ROI / Mask Tab
% =========================================================================

function prepareROITab(fig)
    app = guidata(fig);
    if isempty(app.data.filteredData), return; end

    if isempty(app.data.baseVesselMap)
        [H, W, ~] = size(app.data.filteredData);

        % --- Attempt 1: Load pre-computed vessel map from Results folder ---
        vMap = tryLoadExternalVesselMap(app, [H, W]);

        if ~isempty(vMap)
            app.data.vesselMapSource = 'external';
            fprintf('  Using external vessel map from Results folder.\n');
        else
            % --- Attempt 2: Fallback — compute TMIP from current data ---
            showProgress(fig, 'Calculating vessel map (TMIP)...');

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

            app.data.vesselMapSource = 'tmip';
            fprintf('  No external vessel map found — using TMIP fallback.\n');
            hideProgress(fig);
        end

        app.data.baseVesselMap = vMap;
        app.data.vesselMap     = vMap;

        % Pre-compute the red overlay template for ROI mask display
        app.data.redOverlayTemplate = cat(3, ...
            ones(size(vMap)), zeros(size(vMap)), zeros(size(vMap)));

        app.ui.ROIContrastSlider.Value = 1;
        app.ui.ROIContrastField.Value  = 1;
        app.ui.ROIThreshSlider.Value   = 0;
        app.ui.ROIThreshField.Value    = 0;

        % Update status label
        if isfield(app.ui, 'maskStatusLabel') && isvalid(app.ui.maskStatusLabel)
            if strcmp(app.data.vesselMapSource, 'external')
                app.ui.maskStatusLabel.Text = 'Vessel map: External (full dataset)';
            else
                app.ui.maskStatusLabel.Text = 'Vessel map: TMIP (single file)';
            end
        end
    end

    app.ui.EnhanceMethodDrop.Value   = 'None';
    app.ui.EnhanceAmountSlider.Value = 0.5;
    app.ui.axHist.YScale             = 'log';
    app.state.isROIPreview           = true;
    app.ui.chkPreviewROI.Value       = 1;

    guidata(fig, app);
    updateHistogram(app);
    displayCurrentFrame(fig);
end

function vMap = tryLoadExternalVesselMap(app, targetSize)
% TRYLOADEXTERNALVESSELMAP  Search for a pre-computed vessel map in the
%   experiment's Results folder and return it resized & normalized to match
%   the currently loaded filtered data.
%
%   SEARCH ORDER (walks up to 3 parent levels):
%       1.  <data_folder>/Results/
%       2.  <data_folder>/../Results/
%       3.  <data_folder>/../../Results/
%
%   FILENAME PATTERNS (tried in order):
%       mean_bmode_vessel_*.mat
%       mean_bmode_*.mat
%
%   INPUTS:
%       app        – guidata struct (needs app.data.params.io.data_folder)
%       targetSize – [H, W] that the output must match (filtered data dims)
%
%   OUTPUT:
%       vMap – Normalized [0,1] vessel map sized [H, W], or [] if nothing found.

    vMap = [];
    TAG  = '[VesselMap]';

    % --- Resolve data folder ---
    if ~isfield(app.data.params, 'io') || ...
       ~isfield(app.data.params.io, 'data_folder')
        fprintf('  %s Skipped: app.data.params.io.data_folder field missing.\n', TAG);
        return;
    end

    dataFolder = app.data.params.io.data_folder;
    if isempty(dataFolder)
        fprintf('  %s Skipped: data_folder is empty.\n', TAG);
        return;
    end

    % Clean trailing separators for reliable fileparts behavior
    dataFolder = regexprep(dataFolder, '[/\\]+$', '');

    if ~isfolder(dataFolder)
        fprintf('  %s Skipped: data_folder is not a valid folder: "%s"\n', TAG, dataFolder);
        return;
    end

    fprintf('  %s Data folder: "%s"\n', TAG, dataFolder);

    % --- Build search directories: walk up to 3 parent levels ---
    searchDirs = {};
    current = dataFolder;
    for level = 0:2
        candidate = fullfile(current, 'Results');
        searchDirs{end+1} = candidate; %#ok<AGROW>
        parent = fileparts(current);
        if strcmp(parent, current), break; end  % reached filesystem root
        current = parent;
    end

    % --- Filename patterns to try (most specific first) ---
    patterns = {'mean_bmode_vessel_*.mat', 'mean_bmode_*.mat'};

    % --- Search ---
    matFile  = '';
    foundDir = '';
    for d = 1:numel(searchDirs)
        if ~isfolder(searchDirs{d})
            fprintf('  %s  Search [%d]: NOT FOUND  "%s"\n', TAG, d, searchDirs{d});
            continue;
        end
        fprintf('  %s  Search [%d]: scanning   "%s"\n', TAG, d, searchDirs{d});

        for p = 1:numel(patterns)
            candidates = dir(fullfile(searchDirs{d}, patterns{p}));
            if ~isempty(candidates)
                % Pick the largest file (likely highest magnification)
                [~, idx] = max([candidates.bytes]);
                matFile  = fullfile(searchDirs{d}, candidates(idx).name);
                foundDir = searchDirs{d};
                fprintf('  %s  MATCH: "%s" (%.1f KB, pattern "%s")\n', ...
                    TAG, candidates(idx).name, candidates(idx).bytes/1024, patterns{p});
                break;
            end
        end
        if ~isempty(matFile), break; end

        % Show what IS in the folder so the user can spot naming issues
        allMats = dir(fullfile(searchDirs{d}, '*.mat'));
        if ~isempty(allMats)
            names = {allMats.name};
            fprintf('  %s  No pattern match. Files present: %s\n', TAG, strjoin(names, ', '));
        else
            fprintf('  %s  Folder exists but contains no .mat files.\n', TAG);
        end
    end

    if isempty(matFile)
        fprintf('  %s No external vessel map found in any search path.\n', TAG);
        return;
    end

    % --- Load and extract the 2D vessel image ---
    try
        fprintf('  %s Loading: "%s"\n', TAG, matFile);
        contents = load(matFile);
        fields   = fieldnames(contents);
        raw      = [];
        usedVar  = '';

        for k = 1:numel(fields)
            v = contents.(fields{k});
            if isnumeric(v) && ismatrix(v) && min(size(v)) > 1
                raw     = double(v);
                usedVar = fields{k};
                break;
            end
        end

        if isempty(raw)
            fprintf('  %s File loaded but no 2D matrix found. Variables: %s\n', ...
                TAG, strjoin(fields, ', '));
            return;
        end

        fprintf('  %s Extracted variable "%s" [%d x %d]\n', ...
            TAG, usedVar, size(raw,1), size(raw,2));

    catch ME
        fprintf('  %s Load failed: %s\n', TAG, ME.message);
        return;
    end

    % --- Resize to match filtered data dimensions ---
    H = targetSize(1);
    W = targetSize(2);
    if size(raw, 1) ~= H || size(raw, 2) ~= W
        fprintf('  %s Resizing from [%d x %d] to [%d x %d]...\n', ...
            TAG, size(raw,1), size(raw,2), H, W);
        raw = imresize(raw, [H W], 'bilinear');
    end

    % --- Normalize: abs -> sqrt compression -> scale to [0,1] ---
    raw = abs(raw);
    raw = raw .^ 0.5;
    mx  = max(raw(:));
    if mx > 0
        vMap = raw / mx;
    else
        vMap = raw;
    end
    fprintf('  %s External vessel map ready. Range: [%.4f, %.4f]\n', ...
        TAG, min(vMap(:)), max(vMap(:)));
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

    method  = app.ui.EnhanceMethodDrop.Value;
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
            if mx > 0, procImg = procImg / mx; end
        case 'Sharpen'
            amount = val * 2;
            procImg = imsharpen(baseImg, 'Radius', 1, 'Amount', amount);
        otherwise
            procImg = baseImg;
    end

    % Apply gamma
    currentGamma = app.ui.ROIContrastSlider.Value;
    procImg = procImg .^ currentGamma;
    mn = min(procImg(:));  mx = max(procImg(:));
    if mx > mn
        procImg = (procImg - mn) / (mx - mn);
    else
        procImg = zeros(size(procImg), 'like', procImg);
    end

    app.data.vesselMap = procImg;
    currThresh = app.ui.ROIThreshSlider.Value;
    app.data.mask = procImg >= currThresh;

    app.state.isROIPreview     = true;
    app.ui.chkPreviewROI.Value = 1;
    guidata(fig, app);

    updateHistogram(app);
    displayCurrentFrame(fig);
end

function onContrastChange(fig, gammaVal)
    app = guidata(fig);
    gammaVal = max(0.1, min(3, gammaVal));
    app.ui.ROIContrastSlider.Value = gammaVal;
    app.ui.ROIContrastField.Value  = gammaVal;
    guidata(fig, app);
    applyVesselEnhancement(fig);
end

function updateHistogram(app)
    if isempty(app.data.vesselMap), return; end
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
    app.ui.ROIThreshField.Value  = val;

    % Move existing threshold line instead of full histogram redraw
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
        uialert(app.fig, 'No mask to save. Adjust the threshold first.', 'Save Mask');
        return;
    end
    [file, path] = uiputfile('*.mat', 'Save Mask');
    if isequal(file, 0), return; end
    mask = app.data.mask; %#ok<NASGU>
    save(fullfile(path, file), 'mask');
    setStatus(app, sprintf('Mask saved: %s', file), 'green');
end

function loadMask(fig)
    app = guidata(fig);
    [file, path] = uigetfile('*.mat', 'Load Mask');
    if isequal(file, 0), return; end
    try
        S = load(fullfile(path, file));
        fn = fieldnames(S);
        mask = S.(fn{1});
        if ~isfield(app.data.params, 'proc'), app.data.params.proc = struct(); end
        app.data.params.proc.maskPath = fullfile(path, file);
        app.data.params.proc.ROIMask = logical(mask);
        % Feed the display layer so the red contour renders
        app.data.mask = logical(mask);
        % Build the red overlay template (same logic as prepareROITab)
        app.data.redOverlayTemplate = cat(3, ones(size(mask)), zeros(size(mask)), zeros(size(mask)));

        guidata(fig, app);
        if isfield(app.ui, 'maskStatusLabel') && isvalid(app.ui.maskStatusLabel)
            app.ui.maskStatusLabel.Text = sprintf('Mask loaded: %s', file);
        end
        setStatus(app, sprintf('Mask loaded: %s', file), 'green');
        % Refresh display so the mask contour appears immediately
        displayCurrentFrame(fig);
    catch ME
        uialert(app.fig, sprintf('Load mask failed: %s', ME.message), 'Error');
    end
end

function runCreateMask(fig)
    app = guidata(fig);
    if isempty(app.data.mask)
        onROIChange(fig, app.ui.ROIThreshSlider.Value);
        app = guidata(fig);
    end
    if isempty(app.data.mask), return; end
    if ~isfield(app.data.params, 'proc'), app.data.params.proc = struct(); end
    app.data.params.proc.ROIMask = app.data.mask;
    app.data.params.proc.enableInteractiveMask = true;
    guidata(fig, app);
    setStatus(app, 'Mask applied to detection.', 'green');
end

function resetMask(fig)
    app = guidata(fig);
    if isfield(app.data.params, 'proc')
        app.data.params.proc.ROIMask = [];
        app.data.params.proc.enableInteractiveMask = false;
    end
    app.data.ROIMaskPreview = [];
    guidata(fig, app);
    onROIPreviewToggle(fig);
    if isfield(app.ui, 'maskStatusLabel') && isvalid(app.ui.maskStatusLabel)
        app.ui.maskStatusLabel.Text = 'No mask loaded';
    end
    setStatus(app, 'Mask cleared.', 'blue');
end

function resetROIPanel(fig)
    app = guidata(fig);

    % --- Clear all vessel map / mask data ---
    app.data.baseVesselMap      = [];
    app.data.vesselMapSource    = 'none';
    app.data.vesselMap          = [];
    app.data.mask               = [];
    app.data.redOverlayTemplate = [];
    if isfield(app.data, 'ROIMaskPreview')
        app.data.ROIMaskPreview = [];
    end
    if isfield(app.data, 'ROIEnhanced')
        app.data.ROIEnhanced = [];
    end
    if isfield(app.data, 'ROIMeanImage')
        app.data.ROIMeanImage = [];
    end

    % --- Clear mask from detection params ---
    if isfield(app.data.params, 'proc')
        app.data.params.proc.ROIMask = [];
        app.data.params.proc.enableInteractiveMask = false;
    end

    % --- Update mask status label in Step 1 panel ---
    if isfield(app.ui, 'maskStatusLabel') && isvalid(app.ui.maskStatusLabel)
        app.ui.maskStatusLabel.Text = 'Status: None';
    end

    guidata(fig, app);

    % --- Re-enter the tab as if arriving for the first time ---
    % baseVesselMap is now empty, so prepareROITab will recalculate
    % the filtered mean, reset all sliders, show the image & histogram
    prepareROITab(fig);
end

% =========================================================================
% SECTION: DCC Visualization helper
% =========================================================================

function filteredData = reconstructDCCImage(app)
% Reconstruct blood-only signal from cached SVD using current DCC band indices.
    if ~isfield(app.data, 'U') || isempty(app.data.U) || ...
       ~isfield(app.data, 'blood_indices') || isempty(app.data.blood_indices)
        filteredData = [];
        return;
    end
    cutoff = [min(app.data.blood_indices), max(app.data.blood_indices)];
    filteredData = reconstruct_SVD_Signal(app.data.U, ...
        app.data.S_diag, app.data.V, app.data.svdDims, cutoff);
end

function str = formatIndicesToRanges(indices)
    if isempty(indices), str = ''; return; end
    indices = sort(unique(indices(:)))';
    breaks = [0, find(diff(indices) > 1), numel(indices)];
    parts = cell(1, numel(breaks) - 1);
    for k = 1:numel(breaks) - 1
        seg = indices(breaks(k)+1:breaks(k+1));
        if numel(seg) == 1
            parts{k} = sprintf('%d', seg);
        else
            parts{k} = sprintf('%d-%d', seg(1), seg(end));
        end
    end
    str = strjoin(parts, ', ');
end
% =========================================================================
% SECTION: Advanced Parameter Modals
% =========================================================================
% All three modals share the same visual grammar:
%   - grouped uipanels, each titled with the pipeline stage it belongs to
%   - bold parameter label, edit control, italic one-line descriptor
%   - Save / Cancel buttons at the bottom
% This makes it obvious which stage each parameter affects and what it does.
%
% Each parameter carries a clear section title, and every field is followed
% by a plain-English descriptor — together with the pipeline-stage label in
% the panel title these form a built-in mini reference for the user.

function openAdvancedKalmanGUI(fig)
    app = guidata(fig);
    p = app.data.params;

    d = uifigure('Name', 'Advanced Kalman / Tracking Settings', ...
        'Position', [200 200 680 680], 'WindowStyle', 'modal');
    g = uigridlayout(d, [4 1]);
    g.RowHeight = {'1x', '1x', '1x', 50};
    g.Padding = [12 12 12 12];
    g.RowSpacing = 10;

    % ---- Section 1: Cost Matrix Weights (Stage 4: Track) ----
    sec1 = uipanel(g, 'Title', ...
        'Cost Matrix Weights  |  Stage 4: Track  (active when "Use advanced cost" is enabled)', ...
        'FontWeight', 'bold');
    g1 = uigridlayout(sec1, [3 3]);
    g1.ColumnWidth = {220, 110, '1x'};
    g1.RowHeight = {'fit','fit','fit'};
    g1.Padding = [8 8 8 8];

    eDirW = makeNumRow(g1, 'Direction penalty weight', ...
        getDefault(p,'track.kalman.direction_penalty_weight',2), ...
        'Linear weight on direction deviation (0 disables).');
    eDirS = makeNumRow(g1, 'Angle penalty slope', ...
        getDefault(p,'track.kalman.angle_penalty_slope',0.3), ...
        'Slope of the angle penalty curve.');
    eBri  = makeNumRow(g1, 'Brightness penalty weight', ...
        getDefault(p,'track.kalman.brightness_penalty_weight',2), ...
        'Weight on amplitude mismatch in cost matrix.');

    % ---- Section 2: Angle Gating (Stage 4: Track) ----
    sec2 = uipanel(g, 'Title', ...
        'Angle Gating  |  Stage 4: Track  (rejects implausible direction changes before assignment)', ...
        'FontWeight', 'bold');
    g2 = uigridlayout(sec2, [3 3]);
    g2.ColumnWidth = {220, 110, '1x'};
    g2.RowHeight = {'fit','fit','fit'};
    g2.Padding = [8 8 8 8];

    eMaxA = makeNumRow(g2, 'Max angle change (deg)', ...
        getDefault(p,'track.kalman.max_angle_change_deg',70), ...
        'Hard ceiling on frame-to-frame direction change.');
    eGate = makeNumRow(g2, 'Gating angle change (deg)', ...
        getDefault(p,'track.kalman.gating_max_angle_change_deg',90), ...
        'Soft gate for pre-filtering candidates before cost computation.');
    eHist = makeNumRow(g2, 'Direction history points', ...
        getDefault(p,'track.kalman.direction_history_points',4), ...
        'Number of past positions used to estimate direction.');

    % ---- Section 3: HK Noise Scaling (Stage 4: Track) ----
    sec3 = uipanel(g, 'Title', ...
        'HK Hierarchical Noise Scaling  |  Stage 4: Track  (used by Kalman_Advanced / Hierarchical-Kalman only)', ...
        'FontWeight', 'bold');
    g3 = uigridlayout(sec3, [2 3]);
    g3.ColumnWidth = {220, 110, '1x'};
    g3.RowHeight = {'fit','fit'};
    g3.Padding = [8 8 8 8];

    eAlpha = makeNumRow(g3, 'HK alpha (process-noise scale)', ...
        getDefault(p,'track.kalman.hk_alpha',0.01), ...
        'Multiplier on process noise in HK variant.');
    eBeta  = makeNumRow(g3, 'HK beta (measurement-noise scale)', ...
        getDefault(p,'track.kalman.hk_beta',0.025), ...
        'Multiplier on measurement noise in HK variant.');

    % ---- Live-update: store handles & wire callbacks for Kalman Gain panel ----
    app = guidata(fig);
    app.ui.HKAlpha = eAlpha;
    app.ui.HKBeta  = eBeta;
    guidata(fig, app);
    eAlpha.ValueChangedFcn = @(~,~) updateKalmanGainSummary(fig);
    eBeta.ValueChangedFcn  = @(~,~) updateKalmanGainSummary(fig);

    % ---- Buttons ----
    gb = uigridlayout(g, [1 3]);
    gb.ColumnWidth = {'1x', 120, 120};
    uilabel(gb, 'Text', '');
    uibutton(gb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) cancelModal());
    uibutton(gb, 'Text', 'Save', 'FontWeight', 'bold', ...
        'BackgroundColor', [0.2 0.6 0.2], 'FontColor', 'w', ...
        'ButtonPushedFcn', @(~,~) doSave());

    function doSave()
        saveParamState(fig);
        app2 = guidata(fig);
        if ~isfield(app2.data.params.track,'kalman'), app2.data.params.track.kalman = struct(); end
        app2.data.params.track.kalman.direction_penalty_weight    = eDirW.Value;
        app2.data.params.track.kalman.angle_penalty_slope         = eDirS.Value;
        app2.data.params.track.kalman.brightness_penalty_weight   = eBri.Value;
        app2.data.params.track.kalman.max_angle_change_deg        = eMaxA.Value;
        app2.data.params.track.kalman.gating_max_angle_change_deg = eGate.Value;
        app2.data.params.track.kalman.direction_history_points    = eHist.Value;
        app2.data.params.track.kalman.hk_alpha                    = eAlpha.Value;
        app2.data.params.track.kalman.hk_beta                     = eBeta.Value;
        if isfield(app2.ui, 'HKAlpha'), app2.ui = rmfield(app2.ui, 'HKAlpha'); end
        if isfield(app2.ui, 'HKBeta'),  app2.ui = rmfield(app2.ui, 'HKBeta');  end
        guidata(fig, app2);
        app3 = guidata(fig); app3.data = clearDownstreamData(app3.data, 4); app3 = manageGUIState(app3, 0); guidata(fig, app3);
        setStatus(app2, 'Advanced Kalman settings updated.', 'green');
        close(d);
    end

    function cancelModal()
        app_c = guidata(fig);
        if isfield(app_c.ui, 'HKAlpha'), app_c.ui = rmfield(app_c.ui, 'HKAlpha'); end
        if isfield(app_c.ui, 'HKBeta'),  app_c.ui = rmfield(app_c.ui, 'HKBeta');  end
        guidata(fig, app_c);
        close(d);
    end

    uiwait(d);
end

function openAdvancedLocGUI(fig)
    app = guidata(fig);
    p = app.data.params;

    d = uifigure('Name', 'Advanced Detection / Localization Settings', ...
        'Position', [200 120 720 740], 'WindowStyle', 'modal');
    g = uigridlayout(d, [4 1]);
    g.RowHeight = {'1x', '1x', '1x', 50};
    g.Padding = [12 12 12 12];
    g.RowSpacing = 10;

    % ---- Section 1: PSF Template (Stage 2: Detect) ----
    sec1 = uipanel(g, 'Title', ...
        'PSF Template  |  Stage 2: Detect  (used by NCC and Neyman-Pearson detectors)', ...
        'FontWeight', 'bold');
    g1 = uigridlayout(sec1, [3 3]);
    g1.ColumnWidth = {200, 220, '1x'};
    g1.RowHeight = {'fit','fit','fit'};
    g1.Padding = [8 8 8 8];

    uilabel(g1, 'Text', 'PSF type', 'FontWeight', 'bold');
    ePSFType = uidropdown(g1, 'Items', {'Gaussian','File','Measured'}, ...
        'Value', getDefault(p,'loc.psf_type','Gaussian'));
    uilabel(g1, 'Text', 'Analytic Gaussian, or load from file/measurement.', 'FontAngle', 'italic');

    uilabel(g1, 'Text', 'PSF size (pixels)', 'FontWeight', 'bold');
    psf_sz = getDefault(p,'loc.psf_size',[5 5]);
    ePSFSize = uieditfield(g1, 'text', 'Value', mat2str(psf_sz));
    uilabel(g1, 'Text', '[rows cols] for the template window (e.g. [5 5]).', 'FontAngle', 'italic');

    uilabel(g1, 'Text', 'PSF file path', 'FontWeight', 'bold');
    ePSFPath = uieditfield(g1, 'text', 'Value', getDefault(p,'loc.psf_file_path',''));
    uilabel(g1, 'Text', 'Only used if type = File.', 'FontAngle', 'italic');

    % ---- Section 2: QC Thresholds (Stage 2 / 3) ----
    sec2 = uipanel(g, 'Title', ...
        'Quality Control Thresholds  |  Stage 2: Detect & Stage 3: Localize', ...
        'FontWeight', 'bold');
    g2 = uigridlayout(sec2, [3 3]);
    g2.ColumnWidth = {200, 220, '1x'};
    g2.RowHeight = {'fit','fit','fit'};
    g2.Padding = [8 8 8 8];

    eROIMax = makeNumRow(g2, 'Max local maxima per ROI', ...
        getDefault(p,'loc.qc_max_roi_maxima',3), ...
        'Reject ROIs with more than this many peaks.');
    eGradSq = makeNumRow(g2, 'Min |grad|^2 for fit', ...
        getDefault(p,'loc.min_gradient_squared',1e-6), ...
        'Reject low-contrast fits (Radial Symmetry).');
    eMinDet = makeNumRow(g2, 'Min Hessian determinant', ...
        getDefault(p,'loc.min_determinant',1e-6), ...
        'Reject saddle points and flat extrema.');

    % ---- Section 3: Gaussian Fit QC (Stage 3: Localize) ----
    sec3 = uipanel(g, 'Title', ...
        'Gaussian Fit Quality  |  Stage 3: Localize  (used by Gaussian Fit and Gaussian Fit Fast only)', ...
        'FontWeight', 'bold');
    g3 = uigridlayout(sec3, [2 3]);
    g3.ColumnWidth = {200, 220, '1x'};
    g3.RowHeight = {'fit','fit'};
    g3.Padding = [8 8 8 8];

    eMinRSq = makeNumRow(g3, 'Min R-squared', ...
        getDefault(p,'loc.min_r_squared',0.3), ...
        'Minimum goodness-of-fit. Reject fits below this quality threshold.');

    % ---- Buttons ----
    gb = uigridlayout(g, [1 3]);
    gb.ColumnWidth = {'1x', 120, 120};
    uilabel(gb, 'Text', '');
    uibutton(gb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) close(d));
    uibutton(gb, 'Text', 'Save', 'FontWeight', 'bold', ...
        'BackgroundColor', [0.2 0.6 0.2], 'FontColor', 'w', ...
        'ButtonPushedFcn', @(~,~) doSave());

    function doSave()
        saveParamState(fig);
        app2 = guidata(fig);
        app2.data.params.loc.psf_type             = ePSFType.Value;
        try
            v = str2num(ePSFSize.Value); %#ok<ST2NM>
            if numel(v) == 2
                app2.data.params.loc.psf_size = v;
            elseif isscalar(v)
                app2.data.params.loc.psf_size = [v v];
            end
        catch
        end
        app2.data.params.loc.qc_max_roi_maxima    = eROIMax.Value;
        app2.data.params.loc.min_gradient_squared  = eGradSq.Value;
        app2.data.params.loc.min_determinant       = eMinDet.Value;
        app2.data.params.loc.min_r_squared         = eMinRSq.Value;
        guidata(fig, app2);
        % Sync the main-tab controls if they exist
        app2 = guidata(fig);
        if isfield(app2.ui, 'GaussMinRSquared') && isvalid(app2.ui.GaussMinRSquared)
            app2.ui.GaussMinRSquared.Value = eMinRSq.Value;
        end
        app3 = guidata(fig); app3.data = clearDownstreamData(app3.data, 2); app3 = manageGUIState(app3, 0); guidata(fig, app3);
        updateLocalizationOptions(fig);
        setStatus(app2, 'Advanced detection/localization settings updated.', 'green');
        close(d);
    end

    uiwait(d);
end

function openAdvancedRenderGUI(fig)
    app = guidata(fig);
    p = app.data.params;
    reg = getAlgorithmRegistry();

    d = uifigure('Name', 'Advanced Post-Processing / Render / Analysis', ...
        'Position', [200 160 740 720], 'WindowStyle', 'modal');
    g = uigridlayout(d, [4 1]);
    g.RowHeight = {'1x', '1x', '1x', 50};
    g.Padding = [12 12 12 12];
    g.RowSpacing = 10;

    % ---- Section 1: Track Smoothing (Stage 5) ----
    sec1 = uipanel(g, 'Title', ...
        'Track Smoothing  |  Stage 5: Post-process', 'FontWeight', 'bold');
    g1 = uigridlayout(sec1, [1 3]);
    g1.ColumnWidth = {200, 220, '1x'};
    g1.RowHeight = {'fit'};
    g1.Padding = [8 8 8 8];

    uilabel(g1, 'Text', 'Smoothing method', 'FontWeight', 'bold');
    smoothing_ids = registryIds(reg.smoothing);
    default_smooth = getDefault(p,'track.smoothing_method','sgolay');
    if ~any(strcmp(smoothing_ids, default_smooth)) && ~isempty(smoothing_ids)
        default_smooth = smoothing_ids{1};
    end
    eSmooth = uidropdown(g1, 'Items', smoothing_ids, 'Value', default_smooth);
    uilabel(g1, 'Text', 'Algorithm used to smooth position traces.', 'FontAngle', 'italic');

    % ---- Section 2: Rendering (Stage 6) ----
    sec2 = uipanel(g, 'Title', 'Rendering  |  Stage 6: Render', 'FontWeight', 'bold');
    g2 = uigridlayout(sec2, [3 3]);
    g2.ColumnWidth = {200, 220, '1x'};
    g2.RowHeight = {'fit','fit','fit'};
    g2.Padding = [8 8 8 8];

    uilabel(g2, 'Text', 'Interpolation method', 'FontWeight', 'bold');
    eInterp = uidropdown(g2, 'Items', {'spline','pchip','linear','makima'}, ...
        'Value', getDefault(p,'render.interpolation_method','spline'));
    uilabel(g2, 'Text', 'Sub-step interpolation between localizations.', 'FontAngle', 'italic');

    eSigma = makeNumRow(g2, 'Gaussian sigma (pixels)', ...
        getDefault(p,'render.gaussian_sigma',0.3), ...
        'Spread of each localization when using Gaussian splatting.');

    eStep = makeNumRow(g2, 'Interpolation step', ...
        getDefault(p,'render.interpolation_step',0.2), ...
        'Sub-step used during interpolation (0.2 ~= 5x sub-pixel density).');

    % ---- Section 3: Analysis ----
    sec3 = uipanel(g, 'Title', ...
        'Analysis  |  Post-rendering statistics (tortuosity, velocity histogram, density)', ...
        'FontWeight', 'bold');
    g3 = uigridlayout(sec3, [3 3]);
    g3.ColumnWidth = {200, 220, '1x'};
    g3.RowHeight = {'fit','fit','fit'};
    g3.Padding = [8 8 8 8];

    tort_bins = getDefault(p, 'analysis.tortuosity_bins', 0:0.05:8);
    if numel(tort_bins) >= 2
        tort_step_default = tort_bins(2) - tort_bins(1);
    else
        tort_step_default = 0.05;
    end
    eTort  = makeNumRow(g3, 'Tortuosity bin step', ...
        tort_step_default, ...
        'Bin width for tortuosity histogram.');
    eVBins = makeNumRow(g3, 'Velocity histogram bins', ...
        getDefault(p,'analysis.velocity_hist_num_bins',60), ...
        'Number of bins in velocity histogram.');
    eDens  = makeNumRow(g3, 'Density grid (mm)', ...
        getDefault(p,'analysis.density_grid_size_mm',0.5), ...
        'Cell size of density map used for statistics.');

    % ---- Buttons ----
    gb = uigridlayout(g, [1 3]);
    gb.ColumnWidth = {'1x', 120, 120};
    uilabel(gb, 'Text', '');
    uibutton(gb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) close(d));
    uibutton(gb, 'Text', 'Save', 'FontWeight', 'bold', ...
        'BackgroundColor', [0.2 0.6 0.2], 'FontColor', 'w', ...
        'ButtonPushedFcn', @(~,~) doSave());

    function doSave()
        saveParamState(fig);
        app2 = guidata(fig);
        app2.data.params.track.smoothing_method        = eSmooth.Value;
        app2.data.params.render.interpolation_method   = eInterp.Value;
        app2.data.params.render.gaussian_sigma         = eSigma.Value;
        app2.data.params.render.interpolation_step     = eStep.Value;
        app2.data.params.analysis.tortuosity_bins      = 0 : eTort.Value : 8;
        app2.data.params.analysis.velocity_hist_num_bins = eVBins.Value;
        app2.data.params.analysis.density_grid_size_mm = eDens.Value;
        guidata(fig, app2);
        app3 = guidata(fig); app3.data = clearDownstreamData(app3.data, 5); app3 = manageGUIState(app3, 0); guidata(fig, app3);
        setStatus(app2, 'Advanced render/analysis settings updated.', 'green');
        close(d);
    end

    uiwait(d);
end

function h = makeNumRow(parent, labelText, defaultVal, helpText)
    uilabel(parent, 'Text', labelText, 'FontWeight', 'bold');
    h = uieditfield(parent, 'numeric', 'Value', defaultVal);
    uilabel(parent, 'Text', helpText, 'FontAngle', 'italic', 'FontColor', [0.35 0.35 0.35]);
end

% =========================================================================
% SECTION: Cost Matrix GUI  (legacy, routed through openAdvancedKalmanGUI)
% =========================================================================

function openCostMatrixGUI(fig)
    openAdvancedKalmanGUI(fig);
end

function openHKConfigGUI(fig)
    openAdvancedKalmanGUI(fig);
end

% =========================================================================
% SECTION: Kalman Gain diagnostics
% =========================================================================
% If tracks already carry a `.kalman_gain` field (written by an upgraded
% tracker), that is used directly. Otherwise a tracker-agnostic proxy is
% computed from the residual between each position and the constant-
% velocity prediction from the two previous steps. The proxy saturates
% in [0,1] and matches the intuitive interpretation of gain: high gain ->
% trust the measurement; low gain -> trust the model.

function updateKalmanGainSummary(fig)
% Computes the theoretical Kalman gain K from the current parameter
% settings and updates the Trust Balance visual indicator.
    app = guidata(fig);
    if ~isfield(app.ui, 'KGainSummary') || ~isvalid(app.ui.KGainSummary), return; end
    if ~isfield(app.ui, 'KGainBarGrid'), return; end

    try
        p = app.data.params;

        % --- Read tracker method from the LIVE dropdown (not stale params) ---
        if isfield(app.ui, 'TrackMethodDropdown') && isvalid(app.ui.TrackMethodDropdown)
            method = app.ui.TrackMethodDropdown.Value;
        else
            method = lower(p.track.method);
        end

        % --- Read current Process Noise from the GUI field (live value) ---
        if isfield(app.ui, 'KalmanNoise') && isvalid(app.ui.KalmanNoise)
            currentProcessNoise = app.ui.KalmanNoise.Value;
        else
            currentProcessNoise = getDefault(p, 'track.kalman.process_noise', 0.1);
        end

        % --- Compute Q and R depending on tracker type ---
        reg = getAlgorithmRegistry();
        idx = find(strcmpi({reg.track.id}, method), 1);
        if isempty(idx)
            isHK = false;
        else
            isHK = reg.track(idx).usesHK;
        end

        if isHK
            % Hierarchical Kalman:  Q = alpha * v_max,  R = beta  (level 1)
            % --- Read live values from GUI fields if Advanced dialog is open ---
            if isfield(app.ui, 'HKAlpha') && isvalid(app.ui.HKAlpha)
                alpha = app.ui.HKAlpha.Value;
            else
                alpha = getDefault(p, 'track.kalman.hk_alpha', 0.01);
            end
            if isfield(app.ui, 'HKBeta') && isvalid(app.ui.HKBeta)
                beta = app.ui.HKBeta.Value;
            else
                beta  = getDefault(p, 'track.kalman.hk_beta',  0.025);
            end
            v_max = getDefault(p, 'track.kalman.hk_v_max', 20);
            Q = alpha * v_max;
            R = beta;
            formulaStr = sprintf( ...
                'K = (a * v_max) / (a * v_max + b)   |   a=%.4f   v_max=%.0f mm/s   b=%.4f', ...
                alpha, v_max, beta);
        else
            % Standard / v2 Kalman:  Q = process_noise,  R = (mean(fwhm)/2.355)^2
            Q = currentProcessNoise;
            fwhm_val = getDefault(p, 'loc.fwhm', [3 3]);
            sigma_loc = mean(fwhm_val) / 2.355;   % localization precision (px)
            R = sigma_loc^2;
            formulaStr = sprintf( ...
                'K = Q / (Q + R)   |   Q=%.4f   R=(FWHM/2.355)^2=%.4f', Q, R);
        end

        % --- Compute K ---
        if (Q + R) < eps
            K = 0.5;
        else
            K = Q / (Q + R);
        end
        K = max(0, min(1, K));  % clamp to [0, 1]

        modelPct = round((1 - K) * 100);
        measPct  = round(K * 100);

        % --- Update the split bar widths ---
        n_model = max(1, round(10 * (1 - K)));
        n_meas  = max(1, round(10 * K));
        app.ui.KGainBarGrid.ColumnWidth = { ...
            sprintf('%dx', n_model), sprintf('%dx', n_meas)};

        % --- Update bar labels ---
        app.ui.KGainModelBar.Text = sprintf('Model %d%%', modelPct);
        app.ui.KGainMeasBar.Text  = sprintf('Loc. %d%%', measPct);

        % --- Color gradient: more extreme → more saturated ---
        if K < 0.5
            app.ui.KGainModelBar.BackgroundColor = [0.15 0.60 0.30];
            app.ui.KGainMeasBar.BackgroundColor  = [0.50 0.65 0.80];
        else
            app.ui.KGainModelBar.BackgroundColor = [0.50 0.75 0.55];
            app.ui.KGainMeasBar.BackgroundColor  = [0.20 0.42 0.82];
        end

        % --- Update summary text ---
        app.ui.KGainSummary.Text = sprintf('Effective Kalman Gain:  K = %.3f', K);

        % --- Update formula text ---
        if isfield(app.ui, 'KGainFormula') && isvalid(app.ui.KGainFormula)
            app.ui.KGainFormula.Text = formulaStr;
        end

    catch ME
        app.ui.KGainSummary.Text = sprintf('Trust computation error: %s', ME.message);
    end
end


% =========================================================================
% SECTION: Tooltip system (built-in mini user-guide)
% =========================================================================

function applyTooltips(app)
    % Attach Tooltip text to every known UI control. Missing controls are
    % silently skipped so the tooltip system is fully decoupled from the
    % exact GUI layout.
    dict = getTooltipDictionary();
    fn = fieldnames(dict);
    for i = 1:numel(fn)
        name = fn{i};
        if isfield(app.ui, name) && ~isempty(app.ui.(name))
            h = app.ui.(name);
            try
                if isprop(h, 'Tooltip')
                    h.Tooltip = dict.(name);
                end
            catch
                % Tooltip attachment errors are cosmetic — ignore.
            end
        end
    end
end

function dict = getTooltipDictionary()
    dict = struct();

    % --- Acquisition / pixel size ---
    dict.TopFPSField    = 'Frame rate of the acquisition (Hz). Used to convert frame-to-frame motion into velocity. [Stage 4 / Stage 6]';
    dict.TopPixelXField = 'Lateral pixel size (mm). Sets the physical scale of the x-axis. [All stages]';
    dict.TopPixelZField = 'Axial pixel size (mm). Sets the physical scale of the z/depth axis. [All stages]';

    % --- Filter (Stage 1) ---
    dict.FilterMethodDropdown = 'Clutter-filter algorithm (Stage 1). SVD/SSM use a single global cutoff; DCC is decomposition-based; SVD_Blockwise adapts the cutoff per spatial block.';
    dict.SVDCutoffStart = 'First singular value to KEEP. Values below this are discarded as tissue clutter. [Stage 1]';
    dict.SVDCutoffEnd   = 'Last singular value to KEEP. Values above this are discarded as noise. [Stage 1]';

    dict.BWThresholdMethod = 'Block-wise threshold selection: MP (Marchenko-Pastur), DopplerGradient, Frequency, or Manual. [Stage 1]';
    dict.BWBlockSize       = 'Spatial block size (mm) for block-wise SVD. [Stage 1]';
    dict.BWOverlapPct      = 'Overlap between blocks (%). Higher overlap = smoother, slower. [Stage 1]';
    dict.BWTissueFreq      = 'Tissue cutoff frequency (Hz), used by the Frequency method (-1 = auto). [Stage 1]';
    dict.BWMPSigma         = 'Std-devs beyond the MP bulk edge to treat as signal. [Stage 1]';
    dict.BWGradientPct     = 'Fractional drop in singular-value slope used as threshold. [Stage 1]';
    dict.BWMinBlood        = 'Minimum number of components to keep per block. [Stage 1]';
    dict.BWMaxTissueFrac   = 'Max fraction of components discarded as tissue. [Stage 1]';
    dict.BWPlotMaps        = 'Plot per-block threshold-selection maps after filtering. [Stage 1]';
    dict.BWManualCutoff    = 'Manual cutoff [start end], used only when Method = Manual. [Stage 1]';

    dict.EnableButterworth = 'Optional temporal Butterworth filter applied after SVD/DCC. [Stage 1]';
    dict.ButterCutoff      = 'Cutoff frequency (Hz). Scalar = low/high-pass; [lo hi] = band-pass. [Stage 1]';
    dict.ButterOrder       = 'Order of the Butterworth filter (typical 2-4). [Stage 1]';
    dict.SpatialMethodDrop = 'Optional spatial filter (Gaussian / Median / None) applied frame-by-frame. [Stage 1]';

    % --- Detection (Stage 2) ---
    dict.DetectMethodDropdown = 'Detection algorithm (Stage 2). Intensity = local maxima; NCC = PSF correlation; Neyman-Pearson = matched-filter test.';
    dict.DetectFWHMField      = 'Expected FWHM of a microbubble PSF, in pixels [rows cols]. [Stage 2]';
    dict.LocThreshField       = 'Intensity threshold for peak acceptance. [Stage 2]';
    dict.LocMaxField          = 'Maximum number of detections per frame. [Stage 2]';
    dict.LocShiftFactor       = 'Max allowed sub-pixel shift from integer pixel peak. Used by ALL localizers. [Stage 3 QC]';
    dict.LocQCDivergence      = 'Reject fits whose sub-pixel center diverges from the candidate peak. Used by ALL localizers. [Stage 3 QC]';
    dict.LocQCRoiMaxima       = 'Reject ROIs containing multiple local maxima. [Stage 2 QC]';
    dict.NP_AlphaField        = 'Neyman-Pearson significance level alpha0. [Stage 2]';
    dict.NCC_ThreshField      = 'NCC correlation threshold for peak acceptance. [Stage 2]';
    dict.PeakContrastField    = 'H-maxima contrast: minimum dip between two peaks to report them as separate bubbles. 0 = strict mode (original), 0.02–0.10 = dense-field mode. [Stage 2]';

    % --- Localization (Stage 3) ---
    dict.LocMethodDropdown = 'Localization algorithm (Stage 3). Radial Symmetry is fastest; Gaussian Fit is most accurate.';
    dict.GaussMinRSquared  = 'Minimum R-squared for an accepted Gaussian fit. Higher = stricter QC, fewer but more accurate localizations. [Stage 3 QC]';
    dict.p_loc_qc_radial   = 'Common QC filters applied to all localization methods (divergence, shift, ROI maxima).';
    dict.p_loc_qc_gauss    = 'Gaussian-specific fitting parameters. Only active when using a Gaussian localizer.';

    % --- Track (Stage 4) ---
    dict.TrackMethodDropdown = 'Tracking algorithm (Stage 4). Kalman variants predict motion; Hungarian/NN are assignment-only.';
    dict.MaxDistField        = 'Maximum link distance between consecutive frames (pixels). [Stage 4]';
    dict.GapFramesField      = 'Maximum missed frames before a track is closed. [Stage 4]';
    dict.MinLengthField      = 'Minimum track length (frames) to keep. [Stage 4]';
    dict.KalmanModelDrop     = 'Kalman motion model: ConstantVelocity or ConstantAcceleration. [Stage 4]';
    dict.KalmanNoise         = 'Kalman process noise scalar. Higher = trust measurements more. [Stage 4]';
    dict.AssignmentDrop      = 'Data-association algorithm inside the tracker. [Stage 4]';
    dict.UseAdvancedCostCheckbox = 'Use direction, brightness and angle-gating in the cost matrix (otherwise pure distance). [Stage 4]';

    dict.TrackQCDirection    = 'Reject tracks whose direction jumps exceed the max angle. [Stage 4 QC]';
    dict.QCMaxAngle          = 'Max allowed direction change per frame (deg). [Stage 4 QC]';
    dict.TrackQCAcceleration = 'Reject tracks with large acceleration outliers. [Stage 4 QC]';
    dict.QCAccelFactor       = 'Reject frames whose acceleration exceeds median by this factor. [Stage 4 QC]';
    dict.TrackQCVD           = 'Reject tracks with velocity outliers vs local median. [Stage 4 QC]';
    dict.QCVDRatio           = 'Tolerance ratio for the velocity QC check. [Stage 4 QC]';

    % --- Post-process / Render (Stage 5 / 6) ---
    dict.EnablePostProcessing  = 'Apply smoothing and interpolation to tracks before rendering. [Stage 5]';
    dict.SmoothField           = 'Smoothing window length (odd). 0 disables smoothing. [Stage 5]';
    dict.DisplayMinLengthField = 'Minimum track length for display/rendering. [Stage 5/6]';
    dict.UpsamplingField       = 'Super-resolution upsampling factor for the render grid. [Stage 6]';
    dict.RenderMethodDrop      = 'Render method: histogram (count) or gaussian_splat (soft blobs). [Stage 6]';

    % --- Crop & workspace ---
    dict.CropBoxField       = 'Crop rectangle [x z w h] in pixels. Leave empty to use the full frame.';
    dict.InteractiveCropBtn = 'Draw a crop rectangle interactively on the display.';
    dict.ApplyCropBtn       = 'Apply the current crop rectangle to raw data (speeds up SVD).';
    dict.LoadCropBtn        = 'Load a previously saved crop rectangle.';

    % --- Buttons ---
    dict.btnUndo        = 'Undo the last parameter change.';
    dict.btnRedo        = 'Redo the last undone change.';
    dict.btnLoadSession = 'Load a previously saved session (data + parameters).';
    dict.btnSaveSession = 'Save the current session (data + parameters) to a .mat file.';

    dict.BtnAdvDetect        = 'Advanced detection / localization parameters (PSF template, QC thresholds). [Stage 2/3]';
    dict.BtnAdvLoc           = 'Advanced detection / localization parameters (PSF template, QC thresholds). [Stage 2/3]';
    dict.BtnAdvKalman        = 'Advanced Kalman settings (cost matrix weights, angle gating, HK noise scaling). [Stage 4]';
    dict.BtnAdvRender        = 'Advanced post-processing / render / analysis settings. [Stage 5/6]';
    dict.BtnConfigCostMatrix = 'Configure the advanced cost matrix (only used when "Use advanced cost" is enabled). [Stage 4]';
    dict.BtnConfigHK         = 'Configure Hierarchical-Kalman (HK) noise-scaling parameters. [Stage 4]';

    % --- Runtime controls ---
    dict.LoadButton            = 'Load raw data from a .mat file.';
    dict.RunFilterButton       = 'Run the selected filter on the raw data (Stage 1).';
    dict.RunDetectButton       = 'Run the detection + localization stages (Stages 2+3).';
    dict.GenerateImagesButton  = 'Render the final super-resolution images from the tracks.';
    dict.PlayPauseButton       = 'Play / pause the frame animation.';
    dict.FrameSlider           = 'Scrub through frames of the current data source.';
    dict.FrameField            = 'Current frame index (synced with the slider).';

    % --- ROI / Mask tab ---
    dict.EnhanceMethodDrop    = 'Vessel-enhancement filter method (e.g. Frangi).';
    dict.EnhanceAmountSlider  = 'Vessel-enhancement strength.';
    dict.ROIContrastSlider    = 'Display gamma for the ROI preview.';
    dict.ROIContrastField     = 'Display gamma (paired with slider).';
    dict.ROIThreshSlider      = 'Binary threshold for the ROI mask preview (0..1, fraction of max).';
    dict.ROIThreshField       = 'Binary threshold for the ROI mask preview (paired with slider).';
    dict.CreateMaskButton     = 'Commit the current ROI preview as the active detection mask.';
    dict.LoadMaskButton       = 'Load a binary mask from a .mat file.';
    dict.ResetMaskButton      = 'Clear the active mask.';
    dict.ResetParamsButton    = 'Reset all parameters to defaults from setDefaultParams.m.';

    % --- Kalman trust balance panel ---
    dict.KGainSummary  = 'Effective Kalman Gain K: 0 = trusts motion model fully, 1 = trusts localizations fully.';
    dict.p_kgain       = 'Shows the theoretical trust split between the motion model and the raw localizations, based on your current noise settings.';
    dict.KGainFormula  = 'The formula used to compute the effective gain from the current process-noise and measurement-noise parameters.';
    dict.KGainModelBar = 'Green bar: fraction of trust placed on the motion-model prediction.';
    dict.KGainMeasBar  = 'Blue bar: fraction of trust placed on the raw localization measurements.';
end

% =========================================================================
% SECTION: Rendering helpers
% =========================================================================

function add_scale_bar(params, imgSize, lw, fs)
    % Draws a 1 mm scale bar in the lower-right corner of the current axes.
    try
        pxX = getDefault(params, 'track.pixel_X_size', ...
                         getDefault(params, 'expParams.pixel_X_size', 0.05));
        upFac = getDefault(params, 'render.upsampling_factor', 1);
        pxRender = pxX / max(upFac, 1);
        if pxRender <= 0 || ~isfinite(pxRender), return; end

        bar_mm = 1.0;
        bar_px = bar_mm / pxRender;

        x0 = imgSize(2) - bar_px - 0.05 * imgSize(2);
        y0 = imgSize(1) - 0.06 * imgSize(1);
        hold on;
        plot([x0, x0 + bar_px], [y0, y0], '-w', 'LineWidth', lw);
        text(x0 + bar_px/2, y0 - 0.02 * imgSize(1), '1 mm', ...
            'Color', 'w', 'HorizontalAlignment', 'center', ...
            'FontSize', fs, 'FontWeight', 'bold');
        hold off;
    catch
        % Scale bar is decorative — never abort rendering on failure.
    end
end

function out = applySpatialFilter(frame, filterParams)
% Apply spatial conditioning to a single frame using setDefaultParams
% spatial_method / spatial_kernel / spatial_sigma1 / spatial_sigma2.
    method = getDefault(filterParams, 'spatial_method', 'None');
    switch lower(char(method))
        case 'gaussian'
            sigma = getDefault(filterParams, 'spatial_sigma1', 1.0);
            out = imgaussfilt(frame, sigma);
        case 'dog'
            s1 = getDefault(filterParams, 'spatial_sigma1', 1.0);
            s2 = getDefault(filterParams, 'spatial_sigma2', 2.0);
            out = imgaussfilt(frame, s1) - imgaussfilt(frame, s2);
        case 'median'
            n = getDefault(filterParams, 'spatial_kernel', 3);
            n = max(1, round(n));
            out = medfilt2(frame, [n n]);
        case 'top-hat'
            n = getDefault(filterParams, 'spatial_kernel', 3);
            n = max(1, round(n));
            se = strel('disk', n);
            out = imtophat(frame, se);
        otherwise
            out = frame;
    end
end

% =========================================================================
% SECTION: Fallback parameter struct
% =========================================================================
% If setDefaultParams fails (e.g. path issues), we build a minimal but
% complete parameter struct here. ensureAllParamFields will then fill in
% anything still missing. The field layout mirrors setDefaultParams.m.

function p = createFallbackParams()
    p = struct();

    % -------- IO --------
    p.io = struct('data_folder', '', 'data_subfolder', '', ...
        'file_pattern', '*.mat', 'save_mat_file', true, ...
        'export_tiff', true, 'export_csv', false, ...
        'save_lightweight', true, 'export_figures', true, ...
        'cropBox', []);

    % -------- Proc / mask --------
    p.proc = struct('enableInteractiveMask', false, ...
        'generateNewMask', false, 'maskPath', '', ...
        'enableInteractiveCrop', false, 'generateNewCrop', false, ...
        'cropPath', '', 'ROIMask', [], ...
        'vesselMask', struct('enable', false, 'method', 'frangi', ...
            'strength', 1.0, 'gamma', 1.0, 'threshold', 0.3));

    % -------- Acquisition --------
    p.acq = struct('framerate', 200);

    % -------- Filter --------
    p.filter = struct();
    p.filter.method = 'svd_filter';
    p.filter.svd_cutoff = [5 100];
    p.filter.enable_butterworth = false;
    p.filter.butter_cutoff = [50 250];
    p.filter.butter_order = 2;
    p.filter.butter_cutoff_norm = 0.05;
    p.filter.spatial_method = 'Gaussian';
    p.filter.spatial_kernel = 3;
    p.filter.spatial_sigma1 = 1.0;
    p.filter.spatial_sigma2 = 2.0;
    p.filter.blockwise = struct( ...
        'threshold_method',    'DopplerGradient', ...
        'block_size_mm',        4.0, ...
        'overlap_pct',          75.0, ...
        'manual_cutoff',        [10 200], ...
        'tissue_freq_hz',       -1, ...
        'mp_deviation_sigma',   2.0, ...
        'gradient_pct',         0.10, ...
        'min_blood_comps',      3, ...
        'max_tissue_frac',      0.60, ...
        'plot_maps',            false);
    p.filter.dcc = struct('tissue_start',1,'tissue_end',10, ...
                          'blood_start',11,'blood_end',150, ...
                          'noise_start',151,'noise_end',500);

    % -------- Localization / detection --------
    p.loc = struct();
    p.loc.method = 'radial';
    p.loc.DetectMethod = 'Intensity';
    p.loc.detection_threshold = 0.5;
    p.loc.max_bubbles_per_frame = 200;
    p.loc.fwhm = [3 3];
    p.loc.NP_alpha0 = 0.01;
    p.loc.crosscor_threshold = 0.7;
    p.loc.psf_type = 'Gaussian';
    p.loc.psf_size = [5 5];
    p.loc.psf_file_path = '';
    p.loc.enable_divergence_check = true;
    p.loc.qc_max_shift_factor = 1;
    p.loc.enable_roi_maxima_check = true;
    p.loc.qc_max_roi_maxima = 3;
    p.loc.min_gradient_squared = 1e-6;
    p.loc.min_determinant = 1e-6;
    p.loc.min_r_squared = 0.3;

    % -------- Tracking --------
    p.track = struct();
    p.track.method = 'kalman';
    p.track.pixel_X_size = 0.05;
    p.track.pixel_Z_size = 0.05;
    p.track.dt = 1 / p.acq.framerate;
    p.track.min_track_length = 8;
    p.track.max_gap_closing_frames = 2;
    p.track.max_linking_distance = 5;
    p.track.use_advanced_cost_matrix = false;
    p.track.smoothing_method = 'sgolay';
    p.track.smoothing_factor = 5;
    p.track.enable_postprocessing = true;
    p.track.display_min_length = 10;

    p.track.kalman = struct( ...
        'motion_model', 'ConstantVelocity', ...
        'process_noise', 0.1, ...
        'assignment_method', 'hungarian', ...
        'use_direction', true, ...
        'use_angle', true, ...
        'use_brightness', true, ...
        'direction_penalty_weight', 2, ...
        'angle_penalty_slope', 0.3, ...
        'brightness_penalty_weight', 2, ...
        'direction_history_points', 4, ...
        'max_angle_change_deg', 70, ...
        'gating_max_angle_change_deg', 90, ...
        'hk_alpha', 0.01, ...
        'hk_beta', 0.025, ...
        'hk_forward_backward', true, ...
        'hk_v_max', 20, ...
        'hk_num_levels', 5, ...
        'hk_spacing_power', 1.0, ...
        'hk_enable_overlap', true, ...
        'hk_overlap_mm_s', 2.0);

    p.track.qc = struct( ...
        'enable_direction_constraint', true, ...
        'max_angle_change_deg', 60, ...
        'enable_acceleration_constraint', true, ...
        'acceleration_C_factor', 3, ...
        'enable_vd_constraint', false, ...
        'max_vd_ratio', 0.5);

    % -------- Render --------
    p.render = struct('method', 'histogram', 'upsampling_factor', 10, ...
        'interpolation_method', 'spline', 'gaussian_sigma', 0.3, ...
        'interpolation_step', 0.2);

    % -------- Analysis --------
    p.analysis = struct( ...
        'tortuosity_bins', 0 : 0.05 : 8, ...
        'velocity_hist_num_bins', 60, ...
        'density_grid_size_mm', 0.5);

    % -------- Exp params (sane defaults) --------
    p.expParams = struct('C', 1540, 'size', [], ...
        'pixel_X_size', 0.05, 'pixel_Z_size', 0.05, ...
        'fovX', NaN, 'fovZ', NaN, 'frequency', NaN, 'lambda', NaN, ...
        'shootingRate', 200, 'flowSpeed', NaN, ...
        'channel_cross_section_mm2', NaN);
end