% =========================================================================
% FILENAME: setDefaultParams.m
% =========================================================================
%
% PURPOSE:
%   Defines and returns a structure containing all default parameters for
%   the ULM framework. This function is designed for easy and quick updates
%   between experiments.
%
%   MODIFIED: This function now also loads experiment-specific files
%   (info.txt, sample data) to calculate derived parameters and set up
%   the mask, making it the single source of truth for all parameters.
%
% HOW TO USE:
%   1. Modify the values in the "QUICK EXPERIMENT SETTINGS" section below
%      for your specific experiment. This is where you set the data path
%      and choose the main processing methods and key tuning parameters.
%   2. The "DETAILED PARAMETER STRUCTURE" section below will automatically
%      use your quick settings. You only need to edit the detailed section
%      if you want to fine-tune less common parameters.
%
% AUTHOR: Grigori Shapiro
% DATE: July 28, 2025

function params = setDefaultParams(isKidneyExperiment, user_data_folder)

    % =========================================================================
    %                  QUICK EXPERIMENT SETTINGS
    % =========================================================================
    % --- Define the most frequently changed parameters here for easy access ---
    
    % --- 1a. Core Settings (Paths & Methods) ---
    if nargin < 2 || isempty(user_data_folder)
        user_data_folder = '';  % <-- Optionally hardcode your path here:
                                % e.g. user_data_folder = 'C:\MyData\Experiment_01';
        
        if isempty(user_data_folder)
            user_data_folder = uigetdir(pwd, 'Select your experiment data folder');
            if user_data_folder == 0
                error('ULM:NoFolderSelected', ...
                    'No data folder selected. Please provide a path or select one via the dialog.');
            end
        end
    end
    
    % Final validation — catches empty string passed from outside,
    % non-existent paths, and typos regardless of how the path was obtained.
    if isempty(user_data_folder) || ~exist(user_data_folder, 'dir')
        error('ULM:InvalidFolder', ...
            'Data folder does not exist or is empty: "%s"\nCheck the path and try again.', ...
            user_data_folder);
    end
    
    user_data_subfolder             = 'Bmode'; % 'Bmode'
    user_clutter_filter_method      = 'svd_filter';    % 'svd_filter', 'svd_ssm', 'dcc_svd', 'svd_blockwise'
    user_cutoff_svd                 = [4, 450];      % SVD Cut off fre
    user_enable_butterworth         = false;           % use butterworth filter

    % --- Block-wise SVD Filter Settings ---
    user_bsvd_threshold_method  = 'DopplerGradient'; % 'DopplerGradient','SSM','Hybrid','Manual'
    user_bsvd_block_size_mm     = 4.0;    % Block size in mm (square). Set [] to auto.
                                          % Auto-rule: smallest 10px-multiple with blk^2 > T.
                                          % Song et al. optimal: ~4mm for 0.05mm/px kidney data.
    user_bsvd_overlap_pct       = 75.0;   % Overlap %. Use 93.75 for publication quality.
                                          % 75% gives ~4x speedup vs 93.75% (paper optimal).
    user_bsvd_manual_cutoff     = [10, 200]; % Only used if ThresholdMethod='Manual'
    user_bsvd_tissue_freq_hz    = -1;     % Tissue Doppler freq threshold [Hz]. -1 = auto.
                                          % Auto: max(5, min(20, framerate/50)).
                                          % Adjust up (e.g. 15) for fast-moving tissue.
    user_bsvd_mp_deviation_sigma = 2.0;   % Marchenko-Pastur sensitivity (high cutoff).
                                          % Higher = keeps fewer components as blood.
    user_bsvd_gradient_pct      = 0.10;   % Cutoff 1A inflection sensitivity.
                                          % Lower = detects inflection earlier.
    user_bsvd_min_blood_comps   = 3;      % Floor: min blood components per block.
    user_bsvd_max_tissue_frac   = 0.60;   % Ceiling: max fraction of K as tissue.
    user_bsvd_plot_maps         = false;  % Show spatial threshold maps after filtering.
    
    user_spatial_method             = 'None';    % 'None', 'Gaussian', 'Median', 'DoG', 'Top-Hat'
    user_spatial_kernel             = 3;             % Kernel size (px)
    user_spatial_sigma1             = 0.5;           % Sigma 1
    user_spatial_sigma2             = 4.0;           % Sigma 2 (for DoG)
    
    user_localization_method        = 'radial';          % 'radial', 'gaussian_fit', 'gaussian_fit_fast'
    user_tracking_method            = 'Kalman';          % 'Hungarian', 'nn', 'Kalman', 'Kalman_Advanced'
    user_rendering_method           = 'histogram';       % 'histogram', 'gaussian'

    % --- 1b. Data Pre-processing (Masking & Cropping) ---
    user_enableInteractiveMask      = false;     % Master switch for interactive mask feature
    user_generateNewMask            = false;     % true: draw a new mask; false: load existing 'Mask.mat'
    user_maskPath                   = fullfile(user_data_folder, 'Results'); % Default path to save/load mask
    
    user_enable_vessel_masking      = false;       % Master switch for algorithmic masking
    user_vmask_method               = 'Top-Hat';   % 'None', 'CLAHE', 'Top-Hat', 'Sharpen'
    user_vmask_strength             = 1;           % 0 to 1
    user_vmask_gamma                = 1;           % Contrast
    user_vmask_threshold            = 0.081;       % Noise Threshold
    
    user_enableInteractiveCrop      = false;     % Master switch for interactive cropping feature
    user_generateNewCrop            = false;     % true: draw a new crop rect; false: load existing 'cropBox.mat'
    user_cropPath                   = fullfile(user_data_folder, 'Results'); % Default path to save/load crop rectangle

    % --- 2. Localization & Detection Settings ---
    user_detection_threshold        = 0.15;      % Sensitivity for bubble detection (0.0 to 1.0).
    user_max_bubbles_per_frame      = 100;       % Max candidate bubbles to process per frame.
    user_fwhm_pixels                = [1, 1]*3;  % Estimated bubble PSF size [x, z] in pixels.
                                                 % Set to NaN to auto-compute from acoustic wavelength: fwhm = (floor(λ/pixel_size)*2+1)
                                                 % Set to [x, z] (pixels) to override with a fixed value. A fixed odd integer (e.g., 3)
                                                 % is often more stable than the wavelength-derived value.
                                                 
    % Advanced Detection Settings (NP & NCC) ---
    user_detect_method              = 'Intensity'; % Options: 'Intensity', 'NP', 'NCC'
    user_np_alpha0                  = 1e-4;        % False alarm rate for NP method (0 < alpha < 0.5)
    user_ncc_threshold              = 0.7;         % Minimum correlation threshold for NCC method
    user_psf_type                   = 'Gaussian';  % Options: 'Gaussian', 'Experimental'
    user_psf_size                   = [5, 5];      % Kernel size for NCC template (MUST be odd numbers)
    user_psf_file_path              = fullfile(user_data_folder, 'Results', 'Experimental_PSF.mat'); % Path to custom PSF

    % --- 3. Core Tracking Settings ---
    user_min_track_length           = 5;     % Minimum points for a valid track.
    user_max_gap_frames             = 1;     % Allowed frames to bridge gaps in tracks.
    user_max_linking_distance_px    = 2;     % For 'hungarian' & 'kalman' trackers only. 
                                             % Set to NaN to auto-compute from flow physics: distance = v_max * dt / pixel_size
                                             % Set to a number (pixels) to use a fixed override. Recommended: 2 px for most setups.
    user_channel_cross_section_mm2 = NaN;    % Cross-sectional area of the main flow channel [mm²].
                                             % Set to NaN to auto-compute: rectangular cross-section = width × 0.3 mm (fixed height).
                                             % width is read from info.txt as mainChannelDiameter [um].
                                             % Set to an explicit value [mm²] to override for any other geometry.
                                             %   Examples:
                                             %     Rectangular 300×300 um:  user_channel_cross_section_mm2 = 0.3 * 0.3;
                                             %     Rectangular 500×300 um:  user_channel_cross_section_mm2 = 0.5 * 0.3;
                                             %     Circular r=150 um:       user_channel_cross_section_mm2 = pi * 0.15^2;

    user_enable_postprocessing      = true;  % Master switch
    user_use_advanced_cost_matrix   = false;  

    % --- 4a. Weighted Cost Matrix Settings ---
    user_direction_weight    = 2;     % W_{dir}    % How much to penalize turns (e.g., 0-5).
    user_angle_penalty_slope = 0.3;   % W_slop     % e.g., 0.5
    user_brightness_weight   = 2;     % W_{int}    % How much to penalize brightness changes (e.g., 0-3).
    user_max_angle_gate_deg  = 70.0;               % soft limit for angle change during tracking.
    user_gating_max_angle_gate_deg = 90;          % Hard limit for angle change during tracking.
    user_track_direction_history_points = 4; % 3-5
    % P_{angle} = angle - user_max_angle_gate_deg; above user_gating_max_angle_gate_deg the cost is inf;
    % P_{intensity} = (currentBrightnesses - avg_track_brightness) / (avg_track_brightness + eps);
    % C_{total}(i,j) = C_{dist} * (1 + W_{dir} * W_slop * P_{angle}) * (1 + W_{int} * P_{intensity})

    % --- 4b. Kalman Filter General Settings ---
    user_assignment_method          = 'hungarian'; % hungarian or nn (Nearest Neighbor)
    user_kalman_motion_model        = 'ConstantVelocity'; % 'ConstantVelocity' or 'ConstantAcceleration'
    user_kalman_process_noise       = 3;                 % Kalman adaptability for standard tracker (e.g., 5-50).
    
    % --- Settings for Advanced/Hierarchical Kalman Tracker ---
    % =========================================================================
    % HIERARCHICAL KALMAN TUNING GUIDE (Process & Measurement Noise)
    % -------------------------------------------------------------------------
    % 1. hk_alpha (Process Noise multiplier)
    %    Formula: process_noise = hk_alpha * v_max_mm_s
    %    Meaning: Controls model flexibility (how much we trust the physical motion model).
    %    * INCREASE (e.g., 0.01 to 0.1) if tracks break on sharp turns (filter is too stiff).
    %    * DECREASE (e.g., 0.001) if tracks jump wildly to random noise (filter is too loose).
    %
    % 2. hk_beta (Measurement Noise multiplier)
    %    Formula: measurement_noise = hk_beta / (2^(level - 1))
    %    Meaning: Controls localization trust (decreases exponentially for faster levels).
    %    * INCREASE if the MB localizations are jittery/noisy (smooths the track).
    %    * DECREASE if localizations are highly precise (forces track to follow points exactly).
    % =========================================================================
    user_kalman_hk_alpha            = 0.01; % Process noise scaling factor (alpha * v_max)
    user_kalman_hk_beta             = 0.025; % Measurement noise scaling factor (beta / level)
    user_kalman_hk_forward_backward = true;  % Enable dual-pass tracking for maximum yield
    user_kalman_hk_v_max            = 20;       % [mm/s] Max velocity to track
    user_kalman_hk_num_levels       = 5;        % Number of levels
    user_kalman_hk_spacing_power    = 1.0;      % 1.0 = linear spacing.
                                                % >1.0 = more levels at slow velocities.
                                                % <1.0 = more levels at fast velocities.
    user_kalman_hk_enable_overlap   = true;    % false = contiguous, no overlap
                                                % true  = add overlap band between adjacent levels
    user_kalman_hk_overlap_mm_s     = 2 ;        % [mm/s] Total overlap band width (only used if enable_overlap = true)
    
    
    % --- 5. Quality Control (QC) Settings ---
    % --- QC for Localization ---
    user_loc_qc_enable_divergence   = true;     % Reject localizations that shift too far from the initial guess.
    user_loc_qc_enable_roi_maxima   = true;     % Reject ROIs with multiple peaks.
    % --- QC for Tracks (Post-Processing) ---
    user_track_qc_enable_direction  = false;     % Enable direction change constraint.
    user_track_qc_max_angle_deg     = 150;     % Max allowed turn in degrees.
    user_track_qc_enable_acceleration = false;   % Enable adaptive acceleration constraint.
    user_track_qc_acceleration_factor = 0.75;    % Tuning factor 'C' for acceleration threshold.
    user_track_qc_enable_vd         = false;     % Enable Velocity Dispersion (straightness) constraint.
    user_track_qc_max_vd_ratio      = 3;      % Max ratio of path length to net displacement.
    
    % --- 6. Rendering & Output Settings ---
    user_smoothing_factor           = 11;               % smoothing window size
    user_track_smoothing_method     = 'sgolay';  % Options: 'rloess', 'sgolay', 'gaussian', 'movmean'
    user_track_interpolation_method = 'spline';  % Options: 'spline', 'pchip', 'linear', 'makima'

    user_upsampling_factor          = 10;               % Final image upsampling factor.
    user_interpolation_step = 1/user_upsampling_factor; % Step for interpolation (0.2 means 5x more points)

    % =========================================================================
    %                DETAILED PARAMETER STRUCTURE
    % =========================================================================
    % This section populates the final 'params' struct from the settings above.

    %% 1. General, I/O, and Pre-processing Parameters
    params.io.data_folder = user_data_folder;
    params.io.data_subfolder = user_data_subfolder;
    params.io.file_pattern = '*.mat';
    params.io.save_mat_file = true;
    params.io.export_tiff = true;
    params.io.export_csv = false;
    params.io.save_lightweight = true;   % Recommended for memory-efficient processing.
    params.io.export_figures = true;     % Automatically save .fig and .png of the histograms

    % --- Pre-processing Masking and Cropping ---
    params.proc.enableInteractiveMask = user_enableInteractiveMask;
    params.proc.generateNewMask = user_generateNewMask;
    params.proc.maskPath = user_maskPath;

    % --- Algorithmic Vessel Masking ---
    params.proc.vesselMask.enable = user_enable_vessel_masking;
    params.proc.vesselMask.method = user_vmask_method;
    params.proc.vesselMask.strength = user_vmask_strength;
    params.proc.vesselMask.gamma = user_vmask_gamma;
    params.proc.vesselMask.threshold = user_vmask_threshold;

    params.proc.enableInteractiveCrop = user_enableInteractiveCrop;
    params.proc.generateNewCrop = user_generateNewCrop;
    params.proc.cropPath = user_cropPath;

    %% 2. Clutter Filtering Parameters
    params.filter.method = user_clutter_filter_method;
    params.filter.svd_cutoff = user_cutoff_svd;
    params.filter.enable_butterworth = user_enable_butterworth;
    params.filter.butter_cutoff = [50 250];
    params.filter.butter_order = 2;
    params.filter.butter_cutoff_norm = 0.05;

    params.filter.spatial_method = user_spatial_method;
    params.filter.spatial_kernel = user_spatial_kernel;
    params.filter.spatial_sigma1 = user_spatial_sigma1;
    params.filter.spatial_sigma2 = user_spatial_sigma2;

    %% 2b. Block-wise SVD Filter Parameters
    params.filter.blockwise.threshold_method   = user_bsvd_threshold_method;
    params.filter.blockwise.block_size_mm      = user_bsvd_block_size_mm;
    params.filter.blockwise.overlap_pct        = user_bsvd_overlap_pct;
    params.filter.blockwise.manual_cutoff      = user_bsvd_manual_cutoff;
    params.filter.blockwise.tissue_freq_hz     = user_bsvd_tissue_freq_hz;
    params.filter.blockwise.mp_deviation_sigma = user_bsvd_mp_deviation_sigma;
    params.filter.blockwise.gradient_pct       = user_bsvd_gradient_pct;
    params.filter.blockwise.min_blood_comps    = user_bsvd_min_blood_comps;
    params.filter.blockwise.max_tissue_frac    = user_bsvd_max_tissue_frac;
    params.filter.blockwise.plot_maps          = user_bsvd_plot_maps;

    %% 3. Bubble Detection & Localization Parameters
    params.loc.method = user_localization_method;
    params.loc.detection_threshold = user_detection_threshold;
    params.loc.max_bubbles_per_frame = user_max_bubbles_per_frame;
    params.loc.fwhm = user_fwhm_pixels;

    % Map Advanced Detection Parameters ---
    params.loc.DetectMethod = user_detect_method;
    params.loc.NP_alpha0 = user_np_alpha0;
    params.loc.crosscor_threshold = user_ncc_threshold;
    params.loc.psf_type = user_psf_type;
    params.loc.psf_size = user_psf_size;
    params.loc.psf_file_path = user_psf_file_path;
    
    % --- Detailed Localization QC ---
    params.loc.enable_divergence_check = user_loc_qc_enable_divergence;
    params.loc.qc_max_shift_factor = 1;
    params.loc.enable_roi_maxima_check = user_loc_qc_enable_roi_maxima;
    params.loc.qc_max_roi_maxima = 3;
    params.loc.min_gradient_squared = 1e-6;
    params.loc.min_determinant = 1e-6;

    params.loc.gauss_fit_box_radius = 2;
    params.loc.min_r_squared = 0.3;

    %% 4. Bubble Tracking Parameters
    params.track.method = user_tracking_method;
    params.track.min_track_length = user_min_track_length;
    params.track.max_gap_closing_frames = user_max_gap_frames;
    params.track.max_linking_distance = user_max_linking_distance_px;
    params.track.use_advanced_cost_matrix = user_use_advanced_cost_matrix;
    
    % --- Detailed Kalman Settings ---
    params.track.kalman.assignment_method = user_assignment_method;
    params.track.kalman.motion_model = user_kalman_motion_model;
    params.track.kalman.process_noise = user_kalman_process_noise;
    % Advanced HK specific params
    params.track.kalman.hk_alpha = user_kalman_hk_alpha;
    params.track.kalman.hk_beta = user_kalman_hk_beta;
    params.track.kalman.hk_forward_backward = user_kalman_hk_forward_backward;
    % NOTE: velocity_levels is generated dynamically in calculateDerivedParams().
    % Do NOT set it here — it would be overwritten.
    params.track.kalman.hk_v_max            = user_kalman_hk_v_max;
    params.track.kalman.hk_num_levels       = user_kalman_hk_num_levels;
    params.track.kalman.hk_spacing_power    = user_kalman_hk_spacing_power;
    params.track.kalman.hk_enable_overlap   = user_kalman_hk_enable_overlap;
    params.track.kalman.hk_overlap_mm_s     = user_kalman_hk_overlap_mm_s;
    
    % --- Detailed Weighted Cost Matrix Settings ---
    params.track.kalman.direction_penalty_weight = user_direction_weight;
    params.track.kalman.angle_penalty_slope = user_angle_penalty_slope;
    params.track.kalman.brightness_penalty_weight = user_brightness_weight;
    params.track.kalman.max_angle_change_deg = user_max_angle_gate_deg;
    params.track.kalman.gating_max_angle_change_deg = user_gating_max_angle_gate_deg;
    params.track.kalman.direction_history_points = user_track_direction_history_points;

    % --- Detailed Track QC (Post-Tracking) ---
    params.track.qc.enable_direction_constraint = user_track_qc_enable_direction;
    params.track.qc.max_angle_change_deg = user_track_qc_max_angle_deg;
    params.track.qc.enable_acceleration_constraint = user_track_qc_enable_acceleration;
    params.track.qc.acceleration_C_factor = user_track_qc_acceleration_factor;
    params.track.qc.enable_vd_constraint = user_track_qc_enable_vd;
    params.track.qc.max_vd_ratio = user_track_qc_max_vd_ratio;
    
    % --- Track Post Processing ---
    params.track.enable_postprocessing = user_enable_postprocessing; % Master switch
    params.track.smoothing_factor = user_smoothing_factor;           % window size
    params.track.smoothing_method = user_track_smoothing_method; 
    params.render.interpolation_method = user_track_interpolation_method;

    %% 5. Super-Resolution Image Reconstruction Parameters
    params.render.method = user_rendering_method;
    params.render.upsampling_factor = user_upsampling_factor;
    params.render.gaussian_sigma = 0.3;
    params.render.interpolation_step = user_interpolation_step; % Step for interpolation (0.2 means 5x more points)

    %% 6. Acquisition Parameters (Default Fallback)
    params.acq.framerate = 200; % Hz

    % =========================================================================
    % 7. ANALYSIS PARAMETERS (Feature Extraction)
    % =========================================================================
    % These parameters control how the "analyze_ULM_Features" script calculates stats.
    
    % --- Tortuosity (Vessel complexity) ---
    % Bins for the histogram: [min : step : max]
    % 1.0 is a straight line. Capillaries usually go up to 2-3.
    params.analysis.tortuosity_bins = 0 : 0.05 : 8.0; 
    
    % --- Velocity ---
    % Number of bins for velocity histograms (e.g., 50 or 100)
    params.analysis.velocity_hist_num_bins = 60; 
    
    % --- Density / Perfusion ---
    % Grid size for calculating vessel density (in mm). 
    % Smaller = higher resolution density map, but noisier.
    params.analysis.density_grid_size_mm = 0.5;

    % =========================================================================
    %      DYNAMIC & DERIVED PARAMETER LOADING
    % =========================================================================
    
    % Load experiment-specific parameters from 'info.txt'
    paramsFile = fullfile(params.io.data_folder, 'info.txt');
    if ~exist(paramsFile, 'file'), error('Parameter file (info.txt) not found at: %s', paramsFile); end
    expParams = getExpParams(paramsFile); % Calls the function from getExpParams.m
    params.expParams = expParams;
    
    % Add user-defined fields that are not parsed from info.txt
    params.expParams.channel_cross_section_mm2 = user_channel_cross_section_mm2; 

    % Determine data dimensions from a sample file
    % NOTE: For kidney data, we might point this to a different subfolder later.
    initial_data_path = fullfile(params.io.data_folder, params.expParams.bubbleType, params.io.data_subfolder);
    data_files = dir(fullfile(initial_data_path, params.io.file_pattern));
    if isempty(data_files), error('No data files found in: %s', initial_data_path); end
    dataStruct = load(fullfile(data_files(1).folder, data_files(1).name));
    dataFields = fieldnames(dataStruct);
    rawData = abs(dataStruct.(dataFields{1}));
    params.expParams.size = size(rawData);
    
    % Calculate all derived physical and acquisition parameters
    params = calculateDerivedParams(params, rawData, isKidneyExperiment);

    printExperimentSummary(params);

end % End of main setDefaultParams function


%% ========================================================================
%  Local Helper Functions
%  ========================================================================

function params = calculateDerivedParams(params, rawData, isKidneyExperiment)
% CALCULATEDERIVED PARAMS  Compute all physically-derived parameters from
%                          experiment metadata and raw data dimensions.
%
% INPUTS:
%   params             - Parameter struct (partially populated by main function)
%   rawData            - Sample raw IQ data array [Nz x Nx x Nt], used only
%                        for dimension inference
%   isKidneyExperiment - (logical) If true, uses in vivo velocity estimate
%                        instead of phantom flow-rate-based estimate
%
% OUTPUTS:
%   params             - Updated parameter struct with all derived fields set
%
% DERIVED FIELDS SET:
%   params.expParams.pixel_X_size     [mm]
%   params.expParams.pixel_Z_size     [mm]
%   params.expParams.C                [m/s] speed of sound
%   params.track.pixel_X_size         [mm]
%   params.track.pixel_Z_size         [mm]
%   params.track.dt                   [s]   inter-frame interval
%   params.track.max_linking_distance [px]  (if set to NaN by user)
%   params.acq.framerate              [Hz]
%   params.loc.fwhm                   [px]  (if set to NaN by user)
%   params.loc.qc_max_roi_maxima            (derived from fwhm)
%   params.loc.MB_image                     PSF template for NCC detection
%   params.track.kalman.velocity_levels     HKT velocity band boundaries
%
% NOTE ON AUTO-COMPUTE vs. MANUAL OVERRIDE:
%   Parameters marked as NaN in the QUICK EXPERIMENT SETTINGS block of
%   setDefaultParams.m are auto-computed here from physics. Parameters set
%   to explicit numeric values are passed through unchanged.

    % =====================================================================
    % SECTION 1: Spatial Calibration
    % Compute pixel sizes from FOV and image dimensions.
    % If one FOV axis is missing from info.txt, it is inferred by assuming
    % square pixels (isotropic sampling).
    % =====================================================================

    if isnan(params.expParams.fovX) && ~isnan(params.expParams.fovZ)
        % Infer lateral FOV from axial FOV assuming isotropic pixel size
        params.expParams.fovX = (params.expParams.fovZ / size(rawData, 1)) * size(rawData, 2);
        fprintf('   [Calibration] fovX inferred from fovZ: %.3f mm\n', params.expParams.fovX);

    elseif isnan(params.expParams.fovZ) && ~isnan(params.expParams.fovX)
        % Infer axial FOV from lateral FOV assuming isotropic pixel size
        params.expParams.fovZ = (params.expParams.fovX / size(rawData, 2)) * size(rawData, 1);
        fprintf('   [Calibration] fovZ inferred from fovX: %.3f mm\n', params.expParams.fovZ);

    elseif isnan(params.expParams.fovX) && isnan(params.expParams.fovZ)
        error('ULM:MissingFOV', ...
            'Both fovX and fovZ are NaN. Check that info.txt contains FOV fields.');
    end

    % Speed of sound (tissue standard)
    params.expParams.C = 1540; % m/s

    % Pixel sizes [mm/pixel]
    params.expParams.pixel_X_size = params.expParams.fovX / params.expParams.size(2);
    params.expParams.pixel_Z_size = params.expParams.fovZ / params.expParams.size(1);

    % Propagate pixel sizes to tracking and acquisition structs
    params.track.pixel_X_size = params.expParams.pixel_X_size;
    params.track.pixel_Z_size = params.expParams.pixel_Z_size;
    params.track.dt           = 1 / params.expParams.shootingRate;
    params.acq.framerate      = params.expParams.shootingRate;

    % Acoustic wavelength [mm], used for PSF-based parameter estimates
    if ~isnan(params.expParams.frequency)
        params.expParams.lambda = (params.expParams.C / 1000) / params.expParams.frequency; % mm
    else
        params.expParams.lambda = NaN;
        fprintf('   [Calibration] Warning: frequency not found in info.txt. Lambda-based auto-compute unavailable.\n');
    end

    fprintf('   [Calibration] Pixel size: X=%.4f mm, Z=%.4f mm | dt=%.4f s\n', ...
        params.expParams.pixel_X_size, params.expParams.pixel_Z_size, params.track.dt);

    % =====================================================================
    % SECTION 2: Maximum Linking Distance (Auto-compute or Manual)
    %
    % The linking distance determines the maximum allowed displacement
    % between a bubble's position in frame t and its candidate in frame t+1.
    %
    % Auto-compute logic:
    %   1. Estimate peak flow velocity for the experiment type.
    %   2. Compute maximum single-frame displacement = v_max * dt.
    %   3. Convert from mm to pixels.
    %
    % Manual override:
    %   Set user_max_linking_distance_px to an explicit value (not NaN)
    %   in setDefaultParams.m to skip this calculation entirely.
    % =====================================================================

    avg_pixel_size_mm = mean([params.track.pixel_X_size, params.track.pixel_Z_size]);

    % --- Always resolve channel cross-section ---
    if isnan(params.expParams.channel_cross_section_mm2)
        channel_height_mm = 0.3;
        channel_width_mm  = params.expParams.mainChannelDiameter / 1000;
        params.expParams.channel_cross_section_mm2 = channel_width_mm * channel_height_mm;
        fprintf('   [Linking] Cross-section auto-computed (rectangular): %.3f x %.3f mm = %.5f mm²\n', ...
            channel_width_mm, channel_height_mm, params.expParams.channel_cross_section_mm2);
    else
        fprintf('   [Linking] Cross-section from user input: %.5f mm²\n', ...
            params.expParams.channel_cross_section_mm2);
    end
    
    % --- Max linking distance: auto-compute or use manual value ---
    if isnan(params.track.max_linking_distance)
    
        if isKidneyExperiment
            v_max_mm_s = 10 * 3;
            fprintf('   [Linking] In vivo mode: v_max = %.1f mm/s\n', v_max_mm_s);
        else
            if ~isnan(params.expParams.flowSpeed)
                flow_rate_mm3_s   = params.expParams.flowSpeed * 1000 / 60;
                avg_velocity_mm_s = flow_rate_mm3_s / params.expParams.channel_cross_section_mm2;
                v_max_mm_s        = avg_velocity_mm_s * 2;
                fprintf('   [Linking] Q=%.3f ml/min, v_avg=%.1f mm/s, v_max=%.1f mm/s\n', ...
                    params.expParams.flowSpeed, avg_velocity_mm_s, v_max_mm_s);
            else
                v_max_mm_s = 20;
                fprintf('   [Linking] Warning: flow speed missing. Using fallback v_max=%.1f mm/s\n', v_max_mm_s);
            end
        end
    
        max_displacement_mm               = v_max_mm_s * params.track.dt;
        params.track.max_linking_distance = ceil(max_displacement_mm / avg_pixel_size_mm);
        fprintf('   [Linking] Auto-computed max linking distance: %d px (%.3f mm/frame)\n', ...
            params.track.max_linking_distance, max_displacement_mm);
    else
        fprintf('   [Linking] Using manual max linking distance: %d px\n', ...
            params.track.max_linking_distance);
    end

    % =====================================================================
    % SECTION 3: FWHM and Localization ROI Size (Auto-compute or Manual)
    %
    % The FWHM controls the ROI window used for sub-pixel localization.
    % It should approximate the bubble's PSF width in the filtered image.
    %
    % Auto-compute logic:
    %   fwhm = floor(lambda / pixel_X_size) * 2 + 1  (odd, in pixels)
    %   This gives a window of roughly one wavelength, which matches the
    %   expected PSF width for a diffraction-limited system.
    %
    % Manual override:
    %   Set user_fwhm_pixels to an explicit [x, z] pair (not NaN) in
    %   setDefaultParams.m. A fixed odd integer (e.g., [3, 3]) is often
    %   more stable than the wavelength-derived value and is recommended
    %   as the starting point.
    % =====================================================================

    if any(isnan(params.loc.fwhm))
        if ~isnan(params.expParams.lambda) && ~isnan(params.expParams.pixel_X_size)
            fwhm_px = floor(params.expParams.lambda / params.expParams.pixel_X_size) * 2 + 1;
            fwhm_px = max(fwhm_px, 3); % Enforce minimum of 3 px
            params.loc.fwhm = [1, 1] * fwhm_px;
            fprintf('   [FWHM] Auto-computed from wavelength: %d px (lambda=%.4f mm, pixel=%.4f mm)\n', ...
                fwhm_px, params.expParams.lambda, params.expParams.pixel_X_size);
        else
            params.loc.fwhm = [3, 3]; % Hard fallback
            fprintf('   [FWHM] Warning: lambda unavailable. Using fallback FWHM = [3, 3] px\n');
        end
    else
        fprintf('   [FWHM] Using manual FWHM: [%d, %d] px\n', params.loc.fwhm(1), params.loc.fwhm(2));
    end

    % ROI maxima QC threshold: derived from FWHM.
    % Allows up to floor(fwhm/2) local peaks in the ROI before rejecting the candidate.
    params.loc.qc_max_roi_maxima = max(1, floor(params.loc.fwhm(1) / 2));

    % =====================================================================
    % SECTION 4: Hierarchical Kalman Tracker (HKT) Velocity Levels
    %
    % Generates the velocity band boundaries for each tracking level.
    % Level boundaries are spaced according to spacing_power:
    %   power = 1.0 -> linear spacing (equal-width bands)
    %   power > 1.0 -> more levels at low velocities
    %   power < 1.0 -> more levels at high velocities
    % =====================================================================

    v_max_global  = params.track.kalman.hk_v_max;
    num_levels    = params.track.kalman.hk_num_levels;
    spacing_power = params.track.kalman.hk_spacing_power;
    enable_overlap = params.track.kalman.hk_enable_overlap;
    overlap_mm_s  = params.track.kalman.hk_overlap_mm_s;

    % Generate normalized boundaries in [0,1], raised to spacing_power, then scaled
    norm_boundaries = linspace(0, 1, num_levels + 1) .^ spacing_power;
    boundaries      = v_max_global * norm_boundaries;
    lower_bounds    = boundaries(1:end-1);
    upper_bounds    = boundaries(2:end);

    % Apply optional overlap between adjacent levels
    if enable_overlap && overlap_mm_s > 0
        half_overlap          = overlap_mm_s / 2;
        lower_bounds(2:end)   = lower_bounds(2:end)   - half_overlap;
        upper_bounds(1:end-1) = upper_bounds(1:end-1) + half_overlap;
        lower_bounds = max(lower_bounds, 0);          % Clamp to valid range
        upper_bounds = min(upper_bounds, v_max_global);
    end

    levels_matrix = round([lower_bounds', upper_bounds'], 1);
    params.track.kalman.velocity_levels = mat2cell(levels_matrix, ones(num_levels, 1), 2);

    overlap_str = 'off';
    if enable_overlap, overlap_str = 'on'; end
    fprintf('   [HKT] Velocity levels (power=%.1f, overlap=%s):\n', spacing_power, overlap_str);
    for lv = 1:num_levels
        fprintf('     Level %d: [%5.1f - %5.1f] mm/s\n', lv, levels_matrix(lv,1), levels_matrix(lv,2));
    end

    % =====================================================================
    % SECTION 5: PSF Template for NCC Detection
    %
    % Generates the reference kernel used by the Normalized Cross-Correlation
    % (NCC) detection method. Two modes are supported:
    %   'Gaussian'     - Analytical 2D Gaussian kernel from FWHM
    %   'Experimental' - Measured PSF loaded from a .mat file
    % =====================================================================

    if strcmp(params.loc.psf_type, 'Gaussian')

        % Convert FWHM to Gaussian sigma: FWHM = 2*sqrt(2*ln(2))*sigma ~ 2.355*sigma
        sigma_x = params.loc.fwhm(1) / 2.355;
        sigma_z = params.loc.fwhm(2) / 2.355;

        % Enforce odd kernel dimensions for symmetric NCC cropping
        sz_x = params.loc.psf_size(1);
        sz_z = params.loc.psf_size(2);
        if mod(sz_x, 2) == 0, sz_x = sz_x + 1; end
        if mod(sz_z, 2) == 0, sz_z = sz_z + 1; end

        % Build spatial grid centered at zero
        [X_grid, Z_grid] = meshgrid(-(sz_x-1)/2 : (sz_x-1)/2, ...
                                    -(sz_z-1)/2 : (sz_z-1)/2);

        % Analytical 2D Gaussian kernel (peak = 1, unnormalized)
        params.loc.MB_image = exp(-(X_grid.^2 / (2*sigma_x^2) + ...
                                    Z_grid.^2 / (2*sigma_z^2)));
        fprintf('   [PSF] Gaussian kernel generated: %dx%d px, sigma=[%.2f, %.2f] px\n', ...
            sz_x, sz_z, sigma_x, sigma_z);

    elseif strcmp(params.loc.psf_type, 'Experimental')

        if ~exist(params.loc.psf_file_path, 'file')
            error('ULM:PSFNotFound', ...
                'Experimental PSF file not found at: %s\nSet psf_type to ''Gaussian'' or provide the correct path.', ...
                params.loc.psf_file_path);
        end

        loaded_data = load(params.loc.psf_file_path);
        fields      = fieldnames(loaded_data);

        if isempty(fields)
            error('ULM:PSFEmpty', 'The Experimental PSF .mat file is empty: %s', params.loc.psf_file_path);
        end

        psf_matrix = loaded_data.(fields{1});
        [sz_z, sz_x] = size(psf_matrix);

        % Sanity check: PSF should be a small kernel, not a full image
        if sz_z < 3 || sz_x < 3 || sz_z > 51 || sz_x > 51
            warning('ULM:PSFSize', ...
                'Experimental PSF dimensions [%d x %d] are unusual for a single bubble PSF.', sz_z, sz_x);
        end

        % Enforce odd dimensions (required for symmetric normxcorr2 cropping)
        if mod(sz_z, 2) == 0
            psf_matrix = psf_matrix(1:end-1, :);
            fprintf('   [PSF] Cropped 1 row to enforce odd axial dimension.\n');
        end
        if mod(sz_x, 2) == 0
            psf_matrix = psf_matrix(:, 1:end-1);
            fprintf('   [PSF] Cropped 1 column to enforce odd lateral dimension.\n');
        end

        % Normalize to peak = 1
        params.loc.MB_image = psf_matrix / max(psf_matrix(:));
        fprintf('   [PSF] Experimental PSF loaded from: %s (%dx%d px)\n', ...
            params.loc.psf_file_path, size(params.loc.MB_image, 2), size(params.loc.MB_image, 1));
    end

end

function printExperimentSummary(params)
% PRINTEXPERIMENTSUMMARY  Print a structured summary of experiment-level
%                         parameters to the MATLAB command window.
%
%   Covers only physical acquisition and geometry parameters — not
%   algorithm settings (SVD cutoffs, tracking weights, etc.).
%   Called automatically at the end of setDefaultParams.

    SEP  = repmat('=', 1, 56);
    SEP2 = repmat('-', 1, 56);

    fprintf('\n%s\n', SEP);
    fprintf('  EXPERIMENT PARAMETER SUMMARY\n');
    fprintf('%s\n', SEP);

    % --- Data Source ---
    fprintf('\n  DATA SOURCE\n%s\n', SEP2);
    fprintf('  %-28s %s\n',   'Data folder:',    params.io.data_folder);
    fprintf('  %-28s %s\n',   'Subfolder:',       params.io.data_subfolder);
    fprintf('  %-28s %s\n',   'Bubble type:',     params.expParams.bubbleType);

    % --- Acquisition ---
    fprintf('\n  ACQUISITION\n%s\n', SEP2);
    fprintf('  %-28s %g Hz\n',    'Frame rate:',       params.acq.framerate);
    fprintf('  %-28s %g MHz\n',   'Frequency:',        params.expParams.frequency);
    fprintf('  %-28s %g m/s\n',   'Speed of sound:',   params.expParams.C);
    if ~isnan(params.expParams.lambda)
        fprintf('  %-28s %.4f mm  (%.1f um)\n', ...
            'Wavelength (lambda):', params.expParams.lambda, params.expParams.lambda * 1000);
    else
        fprintf('  %-28s N/A\n', 'Wavelength (lambda):');
    end

    % --- Field of View & Spatial Calibration ---
    fprintf('\n  SPATIAL CALIBRATION\n%s\n', SEP2);
    fprintf('  %-28s [%g x %g] pixels  (Z x X)\n', ...
        'Image size:',  params.expParams.size(1), params.expParams.size(2));
    fprintf('  %-28s %.4f mm\n',  'FOV axial  (Z):',   params.expParams.fovZ);
    fprintf('  %-28s %.4f mm\n',  'FOV lateral (X):',  params.expParams.fovX);
    fprintf('  %-28s %.4f mm/px\n', 'Pixel size Z:',   params.expParams.pixel_Z_size);
    fprintf('  %-28s %.4f mm/px\n', 'Pixel size X:',   params.expParams.pixel_X_size);
    fprintf('  %-28s %.4f s\n',   'Inter-frame dt:',   params.track.dt);

    % --- Flow & Geometry ---
    fprintf('\n  FLOW & CHANNEL GEOMETRY\n%s\n', SEP2);
    fprintf('  %-28s %g ml/min\n', 'Flow rate:',        params.expParams.flowSpeed);
    fprintf('  %-28s %g um\n',     'Main channel width:', params.expParams.mainChannelDiameter);
    if ~isnan(params.expParams.secondaryChannelDiameter)
        fprintf('  %-28s %g um\n', 'Secondary channel:', params.expParams.secondaryChannelDiameter);
    end
    fprintf('  %-28s %g deg\n',   'Inflow angle:',      params.expParams.angle);
    if isnan(params.expParams.channel_cross_section_mm2)
        fprintf('  %-28s N/A (not yet computed)\n', 'Channel cross-section:');
    else
        fprintf('  %-28s %.5f mm²\n', 'Channel cross-section:', ...
            params.expParams.channel_cross_section_mm2);
    end

    % --- Derived Tracking Parameters ---
    fprintf('\n  DERIVED TRACKING VALUES\n%s\n', SEP2);
    fprintf('  %-28s %d px\n',    'Max linking distance:', params.track.max_linking_distance);
    fprintf('  %-28s [%g, %g] px\n', 'FWHM [X, Z]:', params.loc.fwhm(1), params.loc.fwhm(2));
    fprintf('  %-28s %d\n',       'QC max ROI maxima:',  params.loc.qc_max_roi_maxima);

    fprintf('\n%s\n\n', SEP);
end