% =========================================================================
% CLASS: ULM_Processor < handle
% =========================================================================
%
% PURPOSE:
%   Central orchestrator for the Ultrasound Localization Microscopy (ULM)
%   pipeline. As a stateful handle class, it encapsulates the dataset,
%   processing parameters, and execution flow — from raw IQ data ingestion
%   to final super-resolved vascular maps.
%
% TYPICAL USAGE:
%   params    = setDefaultParams(false);         % Load configuration
%   processor = ULM_Processor(params);           % Initialize
%   processor.Step0_InitializeBoundaries();      % Interactive crop/mask
%   processor.run_localization_and_tracking_loop(); % Batch processing
%   processor.run_streaming_reconstruction(10);  % Render maps (minLen=10)
%   processor.saveResults();                     % Export everything
%
% UPDATING PARAMETERS AFTER CONSTRUCTION:
%   Since this is a handle class, params can be modified directly without
%   re-constructing the processor (which would create duplicate folders):
%       processor.params.track.method = 'kalman_v2';
%       processor.params.track.kalman.process_noise = 5;
%   For changes that affect the results directory path (upsampling_factor,
%   tracking method), call updateParams() to regenerate the path:
%       processor.updateParams(modified_params);
%
% =========================================================================
%                           PROPERTIES
% =========================================================================
%
%   params              (Struct) Master configuration struct from
%                        setDefaultParams(). Controls every algorithmic
%                        choice in the pipeline. Key sub-structs:
%                          .io       - I/O paths, file patterns, export flags
%                          .proc     - Crop, mask, and vessel masking config
%                          .filter   - SVD cutoff, spatial filter, Butterworth
%                          .loc      - Detection method, FWHM, PSF template
%                          .track    - Tracking method, linking distance,
%                                      gap closing, Kalman settings, QC
%                          .render   - Upsampling factor, render method, sigma
%                          .acq      - Frame rate
%                          .expParams - Pixel sizes, FOV, frame count, etc.
%
%   data_files          (Struct Array) Output of dir() listing all input
%                        .mat data buffers found in data_dir.
%
%   data_dir            (String) Full path to the active input data folder.
%                        Built as: data_folder / bubbleType / [subfolder].
%
%   results_dir         (String) Timestamped output directory for this run.
%                        Structure: Results/res_<upsample>/<method>/<timestamp>/
%
%   general_results_dir (String) Parent directory for all results runs.
%
%   tracks_dir          (String) Temporary folder where per-buffer track
%                        .mat files are saved during batch processing and
%                        streamed from during reconstruction.
%
%   bubbleTracks        (Struct Array) Final aggregated tracks from all
%                        buffers. Populated by saveResults().
%
%   densityMap          (Matrix) Final super-resolved density map. Each
%                        pixel counts how many track points passed through
%                        it. Populated by run_streaming_reconstruction().
%
%   velocityMap         (Matrix) Final super-resolved mean velocity map.
%                        Each pixel is the weighted average of velocities
%                        from all tracks that passed through it. Populated
%                        by run_streaming_reconstruction().
%
%   metadata            (Struct) Runtime metadata with sub-fields:
%                          .cropRect         - [x, y, w, h] crop rectangle
%                          .roiMask          - Logical mask for ROI
%                          .processing_times - Timing for each pipeline stage
%
% =========================================================================
%                        PUBLIC METHODS
% =========================================================================
%
%   obj = ULM_Processor(params)
%       Constructor. Validates the data folder, builds the directory
%       structure (input path, tracks folder, timestamped results folder),
%       locates all .mat data files, and cleans any previous tracks.
%       INPUT:  params - (Struct) from setDefaultParams().
%       OUTPUT: obj    - Initialized ULM_Processor handle.
%
%   updateParams(newParams)
%       Replaces obj.params and regenerates the results directory path
%       with a fresh timestamp. Use this instead of re-constructing the
%       processor when changing parameters that affect the results path
%       (e.g., track.method, render.upsampling_factor). Does NOT re-scan
%       data files or re-create the tracks folder.
%       INPUT: newParams - (Struct) Updated parameter struct.
%
%   set_data_dir(new_dir)
%       Updates obj.data_dir to a new folder path and automatically
%       reloads the file listing via reload_data_files(). Use when
%       switching to a different data folder without re-constructing.
%       INPUT: new_dir - (String) Path to an existing data folder.
%
%   set_data_files(new_files)
%       Directly replaces obj.data_files with a custom struct array.
%       Expected format is the output of dir() (fields: name, folder, etc.).
%       Use for manual control over which files to process.
%       INPUT: new_files - (Struct Array) File listing to use.
%
%   reload_data_files(file_pattern)
%       Re-scans obj.data_dir for files matching the given pattern and
%       updates obj.data_files. Called automatically by set_data_dir().
%       If no pattern is provided, uses params.io.file_pattern.
%       INPUT: file_pattern - (String, optional) Glob pattern, e.g. '*.mat'.
%
%   Step0_InitializeBoundaries()
%       Interactive or file-based setup of spatial constraints before
%       batch processing. Handles two independent features:
%       - Crop: Draws/loads a rectangle to spatially crop every buffer.
%       - Mask: Draws/loads a freehand ROI mask for detection filtering.
%       Both use a fast SVD-filtered preview from the first data file.
%       Stores results in obj.metadata.cropRect and obj.metadata.roiMask.
%
%   [filteredData, currentMask] = Step1_Filter(rawData, currentMask)
%       Full clutter filtering and spatial conditioning pipeline for one
%       buffer. Steps: adaptive crop -> SVD/DCC filter -> envelope
%       detection -> spatial filter -> adaptive mask application.
%       INPUT:  rawData     - [H x W x T] raw IQ data cube.
%               currentMask - Logical mask (may be resized internally).
%       OUTPUT: filteredData - [H x W x T] filtered data.
%               currentMask  - Possibly resized mask.
%
%   run_localization_and_tracking_loop()
%       Main heavy-lifting loop. Iterates through all data buffers:
%       For each buffer: loads raw data -> Step1_Filter -> detect &
%       localize -> track -> save tracks to disk. Reports per-step
%       timing. Appends localizations to a cumulative file.
%       Populates obj.metadata.processing_times.localization_and_tracking.
%
%   run_streaming_reconstruction(minLength)
%       Memory-efficient map rendering. Streams track files from disk
%       one at a time. Each file's tracks are either post-processed
%       (smoothed, interpolated, QC'd) or used raw, then split into:
%       - Density tracks: ALL tracks, using gap-interpolated paths from
%         Kalman v2 when available.
%       - Velocity tracks: Only tracks with has_reliable_velocity=true,
%         preventing zero-velocity short tracks from corrupting the map.
%       INPUT:  minLength - (scalar, default 1) Minimum track length.
%       OUTPUT: Populates obj.densityMap and obj.velocityMap.
%
%   processed_tracks = postProcessTracks(raw_tracks, min_track_length)
%       Track quality control, smoothing, and interpolation. Stages:
%       1. Length filter (reject shorter than min_track_length).
%       2. Tortuosity QC (reject if path_length/displacement > threshold).
%       3. Smoothing (Savitzky-Golay, rloess, movmean, etc.).
%       4. Spline interpolation (upsample path for sub-pixel rendering).
%       5. Velocity recalculation on the smoothed path.
%       6. Propagate has_reliable_velocity from the original track.
%       Runs in parallel (parfor) for efficiency.
%       INPUT:  raw_tracks       - (Struct Array) tracks from tracker.
%               min_track_length - (scalar) minimum length to keep.
%       OUTPUT: processed_tracks - (Struct Array) filtered & smoothed.
%
%   saveResults()
%       Final aggregation and export. Steps:
%       1. Saves params.mat to results_dir.
%       2. Loads all track files, sorts numerically, applies frame offsets.
%       3. Builds a TrackTable with: TrackID, Length, MeanVelocity_mm_s,
%          HasReliableVelocity, StartFrame, EndFrame, Duration_s.
%       4. Runs analyze_ULM_Features for statistical analysis.
%       5. Generates and saves velocity, length, and tortuosity histograms.
%       6. Exports density_map.tif and velocity_map.tif.
%
% =========================================================================
%                    HELPER METHODS (public access)
% =========================================================================
%
%   filteredData = filter_data(rawData, indent_prefix)
%       Dispatches to the configured clutter filter: svd_filter, svd_ssm,
%       dcc_svd, or svd_blockwise. Optionally applies Butterworth bandpass.
%
%   candidateBubbles = detect_candidates(filteredData)
%       Dispatches to the configured bubble detector: Intensity (classic
%       local-maxima), NP (Neyman-Pearson), or NCC (normalized cross-
%       correlation with PSF template). Passes ROI mask to NP/NCC so
%       noise statistics are computed on unmasked data.
%
%   localizations = localize_bubbles(filteredData, indent_prefix)
%       Runs detection then sub-pixel localization. Dispatches to radial
%       symmetry, gaussian_fit, or gaussian_fit_fast.
%
%   tracks = track_bubbles(localizations, indent_prefix)
%       Dispatches to the configured tracker: hungarian, nn, kalman,
%       kalman_v2, or kalman_advanced. Returns a struct array of tracks.
%
%   [density, velocity, weights] = render_tracks(tracks, sRes_dims)
%       Dispatches to histogram or Gaussian renderer. Returns three
%       accumulation maps: density count, velocity sum, velocity weight.
%
%   [density_tracks, vel_tracks] = splitTracksForRendering(all_tracks)
%       Separates tracks for density vs velocity map rendering:
%       - Density tracks: all tracks; swaps .path with .path_interpolated
%         (Kalman v2 gap predictions) when available.
%       - Velocity tracks: only tracks with has_reliable_velocity=true.
%
%   cropBox = getROICrop(dataFrame)
%       Opens an interactive figure for rectangle crop selection.
%       Saves the crop box as a timestamped .mat file.
%
%   mask = getROIMask(dataFrame)
%       Opens an interactive figure for freehand ROI mask drawing.
%       Saves the mask as a timestamped .mat file.
%
%   validLocs = validate_localizations(localizations)
%       Removes localizations outside the combined vessel mask and ROI
%       mask. Used as an optional safety check.
%
%   par_save(fname, tracks)
%       Saves a tracks variable to a .mat file with retry logic (3
%       attempts, 1-second delays) for network drive resilience.
%
% =========================================================================
%                    PRIVATE METHODS
% =========================================================================
%
%   previewImg = generateSVDPreview(~, cropRect)
%       Loads the first data file, takes up to 200 frames, optionally
%       applies a crop rectangle, then runs a fast SVD filter to produce
%       a clean mean image for interactive crop/mask drawing.
%
%   generateAnalysisFigures(TrackTable, TrackPaths)
%       Runs analyze_ULM_Features on the aggregated track data and
%       generates velocity, track-length, and tortuosity histograms.
%       Saves figures as .png and .fig to results_dir.
%
% DEPENDENCIES:
%   - setDefaultParams.m (parameter configuration)
%   - SVD_filter / SVD_SSM / DCC_SVD / SVD_blockwise (clutter filters)
%   - detectBubbles / detectBubbles_NP / detectBubbles_NCC (detectors)
%   - localizeRadialSymmetry / fit2DGaussian / fit2DGaussian_Fast
%   - trackHungarian / trackNearestNeighbor / trackKalman / trackKalman_v2
%     / trackKalman_Advanced
%   - renderHistogram / renderGaussian
%   - analyze_ULM_Features (optional, for statistics)
%   - Butterworth_bandpass_filter (optional)
%
% AUTHOR: Grigori Shapiro
% =========================================================================

classdef ULM_Processor < handle

    properties
        params              % Master configuration struct
        data_files          % Directory listing of input data buffers
        data_dir            % Path to the active input data folder
        results_dir         % Timestamped output directory for this run
        general_results_dir % Base directory for all results
        tracks_dir          % Temporary folder for per-buffer track files
        bubbleTracks        % Aggregated tracks from all buffers
        densityMap          % Final super-resolved density map
        velocityMap         % Final super-resolved velocity map
        metadata            % Crop coordinates, ROI masks, processing times
    end

    methods

        %% ================================================================
        %  CONSTRUCTOR
        %  ================================================================
        function obj = ULM_Processor(params)
        % ULM_PROCESSOR  Initialize the processor, validate paths, and
        %   set up the directory structure for input data and outputs.

            obj.params = params;
            obj.metadata.processing_times = struct();

            % Validate the base data folder
            if isempty(obj.params.io.data_folder) || ~isfolder(obj.params.io.data_folder)
                error('ULM:InvalidFolder', ...
                    'Data folder does not exist: "%s"', obj.params.io.data_folder);
            end

            % Build the input data path: base / bubbleType / [subfolder]
            base_folder = obj.params.io.data_folder;
            input_data_folder = fullfile(base_folder, params.expParams.bubbleType);

            if isfield(obj.params.io, 'data_subfolder') && ~isempty(obj.params.io.data_subfolder)
                input_data_folder = fullfile(input_data_folder, obj.params.io.data_subfolder);
                fprintf('   - Data subfolder: "%s"\n', obj.params.io.data_subfolder);
            end
            obj.data_dir = input_data_folder;

            % Build output paths
            obj.tracks_dir = fullfile(base_folder, 'Tracks');
            obj.general_results_dir = fullfile(base_folder, 'Results');
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
            obj.results_dir = fullfile(base_folder, 'Results', ...
                ['res ' num2str(obj.params.render.upsampling_factor)], ...
                obj.params.track.method, timestamp);

            fprintf('   - Directory structure:\n');
            fprintf('       Input:   %s\n', input_data_folder);
            fprintf('       Tracks:  %s\n', obj.tracks_dir);
            fprintf('       Results: %s\n', obj.results_dir);

            % Locate data files
            obj.data_files = dir(fullfile(input_data_folder, obj.params.io.file_pattern));
            if isempty(obj.data_files)
                error('ULM:NoDataFiles', ...
                    'No files matching "%s" in "%s".', ...
                    obj.params.io.file_pattern, input_data_folder);
            end

            % Prepare the tracks directory (clean if exists, create if not)
            if exist(obj.tracks_dir, 'dir')
                delete(fullfile(obj.tracks_dir, '*.*'));
            else
                mkdir(obj.tracks_dir);
            end
            if ~exist(obj.results_dir, 'dir'), mkdir(obj.results_dir); end

            fprintf('   - Found %d data files. Processor ready.\n', length(obj.data_files));
        end

        function obj = set_data_dir(obj, new_dir)
            % Update the path to the active input data folder
            if ~ischar(new_dir) && ~isstring(new_dir)
                error('ULM:InvalidInput', 'Directory path must be a string or character array.');
            end
            if ~isfolder(new_dir)
                error('ULM:InvalidFolder', 'Data directory does not exist: "%s"', new_dir);
            end
            
            obj.data_dir = new_dir;
            obj.reload_data_files;
            fprintf('   - Updated data_dir to: "%s"\n', obj.data_dir);
        end
    
        function obj = set_data_files(obj, new_files)
            % Update the directory listing of input data buffers
            % Expected input is a struct array like the one returned by the 'dir' function
            if ~isstruct(new_files) && ~isempty(new_files)
                error('ULM:InvalidInput', 'data_files must be a struct array (output of dir()).');
            end
            
            obj.data_files = new_files;
            fprintf('   - Updated data_files. Total files: %d\n', length(obj.data_files));
        end
        
        function obj = reload_data_files(obj, file_pattern)
            % Helper function: Reloads data_files based on the current data_dir
            if nargin < 2
                file_pattern = obj.params.io.file_pattern;
            end
            
            obj.data_files = dir(fullfile(obj.data_dir, file_pattern));
            
            if isempty(obj.data_files)
                warning('ULM:NoDataFiles', 'No files matching "%s" in "%s".', ...
                    file_pattern, obj.data_dir);
            else
                fprintf('   - Reloaded data_files. Found %d files.\n', length(obj.data_files));
            end
        end

        %% ================================================================
        %  UPDATE PARAMS (avoids re-constructing & double directories)
        %  ================================================================
        function updateParams(obj, newParams)
        % UPDATEPARAMS  Replace the parameter struct and regenerate paths.
        %   Use this instead of re-constructing the processor when you need
        %   to change parameters that affect the results directory path
        %   (e.g., tracking method, upsampling factor).
        %
        %   Usage:
        %       params.track.method = 'kalman_v2';
        %       processor.updateParams(params);

            obj.params = newParams;

            % Regenerate the results directory with a new timestamp
            base_folder = obj.params.io.data_folder;
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
            obj.results_dir = fullfile(base_folder, 'Results', ...
                ['res ' num2str(obj.params.render.upsampling_factor)], ...
                obj.params.track.method, timestamp);

            if ~exist(obj.results_dir, 'dir'), mkdir(obj.results_dir); end
            fprintf('   - Parameters updated. New results dir: %s\n', obj.results_dir);
        end

        %% ================================================================
        %  STEP 0: INITIALIZE SPATIAL BOUNDARIES (CROP & MASK)
        %  ================================================================
        function Step0_InitializeBoundaries(obj)
            fprintf('\n   - [STEP 0] Initializing Spatial Boundaries...\n');

            % --- CROP ---
            if obj.params.proc.enableInteractiveCrop
                if obj.params.proc.generateNewCrop
                    fprintf('     -> Generating SVD preview for interactive crop...\n');
                    previewImg = obj.generateSVDPreview([], []);
                    obj.metadata.cropRect = obj.getROICrop(previewImg);
                else
                    fprintf('     -> Loading existing crop box...\n');
                    [f, p] = uigetfile('*.mat', 'Select Crop Box File', obj.params.proc.cropPath);
                    if f ~= 0
                        load(fullfile(p, f), 'cropBox');
                        obj.metadata.cropRect = cropBox;
                        fprintf('     -> Crop box: [%s]\n', num2str(cropBox));
                    else
                        warning('No crop file selected. Cropping disabled.');
                        obj.params.proc.enableInteractiveCrop = false;
                    end
                end
            end

            % --- MASK ---
            if obj.params.proc.enableInteractiveMask
                if obj.params.proc.generateNewMask
                    fprintf('     -> Generating SVD preview for interactive mask...\n');
                    cropRect = [];
                    if obj.params.proc.enableInteractiveCrop && isfield(obj.metadata, 'cropRect')
                        cropRect = obj.metadata.cropRect;
                    end
                    previewImg = obj.generateSVDPreview([], cropRect);
                    obj.metadata.roiMask = obj.getROIMask(previewImg);
                else
                    fprintf('     -> Loading existing mask...\n');
                    [f, p] = uigetfile('*.mat', 'Select ROI Mask File', obj.params.proc.maskPath);
                    if f ~= 0
                        loaded_mask = load(fullfile(p, f));
                        maskFields = fieldnames(loaded_mask);
                        obj.metadata.roiMask = logical(loaded_mask.(maskFields{1}));
                        fprintf('     -> Mask loaded successfully.\n');
                    else
                        warning('No mask file selected. Masking disabled.');
                        obj.params.proc.enableInteractiveMask = false;
                    end
                end
            end
        end

        %% ================================================================
        %  STEP 1: CLUTTER FILTERING & SPATIAL CONDITIONING
        %  ================================================================
        function [filteredData, currentMask] = Step1_Filter(obj, rawData, currentMask)
            [H, W, ~] = size(rawData);

            % Adaptive crop (clamp to actual data dimensions)
            if obj.params.proc.enableInteractiveCrop && isfield(obj.metadata, 'cropRect')
                rect = obj.metadata.cropRect;
                x_start = max(1, rect(1));
                y_start = max(1, rect(2));
                x_end   = min(W, x_start + rect(3) - 1);
                y_end   = min(H, y_start + rect(4) - 1);
                rawData = rawData(y_start:y_end, x_start:x_end, :);
                [H, W, ~] = size(rawData);
            end

            % SVD / DCC clutter filter
            filteredData = obj.filter_data(rawData, '       ');

            % Envelope detection
            filteredData = abs(filteredData);

            % Spatial conditioning
            if ~strcmp(obj.params.filter.spatial_method, 'None')
                filteredData = obj.apply_spatial_filter(filteredData, obj.params.filter);
            end

            % Adaptive masking: resize mask to match (possibly cropped) data.
            % For Intensity detection, apply mask directly to pixel values.
            % For NP/NCC, the mask is passed to the detector separately so
            % noise statistics are computed on the unmasked data.
            if obj.params.proc.enableInteractiveMask && ~isempty(currentMask)
                [mH, mW] = size(currentMask);
                if mH ~= H || mW ~= W
                    fprintf('      -> Resizing mask from %dx%d to %dx%d.\n', mH, mW, H, W);
                    currentMask = imresize(currentMask, [H, W], 'nearest');
                end
                if strcmp(obj.params.loc.DetectMethod, 'Intensity')
                    filteredData = filteredData .* currentMask;
                end
            end
        end

        function outData = apply_spatial_filter(~, inData, filterParams)
        % APPLY_SPATIAL_FILTER  Frame-by-frame spatial conditioning.

            [~, ~, T] = size(inData);
            outData = zeros(size(inData), 'like', inData);

            method = filterParams.spatial_method;
            kSize  = max(1, round(filterParams.spatial_kernel));
            sig1   = filterParams.spatial_sigma1;
            sig2   = filterParams.spatial_sigma2;

            % Enforce odd kernel for Gaussian/Median
            if ismember(method, {'Gaussian', 'Median'}) && mod(kSize, 2) == 0
                kSize = kSize + 1;
            end
            if strcmp(method, 'Top-Hat'), se = strel('disk', kSize); end

            for t = 1:T
                frame = inData(:,:,t);
                switch method
                    case 'Gaussian', outData(:,:,t) = imgaussfilt(frame, sig1, 'FilterSize', kSize);
                    case 'Median',   outData(:,:,t) = medfilt2(frame, [kSize kSize]);
                    case 'DoG',      outData(:,:,t) = max(0, imgaussfilt(frame, sig1) - imgaussfilt(frame, sig2));
                    case 'Top-Hat',  outData(:,:,t) = imtophat(frame, se);
                end
            end
        end

        %% ================================================================
        %  MAIN BATCH PROCESSING LOOP
        %  ================================================================
        function run_localization_and_tracking_loop(obj)
            num_buffers = length(obj.data_files);
            stage1_timer = tic;

            % Retrieve the mask if it was set in Step 0
            currentMask = [];
            if isfield(obj.metadata, 'roiMask')
                currentMask = obj.metadata.roiMask;
            end

            for i = 1:num_buffers
                fprintf('\n   === Buffer %d/%d: %s ===\n', i, num_buffers, obj.data_files(i).name);
                buffer_timer = tic;

                % Load raw data
                dataStruct = load(fullfile(obj.data_files(i).folder, obj.data_files(i).name));
                fNames = fieldnames(dataStruct);
                rawData = dataStruct.(fNames{1});

                % Step 1: Filter
                fprintf('     [1/3] Filtering...\n');
                filter_t = tic;
                [filteredData, currentMask] = obj.Step1_Filter(rawData, currentMask);
                if ~isempty(currentMask), obj.metadata.roiMask = currentMask; end
                fprintf('           Done in %.1f s.\n', toc(filter_t));

                % Step 2: Detect & Localize
                fprintf('     [2/3] Localizing...\n');
                loc_t = tic;
                localizations = obj.localize_bubbles(filteredData, '       ');
                fprintf('           %d localizations in %.1f s.\n', height(localizations), toc(loc_t));

                % Offset frame numbers for global time continuity
                localizations.Frame = localizations.Frame + (i-1) * obj.params.expParams.size(3);

                % Append localizations to cumulative file
                resultsPath = fullfile(obj.results_dir, 'LocalizationsTable.mat');
                if exist(resultsPath, 'file')
                    data = load(resultsPath);
                    bubbleLocalizations = [data.bubbleLocalizations; localizations];
                else
                    bubbleLocalizations = localizations;
                end
                save(resultsPath, 'bubbleLocalizations');

                % Step 3: Track
                fprintf('     [3/3] Tracking...\n');
                track_t = tic;
                buffer_tracks = obj.track_bubbles(localizations, '       ');
                fprintf('           %d tracks in %.1f s.\n', length(buffer_tracks), toc(track_t));

                % Save tracks to disk
                track_filename = fullfile(obj.tracks_dir, sprintf('tracks_%03d.mat', i));
                obj.par_save(track_filename, buffer_tracks);

                fprintf('   === Buffer %d/%d complete: %.1f s ===\n', i, num_buffers, toc(buffer_timer));
            end

            elapsed = toc(stage1_timer);
            obj.metadata.processing_times.localization_and_tracking = elapsed;
            fprintf('\n   Total batch processing: %.1f s (%.1f min).\n', elapsed, elapsed/60);
        end

        %% ================================================================
        %  STREAMING RECONSTRUCTION (DENSITY & VELOCITY MAPS)
        %  ================================================================
        function run_streaming_reconstruction(obj, minLength)
        % RUN_STREAMING_RECONSTRUCTION  Builds super-resolution maps by
        %   streaming track files from disk one at a time.
        %
        %   Handles two rendering paths:
        %   - Post-processing ON:  Tracks are smoothed, interpolated, and
        %     QC-filtered. Velocity reliability is propagated from the
        %     tracker and respected during rendering.
        %   - Post-processing OFF: Raw tracks are used directly. Gap-
        %     interpolated paths (from Kalman v2) are used for density
        %     rendering. Tracks flagged as velocity-unreliable contribute
        %     to density but not to the velocity map.

            if nargin < 2, minLength = 1; end

            fprintf('\n   - Streaming reconstruction (minLength = %d)...\n', minLength);
            reconstruction_timer = tic;

            track_files = dir(fullfile(obj.tracks_dir, 'tracks_*.mat'));
            if isempty(track_files)
                error('ULM:NoTrackFiles', ...
                    'No track files in %s. Run batch processing first.', obj.tracks_dir);
            end

            % Determine canvas dimensions (account for crop)
            H = obj.params.expParams.size(1);
            W = obj.params.expParams.size(2);
            if obj.params.proc.enableInteractiveCrop && isfield(obj.metadata, 'cropRect')
                rect = obj.metadata.cropRect;
                H = rect(4);
                W = rect(3);
            end
            sRes_dims = [H, W] * obj.params.render.upsampling_factor;

            % Accumulation canvases
            final_density_map         = zeros(sRes_dims);
            final_velocity_sum_map    = zeros(sRes_dims);
            final_velocity_weight_map = zeros(sRes_dims);

            h_wait = waitbar(0, 'Streaming Reconstruction...', 'Name', 'ULM Post-Processing');

            for i = 1:length(track_files)
                file_timer = tic;
                waitbar(i / length(track_files), h_wait, ...
                    sprintf('Processing file %d / %d', i, length(track_files)));

                loaded_data = load(fullfile(track_files(i).folder, track_files(i).name));

                % --- Prepare tracks for rendering ---
                if obj.params.track.enable_postprocessing
                    % Post-processing: smooth, interpolate, QC filter
                    all_tracks = obj.postProcessTracks(loaded_data.tracks, minLength);
                    fprintf('       [PostProcess] File %d/%d: %d -> %d tracks (%.1f s)\n', ...
                        i, length(track_files), length(loaded_data.tracks), ...
                        length(all_tracks), toc(file_timer));
                else
                    % No post-processing: use raw tracks with length filter
                    all_tracks = loaded_data.tracks;
                    if ~isempty(all_tracks)
                        keep = [all_tracks.length] >= minLength;
                        all_tracks = all_tracks(keep);
                    end
                    fprintf('       [Raw] File %d/%d: %d tracks after length filter (%.1f s)\n', ...
                        i, length(track_files), length(all_tracks), toc(file_timer));
                end

                if isempty(all_tracks), continue; end

                % --- Split into density tracks and velocity tracks ---
                % Density: uses gap-interpolated paths when available (fills
                %   holes from Kalman gap predictions). ALL tracks contribute.
                % Velocity: only tracks with reliable velocity contribute,
                %   preventing zero-velocity short tracks from dragging down
                %   the weighted mean.
                [density_tracks, velocity_tracks] = obj.splitTracksForRendering(all_tracks);

                % Render density (from all tracks, with interpolated paths)
                [density_chunk, ~, ~] = obj.render_tracks(density_tracks, sRes_dims);
                final_density_map = final_density_map + density_chunk;

                % Render velocity (only from reliable tracks)
                if ~isempty(velocity_tracks)
                    [~, vel_sum_chunk, vel_weight_chunk] = obj.render_tracks(velocity_tracks, sRes_dims);
                    final_velocity_sum_map    = final_velocity_sum_map    + vel_sum_chunk;
                    final_velocity_weight_map = final_velocity_weight_map + vel_weight_chunk;
                end
            end
            close(h_wait);

            % Finalize maps
            obj.densityMap  = final_density_map;
            obj.velocityMap = zeros(sRes_dims);
            valid_px = final_velocity_weight_map > 1e-6;
            obj.velocityMap(valid_px) = final_velocity_sum_map(valid_px) ./ final_velocity_weight_map(valid_px);

            elapsed = toc(reconstruction_timer);
            obj.metadata.processing_times.reconstruction = elapsed;
            fprintf('   - Reconstruction complete: %.1f s.\n', elapsed);
        end

        %% ================================================================
        %  POST-PROCESS TRACKS (SMOOTH, INTERPOLATE, QC)
        %  ================================================================
        function processed_tracks = postProcessTracks(obj, raw_tracks, min_track_length)
        % POSTPROCESSTRACKS  Filters, smooths, and interpolates raw tracks.
        %
        %   Workflow:
        %     1. Length filter: reject tracks shorter than min_track_length.
        %     2. QC jitter filter: reject tracks with excessive tortuosity
        %        (velocity dispersion ratio > max_vd_ratio).
        %     3. Smoothing: Savitzky-Golay, rloess, or other method.
        %     4. Interpolation: upsample path via spline/pchip/linear.
        %     5. Velocity recalculation on the smoothed path.
        %     6. Propagate velocity reliability flag from the tracker.

            if isempty(raw_tracks)
                processed_tracks = [];
                return;
            end

            % Cache parameters for parfor
            smoothing_factor = obj.params.track.smoothing_factor;
            smoothing_method = obj.params.track.smoothing_method;
            interp_step      = obj.params.render.interpolation_step;
            interp_method    = obj.params.render.interpolation_method;
            enable_vd_qc     = obj.params.track.qc.enable_vd_constraint;
            max_vd_ratio     = obj.params.track.qc.max_vd_ratio;
            pixel_X          = obj.params.track.pixel_X_size;
            pixel_Z          = obj.params.track.pixel_Z_size;
            dt               = obj.params.track.dt;
            sgolay_poly_order = 3;

            processed_tracks_cell = cell(1, length(raw_tracks));

            parfor i = 1:length(raw_tracks)
                track = raw_tracks(i);

                % --- Stage 1: Length filter ---
                if track.length < min_track_length
                    processed_tracks_cell{i} = [];
                    continue;
                end

                % --- Stage 2: QC jitter filter (tortuosity) ---
                if enable_vd_qc
                    path_qc = track.path;
                    net_disp = norm(path_qc(end,:) - path_qc(1,:));
                    if net_disp < 1e-6
                        processed_tracks_cell{i} = [];
                        continue;
                    end
                    total_len = sum(sqrt(sum(diff(path_qc).^2, 2)));
                    if (total_len / net_disp) > max_vd_ratio
                        processed_tracks_cell{i} = [];
                        continue;
                    end
                end

                % --- Stage 3: Smoothing ---
                path_to_process = track.path;
                win = min(smoothing_factor, track.length);

                switch lower(smoothing_method)
                    case 'sgolay'
                        if mod(win, 2) == 0, win = win - 1; end
                        poly_ord = min(sgolay_poly_order, win - 1);
                        if track.length > win && poly_ord > 0
                            path_to_process = [ ...
                                sgolayfilt(double(track.path(:,1)), poly_ord, win), ...
                                sgolayfilt(double(track.path(:,2)), poly_ord, win)];
                        end
                    case {'rloess', 'loess', 'gaussian', 'movmean', 'lowess'}
                        if win >= 3
                            path_to_process = [ ...
                                smoothdata(track.path(:,1), smoothing_method, win), ...
                                smoothdata(track.path(:,2), smoothing_method, win)];
                        end
                end

                % --- Stage 4: Interpolation (upsampling) ---
                orig_idx   = 1:track.length;
                interp_idx = 1:interp_step:track.length;

                x_interp = interp1(orig_idx, path_to_process(:,1), interp_idx, interp_method);
                y_interp = interp1(orig_idx, path_to_process(:,2), interp_idx, interp_method);
                f_interp = interp1(orig_idx, track.frames,         interp_idx, 'linear');

                x_interp = fillmissing(x_interp, 'linear', 'EndValues', 'nearest');
                y_interp = fillmissing(y_interp, 'linear', 'EndValues', 'nearest');

                new_path   = [x_interp', y_interp'];
                new_length = size(new_path, 1);

                % --- Stage 5: Velocity recalculation ---
                if new_length > 1
                    dx = diff(new_path(:,1));
                    dy = diff(new_path(:,2));
                    dt_f = diff(f_interp');

                    disp_mm = sqrt((dx * pixel_X).^2 + (dy * pixel_Z).^2);
                    dt_sec  = dt_f * dt;

                    vel = zeros(size(dt_sec));
                    ok  = dt_sec > 1e-9;
                    vel(ok) = disp_mm(ok) ./ dt_sec(ok);
                    vel = [vel; vel(end)]; %#ok<AGROW>

                    mean_v = mean(vel, 'omitnan');
                    if isnan(mean_v), mean_v = 0; end
                    avg_vel = repmat(mean_v, new_length, 1);
                else
                    vel     = 0;
                    avg_vel = 0;
                end

                % --- Stage 6: Propagate velocity reliability ---
                % Smoothing/interpolating a handful of noisy points does
                % not make velocity trustworthy. Preserve the tracker's flag.
                if isfield(track, 'has_reliable_velocity')
                    reliable = track.has_reliable_velocity;
                else
                    reliable = true;  % v1 tracks have no flag
                end

                % Mean intensity (for export metadata)
                val_intensity = 0;
                if isfield(track, 'localizations') && istable(track.localizations) ...
                        && ismember('Intensity', track.localizations.Properties.VariableNames)
                    val_intensity = mean(track.localizations.Intensity);
                end

                % Package output
                processed_tracks_cell{i} = struct( ...
                    'id',                     track.id, ...
                    'path',                   new_path, ...
                    'frames',                 f_interp', ...
                    'length',                 new_length, ...
                    'velocities_mm_s',        vel, ...
                    'average_velocity_mm_s',  avg_vel, ...
                    'has_reliable_velocity',  reliable, ...
                    'mean_intensity',         val_intensity, ...
                    'localizations',          []);
            end

            processed_tracks = [processed_tracks_cell{:}];
        end

        %% ================================================================
        %  SAVE RESULTS (AGGREGATION, STATISTICS, EXPORT)
        %  ================================================================
        function saveResults(obj)
            fprintf('\n   - Saving results to: %s\n', obj.results_dir);
            save_timer = tic;

            % Save parameters
            params = obj.params; %#ok<PROPLC>
            save(fullfile(obj.results_dir, 'params.mat'), 'params');

            % Locate and sort track files numerically
            trackFiles = dir(fullfile(obj.tracks_dir, '*tracks*.mat'));
            if isempty(trackFiles)
                warning('No track files found. Skipping aggregation.');
                return;
            end

            fileNumbers = zeros(length(trackFiles), 1);
            for i = 1:length(trackFiles)
                numStr = regexp(trackFiles(i).name, '\d+', 'match');
                if ~isempty(numStr)
                    fileNumbers(i) = str2double(numStr{end});
                end
            end
            [~, sortedIdx] = sort(fileNumbers);
            trackFiles = trackFiles(sortedIdx);

            fprintf('     Aggregating %d track files...\n', length(trackFiles));

            % Batch size for frame offset calculation
            if isfield(obj.params.expParams, 'size') && length(obj.params.expParams.size) >= 3
                frames_per_batch = obj.params.expParams.size(3);
            else
                warning('Batch size not in params.expParams.size(3). Frame offsets will be zero.');
                frames_per_batch = 0;
            end

            % Initialize storage cells
            all_IDs        = {};
            all_Lengths    = {};
            all_MeanVel    = {};
            all_Reliable   = {};
            all_StartFrame = {};
            all_EndFrame   = {};
            all_Duration   = {};
            all_Paths      = {};
            total_tracks   = 0;

            for k = 1:length(trackFiles)
                filePath = fullfile(trackFiles(k).folder, trackFiles(k).name);
                current_frame_offset = (k - 1) * frames_per_batch;

                try
                    loadedData = load(filePath);

                    % Robust field detection
                    if isfield(loadedData, 'tracks')
                        bt = loadedData.tracks;
                    elseif isfield(loadedData, 'finalTracks')
                        bt = loadedData.finalTracks;
                    else
                        vars = fieldnames(loadedData);
                        bt = loadedData.(vars{1});
                    end

                    if isempty(bt) || ~isstruct(bt), continue; end

                    n = length(bt);
                    total_tracks = total_tracks + n;

                    % Scalar fields
                    all_IDs{end+1}     = [bt.id]';                                  %#ok<AGROW>
                    all_Lengths{end+1} = [bt.length]';                              %#ok<AGROW>

                    % Paths
                    all_Paths{end+1} = {bt.path}';                                  %#ok<AGROW>

                    % Mean velocity (first element of the vector)
                    vc = {bt.average_velocity_mm_s};
                    all_MeanVel{end+1} = cellfun(@(v) v(1), vc, 'UniformOutput', true)'; %#ok<AGROW>

                    % Velocity reliability flag (backward-compatible)
                    if isfield(bt, 'has_reliable_velocity')
                        all_Reliable{end+1} = [bt.has_reliable_velocity]';          %#ok<AGROW>
                    else
                        all_Reliable{end+1} = true(n, 1);                           %#ok<AGROW>
                    end

                    % Frame timing with global offset
                    fc = {bt.frames};
                    raw_start = cellfun(@(f) double(f(1)),   fc, 'UniformOutput', true)';
                    raw_end   = cellfun(@(f) double(f(end)), fc, 'UniformOutput', true)';
                    all_StartFrame{end+1} = raw_start + current_frame_offset;       %#ok<AGROW>
                    all_EndFrame{end+1}   = raw_end   + current_frame_offset;       %#ok<AGROW>
                    all_Duration{end+1}   = (raw_end - raw_start) * obj.params.track.dt; %#ok<AGROW>

                catch ME
                    warning('Error reading %s: %s', trackFiles(k).name, ME.message);
                end
            end

            % Build and save the master table
            if total_tracks > 0
                TrackTable = table( ...
                    vertcat(all_IDs{:}), ...
                    vertcat(all_Lengths{:}), ...
                    vertcat(all_MeanVel{:}), ...
                    vertcat(all_Reliable{:}), ...
                    vertcat(all_StartFrame{:}), ...
                    vertcat(all_EndFrame{:}), ...
                    vertcat(all_Duration{:}), ...
                    'VariableNames', {'TrackID', 'Length', 'MeanVelocity_mm_s', ...
                                      'HasReliableVelocity', ...
                                      'StartFrame', 'EndFrame', 'Duration_s'});

                TrackPaths = vertcat(all_Paths{:});

                if isfield(obj.params.io, 'save_lightweight') && obj.params.io.save_lightweight
                    save(fullfile(obj.results_dir, 'TrackTable.mat'), ...
                         'TrackTable', 'TrackPaths', 'params');
                    fprintf('     Saved lightweight data (%d tracks).\n', total_tracks);
                end

                % Feature analysis and histograms
                obj.generateAnalysisFigures(TrackTable, TrackPaths);
            else
                fprintf('     No tracks found to aggregate.\n');
            end

            % TIFF exports
            if obj.params.io.export_tiff
                if ~isempty(obj.densityMap)
                    normDensity = obj.densityMap / (max(obj.densityMap(:)) + eps);
                    imwrite(normDensity, fullfile(obj.results_dir, 'density_map.tif'));
                    fprintf('     Exported density_map.tif\n');
                end
                if ~isempty(obj.velocityMap)
                    imwrite(mat2gray(obj.velocityMap), fullfile(obj.results_dir, 'velocity_map.tif'));
                    fprintf('     Exported velocity_map.tif\n');
                end
            end

            fprintf('   - Save complete (%.1f s).\n', toc(save_timer));
        end

    end % public methods

    %% ====================================================================
    %  HELPER METHODS
    %  ====================================================================
    methods (Access = public)

        function filteredData = filter_data(obj, rawData, indent_prefix)
        % FILTER_DATA  Dispatches to the configured clutter filter.

            if nargin < 3, indent_prefix = '       '; end

            switch obj.params.filter.method
                case 'svd_filter'
                    filteredData = SVD_filter(rawData, obj.params.filter.svd_cutoff, ...
                                             'IndentPrefix', indent_prefix);
                case 'svd_ssm'
                    filteredData = SVD_SSM(rawData, 'IndentPrefix', indent_prefix);
                case 'dcc_svd'
                    filteredData = DCC_SVD(rawData, obj.params.acq.framerate, ...
                                          'ReconstructionMode', 'blood', ...
                                          'IndentPrefix', indent_prefix);
                case 'svd_blockwise'
                    bw = obj.params.filter.blockwise;
                    nv_args = { ...
                        'ThresholdMethod',       bw.threshold_method, ...
                        'OverlapPct',            bw.overlap_pct, ...
                        'ManualCutoff',          bw.manual_cutoff, ...
                        'TissueFreqThreshHz',    bw.tissue_freq_hz, ...
                        'MPDeviationSigma',      bw.mp_deviation_sigma, ...
                        'GradientInflectionPct', bw.gradient_pct, ...
                        'MinBloodComponents',    bw.min_blood_comps, ...
                        'MaxTissueFraction',     bw.max_tissue_frac, ...
                        'PlotThresholdMaps',     bw.plot_maps, ...
                        'IndentPrefix',          indent_prefix};
                    if ~isempty(bw.block_size_mm) && ~isnan(bw.block_size_mm)
                        nv_args = [nv_args, {'BlockSizeMm', [bw.block_size_mm, bw.block_size_mm]}];
                    end
                    [filteredData, bw_diag] = SVD_blockwise(rawData, obj.params, nv_args{:});
                    obj.params.filter.blockwise.last_diagnostics = bw_diag;
                otherwise
                    error('Unknown filter method: %s', obj.params.filter.method);
            end

            % Optional Butterworth bandpass
            if isfield(obj.params.filter, 'enable_butterworth') && obj.params.filter.enable_butterworth
                fprintf('%s -> Butterworth bandpass...\n', indent_prefix);
                filteredData = Butterworth_bandpass_filter(filteredData, ...
                    obj.params.filter.butter_cutoff, ...
                    obj.params.acq.framerate, ...
                    obj.params.filter.butter_order);
            end
        end

        function candidateBubbles = detect_candidates(obj, filteredData)
        % DETECT_CANDIDATES  Dispatches to Intensity, NP, or NCC detector.

            roiMask = [];
            if obj.params.proc.enableInteractiveMask && isfield(obj.metadata, 'roiMask')
                roiMask = obj.metadata.roiMask;
            end

            switch upper(obj.params.loc.DetectMethod)
                case 'INTENSITY'
                    candidateBubbles = detectBubbles(filteredData, obj.params.loc);
                case 'NP'
                    candidateBubbles = detectBubbles_NP(filteredData, obj.params.loc, roiMask);
                case 'NCC'
                    if ~isfield(obj.params.loc, 'MB_image') || isempty(obj.params.loc.MB_image)
                        error('ULM:NoPSF', 'NCC requires params.loc.MB_image.');
                    end
                    candidateBubbles = detectBubbles_NCC(filteredData, obj.params.loc, roiMask);
                otherwise
                    error('ULM:BadDetector', 'Unknown DetectMethod: "%s".', obj.params.loc.DetectMethod);
            end
        end

        function localizations = localize_bubbles(obj, filteredData, indent_prefix)
        % LOCALIZE_BUBBLES  Detection + sub-pixel localization.

            if nargin < 3, indent_prefix = '       '; end
            candidateBubbles = obj.detect_candidates(filteredData);
            fprintf('%s-> %s: %d candidates.\n', indent_prefix, obj.params.loc.DetectMethod, height(candidateBubbles));

            switch obj.params.loc.method
                case 'radial'
                    localizations = localizeRadialSymmetry(filteredData, candidateBubbles, obj.params.loc, indent_prefix);
                case 'gaussian_fit'
                    localizations = fit2DGaussian(filteredData, candidateBubbles, obj.params.loc, indent_prefix);
                case 'gaussian_fit_fast'
                    localizations = fit2DGaussian_Fast(filteredData, candidateBubbles, obj.params.loc, indent_prefix);
                otherwise
                    error('Unknown localization method: %s', obj.params.loc.method);
            end
        end

        function tracks = track_bubbles(obj, localizations, indent_prefix)
        % TRACK_BUBBLES  Dispatches to the configured tracking algorithm.

            if nargin < 3, indent_prefix = '       '; end
            switch lower(obj.params.track.method)
                case 'hungarian'
                    tracks = trackHungarian(localizations, obj.params, indent_prefix);
                case 'nn'
                    tracks = trackNearestNeighbor(localizations, obj.params, indent_prefix);
                case 'kalman'
                    tracks = trackKalman(localizations, obj.params, indent_prefix);
                case 'kalman_advanced'
                    tracks = trackKalman_Advanced(localizations, obj.params, indent_prefix);
                otherwise
                    error('Unknown tracking method: %s', obj.params.track.method);
            end
        end

        function [density, velocity, weights] = render_tracks(obj, tracks, sRes_dims)
        % RENDER_TRACKS  Dispatches to histogram or Gaussian renderer.

            switch obj.params.render.method
                case 'histogram'
                    [density, velocity, weights] = renderHistogram(tracks, sRes_dims, obj.params.render);
                case 'gaussian'
                    [density, velocity, weights] = renderGaussian(tracks, sRes_dims, obj.params.render);
                otherwise
                    error('Unknown render method: %s', obj.params.render.method);
            end
        end

        function [density_tracks, velocity_tracks] = splitTracksForRendering(~, all_tracks)
        % SPLITTRACKSFORRENDERING  Separates tracks for density vs velocity.
        %
        %   Density tracks: ALL tracks. If path_interpolated is available
        %     (from Kalman v2 gap prediction), it replaces .path so that
        %     gap frames contribute to the density map.
        %   Velocity tracks: Only tracks with has_reliable_velocity = true.
        %     This prevents short tracks (with zeroed velocity) from
        %     dragging down the weighted mean velocity in the map.

            density_tracks = all_tracks;

            % Swap in interpolated paths for density rendering
            if isfield(all_tracks, 'path_interpolated')
                for k = 1:length(density_tracks)
                    if ~isempty(density_tracks(k).path_interpolated)
                        n_interp = size(density_tracks(k).path_interpolated, 1);
                        density_tracks(k).path = density_tracks(k).path_interpolated;
                        % Extend the velocity vector to match the longer path
                        density_tracks(k).average_velocity_mm_s = ...
                            repmat(density_tracks(k).average_velocity_mm_s(1), n_interp, 1);
                    end
                end
            end

            % Filter velocity tracks by reliability flag
            if isfield(all_tracks, 'has_reliable_velocity')
                rel_mask = [all_tracks.has_reliable_velocity];
                velocity_tracks = all_tracks(rel_mask);
            else
                velocity_tracks = all_tracks;
            end
        end

        function cropBox = getROICrop(obj, dataFrame)
        % GETROICROP  Interactive rectangle selection for spatial cropping.

            fig = figure; imshow(dataFrame, []); title('Draw Crop Region');
            roi = drawrectangle;
            cropBox = round(roi.Position);
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            filename = sprintf('CropBox_%s.mat', timestamp);
            save(fullfile(obj.params.proc.cropPath, filename), 'cropBox');
            fprintf('     -> Crop saved: %s\n', filename);
            close(fig);
        end

        function mask = getROIMask(obj, dataFrame)
        % GETROIMASK  Interactive freehand selection for ROI masking.

            fig = figure; imshow(dataFrame, []); title('Draw ROI');
            roi_2 = images.roi.AssistedFreehand; draw(roi_2);
            mask = createMask(roi_2, size(dataFrame,1), size(dataFrame,2));
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            filename = sprintf('Mask_%s.mat', timestamp);
            save(fullfile(obj.params.proc.maskPath, filename), 'mask');
            fprintf('     -> Mask saved: %s\n', filename);
            close(fig);
        end

        function validLocs = validate_localizations(obj, localizations)
        % VALIDATE_LOCALIZATIONS  Removes points outside active masks.

            if isempty(localizations)
                validLocs = localizations;
                return;
            end

            H = obj.params.expParams.size(1);
            W = obj.params.expParams.size(2);
            combinedMask = true(H, W);
            maskActive = false;

            if isfield(obj.params.proc, 'vesselMask') && ...
               isfield(obj.params.proc.vesselMask, 'enable') && obj.params.proc.vesselMask.enable && ...
               isfield(obj.params.proc.vesselMask, 'vesselMask') && ~isempty(obj.params.proc.vesselMask.vesselMask)
                combinedMask = combinedMask & logical(obj.params.proc.vesselMask.vesselMask);
                maskActive = true;
            end

            if isfield(obj.params.proc, 'enableInteractiveMask') && obj.params.proc.enableInteractiveMask && ...
               isfield(obj.metadata, 'roiMask') && ~isempty(obj.metadata.roiMask)
                combinedMask = combinedMask & logical(obj.metadata.roiMask);
                maskActive = true;
            end

            if ~maskActive
                validLocs = localizations;
                return;
            end

            x_idx = round(localizations.X);
            y_idx = round(localizations.Y);
            in_bounds = (x_idx >= 1 & x_idx <= W & y_idx >= 1 & y_idx <= H);
            keep_mask = false(size(in_bounds));
            if any(in_bounds)
                indices = sub2ind([H, W], y_idx(in_bounds), x_idx(in_bounds));
                keep_mask(in_bounds) = combinedMask(indices);
            end

            validLocs = localizations(keep_mask, :);
            n_removed = height(localizations) - height(validLocs);
            if n_removed > 0
                fprintf('      -> Mask check: removed %d out-of-mask localizations.\n', n_removed);
            end
        end

        function par_save(~, fname, tracks)
        % PAR_SAVE  Save with retry for network drive resilience.

            for attempt = 1:3
                try
                    save(fname, 'tracks');
                    return;
                catch ME
                    if attempt == 3, rethrow(ME); end
                    warning('Save failed (%d/3), retrying: %s', attempt, fname);
                    pause(1);
                end
            end
        end

    end % helper methods

    %% ====================================================================
    %  PRIVATE METHODS
    %  ====================================================================
    methods (Access = private)

        function previewImg = generateSVDPreview(obj, ~, cropRect)
        % GENERATESVDPREVIEW  Produces a quick SVD-filtered mean image
        %   from the first data file, used for interactive crop/mask.

            dataStruct  = load(fullfile(obj.data_files(1).folder, obj.data_files(1).name));
            fNames      = fieldnames(dataStruct);
            previewData = dataStruct.(fNames{1});
            previewData = previewData(:,:,1:min(200, size(previewData,3)));

            % Apply crop if specified
            if nargin >= 3 && ~isempty(cropRect)
                rect = cropRect;
                previewData = previewData(rect(2):rect(2)+rect(4)-1, ...
                                          rect(1):rect(1)+rect(3)-1, :);
            end

            % Fast SVD filter for visualization
            [U, S, V] = svd(reshape(previewData, [], size(previewData,3)), 'econ');
            cutoff = obj.params.filter.svd_cutoff;
            U(:, 1:max(1, cutoff(1)-1)) = 0;
            U(:, min(size(U,2), cutoff(2)+1):end) = 0;
            filtered = reshape(U * S * V', size(previewData));
            previewImg = mean(abs(filtered), 3);
        end

        function generateAnalysisFigures(obj, TrackTable, TrackPaths)
        % GENERATEANALYSISFIGURES  Runs feature analysis and exports plots.

            try
                minimalTracks = struct('path', TrackPaths, ...
                    'velocities_mm_s', num2cell(TrackTable.MeanVelocity_mm_s), ...
                    'length', num2cell(TrackTable.Length));

                analysisStats = analyze_ULM_Features(minimalTracks, obj.params, obj.results_dir);

                if isfield(obj.params.io, 'export_figures') && obj.params.io.export_figures && ~isempty(analysisStats)
                    fprintf('     Generating histograms...\n');

                    % Velocity histogram
                    f = figure('Visible', 'off');
                    histogram(TrackTable.MeanVelocity_mm_s, 100, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
                    xlabel('Velocity [mm/s]'); ylabel('Count');
                    title(sprintf('Velocity Distribution (Mean: %.2f mm/s)', analysisStats.Velocity.Mean));
                    grid on;
                    saveas(f, fullfile(obj.results_dir, 'Hist_Velocity.png'));
                    saveas(f, fullfile(obj.results_dir, 'Hist_Velocity.fig'));
                    close(f);

                    % Track length histogram
                    f = figure('Visible', 'off');
                    histogram(TrackTable.Length, 'BinMethod', 'integers', 'FaceColor', [0.8 0.4 0.2]);
                    xlabel('Track Length [points]'); ylabel('Count');
                    title(sprintf('Track Length Distribution (Mean: %.1f)', mean(TrackTable.Length)));
                    grid on;
                    saveas(f, fullfile(obj.results_dir, 'Hist_Length.png'));
                    saveas(f, fullfile(obj.results_dir, 'Hist_Length.fig'));
                    close(f);

                    % Tortuosity histogram (if available)
                    if isfield(analysisStats, 'Tortuosity')
                        f = figure('Visible', 'off');
                        bar(analysisStats.Tortuosity.HistEdges(1:end-1), analysisStats.Tortuosity.Histogram, ...
                            'FaceColor', [0.4 0.8 0.4]);
                        xlabel('Tortuosity Index'); ylabel('Count');
                        title(sprintf('Tortuosity (Median: %.3f)', analysisStats.Tortuosity.Median));
                        grid on;
                        saveas(f, fullfile(obj.results_dir, 'Hist_Tortuosity.png'));
                        saveas(f, fullfile(obj.results_dir, 'Hist_Tortuosity.fig'));
                        close(f);
                    end

                    fprintf('     Histograms saved.\n');
                end
            catch ME
                warning(ME.identifier, 'Feature analysis failed: %s', ME.message);
            end
        end

    end % private methods

end % classdef