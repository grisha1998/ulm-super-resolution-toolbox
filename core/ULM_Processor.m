% =========================================================================
% CLASS: ULM_Processor < handle
% =========================================================================
%
% PURPOSE:
%   This class is the central orchestrator (the "Heart") of the entire 
%   Ultrasound Localization Microscopy (ULM) pipeline. As a stateful `handle` 
%   class, it encapsulates the complete dataset, the processing parameters, 
%   and the execution flow—from raw ultrasound data ingestion to the final 
%   super-resolved vascular maps. It seamlessly dispatches data to specific 
%   algorithmic functions based on the user's configuration.
%
% CORE PIPELINE WORKFLOW:
%   1. Initialization (Constructor & Step 0): Sets up directory structures, 
%      locates data batches, and allows the user to interactively define 
%      spatial constraints (Crop Box & ROI Mask) using a fast SVD preview.
%   2. Batch Processing (run_localization_and_tracking_loop):
%      Iterates through large datasets buffer-by-buffer to minimize memory overhead:
%      - Filter (Step 1): Applies spatiotemporal clutter filtering (SVD/DCC), 
%        envelope detection, spatial smoothing, and adaptive masking.
%      - Detect & Localize: Dispatches to Intensity, NP, or NCC detection, 
%        followed by sub-pixel localization (Radial Symmetry or 2D Gaussian).
%      - Track: Links localizations into trajectories (Hungarian, NN, or Kalman).
%      - Save Intermediate: Flushes raw tracks for the current buffer to disk.
%   3. Post-Processing & Reconstruction (run_streaming_reconstruction):
%      Streams intermediate track files from disk. Tracks undergo rigorous 
%      Quality Control (e.g., Velocity Dispersion/Tortuosity filtering), 
%      Savitzky-Golay/Loess smoothing, and Spline interpolation. The refined 
%      tracks are then rendered into continuous Density and Velocity maps.
%   4. Finalization (saveResults): Aggregates all batch tracks, calculates 
%      final kinematic statistics (mean velocities, lengths, tortuosity), 
%      generates distribution histograms, and exports the high-resolution TIFFs.
%
% PROPERTIES:
%   --- Configuration ---
%   params              : (Struct) Master configuration object dictating all algorithmic choices.
%   --- Path Management ---
%   data_files          : (Struct Array) Directory listing of all input raw data buffers.
%   data_dir            : (String) Path to the active input data.
%   results_dir         : (String) Timestamped directory for current run's outputs.
%   general_results_dir : (String) Base directory for all results.
%   tracks_dir          : (String) Temporary folder for batch-processed tracks.
%   --- State & Results ---
%   bubbleTracks        : (Struct Array) The final, aggregated list of all valid tracks.
%   densityMap          : (Matrix) The final super-resolved microbubble hit-count image.
%   velocityMap         : (Matrix) The final super-resolved mean velocity image.
%   --- Metadata ---
%   metadata            : (Struct) Stores interactive crop coordinates, ROI masks, and processing times.
%
% KEY PUBLIC METHODS:
%   ULM_Processor(params)                : Constructor. Validates paths and preps the filesystem.
%   Step0_InitializeBoundaries()         : Generates fast SVD previews to define Crop and Mask.
%   Step1_Filter(...)                    : Applies adaptive cropping, clutter filtering, and masking.
%   run_localization_and_tracking_loop() : The main heavy-lifting loop processing all data batches.
%   run_streaming_reconstruction(...)    : Memory-efficient map rendering via advanced track post-processing.
%   postProcessTracks(...)               : Critical track filtering, smoothing, and spline interpolation.
%   saveResults()                        : Master aggregation, statistical plotting, and TIFF export.
%
% INTERNAL HELPER METHODS:
%   filter_data, detect_candidates, localize_bubbles, track_bubbles, 
%   render_tracks, getROICrop, getROIMask, validate_localizations, par_save.
%
% AUTHOR: Grigori Shapiro
% =========================================================================

classdef ULM_Processor < handle
    
    properties
        params              % Configuration parameters
        
        % Paths
        data_files          % List of input data files
        data_dir            % Main data directory
        results_dir         % Main output directory
        general_results_dir % General output directory
        tracks_dir          % Directory for intermediate track files
        
        % Final Results
        bubbleTracks        % Aggregated tracks from all buffers
        densityMap          % Final super-resolved density map
        velocityMap         % Final super-resolved velocity map
        
        % Metadata
        metadata
    end
    
    methods
        function obj = ULM_Processor(params)
            % --- Constructor ---
            % Initializes the processor, sets up paths, and finds data files.
            obj.params = params;
            
            % Validate data path
            if isempty(obj.params.io.data_folder) || ~isfolder(obj.params.io.data_folder)
                error('Base data folder is not specified or does not exist. Please set params.io.data_folder.');
            end
            
            % --- Setup directory structure based on user request ---
            base_folder = obj.params.io.data_folder;
            
            % Start with the default input folder (based on bubble type)
            input_data_folder = fullfile(base_folder, params.expParams.bubbleType);
            
            % Check for the optional data_subfolder and append it to the path if it exists and is not empty.
            if isfield(obj.params.io, 'data_subfolder') && ~isempty(obj.params.io.data_subfolder)
                input_data_folder = fullfile(input_data_folder, obj.params.io.data_subfolder);
                obj.data_dir = input_data_folder;
                fprintf('   - Optional data subfolder specified: "%s"\n', obj.params.io.data_subfolder);
            end

            obj.tracks_dir = fullfile(base_folder, 'Tracks');
            timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss')); % Create a timestamp string, e.g., "2025-07-29_10-40-17"
            obj.results_dir = fullfile(base_folder, 'Results', 'res ' + string(obj.params.render.upsampling_factor), obj.params.track.method, timestamp); % Add the timestamp as the final subfolder in the results path
            obj.general_results_dir = fullfile(base_folder, 'Results');

            fprintf('   - Setting up directory structure:\n');
            fprintf('       - Input Data: %s\n', input_data_folder);
            fprintf('       - Tracks Out: %s\n', obj.tracks_dir);
            fprintf('       - Results Out: %s\n', obj.results_dir);

            % Find all data files in the specified input folder
            obj.data_files = dir(fullfile(input_data_folder, obj.params.io.file_pattern));
            if isempty(obj.data_files)
                error('No data files found in "%s" matching the pattern "%s".', input_data_folder, obj.params.io.file_pattern);
            end
            
            % Create output directories if they don't exist
            if exist(obj.tracks_dir, 'dir')
                % If it exists, delete all files inside it.
                % Note: This will not delete subdirectories.
                delete(fullfile(obj.tracks_dir, '*.*'));
            else
                % If it does not exist, create it.
                mkdir(obj.tracks_dir);
            end
            if ~exist(obj.results_dir, 'dir'), mkdir(obj.results_dir); end
            
            obj.metadata.processing_times = struct();
            fprintf('   - Found %d data files to process.\n', length(obj.data_files));
            fprintf('   - Processor initialized successfully.\n');

        end
        
        function Step0_InitializeBoundaries(obj)
            % STEP 0: Interactively define or load Crop and Mask before the main loop
            fprintf('\n   - [STEP 0] Initializing Spatial Boundaries (Crop & Mask)...\n');
            
            % --- 1. HANDLE CROP ---
            if obj.params.proc.enableInteractiveCrop
                if obj.params.proc.generateNewCrop
                    fprintf('     -> Generating fast preview for Spatial Crop...\n');
                    % Load small subset for fast preview
                    dataStruct = load(fullfile(obj.data_files(1).folder, obj.data_files(1).name));
                    fNames = fieldnames(dataStruct); previewData = dataStruct.(fNames{1});
                    previewData = previewData(:,:,1:min(200, size(previewData,3)));
                    
                    % Fast SVD for clear visualization
                    [U,S,V] = svd(reshape(previewData, [], size(previewData,3)), 'econ');
                    cutoff = obj.params.filter.svd_cutoff;
                    U(:,1:max(1, cutoff(1)-1)) = 0; U(:,min(size(U,2), cutoff(2)+1):end) = 0;
                    filtered = reshape(U*S*V', size(previewData));
                    previewImg = mean(abs(filtered), 3);
                    
                    obj.metadata.cropRect = obj.getROICrop(previewImg);
                else
                    fprintf('     -> Select existing Crop Box file...\n');
                    [f, p] = uigetfile('*.mat', 'Select Crop Box File', obj.params.proc.cropPath);
                    if f ~= 0
                        load(fullfile(p, f), 'crop_box'); 
                        obj.metadata.cropRect = crop_box;
                        fprintf('     -> Loaded Crop Box: [%s]\n', num2str(crop_box));
                    else
                        warning('No crop file selected. Cropping disabled.');
                        obj.params.proc.enableInteractiveCrop = false;
                    end
                end
            end
            
            % --- 2. HANDLE MASK ---
            if obj.params.proc.enableInteractiveMask
                if obj.params.proc.generateNewMask
                    fprintf('     -> Generating fast preview for ROI Mask...\n');
                    % If crop was defined, we must apply it to the preview so the mask aligns
                    dataStruct = load(fullfile(obj.data_files(1).folder, obj.data_files(1).name));
                    fNames = fieldnames(dataStruct); previewData = dataStruct.(fNames{1});
                    previewData = previewData(:,:,1:min(200, size(previewData,3)));
                    
                    if obj.params.proc.enableInteractiveCrop && isfield(obj.metadata, 'cropRect')
                        rect = obj.metadata.cropRect;
                        previewData = previewData(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1, :);
                    end
                    
                    % Fast SVD
                    [U,S,V] = svd(reshape(previewData, [], size(previewData,3)), 'econ');
                    cutoff = obj.params.filter.svd_cutoff;
                    U(:,1:max(1, cutoff(1)-1)) = 0; U(:,min(size(U,2), cutoff(2)+1):end) = 0;
                    filtered = reshape(U*S*V', size(previewData));
                    previewImg = mean(abs(filtered), 3);
                    
                    obj.metadata.roiMask = obj.getROIMask(previewImg);
                else
                    fprintf('     -> Select existing Mask file...\n');
                    [f, p] = uigetfile('*.mat', 'Select ROI Mask File', obj.params.proc.maskPath);
                    if f ~= 0
                        loaded_mask = load(fullfile(p, f));
                        maskFields = fieldnames(loaded_mask);
                        obj.metadata.roiMask = logical(loaded_mask.(maskFields{1}));
                        fprintf('     -> Loaded Mask successfully.\n');
                    else
                        warning('No mask file selected. Masking disabled.');
                        obj.params.proc.enableInteractiveMask = false;
                    end
                end
            end
        end

        function [filteredData, currentMask] = Step1_Filter(obj, rawData, currentMask)
            % STEP 1: Full Clutter and Spatial Filtering Pipeline with Adaptive Dimensions
            [H, W, T] = size(rawData);
            
            % --- 1. Adaptive Crop ---
            if obj.params.proc.enableInteractiveCrop && isfield(obj.metadata, 'cropRect')
                rect = obj.metadata.cropRect;
                % Adaptive Clamping to handle registered data dimension changes
                x_start = max(1, rect(1));
                y_start = max(1, rect(2));
                x_end = min(W, x_start + rect(3) - 1);
                y_end = min(H, y_start + rect(4) - 1);
                
                rawData = rawData(y_start:y_end, x_start:x_end, :);
                [H, W, T] = size(rawData); % Update dimensions
            end
            
            % --- 2. SVD / DCC Clutter Filter ---
            filteredData = obj.filter_data(rawData, '       ');
            
            % --- 3. Envelope Detection ---
            filteredData = abs(filteredData);
            
            % --- 4. Spatial Conditioning ---
            if ~strcmp(obj.params.filter.spatial_method, 'None')
                filteredData = obj.apply_spatial_filter(filteredData, obj.params.filter);
            end
            
            % --- 5. Adaptive Masking ---
            % CHANGED: Apply mask here ONLY for 'Intensity' detection method.
            % For 'NP' and 'NCC', the mask is passed directly to the detect function
            % and applied internally AFTER statistics computation (critical for NP
            % correctness — masking before median/MAD biases the noise estimate).
            if obj.params.proc.enableInteractiveMask && ~isempty(currentMask)
                [mH, mW] = size(currentMask);
                if mH ~= H || mW ~= W
                    fprintf('      -> [Adaptive] Resizing mask from %dx%d to match data %dx%d.\n', mH, mW, H, W);
                    currentMask = imresize(currentMask, [H, W], 'nearest');
                end
                % Only apply mask to filteredData here for Intensity-based detection.
                % NP and NCC will receive the mask as a separate argument.
                if strcmp(obj.params.loc.DetectMethod, 'Intensity')
                    filteredData = filteredData .* currentMask;
                end
            end
        end

        function outData = apply_spatial_filter(obj, inData, filterParams)
            % Applies the requested spatial filter
            [H, W, T] = size(inData);
            outData = zeros(size(inData), 'like', inData);
            
            method = filterParams.spatial_method;
            kSize = max(1, round(filterParams.spatial_kernel));
            sig1 = filterParams.spatial_sigma1;
            sig2 = filterParams.spatial_sigma2;
            
            if (strcmp(method, 'Gaussian') || strcmp(method, 'Median')) && mod(kSize, 2) == 0
                kSize = kSize + 1; 
            end
            if strcmp(method, 'Top-Hat'), se = strel('disk', kSize); end
            
            for t = 1:T
                frame = inData(:,:,t);
                switch method
                    case 'Gaussian'
                        outData(:,:,t) = imgaussfilt(frame, sig1, 'FilterSize', kSize);
                    case 'Median'
                        outData(:,:,t) = medfilt2(frame, [kSize kSize]);
                    case 'DoG'
                        g1 = imgaussfilt(frame, sig1); g2 = imgaussfilt(frame, sig2);
                        outData(:,:,t) = max(0, g1 - g2);
                    case 'Top-Hat'
                        outData(:,:,t) = imtophat(frame, se);
                end
            end
        end

        function run_localization_and_tracking_loop(obj)
            % --- REFACTORED MAIN LOOP ---
            num_buffers = length(obj.data_files);
            stage1_timer = tic;
            
            % Grab the initial mask if it was loaded/created in Step 0
            currentMask = [];
            if isfield(obj.metadata, 'roiMask')
                currentMask = obj.metadata.roiMask;
            end
            
            for i = 1:num_buffers
                fprintf('\n   - Processing buffer %d/%d: %s\n', i, num_buffers, obj.data_files(i).name);
                buffer_timer = tic;
                
                dataStruct = load(fullfile(obj.data_files(i).folder, obj.data_files(i).name));
                fNames = fieldnames(dataStruct); rawData = dataStruct.(fNames{1});
                
                % --- STEP 1: Filter (Adaptive dimensions applied here) ---
                fprintf('     -> [Step 1] Filtering & Formatting...\n');
                [filteredData, currentMask] = obj.Step1_Filter(rawData, currentMask);
                
                % Update metadata mask to the dynamically sized one
                if ~isempty(currentMask), obj.metadata.roiMask = currentMask; end
                
                % --- STEP 2: Detect & Localize ---
                fprintf('     -> [Step 2] Localizing...\n');
                localizations = obj.localize_bubbles(filteredData, '       ');
                fprintf('        Found %d valid localizations.\n', height(localizations));
                
                % Offset frames for continuous time tracking
                localizations.Frame = localizations.Frame + (i-1)*obj.params.expParams.size(3);
                
                % Save localizations
                resultsPath = fullfile(obj.results_dir, 'LocalizationsTable.mat');
                if exist(resultsPath, 'file')
                    data = load(resultsPath);
                    bubbleLocalizations = [data.bubbleLocalizations; localizations];
                    save(resultsPath, 'bubbleLocalizations');
                else
                    bubbleLocalizations = localizations;
                    save(resultsPath, 'bubbleLocalizations');
                end
                
                % --- STEP 3: Track ---
                fprintf('     -> [Step 3] Tracking...\n');
                buffer_tracks = obj.track_bubbles(localizations, '       ');
                
                % Save Tracks
                track_filename = fullfile(obj.tracks_dir, sprintf('tracks_%03d.mat', i));
                obj.par_save(track_filename, buffer_tracks);
                
                fprintf('   - Buffer %d/%d finished in %.2f s.\n', i, num_buffers, toc(buffer_timer));
            end
            
            obj.metadata.processing_times.localization_and_tracking = toc(stage1_timer);
            fprintf('\n   - Total tracking time: %.2f minutes.\n', obj.metadata.processing_times.localization_and_tracking / 60);
        end

        function run_streaming_reconstruction(obj, minLength)
            % --- STAGE 2 & 3 COMBINED: Reconstruct maps by streaming from track files ---
            % This method avoids high memory usage by loading, rendering, and
            % discarding tracks from one file at a time.
        
            if nargin < 2, minLength = 1; end % Default to using all tracks
        
            fprintf('\n   - Starting streaming reconstruction (minLength = %d)...\n', minLength);
            reconstruction_timer = tic;
        
            track_files = dir(fullfile(obj.tracks_dir, 'tracks_*.mat'));
            if isempty(track_files)
                error('No intermediate track files found in %s. Stage 1 may have failed.', obj.tracks_dir);
            end
        
            % --- Initialize empty canvases for the final maps ---
            H = obj.params.expParams.size(1);
            W = obj.params.expParams.size(2);
            sRes_dims = [H, W] * obj.params.render.upsampling_factor;
        
            % These maps will be cumulatively added to in the loop
            final_density_map = zeros(sRes_dims);
            final_velocity_sum_map = zeros(sRes_dims);
            final_velocity_weight_map = zeros(sRes_dims);
        
            h_wait = waitbar(0, 'Streaming Reconstruction...', 'Name', 'ULM Post-Processing');
        
            % --- Loop through each intermediate track file ---
            for i = 1:length(track_files)
                waitbar(i/length(track_files), h_wait, sprintf('Processing track file %d of %d', i, length(track_files)));
                
                % Load tracks from one buffer
                loaded_data = load(fullfile(track_files(i).folder, track_files(i).name));
                
                % Instead of using raw tracks, process them first.
                if obj.params.track.enable_postprocessing 
                    tracks_to_render = obj.postProcessTracks(loaded_data.tracks, minLength);
                    fprintf('       - [PostProcess] File %d/%d: Kept %d of %d tracks for rendering (minLength >= %d).\n', ...
                     i, length(track_files), length(tracks_to_render), length(loaded_data.tracks), minLength);
                else
                    tracks_to_render = loaded_data.tracks;
                    if ~isempty(tracks_to_render)
                        all_lengths = [tracks_to_render.length];
                        keep_mask = all_lengths >= minLength;
                        tracks_to_render = tracks_to_render(keep_mask);
                    end
                end
                
                if isempty(tracks_to_render), continue; end
        
                % Render this small batch of tracks
                [density_chunk, vel_sum_chunk, vel_weight_chunk] = ...
                    obj.render_tracks(tracks_to_render, sRes_dims);
        
                % Add the results to the final canvases
                final_density_map = final_density_map + density_chunk;
                final_velocity_sum_map = final_velocity_sum_map + vel_sum_chunk;
                final_velocity_weight_map = final_velocity_weight_map + vel_weight_chunk;
            end
            close(h_wait);
        
            % --- Finalize and store the completed maps ---
            obj.densityMap = final_density_map;
        
            obj.velocityMap = zeros(sRes_dims);
            valid_pixels = final_velocity_weight_map > 1e-6;
            obj.velocityMap(valid_pixels) = final_velocity_sum_map(valid_pixels) ./ final_velocity_weight_map(valid_pixels);
        
            obj.metadata.processing_times.reconstruction = toc(reconstruction_timer);
            fprintf('   - Total STAGE 2 & 3 Streaming map reconstruction finished in %.2f s.\n', obj.metadata.processing_times.reconstruction);
        end

        function processed_tracks = postProcessTracks(obj, raw_tracks, min_track_length)
            % =========================================================================
            % PURPOSE:
            %   Filters, smooths, and interpolates raw tracks to prepare them
            %   for final rendering. This is the last step of track
            %   post-processing.
            %
            % WORKFLOW:
            %   1.  (Filter) Basic Length Filter: Removes tracks shorter
            %       than 'min_track_length'.
            %   2.  (Filter) QC Jitter Filter: (Optional) Calculates the
            %       Velocity Dispersion (VD) ratio (also known as Tortuosity).
            %       It rejects tracks that are too "jittery" or non-directional.
            %       This is controlled by 'params.track.qc.enable_vd_constraint'.
            %   3.  (Smooth) Savitzky-Golay Filter: Applies an 'sgolayfilt'
            %       to the path to remove high-frequency localization
            %       jitter while preserving the track's true shape better
            %       than a simple moving average.
            %   4.  (Interpolate) spline Interpolation: Upsamples the track
            %       path using a shape-preserving interpolator ('spline') to
            %       create a smoother, denser path for rendering.
            %   5.  (Recalculate) Velocity Calculation: Re-computes velocities
            %       based on the new smoothed and interpolated path.
            %
            % INPUTS:
            %   obj:                The ULM_Processor object instance.
            %   raw_tracks:         (struct array) The array of tracks to process.
            %   min_track_length:   (scalar) The minimum number of points a
            %                       track must have to be kept.
            %
            % OUTPUT:
            %   processed_tracks:   (struct array) A new array containing only
            %                       the valid, processed tracks.
            % =========================================================================
            
            if isempty(raw_tracks)
                processed_tracks = [];
                return;
            end
            
            % --- Get parameters from the object properties ---
            smoothing_factor = obj.params.track.smoothing_factor; 
            smoothing_method = obj.params.track.smoothing_method; % e.g., 'rloess', 'sgolay'
            
            interp_step = obj.params.render.interpolation_step; 
            interp_method = obj.params.render.interpolation_method; % e.g., 'spline', 'pchip'
            
            % --- Quality Control (QC) Parameters ---
            enable_vd_qc = obj.params.track.qc.enable_vd_constraint; 
            max_vd_ratio = obj.params.track.qc.max_vd_ratio; 
            
            % --- Savitzky-Golay Filter Settings (Specific for 'sgolay') ---
            sgolay_poly_order = 3; 
                                                               
            processed_tracks_cell = cell(1, length(raw_tracks));
            
            % Use a parallel loop for efficiency
            parfor i = 1:length(raw_tracks)
                track = raw_tracks(i);
                
                % --- STAGE 1: Basic Length Filter ---
                if track.length < min_track_length
                    processed_tracks_cell{i} = [];
                    continue; 
                end
                
                % --- STAGE 2: QC Jitter Filter (Conditional) ---
                is_valid_track = true; 
                
                if enable_vd_qc 
                    path = track.path;
                    net_displacement = norm(path(end,:) - path(1,:));
                    
                    if net_displacement < 1e-6 
                        is_valid_track = false; 
                    else
                        total_path_length = sum(sqrt(sum(diff(path).^2, 2)));
                        vd_ratio = total_path_length / net_displacement;
                        
                        if vd_ratio > max_vd_ratio 
                            is_valid_track = false; 
                        end
                    end
                end
                
                if ~is_valid_track
                    processed_tracks_cell{i} = [];
                    continue; 
                end
                
                % --- STAGE 3: Smoothing (Variable Method) ---
                path_to_process = track.path;
                
                % Check if track is long enough for the requested window
                actual_smoothing_window = min(smoothing_factor, track.length);
                
                switch lower(smoothing_method)
                    case 'sgolay'
                        % Savitzky-Golay Logic (Original)
                        % Requires odd window size
                        if mod(actual_smoothing_window, 2) == 0 
                            actual_smoothing_window = actual_smoothing_window - 1;
                        end
                        
                        actual_poly_order = min(sgolay_poly_order, actual_smoothing_window - 1);
                        
                        if track.length > actual_smoothing_window && actual_poly_order > 0
                            y_smooth = sgolayfilt(double(track.path(:,2)), actual_poly_order, actual_smoothing_window);
                            x_smooth = sgolayfilt(double(track.path(:,1)), actual_poly_order, actual_smoothing_window);
                            path_to_process = [x_smooth, y_smooth];
                        end
                        
                    case {'rloess', 'loess', 'gaussian', 'movmean', 'lowess'}
                        % Modern Matlab Smoothing (Supports 'rloess' - Robust Local Regression)
                        % smoothdata is robust and handles window sizes automatically better than older functions
                        
                        % Ensure window is at least 3 for meaningful smoothing, else skip
                        if actual_smoothing_window >= 3
                            x_smooth = smoothdata(track.path(:,1), smoothing_method, actual_smoothing_window);
                            y_smooth = smoothdata(track.path(:,2), smoothing_method, actual_smoothing_window);
                            path_to_process = [x_smooth, y_smooth];
                        end
                        
                    otherwise
                        % 'none' or unrecognized - do nothing, keep raw path
                end
        
                % --- STAGE 4: Interpolation (Upsampling) ---
                original_indices = 1:track.length;
                interp_indices = 1:interp_step:track.length;
                
                % Use the user-selected interpolation method (e.g., 'spline', 'pchip', 'linear')
                y_interp = interp1(original_indices, path_to_process(:,2), interp_indices, interp_method); 
                x_interp = interp1(original_indices, path_to_process(:,1), interp_indices, interp_method); 
                frames_interp = interp1(original_indices, track.frames, interp_indices, 'linear');
                
                % Clean up any potential NaNs from interpolation
                x_interp = fillmissing(x_interp, 'linear', 'EndValues', 'nearest');
                y_interp = fillmissing(y_interp, 'linear', 'EndValues', 'nearest');
                
                new_path = [x_interp', y_interp'];
                new_length = size(new_path, 1);
                
                % --- STAGE 5: Recalculate Velocities ---
                if new_length > 1
                    dx_pixels = diff(new_path(:,1));
                    dy_pixels = diff(new_path(:,2));
                    dt_frames = diff(frames_interp');
                    
                    displacement_mm = sqrt((dx_pixels * obj.params.track.pixel_X_size).^2 + (dy_pixels * obj.params.track.pixel_Z_size).^2);
                    dt_sec = dt_frames * obj.params.track.dt;
                    
                    velocities_mm_s = zeros(size(dt_sec));
                    valid_dt = dt_sec > 1e-9; 
                    velocities_mm_s(valid_dt) = displacement_mm(valid_dt) ./ dt_sec(valid_dt);
                    velocities_mm_s = [velocities_mm_s; velocities_mm_s(end)]; 
                    
                    mean_vel = mean(velocities_mm_s, 'omitnan');
                    if isnan(mean_vel), mean_vel = 0; end
                    avg_vel_vector = repmat(mean_vel, new_length, 1);
                else
                    velocities_mm_s = 0;
                    avg_vel_vector = 0;
                end
                
                val_intensity = 0;
                if isfield(track, 'localizations') && istable(track.localizations) ...
                        && ismember('Intensity', track.localizations.Properties.VariableNames)
                    val_intensity = mean(track.localizations.Intensity);
                end

                % Package the final processed track
                processed_tracks_cell{i} = struct(...
                    'id', track.id, ...
                    'path', new_path, ...
                    'frames', frames_interp', ...
                    'length', new_length, ...
                    'velocities_mm_s', velocities_mm_s, ...
                    'average_velocity_mm_s', avg_vel_vector, ...
                    'mean_intensity', val_intensity, ...
                    'localizations', [] ... 
                );
            end
            
            processed_tracks = [processed_tracks_cell{:}];
        end

        function saveResults(obj)
            % =====================================================================
            % SAVE RESULTS (AGGREGATED & OPTIMIZED)
            % This function aggregates track files from the disk, performs 
            % final analysis (Saturation, Features), saves a lightweight summary,
            % and exports visualizations.
            % =====================================================================
            
            fprintf('   - Saving results summary to: %s\n', obj.results_dir);
            
            % --- 1. Save Parameters ---
            params = obj.params;
            save(fullfile(obj.results_dir, 'params.mat'), 'params');
        
            % --- 2. Aggregate Tracks from Disk (Fast Vectorized Mode) ---
            % We scan the output folder for batch track files
            trackFiles = dir(fullfile(obj.tracks_dir, '*tracks*.mat'));
            
            if isempty(trackFiles)
                warning('No track files found in %s. Skipping aggregation.', obj.tracks_dir);
                return;
            end
            
            % --- STEP A: Sort Files Numerically (Critical!) ---
            % Files usually come as 'tracks_1.mat', 'tracks_10.mat', 'tracks_2.mat'.
            % We must sort them by the number in the filename to ensure correct time order.
            fileNumbers = zeros(length(trackFiles), 1);
            for i = 1:length(trackFiles)
                % Extract numbers using RegEx
                numStr = regexp(trackFiles(i).name, '\d+', 'match'); 
                if ~isempty(numStr)
                    % Usually the last number is the batch index
                    fileNumbers(i) = str2double(numStr{end}); 
                end
            end
            [~, sortedIdx] = sort(fileNumbers);
            trackFiles = trackFiles(sortedIdx); % Reorder the structure
            
            fprintf('     - Aggregating data from %d track files on disk (Sorted)...\n', length(trackFiles));
            
            % --- STEP B: Get Batch Size for Offset Calculation ---
            % User defined: expParams.size = [Height, Width, FramesPerBatch]
            if isfield(obj.params.expParams, 'size') && length(obj.params.expParams.size) >= 3
                frames_per_batch = obj.params.expParams.size(3);
            else
                warning('Batch size not found in params.expParams.size(3). Assuming 0 offset (Time will be wrong!).');
                frames_per_batch = 0; 
            end
            
            % Initialize storage
            all_IDs = {};
            all_Lengths = {};
            all_MeanVel = {};
            all_StartFrame = {};
            all_EndFrame = {};  
            all_Duration = {};   
            all_Paths = {}; 
            
            total_tracks_count = 0;
            
            % Loop over sorted batch files
            for k = 1:length(trackFiles)
                filePath = fullfile(trackFiles(k).folder, trackFiles(k).name);
                
                % Calculate the Frame Offset for this batch
                % If k=1 (Batch 1), offset is 0. If k=2, offset is 1*batch_size, etc.
                % Note: This assumes files are 1, 2, 3... sequential.
                current_frame_offset = (k - 1) * frames_per_batch;
                
                try
                    loadedData = load(filePath);
                    
                    % --- Robust Variable Finding ---
                    if isfield(loadedData, 'tracks')
                        batchTracks = loadedData.tracks;
                    elseif isfield(loadedData, 'finalTracks')
                        batchTracks = loadedData.finalTracks;
                    else
                        vars = fieldnames(loadedData);
                        batchTracks = loadedData.(vars{1});
                    end
                    
                    if isempty(batchTracks), continue; end
                    
                    % Ensure it is a struct array
                    if ~isstruct(batchTracks)
                        continue;
                    end
                    
                    numBatch = length(batchTracks);
                    total_tracks_count = total_tracks_count + numBatch;
                    
                    % --- VECTORIZED EXTRACTION (Fast) ---
                    
                    % 1. IDs and Lengths (scalars)
                    all_IDs{end+1} = [batchTracks.id]';
                    all_Lengths{end+1} = [batchTracks.length]';
                    
                    % 2. Paths (Cell array of matrices)
                    all_Paths{end+1} = {batchTracks.path}';
                    
                    % 3. Mean Velocity (Take 1st element)
                    vel_cells = {batchTracks.average_velocity_mm_s};
                    vel_vals = cellfun(@(v) v(1), vel_cells, 'UniformOutput', true)'; 
                    all_MeanVel{end+1} = vel_vals;
                    
                    % 4. Time / Frames (Start, End) with OFFSET CORRECTION
                    frame_cells = {batchTracks.frames};
                    
                    % Extract raw start/end from the batch-relative data
                    raw_start = cellfun(@(f) double(f(1)), frame_cells, 'UniformOutput', true)';
                    raw_end   = cellfun(@(f) double(f(end)), frame_cells, 'UniformOutput', true)';
                    
                    % Apply Offset: Transform to Global Experiment Time
                    global_start = raw_start + current_frame_offset;
                    global_end   = raw_end   + current_frame_offset;
                    
                    all_StartFrame{end+1} = global_start;
                    all_EndFrame{end+1}   = global_end;
                    
                    % Duration is independent of offset (End - Start is constant)
                    all_Duration{end+1}   = (raw_end - raw_start) * obj.params.track.dt;
        
                catch ME
                    warning('Error reading batch file %s: %s', trackFiles(k).name, ME.message);
                end
            end
            
            % --- 3. Build the Master Table & Save ---
            if total_tracks_count > 0
                % Concatenate all chunks
                TrackTable = table(...
                    vertcat(all_IDs{:}), ...
                    vertcat(all_Lengths{:}), ...
                    vertcat(all_MeanVel{:}), ...
                    vertcat(all_StartFrame{:}), ...
                    vertcat(all_EndFrame{:}), ...
                    vertcat(all_Duration{:}), ...
                    'VariableNames', {'TrackID', 'Length', 'MeanVelocity_mm_s', 'StartFrame', 'EndFrame', 'Duration_s'});
                
                % Combine Paths separately (Cell Array)
                TrackPaths = vertcat(all_Paths{:}); 
                
                % Check config for lightweight saving
                do_save_light = isfield(obj.params.io, 'save_lightweight') && obj.params.io.save_lightweight;
                
                if do_save_light
                    save(fullfile(obj.results_dir, 'TrackTable.mat'), ...
                         'TrackTable', 'TrackPaths', 'params');
                    fprintf('     - Saved aggregated lightweight data (%d total tracks).\n', total_tracks_count);
                end
                
                % --- 4. Run Analysis (Saturation & Features) ---
                
                % A. Saturation Graphs (Kinetics)
                try
                    fprintf('     - Generating Saturation Graphs...\n');
                    %analyze_Saturation_Kinetics(fullfile(obj.results_dir, 'TrackTable.mat'), obj.results_dir);
                catch ME
                    warning(ME.identifier, 'Saturation analysis failed: %s', ME.message);
                end
                
                % B. Deep Feature Extraction (ULM Feature Metrics)
                try
                    % Reconstruct minimal struct for compatibility
                    minimalTracks = struct('path', TrackPaths, ...
                                           'velocities_mm_s', num2cell(TrackTable.MeanVelocity_mm_s), ... 
                                           'length', num2cell(TrackTable.Length));
                    
                    % Run the math analysis
                    analysisStats = analyze_ULM_Features(minimalTracks, obj.params, obj.results_dir);

                    % --- GENERATE FIGURES (The part we filled in) ---
                    if isfield(obj.params.io, 'export_figures') && obj.params.io.export_figures && ~isempty(analysisStats)
                         fprintf('     - Generating Statistics Histograms...\n');
                         
                         % -- Figure 1: Velocity Histogram --
                         f_vel = figure('Visible', 'off');
                         histogram(TrackTable.MeanVelocity_mm_s, 100, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
                         xlabel('Velocity [mm/s]'); ylabel('Count');
                         title(sprintf('Velocity Distribution (Mean: %.2f mm/s)', analysisStats.Velocity.Mean));
                         grid on;
                         saveas(f_vel, fullfile(obj.results_dir, 'Hist_Velocity.png'));
                         saveas(f_vel, fullfile(obj.results_dir, 'Hist_Velocity.fig'));
                         close(f_vel);
                         
                         % -- Figure 2: Track Length Histogram --
                         f_len = figure('Visible', 'off');
                         histogram(TrackTable.Length, 'BinMethod', 'integers', 'FaceColor', [0.8 0.4 0.2]);
                         xlabel('Track Length [points]'); ylabel('Count');
                         title(sprintf('Track Length Distribution (Mean: %.1f)', mean(TrackTable.Length)));
                         grid on;
                         saveas(f_len, fullfile(obj.results_dir, 'Hist_Length.png'));
                         saveas(f_len, fullfile(obj.results_dir, 'Hist_Length.fig'));
                         close(f_len);

                         % -- Figure 3: Tortuosity (if calculated) --
                         if isfield(analysisStats, 'Tortuosity')
                             f_tort = figure('Visible', 'off');
                             bar(analysisStats.Tortuosity.HistEdges(1:end-1), analysisStats.Tortuosity.Histogram, 'FaceColor', [0.4 0.8 0.4]);
                             xlabel('Tortuosity Index'); ylabel('Count');
                             title(sprintf('Tortuosity (Median: %.3f)', analysisStats.Tortuosity.Median));
                             grid on;
                             saveas(f_tort, fullfile(obj.results_dir, 'Hist_Tortuosity.png'));
                             saveas(f_tort, fullfile(obj.results_dir, 'Hist_Tortuosity.fig'));
                             close(f_tort);
                         end
                         
                         fprintf('     - Saved statistic figures to disk.\n');
                    end
                catch ME
                     warning(ME.identifier, 'Feature analysis failed: %s', ME.message);
                end
                
            else
                fprintf('     - No tracks found to aggregate.\n');
            end

            % --- 6. TIFF Exports (Density Maps) ---
            if obj.params.io.export_tiff 
                if ~isempty(obj.densityMap)
                    normDensity = obj.densityMap / (max(obj.densityMap(:)) + eps);
                    imwrite(normDensity, fullfile(obj.results_dir, 'density_map.tif'));
                    fprintf('     - Exported density map TIFF.\n');
                end
                if ~isempty(obj.velocityMap)
                     % Save normalized velocity map (simple grayscale visualization)
                     % For better visualization, consider saving the raw map or using a colormap function externally
                     imwrite(mat2gray(obj.velocityMap), fullfile(obj.results_dir, 'velocity_map.tif'));
                end
            end
            
            fprintf('   - Save complete.\n');
        end


    end

    methods (Access = public)
        % --- Private helper methods for processing steps ---
        
        function filteredData = filter_data(obj, rawData, indent_prefix)
                if nargin < 3, indent_prefix = '       '; end % Default indent
                
                % --- Step 1: Apply the primary SVD-based clutter filter ---
                switch obj.params.filter.method
                   case 'svd_filter'
                       % Call the basic SVD filter with a manual cutoff
                       filteredData = SVD_filter(rawData, obj.params.filter.svd_cutoff, ...
                                                 'IndentPrefix', indent_prefix);
            
                   case 'svd_ssm'
                       % Call the adaptive SVD filter using Spatial Similarity Matrix
                       filteredData = SVD_SSM(rawData, 'IndentPrefix', indent_prefix);
            
                   case 'dcc_svd'
                       % Call the adaptive SVD filter using feature-based clustering
                       framerate = obj.params.acq.framerate;
                       filteredData = DCC_SVD(rawData, framerate, ...
                                              'ReconstructionMode', 'blood', ...
                                              'IndentPrefix', indent_prefix);

                   case 'svd_blockwise'
                       bw = obj.params.filter.blockwise;   % shorthand
                    
                       % Build Name-Value pair list dynamically
                       nv_args = { ...
                           'ThresholdMethod',        bw.threshold_method, ...
                           'OverlapPct',             bw.overlap_pct, ...
                           'ManualCutoff',           bw.manual_cutoff, ...
                           'TissueFreqThreshHz',     bw.tissue_freq_hz, ...
                           'MPDeviationSigma',       bw.mp_deviation_sigma, ...
                           'GradientInflectionPct',  bw.gradient_pct, ...
                           'MinBloodComponents',     bw.min_blood_comps, ...
                           'MaxTissueFraction',      bw.max_tissue_frac, ...
                           'PlotThresholdMaps',      bw.plot_maps, ...
                           'IndentPrefix',           indent_prefix ...
                       };
                    
                       % Prefer physical mm spec if pixel size is available
                       if ~isempty(bw.block_size_mm) && ~isnan(bw.block_size_mm)
                           nv_args = [nv_args, {'BlockSizeMm', [bw.block_size_mm, bw.block_size_mm]}];
                       end
                    
                       [filteredData, bw_diag] = SVD_blockwise(rawData, obj.params, nv_args{:});
                    
                       % Optionally store diagnostics for later inspection
                       obj.params.filter.blockwise.last_diagnostics = bw_diag;

                   otherwise
                       error('Unknown clutter filter method: %s', obj.params.filter.method);
                end
                
                % --- Step 2 (Optional): Apply a Butterworth bandpass filter on the result ---
                % This step is controlled by a single flag for all methods.
                if isfield(obj.params.filter, 'enable_butterworth') && obj.params.filter.enable_butterworth
                    fprintf('%s -> Applying Butterworth bandpass filter...\n', indent_prefix);
            
                    % Get filter parameters from the object properties
                    cutoffFreq = obj.params.filter.butter_cutoff;   % e.g., [50 250] in Hz
                    samplingFreq = obj.params.acq.framerate;        % Acquisition framerate in Hz
                    filterOrder = obj.params.filter.butter_order;   % e.g., 2 for a second-order filter
            
                    % Call the dedicated Butterworth filter function
                    filteredData = Butterworth_bandpass_filter(filteredData, cutoffFreq, samplingFreq, filterOrder);
                end
        end

        function candidateBubbles = detect_candidates(obj, filteredData)
            % Dispatches to the correct detection function based on
            % params.loc.DetectMethod ('Intensity', 'NP', or 'NCC').
            % The roiMask is retrieved from obj.metadata and passed to NP/NCC
            % so they can apply it AFTER their internal statistics computation.
        
            % Retrieve mask if available (may be empty if masking is disabled)
            roiMask = [];
            if obj.params.proc.enableInteractiveMask && isfield(obj.metadata, 'roiMask')
                roiMask = obj.metadata.roiMask;
            end
        
            switch upper(obj.params.loc.DetectMethod)
        
                case 'INTENSITY'
                    % Classic intensity-based local-maxima detection.
                    % The mask was already applied to filteredData in Step1_Filter.
                    candidateBubbles = detectBubbles(filteredData, obj.params.loc);
        
                case 'NP'
                    % Neyman-Pearson criterion (Corazza et al., 2023).
                    % roiMask is applied internally AFTER per-pixel noise estimation
                    % so that the median/MAD are unbiased.
                    candidateBubbles = detectBubbles_NP(filteredData, obj.params.loc, roiMask);
        
                case 'NCC'
                    % Normalised Cross-Correlation with PSF template.
                    % Requires params.loc.MB_image to be pre-computed (done in
                    % setDefaultParams > calculateDerivedParams).
                    if ~isfield(obj.params.loc, 'MB_image') || isempty(obj.params.loc.MB_image)
                        error(['ULM_Processor:detect_candidates — NCC method selected but ' ...
                               'params.loc.MB_image is empty. Check that calculateDerivedParams ' ...
                               'ran successfully and that params.loc.psf_type is set.']);
                    end
                    candidateBubbles = detectBubbles_NCC(filteredData, obj.params.loc, roiMask);
        
                otherwise
                    error(['ULM_Processor:detect_candidates — Unknown DetectMethod: "%s". ' ...
                          'Valid options are ''Intensity'', ''NP'', ''NCC''.', ...
                          obj.params.loc.DetectMethod]);
            end
        end
        
        function localizations = localize_bubbles(obj, filteredData, indent_prefix)
            if nargin < 3, indent_prefix = '       '; end % Default indent
            candidateBubbles = obj.detect_candidates(filteredData);
            fprintf('%s-> Detection: %s found %d candidates.\n', ...
                    indent_prefix, obj.params.loc.DetectMethod, height(candidateBubbles));
        
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
            if nargin < 3, indent_prefix =  '       '; end % Default indent
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
            switch obj.params.render.method
                case 'histogram'
                    [density, velocity, weights] = renderHistogram(tracks, sRes_dims, obj.params.render);
                case 'gaussian'
                    [density, velocity, weights] = renderGaussian(tracks, sRes_dims, obj.params.render);
                otherwise
                    error('Unknown reconstruction method: %s', obj.params.render.method);
            end
        end

        function cropBox = getROICrop(obj, dataFrame)
            fig = figure; imshow(dataFrame, []); title('Draw Crop Region');
            roi = drawrectangle;
            cropBox = round(roi.Position); % [xmin, ymin, width, height]
            
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            filename = sprintf('CropBox_%s.mat', timestamp);
            save(fullfile(obj.params.proc.cropPath, filename), 'cropBox');
            fprintf('     -> Crop saved to: %s\n', filename);
            close(fig);
        end

        function mask = getROIMask(obj, dataFrame)
            fig = figure; imshow(dataFrame, []); title('Draw Region of Interest (ROI)');
            roi_2 = images.roi.AssistedFreehand; draw(roi_2);
            mask = createMask(roi_2, size(dataFrame,1), size(dataFrame,2));
            
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            filename = sprintf('Mask_%s.mat', timestamp);
            save(fullfile(obj.params.proc.maskPath, filename), 'mask');
            fprintf('     -> Mask saved to: %s\n', filename);
            close(fig);
        end

        function validLocs = validate_localizations(obj, localizations)
            % Removes localizations that fall outside the active masks.
            
            if isempty(localizations)
                validLocs = localizations;
                return;
            end

            % 1. Combine all active masks into one logical mask
            % We assume the image size matches the experiment params
            H = obj.params.expParams.size(1);
            W = obj.params.expParams.size(2);
            combinedMask = true(H, W); % Start with all valid
            maskActive = false;

            % A. Check Vessel Mask
            if isfield(obj.params.proc, 'vesselMask') && ...
               isfield(obj.params.proc.vesselMask, 'enable') && ...
               obj.params.proc.vesselMask.enable && ...
               isfield(obj.params.proc.vesselMask, 'vesselMask') && ...
               ~isempty(obj.params.proc.vesselMask.vesselMask)
           
                combinedMask = combinedMask & logical(obj.params.proc.vesselMask.vesselMask);
                maskActive = true;
            end

            % B. Check Interactive ROI Mask
            if isfield(obj.params.proc, 'enableInteractiveMask') && ...
               obj.params.proc.enableInteractiveMask && ...
               isfield(obj.metadata, 'roiMask') && ...
               ~isempty(obj.metadata.roiMask)
           
                combinedMask = combinedMask & logical(obj.metadata.roiMask);
                maskActive = true;
            end
            
            % If no masks are active, return original data
            if ~maskActive
                validLocs = localizations;
                return;
            end

            % 2. Check points against the combined mask
            x_idx = round(localizations.X);
            y_idx = round(localizations.Y);

            % Ensure indices are within image bounds first
            in_bounds = (x_idx >= 1 & x_idx <= W & y_idx >= 1 & y_idx <= H);
            
            % Check if the pixel is true in the mask
            keep_mask = false(size(in_bounds));
            
            % Only check mask value for points that are within image bounds
            if any(in_bounds)
                % sub2ind allows us to check specific (y,x) coordinates in the mask matrix
                indices = sub2ind([H, W], y_idx(in_bounds), x_idx(in_bounds));
                keep_mask(in_bounds) = combinedMask(indices);
            end

            % 3. Filter the table
            validLocs = localizations(keep_mask, :);
            
            num_removed = height(localizations) - height(validLocs);
            if num_removed > 0
                fprintf('      -> Mask Safety Check: Removed %d points that drifted outside the mask.\n', num_removed);
            end
        end
        
        function par_save(~, fname, tracks)
            % Helper function for saving inside a parfor loop
            % Includes a retry mechanism for network drive latency/locks
            maxRetries = 3;
            for attempt = 1:maxRetries
                try
                    save(fname, 'tracks');
                    break; % Success, exit the loop
                catch ME
                    if attempt == maxRetries
                        rethrow(ME); % Fail completely on the last attempt
                    end
                    warning('Failed to save %s. Retrying in 1 second... (%d/%d)', fname, attempt, maxRetries);
                    pause(1); % Wait 1 second before retrying
                end
            end
        end

    end
end