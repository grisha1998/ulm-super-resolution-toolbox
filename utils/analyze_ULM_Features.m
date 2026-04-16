function analysisStats = analyze_ULM_Features(bubbleTracks, params, resultsDir)
% =========================================================================
% FUNCTION: analyze_ULM_Features
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Extracts quantitative, physically meaningful biomarkers from raw ULM 
%   microbubble tracking data. This function bridges the gap between raw 
%   trajectories and clinical/ULM Feature Metrics. 
%   - Advantages: It is highly versatile, operating seamlessly either as an 
%     integrated step within the automated `ULM_Processor` pipeline or as a 
%     standalone interactive GUI script. It computes comprehensive statistics 
%     (Velocity, Tortuosity, Vessel Density, and Intensity) and automatically 
%     generates publication-ready statistical histograms.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Mode Selection: Detects if input arguments were provided. If not, it 
%      launches an interactive UI prompt to load a `TrackTable.mat` file 
%      (Standalone Mode).
%   2. Parameter Extraction: Robustly extracts the physical pixel size (in mm) 
%      from the provided parameters to ensure spatial metrics are physically accurate.
%   3. Kinematic Loop: Iterates through all tracks to calculate:
%      a. Mean velocities and intensities.
%      b. Tortuosity Index: The ratio of the actual arc length (sum of step-by-step 
%         Euclidean distances) to the chord length (straight-line distance 
%         from start to end). A perfectly straight vessel yields 1.0.
%   4. Statistical Aggregation: Compiles histograms and computes robust central 
%      tendencies (Mean, Median, Standard Deviation) for all features.
%   5. Perfusion / Density Analysis: Overlays a grid (default 0.5 mm step) 
%      over the FOV to calculate total vessel length density (mm/mm^2) and 
%      grid occupancy percentage.
%   6. Export & Visualization: Saves the computed statistics to a `.mat` file 
%      and conditionally exports visual histograms for velocity, length, and tortuosity.
%
% SYNTAX OPTIONS:
%   % Standalone Interactive Mode
%   analysisStats = analyze_ULM_Features() 
%
%   % Pipeline Mode
%   analysisStats = analyze_ULM_Features(bubbleTracks, params, resultsDir)
%
% EXAMPLES:
%   % Execute within a pipeline and save results to a specific folder:
%   stats = analyze_ULM_Features(finalTracks, pipelineParams, 'C:\Results\Run1\');
%
% INPUTS:
%   bubbleTracks - (Type: Struct Array) Array of valid track structures containing 
%                  '.path', '.velocities_mm_s', '.length', and optionally '.mean_intensity'.
%   params       - (Type: Struct) Configuration parameters. Crucially relies on 
%                  'params.expParams.pixel_size_mm' or 'pixel_X_size' for scaling.
%   resultsDir   - (Type: String) Path to the directory where the output '.mat' 
%                  and figures will be saved.
%
% OUTPUTS:
%   analysisStats - (Type: Struct) Comprehensive results including:
%                   .Velocity (Mean, Median, Std, HistCounts)
%                   .Tortuosity (Mean, Median, Histogram)
%                   .Density (OccupancyPercent, VesselLengthDensityGlobal)
%
% AUTHOR: Grigori Shapiro
% =========================================================================

    % --- 0. Standalone Mode Handling ---
    isStandalone = false;
    TrackTable = []; % Initialize for later plotting if needed
    
    if nargin < 1
        fprintf('   > Standalone Mode: Please select the TrackTable file...\n');
        
        % 1. Select File
        [fname, fpath] = uigetfile('*TrackTable.mat', 'Select Track Data File');
        if isequal(fname,0)
            analysisStats = [];
            return; 
        end
        inputFile = fullfile(fpath, fname);
        
        % 2. Load Data
        fprintf('     Loading: %s ...\n', fname);
        data = load(inputFile);
        
        if ~isfield(data, 'TrackTable') || ~isfield(data, 'TrackPaths') || ~isfield(data, 'params')
            error('Invalid file. Must contain: TrackTable, TrackPaths, and params.');
        end
        
        TrackTable = data.TrackTable;
        TrackPaths = data.TrackPaths;
        params = data.params;
        
        % 3. Transform Data to Standard Format
        % We create a minimal structure array compatible with the existing logic
        bubbleTracks = struct('path', TrackPaths, ...
                              'velocities_mm_s', num2cell(TrackTable.MeanVelocity_mm_s), ... 
                              'length', num2cell(TrackTable.Length));
        
        % 4. Select Results Directory
        selDir = uigetdir(fpath, 'Select Output Folder for Results');
        if isequal(selDir, 0)
            resultsDir = fpath; % Default to source folder if cancelled
        else
            resultsDir = selDir;
        end
        
        isStandalone = true;
    end

    fprintf('   > Starting Deep Feature Extraction (ULM Feature Metrics)...\n');

    % --- 1. Validation ---
    if isempty(bubbleTracks)
        warning('No tracks available for analysis.');
        analysisStats = [];
        return;
    end

    % --- 2. ROBUST PIXEL SIZE EXTRACTION ---
    % Critical: We need to know the pixel size to convert pixels to mm.
    if isfield(params, 'expParams') && isfield(params.expParams, 'pixel_size_mm')
        pixel_size = params.expParams.pixel_size_mm;
    elseif isfield(params, 'expParams') && isfield(params.expParams, 'pixel_X_size')
        pixel_size = mean([params.expParams.pixel_X_size, params.expParams.pixel_Z_size]);
    elseif isfield(params, 'track') && isfield(params.track, 'pixel_X_size')
        pixel_size = params.track.pixel_X_size;
    else
        warning('Pixel size not found in params. Assuming 0.1 mm.');
        pixel_size = 0.1; 
    end

    % --- 3. PRE-ALLOCATION & DATA EXTRACTION LOOP ---
    numTracks = length(bubbleTracks);
    
    all_velocities = [];        
    track_lengths_mm = zeros(numTracks, 1);
    tortuosity_index = zeros(numTracks, 1);
    mean_intensities = zeros(numTracks, 1);
    
    for i = 1:numTracks
        track = bubbleTracks(i);
        path = track.path; % [x, z] in pixels
        
        % --- A. Velocity Collection ---
        if ~isempty(track.velocities_mm_s)
            all_velocities = [all_velocities; track.velocities_mm_s]; %#ok<AGROW>
        end
        
        % --- B. Intensity (Robust Check) ---
        if isfield(track, 'mean_intensity') && ~isempty(track.mean_intensity)
             mean_intensities(i) = track.mean_intensity;
        elseif isfield(track, 'localizations') && istable(track.localizations) ...
                && ismember('Intensity', track.localizations.Properties.VariableNames)
             mean_intensities(i) = mean(track.localizations.Intensity);
        else
             mean_intensities(i) = 0; 
        end
        
        % --- C. Tortuosity Calculation ---
        % 1. Arc Length (Sum of all steps)
        steps = sqrt(sum(diff(path, 1, 1).^2, 2));
        arc_length_px = sum(steps);
        track_lengths_mm(i) = arc_length_px * pixel_size;
        
        % 2. Chord Length (Euclidean distance start-to-end)
        chord_length_px = sqrt(sum((path(end,:) - path(1,:)).^2));
        
        % 3. Calculate Ratio
        if chord_length_px > 1e-6 
            tortuosity_index(i) = arc_length_px / chord_length_px;
        else
            tortuosity_index(i) = 1.0; 
        end
    end

    % --- 4. STATISTICAL ANALYSIS ---
    analysisStats = struct();
    analysisStats.generated_at = datetime('now');
    analysisStats.NumTracks = numTracks;
    analysisStats.PixelSizeUsed = pixel_size;

    % -- Histogram Parameters --
    if isfield(params, 'analysis')
        v_bins = params.analysis.velocity_hist_num_bins;
        t_edges = params.analysis.tortuosity_bins;
        if isscalar(t_edges), t_edges = linspace(1, 5, t_edges); end
    else
        v_bins = 50; 
        t_edges = 1.0 : 0.05 : 4.0; 
    end

    % --- A. Velocity Statistics ---
    if ~isempty(all_velocities)
        analysisStats.Velocity.Mean = mean(all_velocities);
        analysisStats.Velocity.Median = median(all_velocities);
        analysisStats.Velocity.Std = std(all_velocities);
        [v_counts, v_edges] = histcounts(all_velocities, v_bins); 
        analysisStats.Velocity.HistCounts = v_counts;
        analysisStats.Velocity.HistEdges = v_edges;
    else
        analysisStats.Velocity.Mean = 0;
        analysisStats.Velocity.HistCounts = [];
        analysisStats.Velocity.HistEdges = [];
    end

    % --- B. Tortuosity Statistics ---
    analysisStats.Tortuosity.Mean = mean(tortuosity_index);
    analysisStats.Tortuosity.Median = median(tortuosity_index);
    [t_counts, ~] = histcounts(tortuosity_index, t_edges);
    analysisStats.Tortuosity.Histogram = t_counts;
    analysisStats.Tortuosity.HistEdges = t_edges;

    % --- C. Intensity Statistics ---
    analysisStats.Intensity.MeanGlobal = mean(mean_intensities);
    
    % --- D. Density / Perfusion Analysis (Grid Based) ---
    if isfield(params, 'analysis') && isfield(params.analysis, 'density_grid_size_mm')
        grid_step_mm = params.analysis.density_grid_size_mm;
    else
        grid_step_mm = 0.5;
    end
    
    if isfield(params.expParams, 'fovX') && ~isnan(params.expParams.fovX)
        fov_width_mm = params.expParams.fovX;
        fov_height_mm = params.expParams.fovZ;
    else
        all_paths = vertcat(bubbleTracks.path);
        fov_width_mm = max(all_paths(:,1)) * pixel_size;
        fov_height_mm = max(all_paths(:,2)) * pixel_size;
    end
    
    x_edges = 0 : grid_step_mm : fov_width_mm;
    z_edges = 0 : grid_step_mm : fov_height_mm;
    
    all_x_mm = []; all_z_mm = [];
    for k = 1:numTracks
        p = bubbleTracks(k).path; 
        all_x_mm = [all_x_mm; p(:,1) * pixel_size]; %#ok<AGROW>
        all_z_mm = [all_z_mm; p(:,2) * pixel_size]; %#ok<AGROW>
    end
    
    if ~isempty(all_x_mm)
        [counts_grid, ~, ~] = histcounts2(all_x_mm, all_z_mm, x_edges, z_edges);
        num_occupied_cells = sum(counts_grid(:) > 0);
        analysisStats.Density.OccupancyPercent = (num_occupied_cells / numel(counts_grid)) * 100;
        non_zero_counts = counts_grid(counts_grid > 0);
        analysisStats.Density.Heterogeneity = std(non_zero_counts) / mean(non_zero_counts);
    else
        analysisStats.Density.OccupancyPercent = 0;
        analysisStats.Density.Heterogeneity = 0;
    end
    
    analysisStats.Density.VesselLengthDensityGlobal = sum(track_lengths_mm) / (fov_width_mm * fov_height_mm);

    % --- 5. SAVE RESULTS ---
    if ~isempty(resultsDir)
        if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
        savePath = fullfile(resultsDir, 'analysisStats.mat');
        save(savePath, 'analysisStats');
    end
    
    % --- 6. STANDALONE PLOTTING ---
    % This block runs ONLY if the script was called without arguments
    if isStandalone && ~isempty(analysisStats) && ~isempty(TrackTable)
         fprintf('     - Generating Statistics Histograms (Standalone Mode)...\n');
         
         % -- Figure 1: Velocity Histogram --
         f_vel = figure('Visible', 'off');
         histogram(TrackTable.MeanVelocity_mm_s, 100, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
         xlabel('Velocity [mm/s]'); ylabel('Count');
         title(sprintf('Velocity Distribution (Mean: %.2f mm/s)', analysisStats.Velocity.Mean));
         grid on;
         saveas(f_vel, fullfile(resultsDir, 'Hist_Velocity.png'));
         saveas(f_vel, fullfile(resultsDir, 'Hist_Velocity.fig'));
         close(f_vel);
         
         % -- Figure 2: Track Length Histogram --
         f_len = figure('Visible', 'off');
         histogram(TrackTable.Length, 'BinMethod', 'integers', 'FaceColor', [0.8 0.4 0.2]);
         xlabel('Track Length [points]'); ylabel('Count');
         title(sprintf('Track Length Distribution (Mean: %.1f)', mean(TrackTable.Length)));
         grid on;
         saveas(f_len, fullfile(resultsDir, 'Hist_Length.png'));
         saveas(f_len, fullfile(resultsDir, 'Hist_Length.fig'));
         close(f_len);
         
         % -- Figure 3: Tortuosity (if calculated) --
         if isfield(analysisStats, 'Tortuosity')
             f_tort = figure('Visible', 'off');
             bar(analysisStats.Tortuosity.HistEdges(1:end-1), analysisStats.Tortuosity.Histogram, 'FaceColor', [0.4 0.8 0.4]);
             xlabel('Tortuosity Index'); ylabel('Count');
             title(sprintf('Tortuosity (Median: %.3f)', analysisStats.Tortuosity.Median));
             grid on;
             saveas(f_tort, fullfile(resultsDir, 'Hist_Tortuosity.png'));
             saveas(f_tort, fullfile(resultsDir, 'Hist_Tortuosity.fig'));
             close(f_tort);
         end
         
         fprintf('     - Saved statistic figures to: %s\n', resultsDir);
    end
    
    % --- 7. CONSOLE REPORT ---
    fprintf('     -------------------------------------------------\n');
    fprintf('     [ULM Feature Metrics REPORT]\n');
    fprintf('     -------------------------------------------------\n');
    fprintf('       |-> Pixel Size Used: %.3f mm\n', pixel_size);
    fprintf('       |-> Tracks Analyzed: %d\n', numTracks);
    fprintf('       |-> Mean Velocity:   %.2f mm/s\n', analysisStats.Velocity.Mean);
    fprintf('       |-> Mean Intensity:  %.2f a.u.\n', analysisStats.Intensity.MeanGlobal);
    fprintf('       |-> Mean Tortuosity: %.4f (1.0 = Straight)\n', analysisStats.Tortuosity.Mean);
    fprintf('       |-> Vessel Density:  %.2f mm/mm^2\n', analysisStats.Density.VesselLengthDensityGlobal);
    fprintf('       |-> Grid Occupancy:  %.2f%% (Grid Step: %.1f mm)\n', ...
        analysisStats.Density.OccupancyPercent, grid_step_mm);
    fprintf('     -------------------------------------------------\n');
    fprintf('     - Metrics saved to: ULM_Analysis_Metrics.mat\n');
end