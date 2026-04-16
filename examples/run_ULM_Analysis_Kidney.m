% =================================================================================
% FILENAME: run_ULM_Analysis_Kidney.m
% =================================================================================
%
% PURPOSE:
%   A unified, self-contained runner script for the Ultrasound Localization 
%   Microscopy (ULM) framework. This script orchestrates the entire ULM 
%   processing pipeline from raw data to final super-resolution maps. It 
%   features a specialized, highly configurable workflow designed specifically 
%   for in-vivo kidney data, handling complex pre-processing tasks such as 
%   respiratory motion detection, frame deletion, and tissue registration.
%
% SCRIPT WORKFLOW:
%   1. CONFIGURATION: 
%      - Sets boolean flags to activate/deactivate the kidney-specific workflow.
%      - Configures fine-grained controls for motion analysis, bad frame removal, 
%        data registration, and overlay displays.
%
%   2. PARAMETER INITIALIZATION: 
%      - Loads default and experiment-specific parameters (via `setDefaultParams`).
%      - Initializes global structures for IO, acquisition, filtering, and rendering.
%
%   3. ULM PROCESSOR INITIALIZATION: 
%      - Instantiates the main `ULM_Processor` object.
%
%   4. OPTIONAL IN-VIVO KIDNEY PRE-PROCESSING:
%      (Executes only if `isKidneyExperiment = true`)
%      - Data Batching: Divides superframes into manageable batches (e.g., 1500 frames).
%      - ROI Definition: Allows interactive creation of a Region of Interest mask.
%      - Motion Detection: Calculates cross-correlation to identify breathing artifacts.
%      - Frame Deletion: Removes corrupted frames (creates `BreathDeletion` dataset).
%      - Registration: Performs rigid motion correction (and optionally non-rigid
%        B-spline correction) to stabilize the tissue (creates `registered` dataset).
%      - Processor Update: Points the `ULM_Processor` to the newly cleaned/registered data.
%
%   5. ULM PIPELINE EXECUTION:
%      - Boundary Initialization: Defines crop and mask boundaries.
%      - B-Mode Generation: Generates and saves mean B-mode images for vessels and tissue.
%      - Localization & Tracking: Executes the core ULM loop (`run_localization_and_tracking_loop`),
%        extracting bubble coordinates and connecting them into tracks.
%
%   6. RENDERING & FILTERING:
%      - Iterates through a defined array of `MIN_TRACK_LENGTH` values.
%      - Runs streaming reconstruction for each track length threshold to filter noise.
%
%   7. FINAL OUTPUT & VISUALIZATION:
%      - Saves numerical results to the general results directory.
%      - Generates high-quality, publication-ready figures for each tracking threshold.
%
% EXPECTED OUTPUTS:
%   Data Files (.mat):
%      - ROI.mat, Mask.mat: Saved region of interest definitions.
%      - correlation_*.mat: Cross-correlation arrays (with and without mask).
%      - indices_to_remove.mat: List of frames flagged as breathing artifacts.
%      - mean_bmode_vessel_*.mat / mean_bmode_tissue_*.mat: Mean intensity matrices.
%      - Final ULM results (Tracks, Coordinates, Density/Velocity Maps).
%
%   Images (.png):
%      - Deleting_breathing_images_robust.png: Motion detection visualization.
%      - Motion_Correction_Summary.png: Translation and rotation tracking over time.
%      - Density_Map_minLen_*.png: Super-resolution bubble density map.
%      - Velocity_Map_Filtered_minLen_*.png: Gaussian-filtered velocity map.
%      - Velocity_Map_Unfiltered_minLen_*.png: Raw velocity map.
%      - Overlay_minLen_*.png: Density map overlaid onto the mean B-mode background.
%      - Combined_VelDens_Map_minLen_*.png: HSV map showing Velocity (Color) and Density (Brightness).
%
% AUTHOR: Grigori Shapiro
% DATE: July 28, 2025
% =================================================================================

%% ========================================================================
%  0. SCRIPT SETUP & PREPARATION
%  ========================================================================

clear; 
close all; 
clc;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

fprintf('=======================================================\n');
fprintf('    MATLAB Ultrasound Localization Microscopy Pipeline\n');
fprintf('=======================================================\n\n');

%% ========================================================================
%  1. USER CONFIGURATION & WORKFLOW CONTROL
%  ========================================================================
fprintf('\n>> STEP 1: CONFIGURATION\n');

% --- Primary Workflow Switch ---
% Set to 'true' to activate all kidney-specific pre-processing steps.
% Set to 'false' to run the standard pipeline (e.g., for phantom data).
isKidneyExperiment = true;  % false, true

% --- Data Folder ---
data_folder = '';
data_subfolder2 = 'Bmode'; % Bmode or PI

% --- Kidney Workflow Fine-Grained Controls ---
% These switches only have an effect if 'isKidneyExperiment' is true.
% They allow you to enable or disable specific, time-consuming steps.
kidney_workflow.run_motion_analysis      = true; % Find and register frames with motion.
kidney_workflow.clear_bad_frames         = true; % Remove the frames identified by motion analysis.
kidney_workflow.use_registered_data      = true;  % Point the processor to motion-corrected data.
kidney_workflow.use_cleaned_data         = true;  % Point the processor to data with bad frames removed.
kidney_workflow.display_kidney_overlay   = true;  % Generate B-mode overlay image at the end.
kidney_workflow.apply_roi_mask           = false; % Apply the ROI mask before processing?

MIN_TRACK_LENGTH = 3; %<-- TUNE THIS PARAMETER

fprintf('   - Kidney-specific workflow is %s.\n', B_SWITCH(isKidneyExperiment));
fprintf('   - ROI Masking is %s.\n', B_SWITCH(kidney_workflow.apply_roi_mask));

%% ========================================================================
%  2. PARAMETER INITIALIZATION
%  ========================================================================
fprintf('\n>> STEP 2: PARAMETER INITIALIZATION\n');
try
    % Call the main parameter function, passing the workflow switches.
    % This function now handles loading defaults, info.txt, sample data,
    % calculating derived parameters, and loading the mask.
    if exist('data_folder', 'var')
        params = setDefaultParams(isKidneyExperiment, data_folder);
    else
        params = setDefaultParams(isKidneyExperiment);
    end
    fprintf('   - Parameters initialized and calculated successfully.\n');
catch ME
    fprintf('   - ERROR: Failed to initialize parameters. %s\n', ME.message);
    return;
end

%% ========================================================================
%  3. ULM PROCESSOR INITIALIZATION
%  ========================================================================
fprintf('\n>> STEP 3: INITIALIZING ULM PROCESSOR\n');
try
    % The processor is initialized with the final, potentially corrected, data path.
    processor = ULM_Processor(params);
catch ME
    fprintf('   - ERROR: Failed to initialize ULM_Processor: %s\n', ME.message);
    return;
end

%% ========================================================================
%  4. OPTIONAL: IN-VIVO KIDNEY PRE-PROCESSING WORKFLOW
%  ========================================================================
if isKidneyExperiment
    fprintf('\n>> STEP 4: RUNNING KIDNEY PRE-PROCESSING WORKFLOW\n');
    
    % Define paths for kidney data processing
    bmode_folder = fullfile(params.io.data_folder, params.expParams.bubbleType, 'Bmode');
    pi_folder = fullfile(params.io.data_folder, params.expParams.bubbleType, 'PI');

    % --- 4a. Motion Analysis and Image Registration ---
    if kidney_workflow.run_motion_analysis
        fprintf('   - Running motion analysis and registration...\n');
        t_motion = tic;
        
        IQfiles_Bmode = dir(fullfile(bmode_folder, '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        IQfiles_PI = dir(fullfile(pi_folder, '*.mat')); IQfiles_PI = sortFilesByNumber(IQfiles_PI);
        
        % Dividing the superframes into smaller batches of 1500 frames - for both Bmode and PI
        fprintf('--- Divides the superframes into batches --- \n\n')
        split_ImageData_tot_Kidney(fullfile(params.io.data_folder, params.expParams.bubbleType), params.expParams.bubbleType, 1500)
        IQfiles_Bmode = dir(fullfile(bmode_folder, '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        IQfiles_PI = dir(fullfile(pi_folder, '*.mat')); IQfiles_PI = sortFilesByNumber(IQfiles_PI);

        %Build a ROI mask  
        fprintf('--- Build a ROI mask --- \n\n')
        [mask, roi] = buildROIMask(IQfiles_Bmode);
        save([processor.general_results_dir filesep 'ROI.mat'], 'roi');
        save([processor.general_results_dir filesep 'Mask.mat'], 'mask');
        
        load([[processor.general_results_dir filesep 'Mask'] '.mat']);
    
        % Correlation calculation with and without mask
        fprintf('--- Correlation calculation with and without mask --- \n\n')
        [correlation_no_mask, correlation_with_mask, fig1, fig2] = calculateAndPlotCorrelation(IQfiles_Bmode, mask); % IQfiles_Bmode
        figure(fig1); title('Original Data - Without Mask');
        figure(fig2); title('Original Data - With Mask');
        % Save the correlation results
        save([processor.general_results_dir filesep 'correlation_no_mask'], 'correlation_no_mask');
        save([processor.general_results_dir filesep 'correlation_with_mask'], 'correlation_with_mask');

        load([[processor.general_results_dir filesep 'correlation_no_mask'] '.mat']);
        load([[processor.general_results_dir filesep 'correlation_with_mask'] '.mat']);
    
        %% Finding the indexes of the Frames that need to be removed based on the cross correlation
        fprintf('--- Finding Indices of Images to Remove based on Cross-Correlation --- \n\n') 
        indices_to_remove = calculate_indices_to_remove(processor.general_results_dir, 19,  params.acq.framerate);
         
        load([[processor.general_results_dir filesep 'indices_to_remove'] '.mat']);
        
        %% Clearing Bad Frames from data
        IQfiles_Bmode = dir(fullfile(bmode_folder, '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        IQfiles_PI = dir(fullfile(pi_folder, '*.mat')); IQfiles_PI = sortFilesByNumber(IQfiles_PI);
        remove_images(IQfiles_Bmode, indices_to_remove);
        remove_images(IQfiles_PI, indices_to_remove);

        fprintf('   - Re-batching cleaned data into uniform files...\n');
        IQfiles_Bmode = dir(fullfile(bmode_folder, 'BreathDeletion', '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        IQfiles_PI = dir(fullfile(pi_folder, 'BreathDeletion', '*.mat')); IQfiles_PI = sortFilesByNumber(IQfiles_PI);
        % Define the number of frames for each new file
        frames_per_new_file = 200;
        % Run the re-batching function for both B-mode and PI data
        rebatch_data_files(IQfiles_Bmode, frames_per_new_file)
        rebatch_data_files(IQfiles_PI, frames_per_new_file)

        % Correlation calculation
        fprintf('--- Correlation calculation for Data after Breath Deletion --- \n\n')
        IQfiles_Bmode = dir(fullfile(bmode_folder, 'BreathDeletion', '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        [correlation_after_BreathDeletion_no_mask, correlation_after_BreathDeletion_with_mask, fig1, fig2] = calculateAndPlotCorrelation(IQfiles_Bmode, mask); % IQfiles_Bmode
        figure(fig1); title('Data after Breath Deletion - Without Mask');
        figure(fig2); title('Data after Breath Deletion - With Mask');
        % Save the correlation results
        save([processor.general_results_dir filesep 'correlation_after_BreathDeletion_no_mask'], 'correlation_after_BreathDeletion_no_mask');
        save([processor.general_results_dir filesep 'correlation_after_BreathDeletion_with_mask'], 'correlation_after_BreathDeletion_with_mask');

        %% register images
        % =================================================================
        % PASS 1: COARSE CORRECTION (Rigid + Drift + New Folder)
        % =================================================================
        % enable_drift   = true  (Fix breathing drift)
        % overwrite_mode = false (Create 'registered' folder)
        % non_rigid_mode = false (Rigid only for stability)
        
        fprintf('\n>>> STARTING PASS 1: Coarse Correction (Original Data) <<<\n');
        
        %Build a ROI mask  
        fprintf('--- Build a ROI mask --- \n\n')
        IQfiles_Bmode = dir(fullfile(bmode_folder, 'BreathDeletion', '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        [mask2, roi2] = buildROIMask(IQfiles_Bmode);
        save([processor.general_results_dir filesep 'ROI2.mat'], 'roi2');
        save([processor.general_results_dir filesep 'Mask2.mat'], 'mask2');

        load([[processor.general_results_dir filesep 'Mask2'] '.mat']);
    
        IQfiles_Bmode = dir(fullfile(bmode_folder, '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        IQfiles_PI    = dir(fullfile(pi_folder, '*.mat'));    IQfiles_PI    = sortFilesByNumber(IQfiles_PI);

        % Correlation calculation with and without mask
        fprintf('--- Correlation calculation with and without mask --- \n\n')
        [correlation_for_register_no_mask, correlation_for_register_with_mask, fig1, fig2] = calculateAndPlotCorrelation(IQfiles_Bmode, mask2); % IQfiles_Bmode
        figure(fig1); title('Correlation for Register - Without Mask');
        figure(fig2); title('Correlation for Register - With Mask');
        % Save the correlation results
        save([processor.general_results_dir filesep 'correlation_for_register_no_mask'], 'correlation_for_register_no_mask');
        save([processor.general_results_dir filesep 'correlation_for_register_with_mask'], 'correlation_for_register_with_mask');
        
        load([[processor.general_results_dir filesep 'correlation_for_register_no_mask'] '.mat']);
        load([[processor.general_results_dir filesep 'correlation_for_register_with_mask'] '.mat']);

        % Execute Pass 1
        run_rigid_motion_correction_workflow(IQfiles_Bmode, IQfiles_PI, indices_to_remove, mask2, correlation_for_register_with_mask, true, false, false);
        
        % % =================================================================
        % % PASS 2: FINE TUNING (Soft Tissue + Overwrite)
        % % =================================================================
        % % enable_drift   = false (Do not re-correct drift, avoid artifacts)
        % % overwrite_mode = true  (Save inside the existing 'registered' folder)
        % % non_rigid_mode = true  (Enable B-Spline for soft kidney pulsing)
        % 
        % fprintf('\n>>> STARTING PASS 2: Fine Tuning (Soft Tissue Correction) <<<\n');
        % 
        % % Define the path to the data created in Pass 1
        % bmode_reg_folder = fullfile(bmode_folder, 'registered');
        % pi_reg_folder    = fullfile(pi_folder, 'registered');
        % 
        % IQfiles_Bmode_Reg = dir(fullfile(bmode_reg_folder, '*.mat')); IQfiles_Bmode_Reg = sortFilesByNumber(IQfiles_Bmode_Reg);
        % IQfiles_PI_Reg    = dir(fullfile(pi_reg_folder, '*.mat'));    IQfiles_PI_Reg    = sortFilesByNumber(IQfiles_PI_Reg);
        % 
        % if isempty(IQfiles_Bmode_Reg)
        %     error('Pass 1 did not seem to create files. Check path.');
        % end
        % 
        % % Execute Pass 2
        % run_rigid_motion_correction_workflow(IQfiles_Bmode_Reg, [], indices_to_remove, mask, correlation_no_mask, false, true, true);
        % 
        % fprintf('   - Double-Pass Motion analysis complete.\n');
        
        if kidney_workflow.use_registered_data
            IQfiles_Bmode = dir(fullfile(bmode_folder, 'registered', '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
            dataStruct = load(fullfile(IQfiles_Bmode(1).folder, IQfiles_Bmode(1).name));
            dataFields = fieldnames(dataStruct);
            ImageData_tot = dataStruct.(dataFields{1});
            params.expParams.size = size(ImageData_tot);
        end
        % implay(ImageData_tot); ceil(max(ImageData_tot,[],'all')/10)

        %% Clearing Bad Frames from registered data
        IQfiles_Bmode = dir(fullfile(bmode_folder, 'registered', '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        IQfiles_PI = dir(fullfile(pi_folder, 'registered', '*.mat')); IQfiles_PI = sortFilesByNumber(IQfiles_PI);
        remove_images(IQfiles_Bmode, indices_to_remove);
        remove_images(IQfiles_PI, indices_to_remove);

        fprintf('   - Re-batching cleaned data into uniform files...\n');
        IQfiles_Bmode = dir(fullfile(bmode_folder, 'registered', 'BreathDeletion', '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        IQfiles_PI = dir(fullfile(pi_folder, 'registered', 'BreathDeletion', '*.mat')); IQfiles_PI = sortFilesByNumber(IQfiles_PI);
        % Define the number of frames for each new file
        frames_per_new_file = 200;
        % Run the re-batching function for both B-mode and PI data
        rebatch_data_files(IQfiles_Bmode, frames_per_new_file)
        rebatch_data_files(IQfiles_PI, frames_per_new_file)
        
        load([[processor.general_results_dir filesep 'Mask2'] '.mat']); % mask

        % Correlation calculation
        fprintf('--- Correlation calculation for Data after Breath Deletion and Motion Correction --- \n\n')
        IQfiles_Bmode = dir(fullfile(bmode_folder, 'registered', 'BreathDeletion', '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
        [correlation_after_BreathDeletion_MotionCorrection_no_mask, correlation_after_BreathDeletion_MotionCorrection_with_mask, fig1, fig2] = calculateAndPlotCorrelation(IQfiles_Bmode, mask2); % IQfiles_Bmode
        figure(fig1); title('Data after Breath Deletion and Motion Correction - Without Mask');
        figure(fig2); title('Data after Breath Deletion and Motion Correction - With Mask');
        % Save the correlation results
        save([processor.general_results_dir filesep 'correlation_after_BreathDeletion_MotionCorrection_no_mask'], 'correlation_after_BreathDeletion_MotionCorrection_no_mask');
        save([processor.general_results_dir filesep 'correlation_after_BreathDeletion_MotionCorrection_with_mask'], 'correlation_after_BreathDeletion_MotionCorrection_with_mask');

        fprintf('   - Motion analysis complete. (Took %.1f min)\n', toc(t_motion)/60);
    else
        fprintf('   - Motion analysis and registration skipped by user.\n');
    end
    
    % --- 4c. Update Data Path for the Main Processor ---
    % Point the main processor to the corrected data folders.
    if kidney_workflow.use_registered_data
        if kidney_workflow.use_cleaned_data
            data_subfolder2 = fullfile(data_subfolder2, 'registered', 'BreathDeletion');
        else
            data_subfolder2 = fullfile(data_subfolder2, 'registered');
        end
    elseif kidney_workflow.use_cleaned_data
        data_subfolder2 = fullfile(data_subfolder2, 'BreathDeletion');
    end
    
    input_data_folder = fullfile(params.io.data_folder, params.expParams.bubbleType, params.io.data_subfolder);
    IQfiles_Bmode = dir(fullfile(input_data_folder, '*.mat')); IQfiles_Bmode = sortFilesByNumber(IQfiles_Bmode);
    dataStruct = load(fullfile(IQfiles_Bmode(1).folder, IQfiles_Bmode(1).name));
    dataFields = fieldnames(dataStruct);
    ImageData_tot = dataStruct.(dataFields{1});
    params.expParams.size = size(ImageData_tot);

    if kidney_workflow.use_registered_data
        % IMPORTANT: Update the params struct to point the processor to the new data location.
        params.io.data_subfolder = data_subfolder2;
        processor = ULM_Processor(params);
    end

    fprintf('   - Main processor will now use data from: %s\n', processor.data_dir);
    
else
    fprintf('\n>> STEP 4: Kidney pre-processing workflow skipped.\n');
    originalFolder = fullfile(params.io.data_folder, params.expParams.bubbleType, params.io.data_subfolder, 'original');
    if ~exist(originalFolder, 'dir')
        % Dividing the superframes into smaller batches of 1500 frames - for both Bmode and PI
        fprintf('--- Divides the superframes into batches --- \n\n')
        %split_ImageData_tot_Kidney(fullfile(params.io.data_folder, params.expParams.bubbleType), params.expParams.bubbleType, 1500)
    end
    processor = ULM_Processor(params);
end

%% ========================================================================
%  5. EXECUTE THE ULM PIPELINE
%  ========================================================================
fprintf('\n>> STEP 5: RUNNING THE MAIN PROCESSING PIPELINE\n');
tic;

% -------------------------------------------------------------
% Interactively define or load the Crop and Mask boundaries 
% based on the settings in info.txt / setDefaultParams
% -------------------------------------------------------------
processor.Step0_InitializeBoundaries();

% Define full file paths including the .mat extension
% (Required for the 'exist' function to work correctly)
file_vessel_path = fullfile(processor.general_results_dir, ['mean_bmode_vessel_', num2str(params.render.upsampling_factor),'.mat']);
file_tissue_path = fullfile(processor.general_results_dir, ['mean_bmode_tissue_', num2str(params.render.upsampling_factor),'.mat']);

% Define full file paths for the PNG images (Same name, different extension)
file_vessel_png = fullfile(processor.general_results_dir, ['mean_bmode_vessel_', num2str(params.render.upsampling_factor),'.png']);
file_tissue_png = fullfile(processor.general_results_dir, ['mean_bmode_tissue_', num2str(params.render.upsampling_factor),'.png']);

% --- Process Vessel ---
if ~exist(file_vessel_path, 'file')
    disp('Generating mean_bmode_vessel...');
    mean_bmode_vessel = generate_mean_bmode_image(processor.data_files, params, 0);
    save(file_vessel_path, 'mean_bmode_vessel');
else
    load(file_vessel_path, 'mean_bmode_vessel');
end

% Save Vessel PNG
% mat2gray scales the matrix values to [0,1] for proper image display
imwrite(mat2gray(mean_bmode_vessel), file_vessel_png); 

% --- Process Tissue ---
if ~exist(file_tissue_path, 'file')
    disp('Generating mean_bmode_tissue...');
    mean_bmode_tissue = generate_mean_bmode_image(processor.data_files, params, 1);
    save(file_tissue_path, 'mean_bmode_tissue');
else
    load(file_tissue_path, 'mean_bmode_tissue');
end

% Save Tissue PNG
imwrite(mat2gray(mean_bmode_tissue), file_tissue_png);
    
processor.run_localization_and_tracking_loop();
totalTime = toc;
fprintf('   - Pipeline execution and track aggregation finished in %.1f minutes.\n', totalTime / 60);

%% ========================================================================
%  6. RENDER FINAL MAPS WITH FILTERING & INPAINTING
%  ========================================================================
fprintf('\n>> STEP 6: RENDERING FINAL IMAGES\n');

MIN_TRACK_LENGTH = [5, 8, 10, 15, 20, 25, 30];  % 3, 5, 8, 10, 15, 18

for i = 1:length(MIN_TRACK_LENGTH)

    processor.run_streaming_reconstruction(MIN_TRACK_LENGTH(i));

% ========================================================================
%  7. SAVE AND DISPLAY RESULTS
%  ========================================================================
fprintf('\n>> STEP 7: SAVING AND VISUALIZING RESULTS\n');
try
    if i == 1
        processor.saveResults();
    end

catch ME
    fprintf('   - An error occurred while saving results: %s\n', ME.message);
end

try
    fprintf('   - Generating and displaying final images...\n');
    
    res = params.render.upsampling_factor;
    res_points       = [1, 3, 5];
    fontsize_points  = [6, 8, 10];
    linewidth_points = [2, 4, 5];
    
    res_clamped = max(min(res, res_points(end)), res_points(1));
    fontsize    = round(interp1(res_points, fontsize_points, res_clamped));
    linewidth   = round(interp1(res_points, linewidth_points, res_clamped), 1);

    % --- 7a. Create and Save Density Map ---
    fig_density = figure('Name', 'Super-Resolution Density Map', 'NumberTitle', 'off');
    density_processed = processor.densityMap.^(1/3);
    imshow(density_processed, [], 'colormap', hot);
    positive_pixels = density_processed(density_processed > 0);
    if ~isempty(positive_pixels)
        clim_upper = prctile(positive_pixels(:), 99.5);
        clim([0 clim_upper]);
    end
    title_str = sprintf('Super-Resolution Density Map (Min Length = %d)', MIN_TRACK_LENGTH(i));
    title(title_str, 'FontSize', 10);
    clb = colorbar;
    clb.Label.String = 'Counts^{1/3}';
    add_scale_bar(params, size(processor.densityMap), linewidth, fontsize);
    filename = sprintf('Density_Map_minLen_%d.png', MIN_TRACK_LENGTH(i));
    export_fig(fig_density, fullfile(processor.results_dir, filename), '-png', '-r300');

    % --- 7b. Create and Save Velocity Map with Gaussian Filter ---
    fig_vel_filt = figure('Name', 'Publication Velocity Map (Filtered)', 'NumberTitle', 'off');
    SR_vel_filtered = imgaussfilt(processor.velocityMap, 0.6);
    cm = colormap([0 0 0; jet(256)]);
    imshow(SR_vel_filtered, [], 'colormap', cm);
    positive_velocities = SR_vel_filtered(SR_vel_filtered > 0);
    if ~isempty(positive_velocities)
        vlim_upper = prctile(positive_velocities(:), 99.5);
        clim([0 vlim_upper]);
    end
    title_str = sprintf('Super-Resolution Velocity Map (Filtered) (Min Length = %d)', MIN_TRACK_LENGTH(i));
    title(title_str, 'FontSize', 10);
    clb = colorbar;
    clb.Label.String = 'Velocity [mm/s]';
    add_scale_bar(params, size(processor.velocityMap), linewidth, fontsize);
    filename = sprintf('Velocity_Map_Filtered_minLen_%d.png', MIN_TRACK_LENGTH(i));
    export_fig(fig_vel_filt, fullfile(processor.results_dir, filename), '-png', '-r300');

    % --- 7c. Create and Save Unfiltered Velocity Map ---
    fig_vel_unfilt = figure('Name', 'Publication Velocity Map (Unfiltered)', 'NumberTitle', 'off');
    imshow(processor.velocityMap, [], 'colormap', cm);
    positive_velocities = processor.velocityMap(processor.velocityMap > 0);
    if ~isempty(positive_velocities)
        vlim_upper = prctile(positive_velocities(:), 99.5);
        clim([0 vlim_upper]);
    end
    title_str = sprintf('Super-Resolution Velocity Map (Unfiltered) (Min Length = %d)', MIN_TRACK_LENGTH(i));
    title(title_str, 'FontSize', 10);
    clb = colorbar;
    clb.Label.String = 'Velocity [mm/s]';
    add_scale_bar(params, size(processor.velocityMap), linewidth, fontsize);
    filename = sprintf('Velocity_Map_Unfiltered_minLen_%d.png', MIN_TRACK_LENGTH(i));
    export_fig(fig_vel_unfilt, fullfile(processor.results_dir, filename), '-png', '-r300');

    % --- 7d. Kidney-Specific Visualization (Advanced Contrast & Colormap) ---
    if isKidneyExperiment && kidney_workflow.display_kidney_overlay
        fprintf('   - Generating overlay image...\n');

        load(file_vessel_path);
        mean_bmode = mean_bmode_vessel;

        if ndims(mean_bmode) == 3, bmode_gray = rgb2gray(mean_bmode); else, bmode_gray = mean_bmode; end
        bmode_adjusted = adapthisteq(mat2gray(bmode_gray));
        bmode_rgb = cat(3, bmode_adjusted, bmode_adjusted, bmode_adjusted);
        density_data = processor.densityMap.^(1/3);
        clim_max = prctile(density_data(:), 99.5);
        density_norm = mat2gray(density_data, [0, clim_max]);
        hot_colormap = hot(256);
        density_rgb = ind2rgb(gray2ind(density_norm, 256), hot_colormap);
        fig_overlay = figure('Name', 'Kidney Overlay (Publication Quality)', 'NumberTitle', 'off');
        imshow(density_rgb);
        hold on;
        h_bmode_image = imshow(bmode_rgb);
        set(h_bmode_image, 'AlphaData', 0.4);
        hold off;
        axis image;
        % Set color limits to match data range (fixes colorbar scaling)
        if exist('clim', 'builtin')
            clim([0, clim_max]);
        else
            caxis([0, clim_max]);
        end
        
        % Setup Colorbar
        clb = colorbar;
        colormap(clb, hot_colormap);
        clb.Limits = [0, clim_max];
        clb.Label.String = 'Counts^{1/3}';
        title_str = sprintf('ULM Density Map with B-Mode Overlay (Min Length = %d)', MIN_TRACK_LENGTH(i));
        title(title_str, 'FontSize', 10);

        filename = sprintf('Overlay_minLen_%d.png', MIN_TRACK_LENGTH(i));
        export_fig(fig_overlay, fullfile(processor.results_dir, filename), '-png', '-r300');

    elseif ~isKidneyExperiment
        fprintf('   - Generating overlay image...\n');

        load(file_tissue_path);
        mean_bmode = mean_bmode_tissue;
        
        % Prepare B-Mode background
        if ndims(mean_bmode) == 3
            bmode_gray = rgb2gray(mean_bmode); 
        else
            bmode_gray = mean_bmode; 
        end
        bmode_adjusted = adapthisteq(mat2gray(bmode_gray));
        bmode_rgb = cat(3, bmode_adjusted, bmode_adjusted, bmode_adjusted);
        
        % Prepare Density Map data
        density_data = processor.densityMap.^(1/3);
        clim_max = prctile(density_data(:), 99.5);
        
        % Create RGB density image
        density_norm = mat2gray(density_data, [0, clim_max]);
        hot_colormap = hot(256);
        density_rgb = ind2rgb(gray2ind(density_norm, 256), hot_colormap);
        
        % Plot figure
        fig_overlay = figure('Name', 'Overlay', 'NumberTitle', 'off');
        imshow(density_rgb);
        hold on;
        h_bmode_image = imshow(bmode_rgb);
        set(h_bmode_image, 'AlphaData', 0.4);
        hold off;
        axis image;
        
        % Set color limits to match data range (fixes colorbar scaling)
        if exist('clim', 'builtin')
            clim([0, clim_max]);
        else
            caxis([0, clim_max]);
        end
        
        % Setup Colorbar
        clb = colorbar;
        colormap(clb, hot_colormap);
        clb.Limits = [0, clim_max];
        clb.Label.String = 'Counts^{1/3}';
        
        % Title and Export
        title_str = sprintf('ULM Density Map with B-Mode Overlay (Min Length = %d)', MIN_TRACK_LENGTH(i));
        title(title_str, 'FontSize', 10);
        
        filename = sprintf('Overlay_minLen_%d.png', MIN_TRACK_LENGTH(i));
        export_fig(fig_overlay, fullfile(processor.results_dir, filename), '-png', '-r300');
    end

    % --- 7e. Create and Save Combined Velocity-Intensity Map (HSV) ---
    fprintf('   - Generating combined Velocity (Color) + Density (Brightness) map...\n');

    fig_combined = figure('Name', 'Combined Vel-Dens Map', 'NumberTitle', 'off');
    % 1. Prepare Velocity Data (Determines the Color/Hue)
    vel_data = processor.velocityMap; 
    pos_vel = vel_data(vel_data > 0);
    if ~isempty(pos_vel)
        clim_vel = prctile(pos_vel(:), 99.5);
    else
        clim_vel = 1;
    end
    vel_norm = mat2gray(vel_data, [0, clim_vel]);

    % 2. Prepare Density Data (Determines the Brightness/Value)
    % This is the first compression (keeps background 0)
    dens_data = processor.densityMap.^(1/3); 
    pos_dens = dens_data(dens_data > 0);
    if ~isempty(pos_dens)
        % Reduced percentile slightly to 99 to saturate the very brightest pixels
        % and allow the dimmer ones to shine more.
        clim_dens = prctile(pos_dens(:), 99.0); 
    else
        clim_dens = 1;
    end
    % Normalize density to [0, 1]
    dens_norm = mat2gray(dens_data, [0, clim_dens]);
    % --- BRIGHTNESS ADJUSTMENT SECTION ---
    % If the image is too dark, lower this gamma value.
    % 1.0 = Linear (Original)
    % 0.5 = Brighter (Square root)
    % 0.3 = Very Bright
    gamma_correction = 0.7; 
    dens_norm = dens_norm .^ gamma_correction;

    % 3. Construct the Image using HSV space
    n_colors = 256;
    cmap = jet(n_colors); 
    vel_idx = gray2ind(vel_norm, n_colors);
    rgb_base = ind2rgb(vel_idx, cmap);
    % Convert RGB to HSV 
    hsv_img = rgb2hsv(rgb_base);
    % Replace the 'Value' channel (brightness) with our boosted Density map
    hsv_img(:,:,3) = dens_norm; 
    % Convert back to RGB for display
    final_combined_img = hsv2rgb(hsv_img);
    % Ensure background is absolute black
    mask_bg = dens_data == 0;
    r_ch = final_combined_img(:,:,1); r_ch(mask_bg) = 0;
    g_ch = final_combined_img(:,:,2); g_ch(mask_bg) = 0;
    b_ch = final_combined_img(:,:,3); b_ch(mask_bg) = 0;
    final_combined_img = cat(3, r_ch, g_ch, b_ch);

    % 4. Display and Save
    imshow(final_combined_img);
    title_str = sprintf('Combined Velocity & Density (Min Length = %d)', MIN_TRACK_LENGTH(i));
    title(title_str, 'FontSize', 10);
    colormap(gca, cmap);
    clim([0 clim_vel]);
    clb = colorbar;
    clb.Label.String = 'Velocity [mm/s]';
    clb.Color = [0 0 0]; % Explicitly set text to Black
    add_scale_bar(params, size(processor.velocityMap), linewidth, fontsize);
    filename = sprintf('Combined_VelDens_Map_minLen_%d.png', MIN_TRACK_LENGTH(i));
    export_fig(fig_combined, fullfile(processor.results_dir, filename), '-png', '-r300');
    
    fprintf('   - Figures generated and saved successfully.\n');
catch ME
     fprintf('\n   - Could not display final images. Error: %s ---\n', ME.message);
end
end

fprintf('\n=======================================================\n');
fprintf('               ULM Analysis Complete\n');
fprintf('=======================================================\n');


%% ========================================================================
%  Local Helper Functions
%  ========================================================================

function add_scale_bar(params, mapSize, linewidth, fontsize)
    % Local helper function to add a 1 mm scale bar to the current axes.
    hold on;
    
    pixel_size_mm = params.expParams.fovX / mapSize(2);
    scalebar_length_pixels = 1 / pixel_size_mm; % Length of 1 mm in pixels
    
    % Position the scale bar at the bottom-left corner
    x_pos = 0.05 * mapSize(2); % 5% from the left edge
    y_pos = 0.95 * mapSize(1); % 5% from the bottom edge
    
    plot([x_pos, x_pos + scalebar_length_pixels], [y_pos, y_pos], 'w-', 'LineWidth', linewidth);
    text(x_pos + scalebar_length_pixels/2, y_pos - fontsize*linewidth, '1 mm', ...
    'Color', 'w', 'FontSize', fontsize, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    hold off;
end

function str = B_SWITCH(bool_val)
    % Helper to convert boolean to 'ON'/'OFF' string for display
    if bool_val, str = 'ON'; else, str = 'OFF'; end
end

function mean_bmode = generate_mean_bmode_image(data_files, params, vessel_or_tissue)
    % Generates a mean B-mode image from all buffers for overlay purposes.
    % Upsample to match the SR image dimensions
    % vessel_or_tissue: 0 - vessel, 1 - tissue
    sRes_dims = params.expParams.size(1:2) * params.render.upsampling_factor;
    
    % Check if files exist
    if isempty(data_files)
        warning('No data files found to generate mean B-mode image.');
        mean_bmode = zeros(sRes_dims);
        return;
    end
    
    num_files = length(data_files);
    accumulated_mean = 0;
    
    % Loop through all files in the directory/list
    for i = 1:num_files
        % Load the current buffer
        dataStruct = load(fullfile(data_files(i).folder, data_files(i).name));
        dataFields = fieldnames(dataStruct); 
        rawData = dataStruct.(dataFields{1});
        
        % Filter data based on vessel or tissue selection
        if vessel_or_tissue == 0
            filteredData = SVD_filter(rawData, params.filter.svd_cutoff);
        elseif vessel_or_tissue == 1
            c1 = params.filter.svd_cutoff;
            filteredData = SVD_filter(rawData, [0 c1(1)]);
        else
            filteredData = SVD_filter(rawData, params.filter.svd_cutoff);
        end

        % Apply Butterworth filter if enabled
        if isfield(params.filter, 'enable_butterworth') && params.filter.enable_butterworth
            filteredData = Butterworth_bandpass_filter(filteredData, params.filter.butter_cutoff, params.acq.framerate, params.filter.butter_order);
        end
        
        % Calculate the mean of the absolute values for the current file and accumulate
        current_mean = double(mean(abs(filteredData), 3));
        accumulated_mean = accumulated_mean + current_mean;
    end
    
    % Average the accumulated means across all files
    final_mean = accumulated_mean / num_files;
    
    % Apply the power factor
    final_mean = final_mean.^(0.5);
    
    % Resize to target dimensions
    mean_bmode = imresize(final_mean, sRes_dims);
end

function [mask, roi] = buildROIMask(IQfiles)
    % Function to build an ROI mask from mean image of ImageData_tot
    % Input:
    %   IQfiles - Struct array containing file information for loading images
    %   savingpath - Path and base name for saving the ROI and mask

    % Load the ImageData_tot variable from the specified file
    dataStruct = load(fullfile(IQfiles(1).folder, IQfiles(1).name));
    dataFields = fieldnames(dataStruct); rawData = dataStruct.(dataFields{1});
    tmp = SVD_filter(rawData, [0 80]);
    tmp = abs(tmp);

    % Extract and process the specified frame
    frame = mean(tmp,3);

    % Display the frame with contrast adjustment
    figure;
    imshow(frame.^(1/2), []);
    title('Mean Image');

    % Create the ROI using Assisted Freehand tool
    roi = images.roi.AssistedFreehand;
    draw(roi);

    % Create a binary mask from the ROI
    mask = createMask(roi);

    % Close the figure
    close(gcf);
end

function [cor_no_mask_Normalized, cor_with_mask_Normalized, fig1, fig2] = calculateAndPlotCorrelation(IQfiles_Bmode, mask)
% CALCULATEANDPLOTCORRELATION Calculates cross-correlation of consecutive B-mode frames.
%
% PURPOSE:
%   This function computes the frame-to-frame cross-correlation to assess
%   tissue stability or motion. It operates in two modes:
%   1. Global Correlation: Correlation of the entire frame.
%   2. Masked Correlation (Optional): Correlation within a specific ROI.
%
% INPUTS:
%   IQfiles_Bmode : (struct array) Struct containing file information for 
%                   B-mode images (standard output from dir() or similar).
%   mask          : (optional, logical/double) A binary mask to apply during 
%                   correlation calculation. If omitted, empty, or if dimensions
%                   do not match the image, masked analysis is skipped.
%
% OUTPUTS:
%   cor_no_mask_Normalized   : (double vector) Normalized correlation coefficients
%                              calculated on the full image frames.
%   cor_with_mask_Normalized : (double vector) Normalized correlation coefficients
%                              calculated within the masked region. Returns [] 
%                              if no mask is used.
%   fig1                     : (figure handle) Handle to the plot of global correlation.
%   fig2                     : (figure handle) Handle to the plot of masked correlation.
%                              Returns [] if no mask is used.
%
% AUTHOR: M.Sc. Thesis Assistant
% DATE: January 2026

    %% 1. Input Handling and Setup
    % Check if mask is provided and not empty
    if nargin < 2 || isempty(mask)
        useMask = false;
        fprintf('Info: No mask provided. Skipping masked correlation analysis.\n');
    else
        useMask = true;
    end

    Nbuffers = numel(IQfiles_Bmode);

    % Initialize output arrays
    cor_no_mask = [];
    cor_no_mask_Normalized = [];
    
    % Initialize masked outputs as empty
    cor = []; 
    cor_with_mask_Normalized = [];
    fig2 = []; 

    % Buffers for normalization factors
    max_cor_buffer = zeros(1, Nbuffers);
    max_cor_no_mask_buffer = zeros(1, Nbuffers);
    
    %% 2. Main Processing Loop
    for hhh = 1:Nbuffers
        fprintf('Processing block %d/%d...\n', hhh, Nbuffers);

        % Load the B-mode image data
        dataStruct = load(fullfile(IQfiles_Bmode(hhh).folder, IQfiles_Bmode(hhh).name));
        dataFields = fieldnames(dataStruct); 
        tmp = abs(dataStruct.(dataFields{1})); % Convert to magnitude if complex

        [imgH, imgW, nFrames] = size(tmp);
        
        % --- VALIDATION: Check Mask Dimensions ---
        % We check this inside the loop because we need the image dimensions.
        % If dimensions mismatch, we disable the mask for the rest of the run.
        if useMask
            [maskH, maskW] = size(mask);
            if imgH ~= maskH || imgW ~= maskW
                warning('Mask dimensions (%dx%d) do not match Image dimensions (%dx%d). Resizing mask to fit.', ...
                        maskH, maskW, imgH, imgW);
                
                % Resize the mask to match the image dimensions
                % 'nearest' is used to avoid interpolation artifacts (gray values) in logical masks
                mask = imresize(mask, [imgH, imgW], 'nearest');
            end
        end

        % Pre-allocate temporary vectors for this buffer
        cor1 = zeros(1, nFrames);
        if useMask
            cor2 = zeros(1, nFrames);
        end

        % --- Frame-by-Frame Correlation ---
        for frame = 1:nFrames
            if frame < nFrames
                % 1. Global Correlation
                cor1(frame) = corr2(tmp(:, :, frame), tmp(:, :, frame+1));
                
                % 2. Masked Correlation (only if mask is valid)
                if useMask
                    cor2(frame) = corr2(tmp(:, :, frame) .* mask, tmp(:, :, frame+1) .* mask);
                end
                
            elseif frame == nFrames
                % Handle the last frame by copying the previous correlation value
                if frame > 1
                    cor1(frame) = cor1(frame - 1);
                    if useMask
                        cor2(frame) = cor2(frame - 1);
                    end
                else
                    % Edge case: single frame buffer
                    cor1(frame) = 1; 
                    if useMask, cor2(frame) = 1; end
                end
            end
        end
        
        % Accumulate raw results
        cor_no_mask = [cor_no_mask, cor1];
        max_cor_no_mask_buffer(hhh) = max(cor1);

        if useMask
            cor = [cor, cor2];
            max_cor_buffer(hhh) = max(cor2);
        end

        % --- Normalization ---
        if max_cor_no_mask_buffer(hhh) ~= 0
            cor_no_mask_Normalized = [cor_no_mask_Normalized, cor1 ./ max_cor_no_mask_buffer(hhh)];
        else
            cor_no_mask_Normalized = [cor_no_mask_Normalized, cor1];
        end

        if useMask
            if max_cor_buffer(hhh) ~= 0
                cor_with_mask_Normalized = [cor_with_mask_Normalized, cor2 ./ max_cor_buffer(hhh)];
            else
                cor_with_mask_Normalized = [cor_with_mask_Normalized, cor2];
            end
        end

    end

    %% 3. Visualization
    % Plot 1: Normalized Correlation - Without Mask
    fig1 = figure('Name', 'Global Cross-Correlation', 'Color', 'w');
    plot(cor_no_mask_Normalized, 'LineWidth', 1.5);
    xlabel('Frame Number');
    ylabel('Normalized Cross Correlation');
    title('Normalized Cross Correlation - Full Field of View');
    grid on;

    % Plot 2: Normalized Correlation - With Mask (Conditional)
    if useMask && ~isempty(cor_with_mask_Normalized)
        fig2 = figure('Name', 'Masked Cross-Correlation', 'Color', 'w');
        plot(cor_with_mask_Normalized, 'LineWidth', 1.5, 'Color', 'r');
        xlabel('Frame Number');
        ylabel('Normalized Cross Correlation');
        title('Normalized Cross Correlation - Masked ROI');
        grid on;
    end
end

function remove_images(IQfiles, indices_to_remove)
    % This function removes images from a list of batches (IQfiles) whose indices
    % are in the list 'indices_to_remove'. The resulting data is saved with
    % the suffix '_BreathDeletion' in a subfolder.
    %
    % Parameters:
    %   IQfiles           - List of files (struct array from dir command)
    %   indices_to_remove - List of global indices to be removed

    Nbuffers = numel(IQfiles);
    if Nbuffers == 0
        return
    end

    % --- Create the output folder ---
    deletion_folder = fullfile(IQfiles(1).folder, 'BreathDeletion');
    if ~exist(deletion_folder, 'dir'), mkdir(deletion_folder); end
    
    % --- Initialize a running offset for the global frame count ---
    global_frame_offset = 0;
    
    for hhh = 1:Nbuffers
        fprintf('Processing batch %d/%d...\n', hhh, Nbuffers);
        
        % --- Load Data ---
        filePath = fullfile(IQfiles(hhh).folder, IQfiles(hhh).name);
        dataStruct = load(filePath);
        dataFields = fieldnames(dataStruct);
        
        % Assuming the data of interest is the first field in the struct
        tmp_data = dataStruct.(dataFields{1});
        
        % Get the number of frames BEFORE removal (Critical for offset calculation)
        num_frames_in_this_batch = size(tmp_data, 3);
        
        % --- Remove Frames ---
        % Pass the current offset to the helper function
        tmp_data_cleaned = remove_frames(tmp_data, indices_to_remove, global_frame_offset);
        
        % --- Save Modified Data ---
        if ~isempty(tmp_data_cleaned)
            original_name = IQfiles(hhh).name;
            new_name = strrep(original_name, '.mat', '_BreathDeletion.mat');
            
            % Using a fixed variable name for saving to maintain consistency
            ImageData_tot = tmp_data_cleaned; 
            
            save(fullfile(deletion_folder, new_name), 'ImageData_tot', '-v7.3');
        end
        
        % Update the offset for the next iteration
        global_frame_offset = global_frame_offset + num_frames_in_this_batch;
    end
    
    fprintf('All batches in the list have been processed successfully.\n');
end

function tmp_cleaned = remove_frames(tmp, indices_to_remove, frame_offset)
    % Helper function to remove frames from a batch of images using a running offset.
    % Parameters:
    %   tmp               - 3D matrix containing the batch of images (HxWxN)
    %   indices_to_remove - List of global indices to be removed
    %   frame_offset      - The number of frames that came before this batch
    %
    % Returns:
    %   tmp_cleaned - 3D matrix with the specified frames removed

    num_frames_in_batch = size(tmp, 3);
    frames_to_keep_mask = true(1, num_frames_in_batch); % Create a logical mask

    % Loop through each frame in the current batch
    for frame_local_idx = 1:num_frames_in_batch
        % --- Correct global index calculation ---
        global_frame_idx = frame_offset + frame_local_idx;

        % Check if this frame's global index should be removed
        if ismember(global_frame_idx, indices_to_remove)
            frames_to_keep_mask(frame_local_idx) = false; % Mark this frame for removal
        end
    end

    % Use the logical mask to select only the frames we want to keep
    tmp_cleaned = tmp(:, :, frames_to_keep_mask);
end

function rebatch_data_files(IQfiles, framesPerFile)
%% rebatch_data_files - Consolidates data from a file list into uniform batches.
%
% =========================================================================
%                               DESCRIPTION
% =========================================================================
% This function processes a specific list of .mat files (IQfiles), consolidates
% their frame data, and rewrites them into new files containing a fixed
% number of frames per file.
%
% KEY BEHAVIORS:
%   1. TRUSTED ORDER: It processes files in the exact order provided in the
%      'IQfiles' input struct. It does NOT perform any internal sorting.
%   2. VARIABLE STANDARDIZATION: The output .mat files will always contain
%      a variable named 'ImageData_tot', regardless of the variable name
%      in the source files.
%   3. MEMORY EFFICIENCY: Uses a streaming buffer approach to handle large
%      datasets without loading everything into RAM.
%   4. CLEANUP: After successfully processing, the original files listed in
%      'IQfiles' are deleted from the disk.
%
% =========================================================================
%                                 INPUTS
% =========================================================================
%   IQfiles       - (struct array) A standard MATLAB file structure (result 
%                   of the dir() command) containing the list of files to 
%                   process, already sorted in the desired order.
%   framesPerFile - (double) The target number of frames for each new 
%                   output file (e.g., 1000).
%
% =========================================================================

    % --- Validation ---
    if isempty(IQfiles)
        fprintf('WARNING: The input file list (IQfiles) is empty. Exiting.\n');
        return;
    end

    % Determine the working folder from the first file in the list
    sourceFolder = IQfiles(1).folder;
    
    fprintf('--- Starting re-batching process ---\n');
    fprintf('    Source folder: %s\n', sourceFolder);
    fprintf('    Input files: %d\n', length(IQfiles));
    fprintf('    Target frames per file: %d\n', framesPerFile);

    % --- Initialization ---
    outputBuffer = [];
    newFileCounter = 1;
    processedFiles = {}; % Track files that were actually loaded
    
    % --- Step 1: Process files in a streaming buffer ---
    for i = 1:length(IQfiles)
        filePath = fullfile(IQfiles(i).folder, IQfiles(i).name);
        
        % Load the data
        try
            dataStruct = load(filePath);
        catch ME
            fprintf('    WARNING: Failed to load %s. Skipping. Error: %s\n', IQfiles(i).name, ME.message);
            continue; 
        end
        
        % Dynamic variable detection
        varNames = fieldnames(dataStruct);
        if isempty(varNames)
            fprintf('    WARNING: File %s contains no variables. Skipping.\n', IQfiles(i).name);
            continue;
        end
        
        % Assume the first variable is the relevant data
        dataVarName = varNames{1}; 
        framesFromThisFile = dataStruct.(dataVarName);
        
        % Validation: Must be non-empty and 3D
        if isempty(framesFromThisFile) || ndims(framesFromThisFile) < 3
             fprintf('    WARNING: Data in %s is empty or not a 3D matrix. Skipping.\n', IQfiles(i).name);
             continue;
        end
        
        % Mark this file as processed (queued for deletion later)
        processedFiles{end+1} = filePath;
        
        % Append frames to the buffer
        if isempty(outputBuffer)
            outputBuffer = framesFromThisFile;
        else
            outputBuffer = cat(3, outputBuffer, framesFromThisFile);
        end
        
        % --- Step 2: Save chunks when buffer is full ---
        while size(outputBuffer, 3) >= framesPerFile
            % Slice the exact amount of frames
            currentChunk = outputBuffer(:, :, 1:framesPerFile);
            
            % Remove those frames from the buffer
            outputBuffer = outputBuffer(:, :, framesPerFile+1:end);
            
            % Generate new file name
            newFileName = sprintf('Rebatched_Data_%03d.mat', newFileCounter);
            newFilePath = fullfile(sourceFolder, newFileName);
            
            % Save as 'ImageData_tot' for consistency with other tools
            ImageData_tot = currentChunk; 
            save(newFilePath, 'ImageData_tot', '-v7.3');
            
            fprintf('        --> Saved %s (%d frames)\n', newFileName, size(ImageData_tot, 3));
            newFileCounter = newFileCounter + 1;
        end
    end
    
    % --- Step 3: Save remaining frames ---
    if ~isempty(outputBuffer)
        newFileName = sprintf('Rebatched_Data_%03d.mat', newFileCounter);
        newFilePath = fullfile(sourceFolder, newFileName);
        
        ImageData_tot = outputBuffer;
        save(newFilePath, 'ImageData_tot', '-v7.3');
        fprintf('        --> Saved final file %s (%d frames)\n', newFileName, size(ImageData_tot, 3));
    end
    
    % --- Step 4: Delete original files ---
    % Only delete files that were successfully added to the buffer
    if ~isempty(processedFiles)
        fprintf('    Deleting %d original files...\n', length(processedFiles));
        for k = 1:length(processedFiles)
            delete(processedFiles{k});
        end
    end
    
    fprintf('--- Re-batching complete. ---\n\n');
end

function indices_to_remove = calculate_indices_to_remove(savingpath, expected_bpm, framerate)
%% calculate_indices_to_remove - Robust Two-Pass Respiratory Motion Detection
%
% =========================================================================
%                               DESCRIPTION
% =========================================================================
% This function automatically identifies and flags frames containing motion 
% artifacts caused by breathing within an ultrasound cross-correlation signal.
% It employs a robust, multi-stage pipeline designed to handle variations
% in breathing rate and amplitude.
%
% =========================================================================
%                           ALGORITHM WORKFLOW
% =========================================================================
% 1.  DATA PREPARATION:
%     - Loads the 1D cross-correlation signal.
%     - Estimates the dominant breathing frequency using FFT to adaptively
%       set parameters for the rest of the pipeline.
%     - Applies symmetric padding to the signal edges to prevent artifacts
%       during filtering.
%
% 2.  SIGNAL TRANSFORMATION (Pan-Tompkins Inspired):
%     - A band-pass filter isolates the primary breathing signal.
%     - The filtered signal is processed to create a smooth "energy envelope"
%       (MWI signal), where each breath event corresponds to a distinct peak.
%       This involves: Derivative -> Squaring -> Moving Window Integration (MWI).
%
% 3.  PEAK DETECTION (Two-Pass Strategy):
%     - PASS 1: Detects strong, high-confidence breath events by finding
%       prominent peaks in the MWI energy signal. This catches the easy cases.
%     - PASS 2: "Rescues" missed breaths (e.g., those that are weak, shallow,
%       or occur close together) by searching for significant local minimums 
%       (valleys) directly in the filtered signal, within the gaps left by Pass 1.
%
% 4.  LOCATION REFINEMENT:
%     - For every detected breath location (from both passes), a local search
%       is performed on the ORIGINAL raw signal. This "snaps" the marker to
%       the precise frame of maximum decorrelation (the bottom of the dip),
%       ensuring high positional accuracy.
%
% 5.  REGION DEFINITION:
%     - Around each refined breath location, a dynamic threshold is used on the
%       MWI energy signal to determine the precise start and end frames of the
%       motion artifact. This defines the full region to be removed.
%
% 6.  OUTPUT & VISUALIZATION:
%     - All frame indices within the defined regions are collected for removal.
%     - A detailed 3-panel figure is generated to visualize the raw signal,
%       the MWI energy with detected peaks, the filtered signal, and the
%       final removal decisions.
%
% =========================================================================
%                                 USAGE
% =========================================================================
%   indices_to_remove = calculate_indices_to_remove(savingpath);
%   indices_to_remove = calculate_indices_to_remove(savingpath, expected_bpm);
%
% =========================================================================
%                                 INPUTS
% =========================================================================
%   savingpath   - (string) Base path for loading 'correlation_no_mask.mat'.
%   expected_bpm - (double, optional) Estimated breathing rate in breaths 
%                  per minute to guide the initial frequency search.
%
% =========================================================================
%                                OUTPUTS
% =========================================================================
%   indices_to_remove - (vector) A sorted, unique list of frame indices
%                       identified for removal from the original data.
%

    %% 0) ======================= TUNABLE PARAMETERS =======================
    % This section contains the key parameters that control the sensitivity 
    % and behavior of the detection algorithm. Adjust these values to fine-tune
    % the performance for different datasets.
    
    % --- Pass 1: MWI Peak Detection ---
    p1_peak_height_coeff = 0.4;  % TUNABLE: Multiplier for the standard deviation in Pass 1 threshold.
                                 % > HIGHER: More conservative, detects only very prominent peaks.
                                 % > LOWER: More sensitive, may detect noise as peaks.
    p1_min_dist_coeff = 0.4;     % TUNABLE: Multiplier for the cycle duration to set minimum peak distance in Pass 1.
                                 % > HIGHER: Enforces a larger gap between detected breaths.
                                 % > LOWER: Allows detection of more closely spaced breaths.

    % --- Pass 2: Filtered Signal Minimum Detection ---
    p2_peak_height_coeff = 1.3;  % TUNABLE: Multiplier for the standard deviation of the filtered signal.
                                 % Defines how "deep" a valley must be to be considered a breath.
                                 % > HIGHER: More conservative, requires a very clear dip.
                                 % > LOWER: More sensitive, will "rescue" weaker breaths.
    p2_min_dist_coeff = 0.5;    % TUNABLE: Multiplier for minimum peak distance in Pass 2. 
                                 % Should be small to allow detection of closely spaced breaths.

    % --- Refinement Step ---
    refinement_window_seconds = 0.3; % TUNABLE: Search radius in seconds (+/-) around a detected peak
                                      % to find the true minimum in the raw signal.
    
    % --- Region Definition Step ---
    region_width_coeff_start = -0.5; % -0.3 
    region_width_coeff_end = +0.3; % +1.2
        % TUNABLE: Multiplier for standard deviation to set the start/end
        % of the removal region. A more negative value (e.g., -0.5)
        % will result in a WIDER removal region around each breath.
        % A value closer to zero (e.g., -0.1) results in a NARROWER region.
    

    %% 1) ======================= LOAD AND PREPARE DATA =======================
    % Load the cross-correlation data
    data  = load([savingpath filesep 'correlation_no_mask'], 'correlation_no_mask');
    cor   = data.correlation_no_mask(:);
    N     = numel(cor);
    fs    = framerate;
    
    % Detrend the signal by removing the mean (DC component)
    cor_d = cor - mean(cor);
    
    %% 2) ======================= ESTIMATE BREATHING RATE & PAD SIGNAL =======================
    % Use Fast Fourier Transform (FFT) to find the dominant frequency in the signal
    Y       = fft(cor_d .* hamming(N)); % Apply Hamming window to reduce spectral leakage
    P2      = abs(Y/N);
    P1      = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    f       = fs * (0:floor(N/2)) / N;
    
    % Define the frequency band to search for breathing
    if nargin>=2 && ~isempty(expected_bpm)
        hz0 = expected_bpm / 60;
        bw  = max(0.2, min(0.5 * hz0, 2));
        band = [max(hz0 - bw, 0.1), min(hz0 + 15 * bw, 6)];
    else
        band = [0.1, 6]; % Default physiological range for rats
    end

    % Find the frequency with the maximum power within the band
    idx      = f >= band(1) & f <= band(2);
    [~, bi]  = max(P1(idx));
    fband    = f(idx);
    freq_dom = fband(bi);
    cycleSec = 1 / freq_dom; % Dominant cycle period in seconds
    fprintf('Detected dominant breathing rate: %.2f Hz (%.1f bpm), cycle duration: ~%.2f s\n', ...
            freq_dom, freq_dom*60, cycleSec);
    
    % Apply symmetric padding to signal edges to prevent artifacts from filtering
    pad_len = round(0.9*cycleSec * fs);
    pad_len = min(pad_len, floor(N/2) - 1);
    if pad_len < 1, pad_len = 1; end
    
    % Pad both the detrended and the original raw signal
    cor_d_padded = [flip(cor_d(1:pad_len)); cor_d; flip(cor_d(end-pad_len+1:end))];
    cor_padded = [flip(cor(1:pad_len)); cor; flip(cor(end-pad_len+1:end))]; 
    N_padded = numel(cor_d_padded);

    %% 2b) ======================= DEBUG: VISUALIZE FFT SPECTRUM & SELECT PEAK =======================
    % This section plots the power spectrum and allows the user to
    % INTERACTIVELY SELECT which of the top 3 distinct peaks to use as the
    % dominant frequency for the rest of the algorithm.
    
    % --- TUNABLE PARAMETER for this debug section ---
    min_freq_separation_hz = 0.1; % Hz. Defines the minimum required frequency separation 
                                  % between two consecutively chosen peaks.
                                  
    fig_debug_fft = figure('Name', 'FFT Spectrum Debug - CHOOSE A PEAK', 'NumberTitle', 'off');
    semilogy(f, P1, '-o', 'MarkerSize', 4, 'DisplayName', 'Power Spectrum (P1)');
    hold on;
    
    % Highlight the search band
    yl = ylim;
    patch([band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], 'y', ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Search Band');
    
    % --- Find and prepare top 3 distinct peaks ---
    P1_in_band = P1(idx);
    f_in_band = f(idx);
    [sorted_powers, sort_indices] = sort(P1_in_band, 'descend');
    
    final_top_indices = []; % This will store the valid indices (within f_in_band)
    
    for i = 1:numel(sort_indices)
        current_local_idx = sort_indices(i);
        current_freq = f_in_band(current_local_idx);
        is_distinct = true;
        
        for j = 1:numel(final_top_indices)
            existing_freq = f_in_band(final_top_indices(j));
            if abs(current_freq - existing_freq) < min_freq_separation_hz
                is_distinct = false;
                break;
            end
        end
        
        if is_distinct
            final_top_indices(end+1) = current_local_idx;
        end
        
        if numel(final_top_indices) >= 3
            break;
        end
    end
    
    % Plot and label the final, distinct peaks
    if ~isempty(final_top_indices)
        for i = 1:numel(final_top_indices)
            local_idx = final_top_indices(i);
            peak_power = P1_in_band(local_idx);
            peak_freq = f_in_band(local_idx);
            plot(peak_freq, peak_power, 'ro', 'MarkerSize', 15, 'LineWidth', 1.5, 'HandleVisibility', 'off');
            text(peak_freq, peak_power, sprintf('  %d', i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
        end
    else
        error('No significant peaks found in the specified frequency band.');
    end
    
    grid on;
    title(sprintf('FFT Power Spectrum (Top %d Distinct Peaks)', numel(final_top_indices)));
    xlabel('Frequency (Hz)');
    ylabel('Power (log scale)');
    legend('show', 'Location', 'northeast');
    xlim([0, band(2) * 1.5]);
    
    drawnow;
    
    % --- NEW: Interactive User Input ---
    % Prompt the user to select one of the labeled peaks.
    prompt_str = sprintf('\nEnter the number of the peak to use (1-%d): ', numel(final_top_indices));
    user_choice = [];
    
    while isempty(user_choice)
        choice = input(prompt_str);
        % Validate the input
        if isnumeric(choice) && isscalar(choice) && (choice >= 1 && choice <= numel(final_top_indices))
            user_choice = round(choice); % Ensure it's an integer
        else
            fprintf('Invalid input. Please enter a number between 1 and %d.\n', numel(final_top_indices));
        end
    end
    
    % --- Update algorithm parameters based on user's choice ---
    chosen_local_idx = final_top_indices(user_choice);
    freq_dom = f_in_band(chosen_local_idx); 
    cycleSec = 1 / freq_dom;              
    
    fprintf('--> User selected peak #%d. Using %.2f Hz (%.1f bpm) for analysis.\n', ...
            user_choice, freq_dom, freq_dom*60);
            
    % Close the debug figure and continue
    if ishandle(fig_debug_fft)
        close(fig_debug_fft);
    end
    % ======================= END OF DEBUG SECTION =======================

    %% 3) ======================= BAND-PASS FILTERING =======================
    % Design a Butterworth filter to isolate the breathing signal
    [b,a] = butter(3, band/(fs/2));
    % Apply zero-phase filtering to avoid shifting the signal in time
    cor_f_padded = filtfilt(b, a, cor_d_padded);
    
    %% 4) ======================= ENERGY ENVELOPE GENERATION =======================
    % This process transforms the filtered signal into a smooth energy envelope (MWI)
    
    % a. Derivative: Emphasizes sharp changes (start of breath)
    der_kernel = [1 2 0 -2 -1] / 8 * fs;
    cor_der_padded = filtfilt(der_kernel, 1, cor_f_padded);
    
    % b. Squaring: Makes all values positive and amplifies high-energy points
    cor_sq_padded  = cor_der_padded .^ 2;
    
    % c. Moving Window Integration (MWI): Smooths the signal into a clean envelope
    mwi_win = round(cycleSec * fs / 3);
    mwi_padded = movmean(cor_sq_padded, mwi_win);
    
    %% 5a) ======================= PASS 1: DETECT PRIMARY BREATHS =======================
    % Find prominent peaks in the MWI energy signal using a conservative threshold.
    thr_global = mean(mwi_padded) + p1_peak_height_coeff * std(mwi_padded);
    minDist_pass1 = round(cycleSec * fs * p1_min_dist_coeff);
    [~, primary_locs] = findpeaks(mwi_padded, 'MinPeakDistance', minDist_pass1, 'MinPeakHeight', thr_global);
    fprintf('Pass 1 (MWI): Found %d primary breath events.\n', numel(primary_locs));
    
    %% 5b) ======================= PASS 2: RESCUE MISSED BREATHS =======================
    % Search for local minimums (valleys) in the filtered signal within the gaps
    % left by Pass 1. This is excellent for finding weaker or closely-spaced breaths.
    secondary_locs = [];

    % Set the threshold for what qualifies as a significant dip in the filtered signal
    thr_filtered_signal = p2_peak_height_coeff * std(cor_f_padded);
    minDist_pass2_global = round(cycleSec * fs * p2_min_dist_coeff); % Renamed to avoid confusion

    % Define search boundaries, including before the first and after the last primary peak
    search_locs = sort([0; primary_locs; N_padded]);

    for i = 1:numel(search_locs) - 1
        start_idx = search_locs(i) + 1;
        end_idx = search_locs(i+1) - 1;
        if start_idx >= end_idx, continue; end

        % --- Adapt MinPeakDistance to the current window size ---
        window_size = end_idx - start_idx + 1;

        % Use the smaller value between the global minimum distance and the window size
        current_min_dist = min(minDist_pass2_global, window_size - 2);

        % If the window is too small to have a peak, skip it
        if current_min_dist < 1
            continue;
        end

        % Search for peaks in the INVERTED filtered signal (i.e., find minimums)
        search_window_f = -cor_f_padded(start_idx:end_idx);

        % Use the new, adaptive 'current_min_dist'
        [~, locs_in_window] = findpeaks(search_window_f, 'MinPeakDistance', current_min_dist, 'MinPeakHeight', thr_filtered_signal);

        if ~isempty(locs_in_window)
            % Convert local window indices back to global padded indices
            actual_locs = start_idx + locs_in_window - 1;
            secondary_locs = [secondary_locs; actual_locs];
        end
    end
    fprintf('Pass 2 (Filtered Signal): Rescued %d additional breath events.\n', numel(secondary_locs));

    % Combine peaks from both passes into an initial list
    secondary_locs = [];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_found_locs_initial = unique(sort([primary_locs; secondary_locs]));
    
    %% 5c) ======================= REFINE PEAK LOCATIONS =======================
    % For each detected location, "snap" it to the true minimum of the RAW
    % cross-correlation signal in a small local window for maximum accuracy.
    refinement_window_samples = round(refinement_window_seconds * fs);
    refined_locs = zeros(size(all_found_locs_initial));
    refined_locs_max_MWI = zeros(size(all_found_locs_initial));
    
    for i = 1:numel(all_found_locs_initial)
        current_loc = all_found_locs_initial(i);
        
        % Define a local search window around the detected peak
        start_search = max(1, current_loc - refinement_window_samples);
        end_search   = min(N_padded, current_loc + refinement_window_samples);
        
        % Find the index of the minimum value in the PADDED filtered signal
        [~, min_idx_local_1] = min(cor_f_padded(start_search:end_search));

        % Find the index of the max value in the PADDED filtered signal
        [~, max_idx_local] = max(mwi_padded(start_search:end_search));
        
        current_loc = start_search + min_idx_local_1 - 1;
        % Define a local search window around the detected peak
        start_search = max(1, current_loc - refinement_window_samples);
        end_search   = min(N_padded, current_loc + refinement_window_samples);
        
        % Find the index of the minimum value in the PADDED ORIGINAL signal
        [~, min_idx_local_2] = min(cor_padded(start_search:end_search));
        
        % Update the location to this precise minimum
        refined_locs(i) = start_search + min_idx_local_2 - 1;

        refined_locs_max_MWI(i) = start_search + max_idx_local - 1;
    end
    
    % Create the final, refined list of unique breath locations
    all_found_locs = unique(refined_locs);
    all_found_locs_MWI = unique(refined_locs_max_MWI);
    fprintf('Refinement: Finalized %d unique breath locations.\n', numel(all_found_locs));
    
    %% 6) ======================= DEFINE REMOVAL REGIONS =======================
    % For each refined breath location, define the start and end of the region to remove.
    halfWin = round(cycleSec * fs / 2);
    regions_padded = []; 
    
    for k = 1:numel(all_found_locs_MWI)
        pk = all_found_locs_MWI(k);
        % Define a large window around the peak to find its boundaries
        ws = max(1, pk - halfWin);
        we = min(N_padded, pk + halfWin);
        local_e = mwi_padded(ws:we);
        pk_relative = pk - ws + 1;
        
        % Find start boundary: look for where the energy rises above a local threshold
        pre_e  = local_e(1:pk_relative); 
        [~, min_index] = min(pre_e); 
        pre_e = pre_e(min_index:end);
        thr_s  = mean(pre_e) + region_width_coeff_start * std(pre_e); 
        rel1 = find(pre_e < thr_s, 1, 'last')+min_index;
        if isempty(rel1) || rel1 == 1, s = ws; else s = ws + rel1 - 1; end
        if s >= pk, s = ws; end
        %s = pk - ceil(0.25 * halfWin);

        % Find end boundary: look for where the energy falls below a local threshold
        post_e = local_e(pk_relative:end); 
        [~, min_index] = min(post_e); 
        post_e = post_e(1:min_index);
        thr_e  = mean(post_e) + region_width_coeff_end * std(post_e); 
        rel2 = find(post_e < thr_e, 1, 'first');
        if isempty(rel2), e = we; elseif rel2 == 1, e = we; else e = (pk - 1) + rel2; end
        %e = pk + ceil(0.25 * halfWin);
        regions_padded(k, :) = [s, e];
    end
    
    %% 7) ======================= COLLECT FINAL INDICES FOR REMOVAL =======================
    % Convert the padded region indices back to the original signal's coordinates
    regions_adjusted = regions_padded - pad_len;
    
    % Collect all individual frame indices from the defined regions
    indices_to_remove = [];
    for k = 1:size(regions_adjusted, 1)
        s = regions_adjusted(k, 1);
        e = regions_adjusted(k, 2);
        % Ensure indices are within the valid bounds of the ORIGINAL signal
        current_indices = max(1, s):min(N, e);
        if ~isempty(current_indices)
            indices_to_remove = [indices_to_remove, current_indices];
        end
    end
    indices_to_remove = unique(indices_to_remove);
    fprintf('Process complete. Identified %d total frames for removal.\n', numel(indices_to_remove));

    %% 8) ======================= VISUALIZATION =======================
    % Prepare signals and locations for plotting by removing the padding
    mwi_plot = mwi_padded(pad_len + 1 : end - pad_len);
    cor_f_plot = cor_f_padded(pad_len + 1 : end - pad_len);
    final_locs_plot = all_found_locs(all_found_locs > pad_len & all_found_locs <= N + pad_len) - pad_len;

    figure('Units','normalized','OuterPosition',[0 0 1 1]);
    
    % Plot 1: Raw signal, final removal frames, and the refined peak locations
    ax1 = subplot(3,1,1);
    plot(ax1, cor, 'b', 'DisplayName', 'Raw Cross-Correlation'); hold(ax1, 'on');
    scatter(ax1, indices_to_remove, cor(indices_to_remove), 10, 'red', 'filled', 'DisplayName', 'Frames to Remove');
    scatter(ax1, final_locs_plot, cor(final_locs_plot), 50, 'k', 'o', 'LineWidth', 1.5, 'DisplayName', 'Refined Breath Minimums');
    title(ax1, 'Raw Signal & Final Removal Decision');
    xlabel(ax1, 'Frame'); ylabel(ax1, 'Cross-Correlation');
    legend(ax1, 'Location', 'best'); grid on;

    % Plot 2: MWI energy, the calculated removal regions, and final peak locations
    ax2 = subplot(3,1,2);
    plot(ax2, mwi_plot,'b','DisplayName','MWI Energy'); hold on;
    yL = ylim(ax2);
    % Dummy patch for a clean legend entry
    patch(ax2, [NaN NaN NaN NaN], [NaN NaN NaN NaN], [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha', 0.4, 'DisplayName', 'Detected Breath Region');
    for k=1:size(regions_adjusted,1)
        s = max(1, regions_adjusted(k,1)); e = min(N, regions_adjusted(k,2));
        patch(ax2, [s e e s], [yL(1) yL(1) yL(2) yL(2)], [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
    end
    scatter(ax2, final_locs_plot, mwi_plot(final_locs_plot), 60, 'k', 'o', 'LineWidth', 1.5, 'DisplayName', 'Final Peak Locations');
    title(ax2, 'MWI Energy and Defined Removal Regions');
    xlabel(ax2, 'Frame'); ylabel(ax2, 'Energy');
    legend(ax2, 'Location', 'best'); grid on;
    
    % Plot 3: Filtered signal with the final, refined peak locations
    ax3 = subplot(3,1,3);
    plot(ax3, cor_f_plot, 'b', 'DisplayName', 'Filtered Signal'); hold(ax3, 'on');
    scatter(ax3, final_locs_plot, cor_f_plot(final_locs_plot), 60, 'k', 'o', 'LineWidth', 1.5, 'DisplayName', 'Final Peak Locations');
    title(ax3, 'Filtered Signal & All Detected Breath Locations');
    xlabel(ax3, 'Frame'); ylabel(ax3, 'Filtered Correlation');
    legend(ax3, 'Location', 'best'); grid on;
    
    % Link all x-axes for synchronized zooming and panning
    linkaxes([ax1, ax2, ax3], 'x');
    
    %% 9) ======================= SAVE OUTPUTS =======================
    save([savingpath filesep 'indices_to_remove'], 'indices_to_remove');
    saveas(gcf, [savingpath filesep 'Deleting_breathing_images_robust.png']);
    savefig(gcf, [savingpath filesep 'Deleting_breathing_images_robust.fig']);
end

function run_rigid_motion_correction_workflow(IQfiles_Bmode, IQfiles_PI, indices_to_remove, ROI_Mask, global_correlation, enable_drift, overwrite_mode, non_rigid_mode)
% RUN_RIGID_MOTION_CORRECTION_WORKFLOW Executes a multi-stage motion correction pipeline for ultrasound data.
%
%   This function performs rigid (translation + rotation) and optional 
%   non-rigid (B-spline) registration on Ultrasound B-Mode and PI 
%   datasets. It is designed to handle "Gated" data, where specific frames 
%   (e.g., during heavy respiration) are excluded from the registration process.
%
%   ALGORITHM OVERVIEW:
%   -------------------
%   1. GATING & CYCLE DETECTION:
%      The function uses 'indices_to_remove' to ignore bad frames. It uses 
%      'regionprops' to identify contiguous blocks of valid frames, treating 
%      each block as a distinct cardiac cycle.
%
%   2. PHASE 1 - RIGID CALCULATION (Intra & Inter):
%      - Intra-Cycle: Stabilizes motion within a single heartbeat.
%      - Master Selection: Identifies the cycle with the highest 'global_correlation'.
%      - Inter-Cycle (Drift): Aligns all other cycles to the Master Cycle 
%        to correct for gradual probe movement or organ shift.
%
%   3. PHASE 2 - TRANSFORMATION & NON-RIGID REFINEMENT:
%      - Applies the calculated rigid transforms.
%      - (Optional) Calculates Non-Rigid B-Spline grids to correct local 
%        tissue deformation.
%      - Crops the image to a "Safe Region" to remove edge artifacts.
%      - Saves the result to disk (overwriting or creating a new folder).
%
%   INPUTS:
%   -------
%   IQfiles_Bmode    : (struct) File list for B-Mode .mat files (dir output).
%   IQfiles_PI       : (struct) File list for Pulse Inversion .mat files. 
%                      Can be empty [] if only processing B-Mode.
%   indices_to_remove: (vector) Global indices of frames to skip (e.g., respiration).
%   ROI_Mask         : (logical matrix) Binary mask indicating the tissue area 
%                      to use for correlation/registration.
%   global_correlation: (vector) Quality metric per frame (0 to 1). Used to 
%                       weight the registration and select the Master Cycle.
%   enable_drift     : (bool, default=true) If true, aligns cycles to each other. 
%                      If false, only stabilizes within the cycle (resets origin).
%   overwrite_mode   : (bool, default=false) If true, overwrites input files. 
%                      If false, saves to a 'registered' subfolder.
%   non_rigid_mode   : (bool, default=false) If true, performs a second pass 
%                      of deformable registration after rigid correction.
%
%   OUTPUTS:
%   --------
%   Files are saved to disk. No variables are returned to the workspace.
%   A 'motion_history' plot is generated in the output folder.
%
%   EXAMPLE:
%       run_rigid_motion_correction_workflow(filesB, filesPI, bad_idxs, ...
%           mask, corr_scores, true, false, true);
%
%   DEPENDENCIES:
%       - calculate_rigid_motion_gated
%       - calculate_inter_cycle_drift
%       - calculate_safe_crop_from_transforms
%       - apply_transforms_to_stack
%       - calculate_non_rigid_motion_gated (if non-rigid enabled)
%
%   Author: Grigori Shapiro
%   Date: 16.1.2026
% -------------------------------------------------------------------------

    %% 0. DEFAULTS & CONFIG
    if nargin < 6, enable_drift = true; end
    if nargin < 7, overwrite_mode = false; end
    if nargin < 8, non_rigid_mode = false; end 

    ENABLE_INTER_CYCLE_CORRECTION = enable_drift; 
    
    fprintf('--- Starting Motion Correction ---\n');
    fprintf('    [Config] Drift: %s | Overwrite: %s | Non-Rigid (Soft): %s\n', ...
            mat2str(ENABLE_INTER_CYCLE_CORRECTION), mat2str(overwrite_mode), mat2str(non_rigid_mode));

    % 1. Initialization
    process_PI = false;
    if nargin > 1 && ~isempty(IQfiles_PI)
        process_PI = true;
        if numel(IQfiles_Bmode) ~= numel(IQfiles_PI)
            warning('Mismatch in number of B-mode and PI files!');
        end
    end

    % --- FOLDER LOGIC ---
    if overwrite_mode
        reg_folder_Bmode = IQfiles_Bmode(1).folder;
        if process_PI, reg_folder_PI = IQfiles_PI(1).folder; end
        fprintf('    [Target] Overwriting files in: %s\n', reg_folder_Bmode);
    else
        reg_folder_Bmode = fullfile(IQfiles_Bmode(1).folder, 'registered');
        if ~exist(reg_folder_Bmode, 'dir'), mkdir(reg_folder_Bmode); end
        if process_PI
            reg_folder_PI = fullfile(IQfiles_PI(1).folder, 'registered');
            if ~exist(reg_folder_PI, 'dir'), mkdir(reg_folder_PI); end
        end
        fprintf('    [Target] Creating new subfolder: %s\n', reg_folder_Bmode);
    end

    Nbuffers = numel(IQfiles_Bmode);
    Batch_Transforms_List = cell(Nbuffers, 1);
    
    % Cycle Registry
    Cycle_Registry = {}; 
    cycle_count = 0;

    % =====================================================================
    % 2. PHASE 1: RIGID CALCULATION
    % =====================================================================
    fprintf('\n>> PHASE 1: Calculating Rigid Motion Models...\n');
    
    global_frame_offset = 0;
    
    for i = 1:Nbuffers
        fprintf('   [Analyze] Batch %d/%d (Start: %s)...\n', i, Nbuffers, datetime("now", "Format", "HH:mm:ss"));
        
        dataStruct = load(fullfile(IQfiles_Bmode(i).folder, IQfiles_Bmode(i).name));
        fieldname = fieldnames(dataStruct);
        batch_Bmode = dataStruct.(fieldname{1});
        [H, W, nFrames] = size(batch_Bmode);
        
        Batch_Transforms_List{i} = zeros(nFrames, 3);
        
        if size(ROI_Mask, 1) ~= H || size(ROI_Mask, 2) ~= W
             resized_Mask = imresize(ROI_Mask, [H, W], 'nearest');
        else
             resized_Mask = ROI_Mask;
        end
        soft_Mask = imgaussfilt(double(resized_Mask), 5.0);
        
        batch_indices_global = (global_frame_offset + 1) : (global_frame_offset + nFrames);
        valid_indices = batch_indices_global(batch_indices_global <= length(global_correlation));
        batch_corr = global_correlation(valid_indices);
        if length(batch_corr) < nFrames, batch_corr(end+1:nFrames) = 0; end
        is_bad_frame_local = ismember(batch_indices_global, indices_to_remove);
        stats = regionprops(~is_bad_frame_local, 'PixelIdxList'); 
        
        % --- Iterate Cycles ---
        for k = 1:length(stats)
            idxList = stats(k).PixelIdxList;
            if length(idxList) > 5
                sub_stack = batch_Bmode(:, :, idxList);
                sub_corr  = batch_corr(idxList);
                
                % A. Calculate Intra-Cycle (Rigid)
                [~, intra_trans, Local_Template] = calculate_rigid_motion_gated(sub_stack, soft_Mask, sub_corr);
                
                % B. Registry (Keeping ONLY the stable baseline template)
                cycle_count = cycle_count + 1;
                Cycle_Registry{cycle_count}.BatchID = i;
                Cycle_Registry{cycle_count}.Indices = idxList;
                Cycle_Registry{cycle_count}.IntraTrans = intra_trans;
                Cycle_Registry{cycle_count}.Template = Local_Template;
                Cycle_Registry{cycle_count}.Score = mean(sub_corr);
            end
        end
        global_frame_offset = global_frame_offset + nFrames;
    end
    
    % 3. SELECTION
    fprintf('\n>> SELECTING BEST CYCLE (GLOBAL MASTER)...\n');
    if cycle_count == 0, warning('No valid cycles found!'); return; end
    
    all_scores = cellfun(@(x) x.Score, Cycle_Registry);
    [max_score, best_idx] = max(all_scores);
    Global_Master_Template = Cycle_Registry{best_idx}.Template;
    
    fprintf('   Best Cycle: #%d (Stability: %.4f)\n', best_idx, max_score);
    
    % =====================================================================
    % 4. PHASE 1.5: GLOBAL DRIFT CORRECTION (STAR TOPOLOGY)
    % =====================================================================
    fprintf('\n>> PHASE 1.5: Aligning Baselines to Global Master...\n');
    
    if ENABLE_INTER_CYCLE_CORRECTION 
        for c = 1:cycle_count
            % The master cycle is our anchor, it doesn't need to move
            if c == best_idx
                continue; 
            end
            
            curr_template = Cycle_Registry{c}.Template;
            
            % Compare THIS cycle's stable baseline directly to the GLOBAL master.
            % Since BOTH are diastolic templates, there is no phase mismatch.
            % This direct comparison mathematically prevents accumulated drift error!
            link_offset = calculate_inter_cycle_drift(curr_template, Global_Master_Template, resized_Mask);
            
            % Apply the absolute offset to the cycle's transforms
            batch_id = Cycle_Registry{c}.BatchID;
            indices = Cycle_Registry{c}.Indices;
            Batch_Transforms_List{batch_id}(indices, :) = Cycle_Registry{c}.IntraTrans + link_offset;
        end
    else
        % If drift correction is disabled, just assign the intra-cycle transforms
        for c = 1:cycle_count
            batch_id = Cycle_Registry{c}.BatchID;
            indices = Cycle_Registry{c}.Indices;
            Batch_Transforms_List{batch_id}(indices, :) = Cycle_Registry{c}.IntraTrans;
        end
    end
    clear Cycle_Registry;
    
    % 5. CROP
    fprintf('\n>> CALCULATING GLOBAL SMART CROP...\n');
    All_Transforms = cat(1, Batch_Transforms_List{:});
    [crop_rect, crops_info] = calculate_safe_crop_from_transforms(All_Transforms, H, W);
    fprintf('   Crop X: %d-%d | Y: %d-%d\n', ...
            crop_rect(1), crop_rect(1)+crop_rect(3), crop_rect(2), crop_rect(2)+crop_rect(4));

    % =====================================================================
    % 6. PHASE 2: APPLICATION + NON-RIGID
    % =====================================================================
    fprintf('\n>> PHASE 2: Applying Transforms, Soft-Correction, and Saving...\n');
    
    motion_history = struct('global_idx', [], 'tx', [], 'ty', [], 'theta', []);
    global_frame_offset = 0;
    
    for i = 1:Nbuffers
        fprintf('   [Save] Batch %d/%d (Start: %s)...\n', i, Nbuffers, datetime("now", "Format", "HH:mm:ss"));
        
        dataStruct = load(fullfile(IQfiles_Bmode(i).folder, IQfiles_Bmode(i).name));
        fieldname = fieldnames(dataStruct);
        batch_Bmode = dataStruct.(fieldname{1});
        [H, W, ~] = size(batch_Bmode);
        
        if process_PI
            dataStructPI = load(fullfile(IQfiles_PI(i).folder, IQfiles_PI(i).name));
            fieldname = fieldnames(dataStructPI);
            batch_PI = dataStructPI.(fieldname{1});
        end
        
        % --- STEP A: Apply Rigid Transforms ---
        transforms = Batch_Transforms_List{i};
        [nFramesBatch, ~] = size(transforms);
        
        % --- ANTI-JITTER SMOOTHING (Per-Cycle Valid Blocks Only!) ---
        % Find the valid (non-breathing) blocks to avoid smoothing across zeros
        batch_indices = (global_frame_offset + 1) : (global_frame_offset + nFramesBatch);
        is_bad = ismember(batch_indices, indices_to_remove);
        stats_smoothing = regionprops(~is_bad, 'PixelIdxList');
        
        % With a framerate of ~200 Hz, an 11-frame window represents ~55 ms.
        % Using a 3rd-order Savitzky-Golay filter preserves the true peak 
        % amplitudes of the cardiac cycle while eliminating high-frequency jitter.
        sgolay_order = 3;
        sgolay_window = 11; % Must be an odd integer
        
        for k = 1:length(stats_smoothing)
            idxList = stats_smoothing(k).PixelIdxList;
            
            % Only smooth if the block is long enough for the filter window
            if length(idxList) >= sgolay_window
                transforms(idxList, 1) = sgolayfilt(transforms(idxList, 1), sgolay_order, sgolay_window); % Tx
                transforms(idxList, 2) = sgolayfilt(transforms(idxList, 2), sgolay_order, sgolay_window); % Ty
                transforms(idxList, 3) = sgolayfilt(transforms(idxList, 3), sgolay_order, sgolay_window); % Theta
            end
        end
        
        fprintf('      Applying Smoothed Rigid Transforms...\n');
        batch_Bmode_reg = apply_transforms_to_stack(batch_Bmode, transforms);
        if process_PI
            batch_PI_reg = apply_transforms_to_stack(batch_PI, transforms); 
        end
        
        % --- STEP B: Non-Rigid (Soft) Refinement ---
        if non_rigid_mode
            fprintf('      Applying Non-Rigid (Soft) Correction...\n');
            
            batch_indices = (global_frame_offset + 1) : (global_frame_offset + nFramesBatch);
            valid_idx = batch_indices(batch_indices <= length(global_correlation));
            batch_corr = global_correlation(valid_idx); 
            if length(batch_corr)<nFramesBatch, batch_corr(end+1:nFramesBatch)=0; end
            
            is_bad = ismember(batch_indices, indices_to_remove);
            stats = regionprops(~is_bad, 'PixelIdxList');
            
            if size(ROI_Mask, 1) ~= H || size(ROI_Mask, 2) ~= W
                 curr_mask = imresize(ROI_Mask, [H, W], 'nearest');
            else
                 curr_mask = ROI_Mask; 
            end
            
            % Dilate mask for non-rigid to cover full range of motion
            max_disp = max(abs(transforms(:,1:2)), [], 1);  % [max_tx, max_ty]
            pad_pixels = ceil(max(max_disp)) + 2;           % +2 safety margin
            se = strel('disk', pad_pixels);
            curr_mask = imdilate(curr_mask, se);

            % Process non-rigid cycles
            for k = 1:length(stats)
                idxList = stats(k).PixelIdxList;
                if length(idxList) > 5
                    sub_stack_B = batch_Bmode_reg(:, :, idxList);
                    sub_corr = batch_corr(idxList);
                    
                    % Calculate B-Spline Grids (Now with progress bar)
                    [corrected_B, grids] = calculate_non_rigid_motion_gated(sub_stack_B, curr_mask, sub_corr);
                    
                    batch_Bmode_reg(:, :, idxList) = corrected_B;
                    
                    if process_PI
                        batch_PI_reg(:, :, idxList) = apply_bspline_to_stack(batch_PI_reg(:, :, idxList), grids);
                    end
                end
            end
        end

        % --- STEP C: Crop ---
        c_xmin = crop_rect(1); c_ymin = crop_rect(2); c_xmax = c_xmin + crop_rect(3); c_ymax = c_ymin + crop_rect(4);
        batch_Bmode_cropped = batch_Bmode_reg(c_ymin:c_ymax, c_xmin:c_xmax, :);
        if process_PI, batch_PI_cropped = batch_PI_reg(c_ymin:c_ymax, c_xmin:c_xmax, :); end
        
        % --- STEP D: Save ---
        [~, nameB, ext] = fileparts(IQfiles_Bmode(i).name);
        
        if overwrite_mode
            save_name_B = [nameB ext];
        else
            save_name_B = [nameB '_registered' ext];
        end
        
        ImageData_tot = batch_Bmode_cropped; 
        save(fullfile(reg_folder_Bmode, save_name_B), 'ImageData_tot', '-v7.3');
        
        if process_PI
            [~, namePI, ext] = fileparts(IQfiles_PI(i).name);
            
            if overwrite_mode
                save_name_PI = [namePI ext];
            else
                save_name_PI = [namePI '_registered' ext];
            end
            
            ImageData_tot = batch_PI_cropped;
            save(fullfile(reg_folder_PI, save_name_PI), 'ImageData_tot', '-v7.3');
        end
        
        % History
        indices = (global_frame_offset + 1) : (global_frame_offset + nFramesBatch);
        motion_history.global_idx = [motion_history.global_idx; indices(:)];
        motion_history.tx = [motion_history.tx; transforms(:,1)];
        motion_history.ty = [motion_history.ty; transforms(:,2)];
        motion_history.theta = [motion_history.theta; transforms(:,3)];
        global_frame_offset = global_frame_offset + nFramesBatch;
    end
    
    visualize_motion_correction(motion_history, reg_folder_Bmode);
    fprintf('--- Workflow Complete. ---\n');
end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [corrected_stack, grids] = calculate_non_rigid_motion_gated(stack, mask, cycle_corr)
    [H, W, nFrames] = size(stack);
    
    % --- Contiguous template (same fix as rigid) ---
    win_len = max(3, min(ceil(nFrames * 0.05), nFrames));
    if nFrames <= win_len
        best_start = 1;
    else
        conv_scores = conv(cycle_corr(:)', ones(1, win_len), 'valid') / win_len;
        [~, best_start] = max(conv_scores);
    end
    best_indices = best_start : (best_start + win_len - 1);
    
    RefFrame_Raw = mean(abs(double(stack(:, :, best_indices))), 3);
    RefFrame_Smooth = imgaussfilt(RefFrame_Raw, 2.0);
    Ref_Norm = mat2gray(RefFrame_Smooth);
    
    corrected_stack = zeros(size(stack), 'like', stack);
    grids = cell(nFrames, 1);
    
    Options.Verbose = false;
    Options.Registration = 'NonRigid';
    Options.Similarity = 'cc';       % More robust for ultrasound
    Options.MaxRef = 1;              % Coarse-to-fine: 2 levels
    Options.Penalty = 1.0;           % Tighter regularization
    Options.Interpolation = 'Cubic'; 
    Options.Spacing = [64 64];       % Starts coarse, refines to ~16

    updater = make_progress_updater(nFrames);
    q = parallel.pool.DataQueue;
    afterEach(q, updater);
    
    fprintf('      Soft-Correction (%d frames): ', nFrames);

    parfor i = 1:nFrames
        Moving_Raw = abs(double(stack(:, :, i)));
        Moving_Smooth = imgaussfilt(Moving_Raw, 2.0);
        Moving_Norm = mat2gray(Moving_Smooth);
        
        [~, grid, spacing, ~, ~, ~] = image_registration(Moving_Norm, Ref_Norm, Options);
        
        % Apply grid to ORIGINAL (unmasked) frame
        corrected_stack(:, :, i) = bspline_transform(grid, stack(:, :, i), spacing, 1);
        grids{i} = struct('grid', grid, 'spacing', spacing);
        
        send(q, i);
    end
    fprintf('\n'); 
end

function corrected_stack = apply_bspline_to_stack(stack, grids)
    [H, W, nFrames] = size(stack);
    corrected_stack = zeros(size(stack), 'like', stack);
    
    % --- PROGRESS BAR SETUP (FIXED) ---
    updater = make_progress_updater(nFrames);
    q = parallel.pool.DataQueue;
    afterEach(q, updater);
    
    fprintf('      Applying to PI (%d frames): ', nFrames);

    parfor i = 1:nFrames
        g_struct = grids{i};
        if ~isempty(g_struct)
            corrected_stack(:, :, i) = bspline_transform(g_struct.grid, stack(:, :, i), g_struct.spacing, 1);
        else
            corrected_stack(:, :, i) = stack(:, :, i);
        end
        send(q, i);
    end
    fprintf('\n');
end

function [corrected_stack, transforms, RefFrame_Raw] = calculate_rigid_motion_gated(stack, mask, cycle_corr)
% CALCULATE_RIGID_MOTION_GATED (Robust Percentile Normalization)
    [H, W, nFrames] = size(stack);
    
    % --- Smart Template Generation ---
    num_for_template = max(3, ceil(nFrames * 0.05));  
    win_len = min(num_for_template, nFrames);
    
    if nFrames <= win_len
        best_start = 1;
    else
        conv_scores = conv(cycle_corr(:)', ones(1, win_len), 'valid') / win_len;
        [~, best_start] = max(conv_scores);
    end
    best_indices = best_start : (best_start + win_len - 1);
    
    anchor_frame = round(mean(best_indices));
    
    % --- Robust Normalization ---
    RefFrame_Raw = mean(abs(double(stack(:, :, best_indices))), 3);
    
    % Use the 99.5th percentile to ignore bright bubbles/noise outliers
    p99 = prctile(RefFrame_Raw(:), 99.5); 
    if p99 == 0, p99 = 1; end 
    
    RefFrame_Norm = RefFrame_Raw / p99;
    RefFrame_Norm(RefFrame_Norm > 1) = 1; % Cap at 1
    RefFrame_Smooth = imgaussfilt(RefFrame_Norm, 1.0);
    
    transforms = zeros(nFrames, 3); 
    [x, y] = meshgrid(1:W, 1:H); x_c = x - W/2; y_c = y - H/2;
    
    opts = optimoptions('lsqnonlin', 'Display', 'off', ...
        'FiniteDifferenceStepSize', 0.5, 'TolFun', 1e-6, 'TolX', 1e-4, 'MaxIter', 80);
    
    fprintf('      Rigid Seq. (%d frames)... ', nFrames);
    
    % 1. FORWARD LOOP
    current_guess = [0, 0, 0]; 
    for i = anchor_frame : nFrames
        MovingFrame_Raw = abs(double(stack(:, :, i)));
        MovingFrame_Norm = MovingFrame_Raw / p99;
        MovingFrame_Norm(MovingFrame_Norm > 1) = 1;
        MovingFrame_Smooth = imgaussfilt(MovingFrame_Norm, 1.0);
        
        cost_fun = @(p) rigid_cost_function_masked(p, MovingFrame_Smooth, RefFrame_Smooth, mask, x_c, y_c);
        
        lb = [-20, -20, -5]; ub = [20, 20, 5];
        try
            p_opt = lsqnonlin(cost_fun, current_guess, lb, ub, opts);
        catch
            p_opt = current_guess; 
        end
        transforms(i, :) = p_opt;
        current_guess = p_opt; 
    end
    
    % 2. BACKWARD LOOP
    current_guess = transforms(anchor_frame, :); 
    for i = (anchor_frame - 1) : -1 : 1
        MovingFrame_Raw = abs(double(stack(:, :, i)));
        MovingFrame_Norm = MovingFrame_Raw / p99;
        MovingFrame_Norm(MovingFrame_Norm > 1) = 1;
        MovingFrame_Smooth = imgaussfilt(MovingFrame_Norm, 1.0);
        
        cost_fun = @(p) rigid_cost_function_masked(p, MovingFrame_Smooth, RefFrame_Smooth, mask, x_c, y_c);
        
        lb = [-20, -20, -5]; ub = [20, 20, 5];
        try
            p_opt = lsqnonlin(cost_fun, current_guess, lb, ub, opts);
        catch
            p_opt = current_guess; 
        end
        transforms(i, :) = p_opt;
        current_guess = p_opt; 
    end
    fprintf('Done.\n');
    corrected_stack = [];
end

function corrected_stack = apply_transforms_to_stack(stack, transforms)
    [H, W, nFrames] = size(stack); corrected_stack = zeros(size(stack), 'like', stack);
    [x, y] = meshgrid(1:W, 1:H); x_c = x - W/2; y_c = y - H/2;
    
    % --- PROGRESS BAR SETUP ---
    updater = make_progress_updater(nFrames);
    q = parallel.pool.DataQueue;
    afterEach(q, updater);
    
    fprintf('      Applying Rigid (%d frames): ', nFrames);

    parfor i = 1:nFrames
        p = transforms(i, :);
        if all(p == 0), corrected_stack(:,:,i) = stack(:,:,i);
        else, corrected_stack(:,:,i) = apply_rigid_transform_fast(double(stack(:,:,i)), p, x_c, y_c); end
        send(q, i);
    end
    fprintf('\n');
end

function updateFn = make_progress_updater(N)
    % This function creates a UNIQUE counter for each loop.
    % It returns a handle to a nested function that has access to 'count'.
    count = 0;
    has_printed = false;
    
    updateFn = @update_inner;
    
    function update_inner(~)
        count = count + 1;
        
        % Update visually every ~2% to avoid console flicker
        if mod(count, ceil(N/50)) == 0 || count == N
            progress = count / N;
            percent = round(progress * 100);
            
            bar_len = 20;
            filled = round(bar_len * progress);
            bar_str = [repmat('=', 1, filled), repmat(' ', 1, bar_len - filled)];
            
            msg = sprintf('[%s] %3d%%', bar_str, percent);
            
            % Backspace only if we have printed before
            if has_printed
                fprintf(repmat('\b', 1, length(msg)));
            end
            fprintf('%s', msg);
            has_printed = true;
        end
    end
end

function offset = calculate_inter_cycle_drift(Local_Template, Global_Master, mask)
% CALCULATE_INTER_CYCLE_DRIFT (Coarse-to-Fine Approach)
    [H, W] = size(Local_Template); 
    [x, y] = meshgrid(1:W, 1:H); x_c = x - W/2; y_c = y - H/2;
    
    Fixed = imgaussfilt(Global_Master, 1.0); 
    Moving = imgaussfilt(Local_Template, 1.0);
    
    % --- STEP 1: Coarse Guess using 2D Cross-Correlation ---
    % CRITICAL FIX: Do NOT multiply by the mask here! The zero-edges of the mask 
    % create artificial high-frequency boundaries that confuse normxcorr2.
    % We perform correlation on the full smoothed images to get the global shift.
    c = normxcorr2(Moving, Fixed);
    [ypeak, xpeak] = find(c == max(c(:)), 1);
    
    ty_guess = ypeak - H;
    tx_guess = xpeak - W;
    
    % TIGHTER SAFETY CHECK: The kidney should not jump by massive amounts 
    % between two stable breathing cycles. If it does, ignore the guess.
    if abs(tx_guess) > 15 || abs(ty_guess) > 15
        tx_guess = 0; ty_guess = 0;
    end
    
    % --- STEP 2: Fine Sub-pixel & Rotation Optimization ---
    opts = optimoptions('lsqnonlin', 'Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-4, 'MaxIter', 100);
    
    % The mask is used ONLY during the fine-tuning step to focus on the kidney
    cost_fun = @(p) rigid_cost_function_masked(p, Moving, Fixed, mask, x_c, y_c);
    
    x0 = [tx_guess, ty_guess, 0]; 
    
    % Tighten the bounds around the coarse guess to prevent the optimizer 
    % from running away due to local speckle noise.
    lb = [tx_guess - 5, ty_guess - 5, -3]; 
    ub = [tx_guess + 5, ty_guess + 5, 3];
    
    try
        offset = lsqnonlin(cost_fun, x0, lb, ub, opts);
    catch
        offset = x0; % Fallback to cross-correlation result if optimization fails
    end
end

function residuals = rigid_cost_function_masked(p, Moving, Ref, mask, x_c, y_c)
    % Apply the current rigid transformation guess
    Warped = apply_rigid_transform_fast(Moving, p, x_c, y_c);
    
    % Calculate the difference and apply the region of interest mask
    diff_img = Warped - Ref; 
    weighted_diff = diff_img .* mask; 
    
    % --- L1 Norm Optimization Trick ---
    % lsqnonlin minimizes the SUM OF SQUARES of the returned array.
    % By returning the square root of the absolute value, lsqnonlin computes 
    % (sqrt(|x|))^2 = |x|, effectively minimizing the L1 norm instead of L2.
    % This makes the registration highly robust to speckle noise outliers.
    % We add a tiny epsilon (1e-8) to prevent division by zero in the 
    % Jacobian calculation when the difference is exactly zero.
    residuals = sqrt(abs(weighted_diff(:)) + 1e-8);
end

function Warped = apply_rigid_transform_fast(Img, p, x_c, y_c)
    tx = p(1); ty = p(2); theta_deg = p(3);
    th = theta_deg * (pi/180); cos_t = cos(th); sin_t = sin(th);
    u = (x_c - tx) * cos_t + (y_c - ty) * sin_t; v = -(x_c - tx) * sin_t + (y_c - ty) * cos_t;
    [H, W] = size(Img); u = u + W/2; v = v + H/2;
    fillVal = mean(Img(:)); Warped = interp2(Img, u, v, 'spline', fillVal);
end

function [crop_rect, info] = calculate_safe_crop_from_transforms(transforms, H, W)
    min_valid_x = 1; max_valid_x = W; min_valid_y = 1; max_valid_y = H;
    xc = W/2; yc = H/2; corners = [1, 1; W, 1; 1, H; W, H];
    for i = 1:size(transforms, 1)
        tx = transforms(i, 1); ty = transforms(i, 2); theta = transforms(i, 3);
        th_rad = theta * (pi/180); c = cos(th_rad); s = sin(th_rad);
        curr_x = zeros(4,1); curr_y = zeros(4,1);
        for c_idx = 1:4
            x0 = corners(c_idx, 1); y0 = corners(c_idx, 2);
            rel_x = x0 - xc; rel_y = y0 - yc;
            curr_x(c_idx) = (rel_x * c - rel_y * s) + xc + tx;
            curr_y(c_idx) = (rel_x * s + rel_y * c) + yc + ty;
        end
        min_valid_x = max(min_valid_x, ceil(max(curr_x([1, 3])))); max_valid_x = min(max_valid_x, floor(min(curr_x([2, 4]))));
        min_valid_y = max(min_valid_y, ceil(max(curr_y([1, 2])))); max_valid_y = min(max_valid_y, floor(min(curr_y([3, 4]))));
    end
    if max_valid_x <= min_valid_x || max_valid_y <= min_valid_y, crop_rect = [1, 1, W-1, H-1];
    else, crop_rect = [min_valid_x, min_valid_y, max_valid_x - min_valid_x, max_valid_y - min_valid_y]; end
    info.left = min_valid_x - 1; info.right = W - max_valid_x; info.top = min_valid_y - 1; info.bottom = H - max_valid_y;
end

function visualize_motion_correction(history, save_folder)
% VISUALIZE_MOTION_CORRECTION Generates a summary plot of the calculated motion.
%
%   Displays translation (Tx, Ty) in pixels and rotation (Theta) in degrees
%   over the entire recording session.
% -------------------------------------------------------------------------

    if isempty(history.global_idx), return; end

    % Sort data by frame index to ensure chronological order
    [sorted_idx, sort_order] = sort(history.global_idx);
    tx = history.tx(sort_order); 
    ty = history.ty(sort_order); 
    theta = history.theta(sort_order);

    % Graphic Settings
    color_tx = [0.00, 0.45, 0.74]; % Blue
    color_ty = [0.49, 0.18, 0.56]; % Purple
    color_rot = [0.47, 0.67, 0.19]; % Green
    axis_color_trans = [0.00, 0.00, 0.50]; 
    axis_color_rot = color_rot;

    fig = figure('Name', 'Motion Correction Summary', 'Position', [100, 100, 1200, 800], 'Color', 'w', 'Visible', 'off');

    % --- Subplot 1: Full Timeline ---
    ax1 = subplot(2,1,1); 
    yyaxis left;
    plot(sorted_idx, tx, '-', 'DisplayName', 'Tx', 'Color', color_tx, 'LineWidth', 1); hold on; 
    plot(sorted_idx, ty, '-', 'DisplayName', 'Ty', 'Color', color_ty, 'LineWidth', 1);
    ylabel('Translation [pix]', 'FontWeight', 'bold'); 
    ax1.YColor = axis_color_trans;
    
    yyaxis right; 
    plot(sorted_idx, theta, '-', 'DisplayName', 'Rot', 'Color', color_rot, 'LineWidth', 1); 
    ylabel('Rotation [deg]', 'FontWeight', 'bold'); 
    ax1.YColor = axis_color_rot;
    
    xlabel('Frame Number');
    title('Global Motion Correction Summary'); 
    legend('Location', 'best'); grid on; xlim([min(sorted_idx), max(sorted_idx)]);

    % --- Subplot 2: Zoomed View (Start of recording) ---
    % Zoom into the first 600 frames (or less if recording is short)
    zoom_end = min(sorted_idx(end), sorted_idx(1) + 600); 
    mask = sorted_idx <= zoom_end;
    
    ax2 = subplot(2,1,2); 
    yyaxis left;
    plot(sorted_idx(mask), tx(mask), '.-', 'Color', color_tx, 'LineWidth', 1.2); hold on; 
    plot(sorted_idx(mask), ty(mask), '.-', 'Color', color_ty, 'LineWidth', 1.2);
    ylabel('Translation [pix]', 'FontWeight', 'bold'); 
    ax2.YColor = axis_color_trans;
    
    yyaxis right; 
    plot(sorted_idx(mask), theta(mask), '.-', 'Color', color_rot, 'LineWidth', 1.2); 
    ylabel('Rotation [deg]', 'FontWeight', 'bold');
    ax2.YColor = axis_color_rot;
    
    xlabel('Frame Number');
    title('Zoomed View (First ~600 Frames)'); grid on;
    
    % Save
    saveas(fig, fullfile(save_folder, 'Motion_Correction_Summary.png'));
    saveas(fig, fullfile(save_folder, 'Motion_Correction_Summary.fig')); 
    close(fig);
end

function IQfiles = sortFilesByNumber(IQfiles)
    % Initialize array for file numbers
    fileNumbers = NaN(length(IQfiles), 1);
    
    for i = 1:length(IQfiles)
        % Get the full name (extension doesn't usually matter for number extraction, 
        % but ignoring it is safer if the extension has numbers like .mp4)
        [~, name, ~] = fileparts(IQfiles(i).name);
        
        % Find ALL consecutive digits in the filename
        % 'match' returns a cell array of strings, e.g., {'2025', '001'}
        numbersFound = regexp(name, '\d+', 'match');
        
        if ~isempty(numbersFound)
            % Assume the batch index is the LAST number in the filename
            % This handles "Data_2025_Set_001" correctly (takes 001)
            lastNumberStr = numbersFound{end};
            fileNumbers(i) = str2double(lastNumberStr);
        else
            % If no number found, assign a value to put it at the end (or start)
            fileNumbers(i) = inf; 
        end
    end
    
    % Sort the files based on the extracted numbers
    % checking for NaNs ensures we don't crash if a file has no numbers
    [~, sortIdx] = sort(fileNumbers);
    IQfiles = IQfiles(sortIdx);
end

