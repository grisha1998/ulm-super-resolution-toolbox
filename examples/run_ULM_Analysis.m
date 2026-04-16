% =================================================================================
% FILENAME: run_ULM_Analysis.m
% =================================================================================
%
% PURPOSE:
%   Main, self-contained runner script for the ULM framework. Orchestrates
%   the full ULM processing pipeline: parameter initialization, localization
%   and tracking, SR reconstruction, and visualization.
%
% SCRIPT WORKFLOW:
%   1.  CONFIGURATION:    Set the data folder path.
%   2.  PARAMETERS:       Load defaults and experiment-specific settings.
%   3.  PROCESSOR INIT:   Initialize the ULM_Processor object.
%   4.  PIPELINE:         Run Step0 boundaries, generate mean B-mode images,
%                         then run the localization & tracking loop.
%   5.  RENDER & SAVE:    Reconstruct SR maps at multiple min-track-length
%                         thresholds, save all figures.
%
% HOW TO USE:
%   1.  Set 'root_data_folder' in Section 1.
%   2.  Ensure all CLASS files and helper functions are on the MATLAB path.
%   3.  Run the script.
%
% AUTHOR: Grigori Shapiro
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
%  1. USER CONFIGURATION
%  ========================================================================
fprintf('\n>> STEP 1: CONFIGURATION\n');

% ---> Set the path to the main data folder of the experiment here.
root_data_folder = '';
fprintf('   - Data folder set to: %s\n', root_data_folder);

%% ========================================================================
%  2. PARAMETER INITIALIZATION
%  ========================================================================
fprintf('\n>> STEP 2: PARAMETER INITIALIZATION\n');
try
    % Pass the data folder to setDefaultParams.
    % If your version of setDefaultParams does not accept a folder argument,
    % use: params = setDefaultParams(); params.io.data_folder = root_data_folder;
    data_folder = root_data_folder;
    if exist('data_folder', 'var')
        params = setDefaultParams(false, data_folder);
    else
        params = setDefaultParams(false);
    end
    fprintf('   - Parameters initialized successfully.\n');
catch ME
    fprintf('   - ERROR: Failed to initialize parameters. %s\n', ME.message);
    return;
end

%% ========================================================================
%  3. ULM PROCESSOR INITIALIZATION
%  ========================================================================
fprintf('\n>> STEP 3: INITIALIZING ULM PROCESSOR\n');
try
    processor = ULM_Processor(params);
catch ME
    fprintf('   - ERROR: Failed to initialize ULM_Processor: %s\n', ME.message);
    return;
end

%% ========================================================================
%  4. EXECUTE THE ULM PIPELINE
%  ========================================================================
fprintf('\n>> STEP 4: RUNNING THE MAIN PROCESSING PIPELINE\n');
tic;

% --- Step 0: Interactively define or load Crop and Mask boundaries ---
processor.Step0_InitializeBoundaries();

% --- Define paths for the mean B-mode images (used for overlay) ---
res = params.render.upsampling_factor;
file_vessel_path = fullfile(processor.general_results_dir, ['mean_bmode_vessel_', num2str(res), '.mat']);
file_tissue_path = fullfile(processor.general_results_dir, ['mean_bmode_tissue_', num2str(res), '.mat']);
file_vessel_png  = fullfile(processor.general_results_dir, ['mean_bmode_vessel_', num2str(res), '.png']);
file_tissue_png  = fullfile(processor.general_results_dir, ['mean_bmode_tissue_', num2str(res), '.png']);

% --- Generate (or load) mean B-mode images ---
if ~exist(file_vessel_path, 'file')
    disp('Generating mean_bmode_vessel...');
    mean_bmode_vessel = generate_mean_bmode_image(processor.data_files, params, 0);
    save(file_vessel_path, 'mean_bmode_vessel');
else
    load(file_vessel_path, 'mean_bmode_vessel');
end
imwrite(mat2gray(mean_bmode_vessel), file_vessel_png);

if ~exist(file_tissue_path, 'file')
    disp('Generating mean_bmode_tissue...');
    mean_bmode_tissue = generate_mean_bmode_image(processor.data_files, params, 1);
    save(file_tissue_path, 'mean_bmode_tissue');
else
    load(file_tissue_path, 'mean_bmode_tissue');
end
imwrite(mat2gray(mean_bmode_tissue), file_tissue_png);

% --- Main localization and tracking loop ---
processor.run_localization_and_tracking_loop();

totalTime = toc;
fprintf('   - Pipeline execution finished in %.1f minutes.\n', totalTime / 60);

%% ========================================================================
%  5. RENDER FINAL MAPS AT MULTIPLE TRACK-LENGTH THRESHOLDS
%  ========================================================================
fprintf('\n>> STEP 5: RENDERING FINAL IMAGES\n');

MIN_TRACK_LENGTH = [5, 8, 10, 15, 20, 25, 30]; % <-- TUNE AS NEEDED

for i = 1:length(MIN_TRACK_LENGTH)

    processor.run_streaming_reconstruction(MIN_TRACK_LENGTH(i));

    %% ====================================================================
    %  6. SAVE AND DISPLAY RESULTS
    %  ====================================================================
    fprintf('\n>> STEP 6: SAVING AND VISUALIZING RESULTS (Min Length = %d)\n', MIN_TRACK_LENGTH(i));

    % Save core results only on the first iteration
    if i == 1
        try
            processor.saveResults();
        catch ME
            fprintf('   - An error occurred while saving results: %s\n', ME.message);
        end
    end

    try
        % --- Determine line/font sizes from upsampling factor ---
        res_points       = [1, 3, 5];
        fontsize_points  = [6, 8, 10];
        linewidth_points = [2, 4, 5];
        res_clamped = max(min(res, res_points(end)), res_points(1));
        fontsize    = round(interp1(res_points, fontsize_points, res_clamped));
        linewidth   = round(interp1(res_points, linewidth_points, res_clamped), 1);

        % --- 6a. Density Map ---
        fig_density = figure('Name', 'Super-Resolution Density Map', 'NumberTitle', 'off');
        density_processed = processor.densityMap .^ (1/3);
        imshow(density_processed, [], 'colormap', hot);
        positive_pixels = density_processed(density_processed > 0);
        if ~isempty(positive_pixels)
            clim([0, prctile(positive_pixels(:), 99.5)]);
        end
        title(sprintf('Super-Resolution Density Map (Min Length = %d)', MIN_TRACK_LENGTH(i)), 'FontSize', 10);
        clb = colorbar; clb.Label.String = 'Counts^{1/3}';
        add_scale_bar(params, size(processor.densityMap), linewidth, fontsize);
        filename = sprintf('Density_Map_minLen_%d.png', MIN_TRACK_LENGTH(i));
        export_fig(fig_density, fullfile(processor.results_dir, filename), '-png', '-r300');

        % --- 6b. Velocity Map (Filtered) ---
        fig_vel_filt = figure('Name', 'Publication Velocity Map (Filtered)', 'NumberTitle', 'off');
        SR_vel_filtered = imgaussfilt(processor.velocityMap, 0.6);
        cm = colormap([0 0 0; jet(256)]);
        imshow(SR_vel_filtered, [], 'colormap', cm);
        pos_vel = SR_vel_filtered(SR_vel_filtered > 0);
        if ~isempty(pos_vel), clim([0, prctile(pos_vel(:), 99.5)]); end
        title(sprintf('Super-Resolution Velocity Map - Filtered (Min Length = %d)', MIN_TRACK_LENGTH(i)), 'FontSize', 10);
        clb = colorbar; clb.Label.String = 'Velocity [mm/s]';
        add_scale_bar(params, size(processor.velocityMap), linewidth, fontsize);
        filename = sprintf('Velocity_Map_Filtered_minLen_%d.png', MIN_TRACK_LENGTH(i));
        export_fig(fig_vel_filt, fullfile(processor.results_dir, filename), '-png', '-r300');

        % --- 6c. Velocity Map (Unfiltered) ---
        fig_vel_unfilt = figure('Name', 'Publication Velocity Map (Unfiltered)', 'NumberTitle', 'off');
        imshow(processor.velocityMap, [], 'colormap', cm);
        pos_vel = processor.velocityMap(processor.velocityMap > 0);
        if ~isempty(pos_vel), clim([0, prctile(pos_vel(:), 99.5)]); end
        title(sprintf('Super-Resolution Velocity Map - Unfiltered (Min Length = %d)', MIN_TRACK_LENGTH(i)), 'FontSize', 10);
        clb = colorbar; clb.Label.String = 'Velocity [mm/s]';
        add_scale_bar(params, size(processor.velocityMap), linewidth, fontsize);
        filename = sprintf('Velocity_Map_Unfiltered_minLen_%d.png', MIN_TRACK_LENGTH(i));
        export_fig(fig_vel_unfilt, fullfile(processor.results_dir, filename), '-png', '-r300');

        % --- 6d. B-Mode Overlay (using tissue image) ---
        fprintf('   - Generating overlay image...\n');
        load(file_tissue_path, 'mean_bmode_tissue');
        mean_bmode = mean_bmode_tissue;

        if ndims(mean_bmode) == 3
            bmode_gray = rgb2gray(mean_bmode);
        else
            bmode_gray = mean_bmode;
        end
        bmode_adjusted = adapthisteq(mat2gray(bmode_gray));
        bmode_rgb = cat(3, bmode_adjusted, bmode_adjusted, bmode_adjusted);

        density_data = processor.densityMap .^ (1/3);
        clim_max = prctile(density_data(:), 99.5);
        density_norm = mat2gray(density_data, [0, clim_max]);
        hot_colormap = hot(256);
        density_rgb = ind2rgb(gray2ind(density_norm, 256), hot_colormap);

        fig_overlay = figure('Name', 'Overlay', 'NumberTitle', 'off');
        imshow(density_rgb);
        hold on;
        h_bmode = imshow(bmode_rgb);
        set(h_bmode, 'AlphaData', 0.4);
        hold off;
        axis image;
        if exist('clim', 'builtin')
            clim([0, clim_max]);
        else
            caxis([0, clim_max]);
        end
        clb = colorbar; colormap(clb, hot_colormap);
        clb.Limits = [0, clim_max]; clb.Label.String = 'Counts^{1/3}';
        title(sprintf('ULM Density Map with B-Mode Overlay (Min Length = %d)', MIN_TRACK_LENGTH(i)), 'FontSize', 10);
        filename = sprintf('Overlay_minLen_%d.png', MIN_TRACK_LENGTH(i));
        export_fig(fig_overlay, fullfile(processor.results_dir, filename), '-png', '-r300');

        % --- 6e. Combined Velocity-Density Map (HSV) ---
        fprintf('   - Generating combined Velocity (Color) + Density (Brightness) map...\n');
        fig_combined = figure('Name', 'Combined Vel-Dens Map', 'NumberTitle', 'off');

        vel_data = processor.velocityMap;
        pos_vel = vel_data(vel_data > 0);
        clim_vel = ~isempty(pos_vel) * prctile(pos_vel(:), 99.5) + isempty(pos_vel) * 1;
        vel_norm = mat2gray(vel_data, [0, clim_vel]);

        dens_data = processor.densityMap .^ (1/3);
        pos_dens = dens_data(dens_data > 0);
        clim_dens = ~isempty(pos_dens) * prctile(pos_dens(:), 99.0) + isempty(pos_dens) * 1;
        dens_norm = mat2gray(dens_data, [0, clim_dens]);
        dens_norm = dens_norm .^ 0.7; % gamma correction (lower = brighter)

        n_colors = 256;
        cmap = jet(n_colors);
        rgb_base = ind2rgb(gray2ind(vel_norm, n_colors), cmap);
        hsv_img = rgb2hsv(rgb_base);
        hsv_img(:,:,3) = dens_norm;
        final_combined_img = hsv2rgb(hsv_img);
        mask_bg = dens_data == 0;
        r_ch = final_combined_img(:,:,1); r_ch(mask_bg) = 0;
        g_ch = final_combined_img(:,:,2); g_ch(mask_bg) = 0;
        b_ch = final_combined_img(:,:,3); b_ch(mask_bg) = 0;
        final_combined_img = cat(3, r_ch, g_ch, b_ch);

        imshow(final_combined_img);
        title(sprintf('Combined Velocity & Density (Min Length = %d)', MIN_TRACK_LENGTH(i)), 'FontSize', 10);
        colormap(gca, cmap); clim([0, clim_vel]);
        clb = colorbar; clb.Label.String = 'Velocity [mm/s]'; clb.Color = [0 0 0];
        add_scale_bar(params, size(processor.velocityMap), linewidth, fontsize);
        filename = sprintf('Combined_VelDens_Map_minLen_%d.png', MIN_TRACK_LENGTH(i));
        export_fig(fig_combined, fullfile(processor.results_dir, filename), '-png', '-r300');

        fprintf('   - Figures generated and saved successfully.\n');
    catch ME
        fprintf('\n   - Could not display final images. Error: %s ---\n', ME.message);
    end

end % end MIN_TRACK_LENGTH loop

fprintf('\n=======================================================\n');
fprintf('               ULM Analysis Complete\n');
fprintf('=======================================================\n');


%% ========================================================================
%  Local Helper Functions
%  ========================================================================

function add_scale_bar(params, mapSize, linewidth, fontsize)
    % Adds a 1 mm scale bar to the current axes.
    hold on;
    pixel_size_mm = params.expParams.fovX / mapSize(2);
    scalebar_length_pixels = 1 / pixel_size_mm;
    x_pos = 0.05 * mapSize(2);
    y_pos = 0.95 * mapSize(1);
    plot([x_pos, x_pos + scalebar_length_pixels], [y_pos, y_pos], 'w-', 'LineWidth', linewidth);
    text(x_pos + scalebar_length_pixels/2, y_pos - fontsize * linewidth, '1 mm', ...
        'Color', 'w', 'FontSize', fontsize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    hold off;
end

function str = B_SWITCH(bool_val)
    % Converts a boolean to 'ON'/'OFF' string.
    if bool_val, str = 'ON'; else, str = 'OFF'; end
end

function mean_bmode = generate_mean_bmode_image(data_files, params, vessel_or_tissue)
    % Generates a mean B-mode image from all data buffers for overlay purposes.
    % vessel_or_tissue: 0 = vessel (SVD high-pass), 1 = tissue (SVD low-pass)
    sRes_dims = params.expParams.size(1:2) * params.render.upsampling_factor;

    if isempty(data_files)
        warning('No data files found to generate mean B-mode image.');
        mean_bmode = zeros(sRes_dims);
        return;
    end

    num_files = length(data_files);
    accumulated_mean = 0;

    for i = 1:num_files
        dataStruct = load(fullfile(data_files(i).folder, data_files(i).name));
        dataFields = fieldnames(dataStruct);
        rawData = dataStruct.(dataFields{1});

        if vessel_or_tissue == 0
            filteredData = SVD_filter(rawData, params.filter.svd_cutoff);
        elseif vessel_or_tissue == 1
            c1 = params.filter.svd_cutoff;
            filteredData = SVD_filter(rawData, [0, c1(1)]);
        else
            filteredData = SVD_filter(rawData, params.filter.svd_cutoff);
        end

        if isfield(params.filter, 'enable_butterworth') && params.filter.enable_butterworth
            filteredData = Butterworth_bandpass_filter(filteredData, params.filter.butter_cutoff, ...
                params.acq.framerate, params.filter.butter_order);
        end

        accumulated_mean = accumulated_mean + double(mean(abs(filteredData), 3));
    end

    final_mean = (accumulated_mean / num_files) .^ 0.5;
    mean_bmode = imresize(final_mean, sRes_dims);
end