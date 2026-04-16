%% ====================================================================================================
%                               PHANTOM VIDEO GENERATION SCRIPT
% ====================================================================================================
% SCRIPT: Bmode_video
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   A comprehensive video rendering engine for Ultrasound Phantom experiments. 
%   Unlike the Kidney script, this handles both structural B-mode imaging 
%   and functional Pulse Inversion (PI) contrast imaging, incorporating 
%   on-the-fly signal processing.
%   - Advantages: Seamlessly integrates SVD clutter filtering and Butterworth 
%     temporal bandpass filtering *directly* into the rendering pipeline. 
%     It applies physically accurate Logarithmic Compression (dB scale) for 
%     B-mode (tissue) and Gamma correction for PI (flow), mimicking the 
%     display pipelines of clinical ultrasound scanners.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Automatic Parameter Parsing: Reads acquisition parameters directly 
%      from a companion `info.txt` file (e.g., PRF/Framerate).
%   2. Modality Detection: Automatically detects and processes both 'Bmode' 
%      and 'PI' subdirectories if present.
%   3. B-Mode Processing Pipeline (Structural):
%      - Extracts the envelope of the RF/IQ signal (`abs()`).
%      - Applies Logarithmic Compression (`logCompress` helper). Since acoustic 
%        impedance spans several orders of magnitude, linear scaling hides 
%        soft tissue details. The signal is converted to decibels (dB), thresholded 
%        by a dynamic range (e.g., 60 dB), and normalized.
%   4. PI Processing Pipeline (Contrast/Flow):
%      - Optionally applies SVD filtering to separate tissue and flow.
%      - Optionally applies a Butterworth bandpass filter to eliminate residual 
%        low-frequency tissue motion or high-frequency electronic noise.
%      - Extracts envelope and applies Gamma correction to enhance contrast agents.
%   5. Rendering: Exports standardized AVI files using the specified framerate.
%
% SYNTAX OPTIONS:
%   Run directly from the MATLAB editor or command window:
%   >> Bmode_video
%
% USER CONFIGURATIONS (Internal Flags):
%   - apply_SVD         : (Logical) Enable/Disable SVD clutter filter.
%   - apply_Butterworth : (Logical) Enable/Disable temporal bandpass filter.
%   - db_range          : (Double) Dynamic range for B-mode (e.g., 50 dB).
%   - db_cutoff         : (Double) Noise floor rejection limit.
%
% EXPECTED DIRECTORY STRUCTURE (INPUTS):
%   [Experiment_Folder]\
%      ├── info.txt    (Must contain parameters, especially PRF)
%      ├── Bmode\      (Folder with structural .mat files)
%      └── PI\         (Folder with contrast/pulse-inversion .mat files)
%
% OUTPUTS:
%   Renders and saves processed `.avi` files within 'Videos' subdirectories 
%   inside both the Bmode and PI folders.
% =========================================================================

%% 1. ENVIRONMENT SETUP
% Clear workspace and command window, close figures
clear;
close all;
clc;

% Add the current directory and all subfolders (including PALA addons) to the MATLAB path
addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('=======================================================\n');
fprintf('         Phantom Video Generation Script\n');
fprintf('=======================================================\n\n');

%% 2. USER CONFIGURATION AND INITIALIZATION

% --- Interactive Folder Selection ---
% Prompt the user to select the root folder containing the experiment data.
fprintf('>> Please select the PHANTOM experiment folder...\n');
default_path = ''; % <-- Set to the folder containing setDefaultParams.m
                   % Leave '' to use current MATLAB path

% If default path doesn't exist, open in current directory
if ~exist(default_path, 'dir'), default_path = ''; end

root_data_folder = uigetdir(default_path, 'Select Phantom Experiment Folder');
if root_data_folder == 0
    error('No folder selected. Script terminated.'); 
end
fprintf('   - Selected Folder: %s\n', root_data_folder);

% --- Load Experiment Parameters ---
% Attempts to load acquisition parameters (e.g., Frame Rate) from 'info.txt'.
try
    paramsFile = fullfile(root_data_folder, 'info.txt');
    if exist(paramsFile, 'file')
        expParams = getExpParams(paramsFile);
        fprintf('   - Loaded parameters from info.txt\n');
    else
        % Fallback defaults if info.txt is missing
        warning('info.txt not found. Using default parameters.');
        expParams.shootingRate = 1000; 
        expParams.bubbleType = ''; 
    end
    
    % Adjust base path if the data is nested inside a subfolder named after the bubble type
    if ~isempty(expParams.bubbleType) && exist(fullfile(root_data_folder, expParams.bubbleType), 'dir')
        base_data_path = fullfile(root_data_folder, expParams.bubbleType);
    else
        base_data_path = root_data_folder;
    end
    
catch ME
    error('Error initializing parameters: %s', ME.message);
end

% --- Visualization Settings ---
% Define global video settings
output_video_fps = expParams.shootingRate;          % Framerate for the output AVI file
gamma_correction_factor = 0.6;  % Gamma correction for linear data (PI/Flow)

% Define specific settings for B-mode visualization (Log Compression)
bmode_settings = struct();
bmode_settings.use_dB = false;   % Enable Log compression (dB scale)
bmode_settings.dB_limit = 45;   % Dynamic range in dB (e.g., 40-60 dB)
bmode_settings.dB_cutoff = 5;   % Noise floor cutoff (lifts the black level)

% --- Workflow Automation ---
% Automatically detect which data folders exist to avoid asking unnecessary questions
workflow_switches = struct();

has_bmode = exist(fullfile(base_data_path, 'Bmode'), 'dir');
has_pi    = exist(fullfile(base_data_path, 'PI'), 'dir');

% Default configuration: Process whatever exists, disable filtering by default for B-mode
workflow_switches.create_bmode_unfiltered = has_bmode;
workflow_switches.create_bmode_filtered   = false; 
workflow_switches.create_pi_unfiltered    = has_pi;
workflow_switches.create_pi_filtered      = false; 

% Display the auto-detected configuration
fprintf('>> Workflow Configuration (Auto-Detected):\n');
fprintf('   - B-mode (Unfiltered): %s\n', B_SWITCH(workflow_switches.create_bmode_unfiltered));
fprintf('   - B-mode (Filtered):   %s\n', B_SWITCH(workflow_switches.create_bmode_filtered));
fprintf('   - PI (Unfiltered):     %s\n', B_SWITCH(workflow_switches.create_pi_unfiltered));
fprintf('   - PI (Filtered):       %s\n\n', B_SWITCH(workflow_switches.create_pi_filtered));

% Allow user to override the auto-configuration
user_override = input('>> Press ENTER to continue, or "c" to configure switches manually: ', 's');
if strcmpi(user_override, 'c')
    workflow_switches.create_bmode_unfiltered = input('   Create B-mode Unfiltered? (1/0): ');
    workflow_switches.create_bmode_filtered   = input('   Create B-mode Filtered (SVD)? (1/0): ');
    workflow_switches.create_pi_unfiltered    = input('   Create PI Unfiltered? (1/0): ');
    workflow_switches.create_pi_filtered      = input('   Create PI Filtered (SVD)? (1/0): ');
end

%% 3. FILTER CONFIGURATION (SVD & BUTTERWORTH)
params = struct();
params.butterworth_enabled = false;

% Only prompt for filter settings if at least one filtered video is requested
if workflow_switches.create_bmode_filtered || workflow_switches.create_pi_filtered
    filter_options = {
        'SVD_filter', 'Efficient Manual SVD (requires range)';
        'SVD_SSM',    'Adaptive SVD using SSM (automatic)';
        'DCC_SVD',    'DCC-SVD Clustering (automatic)'
    };
    
    fprintf('\n>> Filtered video requested. Please select a clutter filter:\n');
    for i = 1:size(filter_options, 1)
        fprintf('   %d. %s\n', i, filter_options{i, 2});
    end
    
    choice = 0;
    while ~ismember(choice, 1:size(filter_options, 1))
        choice = input('Enter your choice: ');
    end
    
    params.filter_method = filter_options{choice, 1};
    
    % If Manual SVD is selected, ask for cutoff ranks
    if strcmpi(params.filter_method, 'SVD_filter')
        params.svd_range = input('   Enter SVD range [start, end] (e.g., [10, 100]): ');
        if isempty(params.svd_range), params.svd_range = [10 100]; end
    end

    % Optional Butterworth Bandpass Filter
    butter_choice = input('>> Apply Butterworth bandpass filter after SVD? (y/n): ', 's');
    if strcmpi(butter_choice, 'y')
        params.butterworth_enabled = true;
        params.butter_order = 2;
        params.butter_cutoff_freq = [50, 250]; % Wide range suitable for phantoms
        fprintf('   - Butterworth filter ENABLED.\n');
    end
    
    % Store framerate for temporal filtering functions
    params.framerate = expParams.shootingRate;

else
    params.filter_method = 'None';
end

%% 4. MAIN EXECUTION LOOP
t_start = tic;
try
    % Create main results directory
    video_root_folder = fullfile(root_data_folder, 'Results', 'Videos');
    
    % Build list of data types to process based on user selection
    data_types_to_process = {};
    if workflow_switches.create_bmode_unfiltered || workflow_switches.create_bmode_filtered
        data_types_to_process{end+1} = 'Bmode';
    end
    if workflow_switches.create_pi_unfiltered || workflow_switches.create_pi_filtered
        data_types_to_process{end+1} = 'PI';
    end

    % --- Iterate through selected data types (Bmode / PI) ---
    for i = 1:numel(data_types_to_process)
        current_data_type = data_types_to_process{i};
        
        % Define input directory
        data_path = fullfile(base_data_path, current_data_type);
        
        fprintf('\n>> Processing %s files from: %s\n', current_data_type, data_path);
        
        % Find all .mat files
        data_files = dir(fullfile(data_path, '*.mat'));
        
        if ~isempty(data_files)
            % Sort files numerically (e.g., file1, file2, file10)
            data_files = sortFilesByNumber(data_files);
            
            % Define output directory for this batch
            proc_state_folder = 'RawData'; 
            final_output_folder = fullfile(video_root_folder, current_data_type, proc_state_folder);
            
            if ~exist(final_output_folder, 'dir')
                mkdir(final_output_folder);
            end

            video_basename = sprintf('SuperFrame_%s', current_data_type);

            % --- Call the Video Generation Function ---
            MakeVideos(data_files, video_basename, final_output_folder, output_video_fps, ...
                       params, gamma_correction_factor, workflow_switches, bmode_settings);
        else
            fprintf('   - WARNING: No %s files found in %s.\n', current_data_type, data_path);
        end
    end

catch ME
    % Error Handling
    fprintf('\n\n--- SCRIPT FAILED ---\n');
    fprintf('ERROR: %s\n', ME.message);
    fprintf('At line %d in file %s\n', ME.stack(1).line, ME.stack(1).name);
end

% Final Timing Report
t_end = toc(t_start);
fprintf('\n-------------------------------------------------------\n');
fprintf('Script finished in %.1f minutes.\n', t_end / 60);
fprintf('Videos saved to: %s\n', video_root_folder);
fprintf('=======================================================\n');


%% ========================================================================
%                            HELPER FUNCTIONS
% ========================================================================

function MakeVideos(IQfiles, baseName, outputFolder, outputFPS, params, gamma_factor, switches, bmode_settings)
% MakeVideos - Core function to process data files and generate AVI videos.
%
% Inputs:
%   IQfiles        - Struct array of .mat files to process.
%   baseName       - Naming prefix for the video files.
%   outputFolder   - Directory where videos will be saved.
%   outputFPS      - Framerate for the video.
%   params         - Filtering parameters.
%   gamma_factor   - Gamma value for linear visualization.
%   switches       - Workflow boolean flags.
%   bmode_settings - Struct containing Log Compression settings (dB).

    % Check if current run is B-mode (to determine visualization method)
    is_bmode_run = contains(baseName, 'Bmode', 'IgnoreCase', true);
    
    % Determine tasks for this specific call based on global switches
    create_unfiltered = (is_bmode_run && switches.create_bmode_unfiltered) || (~is_bmode_run && switches.create_pi_unfiltered);
    create_filtered = (is_bmode_run && switches.create_bmode_filtered) || (~is_bmode_run && switches.create_pi_filtered);

    % --- 1. Prepare Output Directories ---
    if create_unfiltered
        outputDirUnfiltered = fullfile(outputFolder, 'Unfiltered');
        if ~exist(outputDirUnfiltered, 'dir'), mkdir(outputDirUnfiltered); end
    end
    
    if create_filtered
        % Create a descriptive folder name based on filter settings
        filter_name_tag = params.filter_method;
        if strcmpi(filter_name_tag, 'SVD_filter')
            filter_name_tag = sprintf('%s_range-%d-%d', filter_name_tag, params.svd_range(1), params.svd_range(2));
        end
        if params.butterworth_enabled
            filter_name_tag = [filter_name_tag, '_BP-filtered'];
        end
        outputDirFiltered = fullfile(outputFolder, filter_name_tag);
        if ~exist(outputDirFiltered, 'dir'), mkdir(outputDirFiltered); end
    end
    
    outputVideo = [];
    outputVideoFilt = [];
    
    % --- 2. Iterate over all .mat files ---
    for n = 1:numel(IQfiles)
        fprintf('     - Processing file %d of %d: %s\n', n, numel(IQfiles), IQfiles(n).name);
        
        % Load Data
        dataStruct = load(fullfile(IQfiles(n).folder, IQfiles(n).name));
        dataFields = fieldnames(dataStruct);
        imageData = dataStruct.(dataFields{1});
        
        % Detect if data is RGB (cannot filter/compress RGB)
        is_rgb_data = (ndims(imageData) == 4 && size(imageData, 3) == 3);
        
        % --- Initialize Video Writers (Only on the first file) ---
        if n == 1 
            if create_unfiltered
                suffix = '';
                if is_bmode_run, suffix = '_dB'; end
                videoName = sprintf('%s_Raw%s_%dFPS.avi', baseName, suffix, outputFPS);
                
                outputVideo = VideoWriter(fullfile(outputDirUnfiltered, videoName));
                outputVideo.Quality = 95;
                outputVideo.FrameRate = outputFPS;
                open(outputVideo);
            end

            if create_filtered
                if is_rgb_data
                     fprintf('   - WARNING: Cannot filter RGB data. Skipping filtered video.\n');
                     create_filtered = false;
                else
                    videoNameFilt = sprintf('%s_filtered_%dFPS.avi', baseName, outputFPS);
                    outputVideoFilt = VideoWriter(fullfile(outputDirFiltered, videoNameFilt));
                    outputVideoFilt.Quality = 95;
                    outputVideoFilt.FrameRate = outputFPS;
                    open(outputVideoFilt);
                end
            end
        end
        
        % --- PROCESSING LOGIC (Identical to Kidney Pipeline) ---
        imageDataFilt = [];
        if create_filtered && ~is_rgb_data
            % NOTE: We pass 'imageData' directly. We DO NOT take abs() here.
            % This preserves phase information for complex SVD processing.
            
            switch params.filter_method
                case 'SVD_filter'
                    % ReconstructionMode='blood' handles complex logic internally
                    imageDataFilt = SVD_filter(imageData, params.svd_range, 'ReconstructionMode', 'blood');
                case 'SVD_SSM'
                    imageDataFilt = SVD_SSM(imageData, 'ReconstructionMode', 'blood');
                case 'DCC_SVD'
                    imageDataFilt = DCC_SVD(imageData, params.framerate, 'ReconstructionMode', 'blood');
            end
            
            if params.butterworth_enabled
                imageDataFilt = Butterworth_bandpass_filter(imageDataFilt, params.butter_cutoff_freq, params.framerate, params.butter_order);
            end
            
            % Remove artifacts (NaN/Inf)
            imageDataFilt(~isfinite(imageDataFilt)) = 0;
        end
        
        % --- VISUALIZATION LOOP ---
        if is_rgb_data
            % Handle RGB data (rare for raw Bmode/PI but possible)
            num_frames = size(imageData, 4);
            for f = 1:num_frames
                 if create_unfiltered
                    writeVideo(outputVideo, im2double(imageData(:,:,:,f)));
                 end
            end
        else
            % Handle Standard Grayscale / Complex Data
            num_frames = size(imageData, 3);
            for f = 1:num_frames
                
                % 1. Visualize Unfiltered Frame
                if create_unfiltered
                    rawFrame = abs(imageData(:, :, f)); % Convert to magnitude for display
                    
                    if is_bmode_run && bmode_settings.use_dB
                        % Apply Log Compression (dB) for B-mode
                        frameToShow = logCompress(rawFrame, bmode_settings.dB_limit, bmode_settings.dB_cutoff);
                    else
                        % Apply Linear/Gamma for PI or standard
                        frameToShow = mat2gray(rawFrame .^ gamma_factor);
                    end
                    writeVideo(outputVideo, frameToShow);
                end
                
                % 2. Visualize Filtered Frame
                if create_filtered && ~isempty(imageDataFilt)
                    % Filtered data (Blood flow) is viewed linearly/Gamma (No Log Compression)
                    filtFrame = abs(imageDataFilt(:, :, f));
                    frameFiltShow = mat2gray(filtFrame .^ gamma_factor);
                    writeVideo(outputVideoFilt, frameFiltShow);
                end
            end
        end
    end
    
    % Close Video Writers
    if ~isempty(outputVideo), close(outputVideo); end
    if ~isempty(outputVideoFilt), close(outputVideoFilt); end
    fprintf('   -> Saved video: %s\n', baseName);
end

% --- Helper: Log Compression for B-mode ---
function imgOut = logCompress(imgIn, db_range, db_cutoff)
    % Converts linear amplitude to dB scale for visualization.
    % db_range: Dynamic range (e.g., 40-60 dB)
    % db_cutoff: Noise floor rejection
    
    imgIn(imgIn == 0) = eps; % Avoid log(0)
    
    mx = max(imgIn(:));
    if mx == 0, mx = 1; end
    
    % Convert to dB normalized to max
    img_dB = 20 * log10(imgIn / mx);
    
    % Clip bottom range
    img_dB(img_dB < -db_range) = -db_range;
    
    % Normalize to [0, 1]
    imgOut = (img_dB + db_range) / db_range;
    
    % Apply cutoff to remove noise floor
    if db_cutoff > 0
        cut_norm = db_cutoff / db_range;
        imgOut(imgOut < cut_norm) = 0;
        % Rescale the remaining dynamic range to 0-1
        imgOut = (imgOut - cut_norm) / (1 - cut_norm);
        imgOut(imgOut < 0) = 0;
    end
end

% --- Helper: Sort Files Numerically ---
function files = sortFilesByNumber(files)
    fileNumbers = zeros(length(files), 1);
    for i = 1:length(files)
        % Extract number from filename using regex
        numStr = regexp(files(i).name, '\d+', 'match');
        if ~isempty(numStr)
            fileNumbers(i) = str2double(numStr{end});
        end
    end
    [~, sortIdx] = sort(fileNumbers);
    files = files(sortIdx);
end

% --- Helper: Boolean to String ---
function str = B_SWITCH(bool_val)
    if bool_val, str = 'ON'; else, str = 'OFF'; end
end