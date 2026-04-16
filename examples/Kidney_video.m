%% ========================================================================
%           MODULAR VIDEO GENERATION SCRIPT (IMPROVED VERSION)
% =========================================================================
% SCRIPT: Kidney_video
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   A highly modular, interactive script designed to generate playback videos 
%   (AVIs) from heavy in-vivo ultrasound datasets (specifically Kidney scans).
%   - Advantages: Eliminates the need to manually plot frames for inspection. 
%     It features an interactive GUI for folder selection, robust numerical 
%     file sorting, and automatically generates both raw (unfiltered) and 
%     processed (filtered) videos side-by-side. This is critical for visually 
%     validating clutter-filter performance in real-time.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Interactive Setup: Prompts the user to select the root experiment 
%      directory via a native UI dialog (`uigetdir`).
%   2. File Management & Sorting (`sortFilesByNumber`): Ultrasound machines 
%      often save sequential files as 'data_1, data_10, data_2'. This script 
%      uses regular expressions to extract the true integer index and sorts 
%      the files chronologically, ensuring smooth video playback.
%   3. Data Parsing: Iterates through large 3D Casorati matrices (.mat files) 
%      representing consecutive buffers of ultrasound data.
%   4. Video Rendering & Gamma Correction:
%      - For complex (IQ) or real data, it computes the signal envelope (absolute value).
%      - Applies non-linear Gamma Correction (`img .^ gamma_factor`) to 
%        compress the dynamic range. This visually enhances weak microbubble 
%        signals while preventing saturation from highly echogenic tissue.
%      - Normalizes frames to [0, 1] using `mat2gray` and writes them 
%        sequentially to the output `.avi` file.
%
% SYNTAX OPTIONS:
%   Run directly from the MATLAB editor or command window:
%   >> Kidney_video
%
% USER CONFIGURATIONS (Hardcoded in script):
%   - create_unfiltered : (Logical) Flag to render the raw data video.
%   - create_filtered   : (Logical) Flag to render the clutter-filtered video.
%   - gamma_factor      : (Double) Typical value: 0.5 to 0.7 for flow visualization.
%
% EXPECTED DIRECTORY STRUCTURE (INPUTS):
%   The selected target folder should contain:
%   [Target_Folder]\
%      ├── file_1.mat  (Containing variables like 'imageData' and 'imageDataFilt')
%      ├── file_2.mat
%
% OUTPUTS:
%   Automatically creates a 'Videos' subdirectory inside the target folder 
%   and saves the compiled `.avi` video files (e.g., 'file_1_unfiltered.avi').
% =========================================================================

%% 1. SCRIPT SETUP
% clear;
% close all;
% clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('=======================================================\n');
fprintf('         Ultrasound Video Generation Script\n');
fprintf('=======================================================\n\n');

%% 2. ======================= USER CONFIGURATION =======================

% --- Improvement: Interactive folder selection with default path ---
fprintf('>> Please select the main experiment folder...\n');
% --- Change: Added default path for opening the window ---
default_path = pwd; %default_path = 'E:\Ilovich Lab\In Vivo Data'; 
if ~exist(default_path, 'dir')
    default_path = ''; % If folder doesn't exist, open in last used
end
root_data_folder = uigetdir(default_path, 'Select the Main Experiment Folder');

if root_data_folder == 0, error('No folder selected. Script terminated.'); end
fprintf('   - Selected Folder: %s\n', root_data_folder);

% --- Improvement: Dynamic subfolder selection including Processed Data ---
target_special_folder = ''; % Flag for special path handling
try
    paramsFile = fullfile(root_data_folder, 'info.txt');
    expParams = getExpParams(paramsFile);
    base_data_path = fullfile(root_data_folder, expParams.bubbleType);
    
    % Find all standard subfolders in Bmode folder
    bmode_root = fullfile(base_data_path, 'Bmode');
    potential_subfolders = dir(bmode_root);
    potential_subfolders = potential_subfolders([potential_subfolders.isdir]);
    potential_subfolders = potential_subfolders(~ismember({potential_subfolders.name},{'.','..'}));
    
    % Create list of options
    folder_names = {potential_subfolders.name};
    
    % --- NEW FEATURE: Check for Registered + Breath Deletion folder ---
    processed_path_relative = fullfile('registered', 'BreathDeletion');
    processed_full_path = fullfile(bmode_root, processed_path_relative);
    
    has_processed_data = exist(processed_full_path, 'dir');
    if has_processed_data
        % Add as the last option
        folder_names{end+1} = '>> PROCESSED DATA (Registered & Breath Deletion)';
    end

    if isempty(folder_names)
        error('No subfolders found in %s', bmode_root);
    end
    
    fprintf('\n>> Please select the data source to process:\n');
    for i = 1:numel(folder_names)
        fprintf('   %d. %s\n', i, folder_names{i});
    end
    
    choice = 0;
    while ~ismember(choice, 1:numel(folder_names))
        choice = input('Enter your choice: ');
    end
    
    selected_name = folder_names{choice};
    
    % Determine the actual subfolder path to use
    if has_processed_data && choice == numel(folder_names)
        % User selected the special processed option
        subfolder_to_process = 'BreathDeletion'; % Naming for output file
        target_special_folder = processed_path_relative; % Path modifier
        fprintf('   - Selected Data Source: Registered & Breath Deleted Data\n\n');
    else
        % User selected a standard subfolder
        subfolder_to_process = selected_name;
        target_special_folder = selected_name;
        fprintf('   - Processing Subfolder: %s\n\n', subfolder_to_process);
    end

catch ME
    error('Could not read experiment parameters or find subfolders. Check info.txt and folder structure.\nERROR: %s', ME.message);
end


% --- Video Settings ---
output_video_fps = expParams.shootingRate;
gamma_correction_factor = 0.5;

% --- Workflow Switches ---
workflow_switches = struct();
workflow_switches.create_bmode_unfiltered = true;
workflow_switches.create_bmode_filtered   = true;
workflow_switches.create_pi_unfiltered    = false;
workflow_switches.create_pi_filtered      = false;

fprintf('>> Workflow Configuration:\n');
fprintf('   - B-mode (Unfiltered): %s\n', B_SWITCH(workflow_switches.create_bmode_unfiltered));
fprintf('   - B-mode (Filtered):   %s\n', B_SWITCH(workflow_switches.create_bmode_filtered));
fprintf('   - PI (Unfiltered):     %s\n', B_SWITCH(workflow_switches.create_pi_unfiltered));
fprintf('   - PI (Filtered):       %s\n\n', B_SWITCH(workflow_switches.create_pi_filtered));


%% 3. ======================= FILTER SELECTION =======================
params = struct();
params.butterworth_enabled = false;

if workflow_switches.create_bmode_filtered || workflow_switches.create_pi_filtered
    filter_options = {
        'SVD_filter', 'Efficient Manual SVD (requires range)';
        'SVD_SSM',    'Adaptive SVD using SSM (automatic)';
        'DCC_SVD',    'DCC-SVD Clustering (automatic)'
    };
    
    fprintf('>> A filtered video is requested. Please select a clutter filter:\n');
    for i = 1:size(filter_options, 1)
        fprintf('   %d. %s\n', i, filter_options{i, 2});
    end
    
    choice = 0;
    while ~ismember(choice, 1:size(filter_options, 1))
        choice = input('Enter your choice: ');
    end
    
    params.filter_method = filter_options{choice, 1};
    
    if strcmpi(params.filter_method, 'SVD_filter')
        params.svd_range = [];
        while isempty(params.svd_range) || numel(params.svd_range) ~= 2 || params.svd_range(1) >= params.svd_range(2)
            try
                params.svd_range = input('Enter the SVD cutoff range for blood signal [start, end] (e.g., [15, 80]): ');
                if numel(params.svd_range) ~= 2
                    fprintf('   ERROR: Please enter exactly two numbers.\n');
                elseif params.svd_range(1) >= params.svd_range(2)
                    fprintf('   ERROR: The start value must be smaller than the end value.\n');
                end
            catch
                fprintf('   ERROR: Invalid input format. Please use brackets, e.g., [15, 80].\n');
                params.svd_range = []; % Reset on error
            end
        end
    end

    butter_choice = '';
    while ~ismember(lower(butter_choice), {'y', 'n'})
        butter_choice = input('>> Apply Butterworth bandpass filter after SVD? (y/n): ', 's');
    end
    if lower(butter_choice) == 'y'
        params.butterworth_enabled = true;
        params.butter_order = 4;
        params.butter_cutoff_freq = [5, 50];
        fprintf('   - Butterworth filter ENABLED.\n');
    end

else
    % --- Correction Here ---
    % Set default value when filtering is not performed
    params.filter_method = 'None';
    fprintf('>> No filtered video requested. Skipping filter selection.\n');
    % --- End Correction ---
end

%% 4. ======================= SCRIPT EXECUTION =======================
t_start = tic;
try
    params.framerate = expParams.shootingRate;
    video_root_folder = fullfile(root_data_folder, 'Results', 'Videos');
    
    % --- Loop structure to process Bmode/PI ---
    data_types_to_process = {};
    if workflow_switches.create_bmode_unfiltered || workflow_switches.create_bmode_filtered
        data_types_to_process{end+1} = 'Bmode';
    end
    if workflow_switches.create_pi_unfiltered || workflow_switches.create_pi_filtered
        data_types_to_process{end+1} = 'PI';
    end

    for i = 1:numel(data_types_to_process)
        current_data_type = data_types_to_process{i};
        
        % Construct the input path
        if strcmp(target_special_folder, 'original')
            target_special_folder = '';
        end
        data_path = fullfile(base_data_path, current_data_type, target_special_folder);
        
        fprintf('\n>> Processing %s files from: %s\n', current_data_type, target_special_folder);
        
        data_files = dir(fullfile(data_path, '*.mat'));
        
        if ~isempty(data_files)
            data_files = sortFilesByNumber(data_files);
            
            % --- NEW: Determine Processing State Folder Name ---
            % Analyze the path to decide the category name
            is_reg = contains(target_special_folder, 'registered', 'IgnoreCase', true);
            is_bd  = contains(target_special_folder, 'BreathDeletion', 'IgnoreCase', true);

            if is_reg && is_bd
                proc_state_folder = 'Registered_BreathDeleted';
            elseif is_reg
                proc_state_folder = 'Registered';
            elseif is_bd
                proc_state_folder = 'BreathDeleted';
            else
                % If it's a standard subfolder (like 'Batch1'), we label it RawData
                % You can append the subfolder name if you want: ['RawData_' subfolder_to_process]
                proc_state_folder = 'RawData'; 
            end

            % --- NEW: Construct the Hierarchical Output Path ---
            % Structure: Results/Videos / DataType / ProcessingState
            final_output_folder = fullfile(video_root_folder, current_data_type, proc_state_folder);
            
            % Check/Create the folder
            if ~exist(final_output_folder, 'dir')
                mkdir(final_output_folder);
            end

            video_basename = sprintf('SuperFrame_%s', current_data_type);

            % Pass the specific final folder to the function
            MakeVideos(data_files, video_basename, final_output_folder, output_video_fps, params, gamma_correction_factor, workflow_switches);
        else
            fprintf('   - WARNING: No %s files found in %s.\n', current_data_type, data_path);
        end
    end

catch ME
    fprintf('\n\n--- SCRIPT FAILED ---\n');
    fprintf('ERROR: %s\n', ME.message);
    fprintf('At line %d in file %s\n', ME.stack(1).line, ME.stack(1).name);
    rethrow(ME); 
end

t_end = toc(t_start);
fprintf('\n-------------------------------------------------------\n');
fprintf('Script finished in %.1f minutes.\n', t_end / 60);
fprintf('=======================================================\n');


%% ========================================================================
%                            HELPER FUNCTIONS
% ========================================================================

function MakeVideos(IQfiles, baseName, outputFolder, outputFPS, params, gamma_factor, switches)
    % Determine data type (Bmode/PI) from basename
    is_bmode_run = contains(baseName, 'Bmode', 'IgnoreCase', true);
    
    % Determine which videos to create based on global switches
    create_unfiltered = (is_bmode_run && switches.create_bmode_unfiltered) || (~is_bmode_run && switches.create_pi_unfiltered);
    create_filtered = (is_bmode_run && switches.create_bmode_filtered) || (~is_bmode_run && switches.create_pi_filtered);

    % --- 1. Define Improved Output Paths ---
    if create_unfiltered
        outputDirUnfiltered = fullfile(outputFolder, 'Unfiltered');
        if ~exist(outputDirUnfiltered, 'dir'), mkdir(outputDirUnfiltered); end
    end
    if create_filtered
        filter_name_tag = params.filter_method;
        if strcmpi(filter_name_tag, 'SVD_filter')
            % --- Correction found in next line ---
            % Used params.svd_range instead of old params.svd_cutoff
            filter_name_tag = sprintf('%s_range-%d-%d', filter_name_tag, params.svd_range(1), params.svd_range(2));
        end
        if params.butterworth_enabled
            filter_name_tag = [filter_name_tag, '_BP-filtered'];
        end
        outputDirFiltered = fullfile(outputFolder, filter_name_tag);
        if ~exist(outputDirFiltered, 'dir'), mkdir(outputDirFiltered); end
    end
    
    % Initialize video objects
    outputVideo = [];
    outputVideoFilt = [];
    
    % --- 2. Loop through files, detect data type and write frames ---
    for n = 1:numel(IQfiles)
        fprintf('     - Processing file %d of %d: %s\n', n, numel(IQfiles), IQfiles(n).name);
        
        dataStruct = load(fullfile(IQfiles(n).folder, IQfiles(n).name));
        dataFields = fieldnames(dataStruct);
        imageData = dataStruct.(dataFields{1});
        
        is_rgb_data = (ndims(imageData) == 4 && size(imageData, 3) == 3);
        
        if n == 1 % Define video writers on first file
            if create_unfiltered
                videoName_unfiltered = sprintf('%s_gamma-%.2f_%dFPS.avi', baseName, gamma_factor, outputFPS);
                outputVideo = VideoWriter(fullfile(outputDirUnfiltered, videoName_unfiltered));
                outputVideo.Quality = 100;
                outputVideo.FrameRate = outputFPS;
                open(outputVideo);
                fprintf('   - Creating Unfiltered video in: %s\n', outputDirUnfiltered);
            end

            if create_filtered
                if is_rgb_data
                    fprintf('   - WARNING: Cannot apply SVD filter to RGB data. Filtered video creation skipped.\n');
                    create_filtered = false;
                else
                    videoName_filtered = sprintf('%s_filtered_%dFPS.avi', baseName, outputFPS);
                    outputVideoFilt = VideoWriter(fullfile(outputDirFiltered, videoName_filtered));
                    outputVideoFilt.Quality = 100;
                    outputVideoFilt.FrameRate = outputFPS;
                    open(outputVideoFilt);
                    fprintf('   - Creating Filtered video in:   %s\n', outputDirFiltered);
                end
            end
        end
        
        imageDataFilt = [];
        if create_filtered && ~is_rgb_data
            switch params.filter_method
                case 'SVD_filter',   imageDataFilt = SVD_filter(imageData, params.svd_range, ReconstructionMode="blood");
                case 'SVD_SSM',      imageDataFilt = SVD_SSM(imageData, ReconstructionMode="blood");
                case 'DCC_SVD',      imageDataFilt = DCC_SVD(imageData, params.framerate, ReconstructionMode="blood");
            end
            
            if params.butterworth_enabled
                fprintf('       - Applying Butterworth bandpass filter...\n');
                imageDataFilt = Butterworth_bandpass_filter(imageDataFilt, params.butter_cutoff_freq, params.framerate, params.butter_order);
            end
            imageDataFilt(~isfinite(imageDataFilt)) = 0;
        end
        
        if is_rgb_data
            num_frames = size(imageData, 4);
            for i = 1:num_frames
                if create_unfiltered
                    frame_enhanced = im2double(imageData(:, :, :, i)).^(gamma_factor);
                    writeVideo(outputVideo, frame_enhanced);
                end
            end
        else % Grayscale data
            num_frames = size(imageData, 3);
            for i = 1:num_frames
                if create_unfiltered
                    frame_unfiltered = mat2gray(abs(imageData(:, :, i)).^(gamma_factor));
                    writeVideo(outputVideo, frame_unfiltered);
                end
                if create_filtered && ~isempty(imageDataFilt)
                    frame_filtered = mat2gray(abs(imageDataFilt(:, :, i)));
                    writeVideo(outputVideoFilt, frame_filtered);
                end
            end
        end
    end
    
    if ~isempty(outputVideo), close(outputVideo); end
    if ~isempty(outputVideoFilt), close(outputVideoFilt); end
    fprintf('   - Video processing complete for %s.\n', baseName);
end

function IQfiles = sortFilesByNumber(IQfiles)
    % This function is fine, no changes needed
    fileNumbers = zeros(length(IQfiles), 1);
    for i = 1:length(IQfiles)
        numStr = regexp(IQfiles(i).name, '\d+', 'match');
        if ~isempty(numStr), fileNumbers(i) = str2double(numStr{end}); end
    end
    [~, sortIdx] = sort(fileNumbers);
    IQfiles = IQfiles(sortIdx);
end

function str = B_SWITCH(bool_val)
    if bool_val, str = 'ON'; else, str = 'OFF'; end
end