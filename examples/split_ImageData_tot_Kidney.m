function split_ImageData_tot_Kidney(data_folder, M_or_N, Sub_Batch_Size)
% =========================================================================
% FUNCTION: split_ImageData_tot_Kidney
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   A robust data management utility designed to partition massive high-frame-rate 
%   ultrasound datasets into manageable sub-batches. 
%   - Advantages: Ultrasound Localization Microscopy (ULM) often produces 
%     gigabytes of data per acquisition, which can easily exceed available RAM 
%     during processing (e.g., SVD filtering). This function automatically 
%     slices these massive 3D matrices along the temporal dimension into 
%     smaller, sequential files. It ensures data integrity by archiving the 
%     original files rather than overwriting them.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Initialization & Scanning: The function accepts a root directory and 
%      scans for `.mat` files in two expected subdirectories: 'Bmode' and 'PI' 
%      (Pulse Inversion).
%   2. Dynamic Size Detection: It loads a single sample file to dynamically 
%      determine the original `Batch_Size` (the length of the 3rd dimension, 
%      representing frames/time).
%   3. Partitioning Math: Calculates the required number of sub-batches based 
%      on the user-defined `Sub_Batch_Size`.
%   4. Slicing & Export (B-mode & PI):
%      - Iterates through every file in the directory.
%      - Slices the 3D matrix along the 3rd dimension: `(:,:, idxStart:idxEnd)`.
%      - Saves the newly created sub-batch with a sequential naming convention.
%   5. Data Archiving (Safety Mechanism): To prevent data loss or infinite 
%      processing loops, the original massive `.mat` file is moved into a 
%      newly created 'original' subfolder.
%
% SYNTAX OPTIONS:
%   split_ImageData_tot_Kidney(data_folder, M_or_N, Sub_Batch_Size)
%
% EXAMPLES:
%   % Split a folder containing batches of 3000 frames into chunks of 1500 frames:
%   split_ImageData_tot_Kidney('C:\Data\Kidney_Scan_01', 'M', 1500);
%
% INPUTS:
%   data_folder    - (Type: String/Char) 
%                    The root directory path containing the data. This folder 
%                    MUST contain 'Bmode' and 'PI' subdirectories.
%                    Example: 'D:\Experiments\Subject_01'
%   M_or_N         - (Type: String/Char) 
%                    A prefix identifier used to name the output files. Useful 
%                    for distinguishing between different subjects, acquisition 
%                    modes, or states.
%                    Example: 'M' or 'Subject1_ConditionA'.
%   Sub_Batch_Size - (Type: Integer) 
%                    The maximum number of temporal frames to include in each 
%                    new sub-batch file. 
%                    Example: 1500 or 2000.
%
% OUTPUTS:
%   None returned to the MATLAB workspace. The function modifies the filesystem 
%   by generating new `.mat` files and moving original files into an 'original' 
%   archive folder.
% =========================================================================

    % Displaying a start message
    fprintf('Starting to split batches in the folder: %s\n', data_folder);
    
    % Get the list of B-mode and PI files
    % Using fullfile is safer for cross-platform compatibility
    IQfiles_Bmode = dir(fullfile(data_folder, 'Bmode', '*.mat'));
    IQfiles_PI    = dir(fullfile(data_folder, 'PI', '*.mat'));
    
    % ---------------------------------------------------------------------
    % Step 1: Initialization and Parameter Determination
    % ---------------------------------------------------------------------
    
    % Check if we have any files at all
    if isempty(IQfiles_Bmode) && isempty(IQfiles_PI)
        fprintf('No .mat files found in Bmode or PI folders. Exiting.\n');
        return;
    end
    
    % Try to load one file to determine dimensions (Batch_Size)
    % We prefer B-mode, but if it's empty, we'll use PI
    if ~isempty(IQfiles_Bmode)
        refFile = IQfiles_Bmode(1);
    else
        refFile = IQfiles_PI(1);
    end
    
    % Load the reference file to get sizes
    fprintf('Loading reference file to determine batch size: %s\n', refFile.name);
    dataStruct = load(fullfile(refFile.folder, refFile.name));
    dataFields = fieldnames(dataStruct);
    ImageData_tot = dataStruct.(dataFields{1});
    
    % Determine the batch size (number of frames in the original file)
    Batch_Size = size(ImageData_tot, 3);
    
    % Calculate the number of sub-batches needed per file
    N_Sub_Batches = ceil(Batch_Size / Sub_Batch_Size);
    
    % ---------------------------------------------------------------------
    % Step 2: Process B-mode files (Only if files exist)
    % ---------------------------------------------------------------------
    if ~isempty(IQfiles_Bmode)
        N_Batches_Bmode = numel(IQfiles_Bmode);
        fprintf('Found %d B-mode files. Processing...\n', N_Batches_Bmode);
        
        for i = 1:N_Batches_Bmode
            fprintf('Processing B-mode batch %d/%d\n', i, N_Batches_Bmode);
            tmp_Bmode = load(fullfile(IQfiles_Bmode(i).folder, IQfiles_Bmode(i).name));
            dataFields = fieldnames(tmp_Bmode);
            
            for j = 1:N_Sub_Batches
                % Determine the frames to extract for the current sub-batch
                idxStart = (j - 1) * Sub_Batch_Size + 1;
                if j * Sub_Batch_Size >= Batch_Size
                    idxEnd = size(tmp_Bmode.(dataFields{1}), 3); % Use actual end
                else
                    idxEnd = j * Sub_Batch_Size;
                end
                
                ImageData_tot = tmp_Bmode.(dataFields{1})(:, :, idxStart:idxEnd);
                
                % Create a unique count identifier for the sub-batch
                count = (i - 1) * N_Sub_Batches + j;
                savingpath = fullfile(data_folder, 'Bmode', [M_or_N 'SuperFrame_Bmode_' num2str(count) '.mat']);
                
                % Ensure directory exists (though it should if files were found)
                if ~exist(fileparts(savingpath), 'dir')
                    mkdir(fileparts(savingpath));
                end
                
                save(savingpath, 'ImageData_tot');
                fprintf('  Saved sub-batch %d as %s\n', count, savingpath);
            end
            
            % Move the original file
            originalFolder = fullfile(IQfiles_Bmode(i).folder, 'original');
            if ~exist(originalFolder, 'dir')
                mkdir(originalFolder);
            end
            movefile(fullfile(IQfiles_Bmode(i).folder, IQfiles_Bmode(i).name), originalFolder);
        end
    else
        fprintf('No B-mode files found. Skipping B-mode processing.\n');
    end
    
    % ---------------------------------------------------------------------
    % Step 3: Process PI files (Only if files exist)
    % ---------------------------------------------------------------------
    if ~isempty(IQfiles_PI)
        N_Batches_PI = numel(IQfiles_PI);
        fprintf('Found %d PI files. Processing...\n', N_Batches_PI);
        
        for i = 1:N_Batches_PI
            fprintf('Processing PI batch %d/%d\n', i, N_Batches_PI);
            tmp_PI = load(fullfile(IQfiles_PI(i).folder, IQfiles_PI(i).name));
            dataFields = fieldnames(tmp_PI);
            
            for j = 1:N_Sub_Batches
                % Determine the frames to extract for the current sub-batch
                idxStart = (j - 1) * Sub_Batch_Size + 1;
                if j * Sub_Batch_Size >= Batch_Size
                    idxEnd = size(tmp_PI.(dataFields{1}), 3);
                else
                    idxEnd = j * Sub_Batch_Size;
                end
                
                ImageData_tot = tmp_PI.(dataFields{1})(:, :, idxStart:idxEnd);
                
                % Create a unique count identifier for the sub-batch
                count = (i - 1) * N_Sub_Batches + j;
                savingpath = fullfile(data_folder, 'PI', [M_or_N '_SuperFrame_PI_' num2str(count) '.mat']);
                
                if ~exist(fileparts(savingpath), 'dir')
                    mkdir(fileparts(savingpath));
                end
                
                save(savingpath, 'ImageData_tot');
                fprintf('  Saved sub-batch %d as %s\n', count, savingpath);
            end
            
            % Move the original file
            originalFolder = fullfile(IQfiles_PI(i).folder, 'original');
            if ~exist(originalFolder, 'dir')
                mkdir(originalFolder);
            end
            movefile(fullfile(IQfiles_PI(i).folder, IQfiles_PI(i).name), originalFolder);
        end
    else
        fprintf('No PI files found. Skipping PI processing.\n');
    end
    
    % Completion message
    fprintf('Completed splitting operations.\n');
end
