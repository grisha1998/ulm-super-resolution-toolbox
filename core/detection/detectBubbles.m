function candidateBubbles = detectBubbles(filteredData, locParams, roiMask)
% detectBubbles - Basic intensity-based detection of microbubble candidates.
%
% DESCRIPTION:
%   This function serves as the initial detection stage for Ultrasound 
%   Localization Microscopy (ULM). It scans the clutter-filtered data to find 
%   local intensity peaks (regional maxima) that potentially represent microbubbles. 
%   It acts as a fast, first-pass filter before applying computationally heavier 
%   sub-pixel localization algorithms.
%
% METHODOLOGY:
%   1. Masking: Applies an optional user-defined spatial mask to restrict 
%      detection to specific regions (e.g., inside vessels), zeroing out noise.
%   2. Normalization: Normalizes each frame's intensity.
%   3. Maxima Detection: Uses 8-connected neighborhood comparisons (`imregionalmax`) 
%      to identify local peaks.
%   4. Thresholding: Discards peaks falling below a user-defined intensity threshold.
%   5. Sorting & Limiting: Sorts the surviving candidates by intensity and retains 
%      only the top N strongest peaks per frame to prevent computational overload.
%
% INPUTS:
%   filteredData - (H x W x T double) 3D matrix of clutter-filtered ultrasound data.
%   locParams    - (Struct) Parameters structure containing:
%                  .detection_threshold (double): Minimum normalized intensity (0 to 1).
%                  .max_bubbles_per_frame (int): Maximum number of candidates to keep per frame.
%   roiMask      - (H x W logical, Optional) Binary spatial mask to confine detection.
%
% OUTPUTS:
%   candidateBubbles - (Table) Integer-pixel locations of detected candidates.
%                      Columns: ['Frame', 'X', 'Y', 'Intensity']
%
% REJECTION CRITERIA (Filtering):
%   - Mask Rejection: Peak falls outside the defined `roiMask`.
%   - Threshold Rejection: Peak intensity is lower than `locParams.detection_threshold`.
%   - Limit Rejection: Peak is valid but not within the top `max_bubbles_per_frame` strongest signals.
%
% AUTHOR: Grigori Shapiro

    if nargin < 3 || isempty(roiMask)
            roiMask = true(size(filteredData, 1), size(filteredData, 2));
    end

    [H, W, T] = size(filteredData);
    
    % Pre-allocate a cell array to store results from each frame
    frame_results = cell(T, 1);
    
    % Process each frame in parallel to speed up detection
    parfor t = 1:T
        frame = filteredData(:,:,t);
        
        % 1. Apply ROI Mask (Zero out noise areas)
        % This ensures no peaks are detected outside the vessels
        frame = frame .* roiMask;
        
        % Normalize frame for consistent thresholding
        frame_max = max(frame(:));
        if frame_max > 0
            frame_norm = frame / frame_max;
        else
            frame_norm = frame;
        end
        
        % Smooth frame to suppress multi-peak artifacts on broad PSFs.
        % Sigma is derived from FWHM; for small PSFs (<=3px) this is effectively a no-op.
        sigma_x = locParams.fwhm(1) / 2.355;
        sigma_z = locParams.fwhm(2) / 2.355;
        frame_detect = imgaussfilt(frame_norm, [sigma_z, sigma_x]);

        % Find regional maxima on the smoothed frame (scale-aware peak locations),
        % but threshold on the ORIGINAL frame_norm (accurate intensities).
        regional_maxima = imregionalmax(frame_detect, 8);
        thresholded_peaks = regional_maxima & (frame_norm > locParams.detection_threshold);
        
        % Get the linear indices, coordinates, and intensity of the detected peaks
        peak_indices = find(thresholded_peaks);
        if isempty(peak_indices)
            continue; % No peaks found in this frame
        end
        
        [peak_y, peak_x] = ind2sub(size(frame), peak_indices);
        peak_intensities = frame(peak_indices);
        
        % Sort peaks by intensity in descending order
        [sorted_intensities, sort_idx] = sort(peak_intensities, 'descend');
        sorted_x = peak_x(sort_idx);
        sorted_y = peak_y(sort_idx);
        
        % Keep only the top N peaks as specified by the user
        num_peaks = min(length(sorted_intensities), locParams.max_bubbles_per_frame);
        
        % Store results for this frame
        frame_results{t} = [repmat(t, num_peaks, 1), sorted_x(1:num_peaks), sorted_y(1:num_peaks), sorted_intensities(1:num_peaks)];
    end
    
    % Concatenate results from all frames into a single matrix
    all_candidates_matrix = cat(1, frame_results{:});
    
    % Convert the final matrix to a MATLAB table for clarity and ease of use
    if isempty(all_candidates_matrix)
        candidateBubbles = table('Size', [0 4], 'VariableTypes', {'double', 'double', 'double', 'double'},...
                               'VariableNames', {'Frame', 'X', 'Y', 'Intensity'});
    else
        candidateBubbles = array2table(all_candidates_matrix,...
            'VariableNames', {'Frame', 'X', 'Y', 'Intensity'});
    end
    %plotAllDetections(candidateBubbles, filteredData);
end


function plotAllDetections(candidateBubbles, filteredData)
% FILENAME: plotAllDetections.m
%
% PURPOSE:
%   Visualizes all detected bubble candidates from all frames on a single
%   background image. The background is an AVERAGE of all frames in the
%   data sequence.
%
% SYNTAX:
%   plotAllDetections(candidateBubbles, filteredData)
%
% INPUTS:
%   candidateBubbles: (table) The output table from detectBubbles.m.
%   filteredData:     (H x W x T, double) The original filtered data, used
%                     to generate the background image.
%    
    if isempty(candidateBubbles)
        disp('No bubbles were detected. Nothing to plot.');
        return;
    end

    % 1. Create a background image by AVERAGING all frames
    % This collapses the time dimension (3rd dim) by taking the mean value.
    backgroundImage = mean(filteredData, 3);
    
    % 2. Display the background image in a new figure
    figure;
    imagesc(backgroundImage);
    colormap('gray');
    axis image;
    hold on;
    
    % 3. Plot all detected bubble locations
    plot(candidateBubbles.X, candidateBubbles.Y, 'r.', 'MarkerSize', 2, 'LineWidth', 1.5);
    
    % 4. Add labels and a title for context
    title('All Detected Candidates Superimposed on Average-Intensity Image');
    xlabel('X Coordinate (pixels)');
    ylabel('Y Coordinate (pixels)');
    legend('Detected Bubbles');
    
    hold off;
    
end