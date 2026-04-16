function [densityMap, velocityMap, velocityCountMap] = renderHistogram(tracks, sRes_dims, renderParams)
% renderHistogram - Fast super-resolution image reconstruction using 2D histogram binning.
%
% DESCRIPTION:
%   This function reconstructs high-resolution density and velocity maps by 
%   quantizing (binning) interpolated microbubble tracks into a discrete grid. 
%   It is highly computationally efficient compared to Gaussian rendering, making 
%   it ideal for rapid prototyping, debugging, and processing massive datasets 
%   where rendering speed is prioritized over sub-pixel smoothness.
%
% METHODOLOGY:
%   1. Upsampling & Quantization: Sub-pixel coordinates from the dense tracks 
%      are scaled up by the `upsampling_factor` and then rounded (`round()`) 
%      to the nearest integer index on the super-resolution grid.
%   2. Redundancy Filtering: Because tracks are densely interpolated, multiple 
%      points in a single track might fall into the exact same integer pixel. 
%      The function uses `unique()` to ensure a track only contributes once to 
%      a specific pixel during a single pass, preventing artificial intensity inflation.
%   3. Velocity Matching: For each unique pixel a track occupies, the algorithm 
%      finds the closest original sub-pixel point and assigns its velocity value 
%      to that pixel.
%   4. Accumulation: The `densityMap` and `velocityCountMap` are incremented 
%      by exactly 1 for each occupied bin. The corresponding velocity is added 
%      to the `velocityMap`.
%
% INPUTS:
%   tracks       - (Struct Array) Filtered and densely interpolated track data. 
%                  Must contain fields: .path (Nx2 coordinates) and 
%                  .average_velocity_mm_s (Nx1 velocity vectors).
%   sRes_dims    - (1x2 double) Target dimensions of the super-resolution grid [Height, Width].
%   renderParams - (Struct) Configuration parameters:
%                  .upsampling_factor (double): The scale factor converting original pixels to SR pixels.
%
% OUTPUTS:
%   densityMap       - (H x W double) Super-resolved map representing raw microbubble counts per pixel.
%   velocityMap      - (H x W double) Super-resolved map of accumulated velocities.
%                      (Note: Must be divided by velocityCountMap later to get mean velocity).
%   velocityCountMap - (H x W double) Map tracking how many discrete velocity measurements fell into each pixel.
%
% NOTE ON USAGE:
%   While much faster, this method introduces quantization errors by discarding 
%   the decimal portion of the sub-pixel coordinates. For final publication 
%   images requiring maximum structural fidelity, `renderGaussian` is recommended.
%
% AUTHOR: Grigori Shapiro

    % Initialize the output maps
    densityMap = zeros(sRes_dims);
    velocityMap = zeros(sRes_dims);
    velocityCountMap = zeros(sRes_dims);
    
    % If there are no tracks to render, return the empty maps
    if isempty(tracks)
        return;
    end
    
    % --- Main Rendering Loop ---
    % This function now assumes 'tracks' are already filtered by length and densely interpolated.
    for i = 1:length(tracks)
        track = tracks(i);
        
        % Directly use the pre-interpolated path from the track struct.
        % Convert path from (X,Y) to (row,col) for matrix indexing and scale up.
        path_sRes = [track.path(:,2), track.path(:,1)] * renderParams.upsampling_factor;
        all_indices = round(path_sRes);
        
        unique_indices = unique(all_indices, 'rows'); % Get unique [row, col] indices

        % Use the pre-calculated average velocity vector.
        velocities_to_render = track.average_velocity_mm_s;
        
        % Loop through EVERY dense point in the track path
        for j = 1:size(unique_indices, 1)
            y_idx = unique_indices(j, 1); % Row index
            x_idx = unique_indices(j, 2); % Column index
            
            % To get the correct velocity, find which point in the original dense
            % path is closest to this unique pixel and use its velocity.
            [~, min_idx] = min(sum((all_indices - unique_indices(j,:)).^2, 2));
            velocity_value = velocities_to_render(min_idx);
    
            % Boundary check 
            if y_idx > 0 && y_idx <= sRes_dims(1) && x_idx > 0 && x_idx <= sRes_dims(2)
                densityMap(y_idx, x_idx) = densityMap(y_idx, x_idx) + 1;
                velocityMap(y_idx, x_idx) = velocityMap(y_idx, x_idx) + velocity_value;
                velocityCountMap(y_idx, x_idx) = velocityCountMap(y_idx, x_idx) + 1;
            end
        end
    end
end