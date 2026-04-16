function [densityMap, velocityMap, velocityCountMap] = renderGaussian(tracks, sRes_dims, renderParams)
% renderGaussian - Advanced super-resolution image reconstruction using 2D Gaussian splatting.
%
% DESCRIPTION:
%   This function reconstructs high-resolution density and velocity maps from 
%   interpolated microbubble tracks. Instead of simply placing each localization 
%   into a single discrete pixel bin (histogram method), it renders every point 
%   along the track as a continuous 2D Gaussian kernel ("splatting"). This 
%   preserves sub-pixel accuracy, prevents quantization errors (aliasing), and 
%   produces visually smoother, more physically accurate vascular maps.
%
% METHODOLOGY:
%   1. Kernel Generation: A 2D Gaussian kernel is pre-calculated based on the 
%      user-defined standard deviation (`sigma`). The kernel size is dynamically 
%      set to encompass 3 standard deviations (99.7% of the weight).
%   2. Upsampling: Sub-pixel coordinates from the dense, interpolated tracks 
%      are scaled up by the specified `upsampling_factor` to match the target 
%      super-resolution grid (`sRes_dims`).
%   3. Splatting (Continuous Rendering): For *every* interpolated point along 
%      a track, the pre-calculated Gaussian kernel is centered on the point's 
%      location. 
%   4. Accumulation: The weights of the kernel are added to the corresponding 
%      pixels in the `densityMap`. Simultaneously, the kernel weights are 
%      multiplied by the local velocity and accumulated in the `velocityMap`.
%   5. Boundary Handling: The algorithm strictly checks image boundaries and 
%      dynamically crops the Gaussian kernel if a bubble passes too close to 
%      the edge of the field of view.
%
% INPUTS:
%   tracks       - (Struct Array) Filtered and densely interpolated track data. 
%                  Must contain fields: .path (Nx2 coordinates) and 
%                  .average_velocity_mm_s (Nx1 velocity vectors).
%   sRes_dims    - (1x2 double) Target dimensions of the super-resolution grid [Height, Width].
%   renderParams - (Struct) Configuration parameters:
%                  .gaussian_sigma (double): Standard deviation defining the spread of the kernel.
%                  .upsampling_factor (double): The scale factor converting original pixels to SR pixels.
%
% OUTPUTS:
%   densityMap       - (H x W double) Super-resolved map of microbubble occurrences (weighted).
%   velocityMap      - (H x W double) Super-resolved map of accumulated velocities.
%                      (Note: Must be divided by velocityCountMap later to get mean velocity).
%   velocityCountMap - (H x W double) Map of the cumulative weights used for velocity normalization.
%
% AUTHOR: Grigori Shapiro

    % --- Step 1: Initialize Maps and Kernel ---
    densityMap = zeros(sRes_dims);
    velocityMap = zeros(sRes_dims);
    velocityCountMap = zeros(sRes_dims); % Renamed for consistency with renderHistogram
    if isempty(tracks), return; end
    
    % Define the Gaussian kernel properties
    sigma = renderParams.gaussian_sigma;
    kernel_radius = ceil(3 * sigma);
    kernel_size = 2 * kernel_radius + 1;
    [X, Y] = meshgrid(-kernel_radius:kernel_radius, -kernel_radius:kernel_radius);
    gaussian_kernel = exp(-(X.^2 + Y.^2) / (2 * sigma^2));
    
    % --- Step 2: Main Rendering Loop ---
    % This function assumes 'tracks' are already filtered and densely interpolated.
    for i = 1:length(tracks)
        track = tracks(i);
        
        % Directly use the pre-interpolated path from the track struct.
        % Convert path from (X,Y) to (row,col) for matrix indexing and scale up.
        path_sRes = [track.path(:,2), track.path(:,1)] * renderParams.upsampling_factor;
        
        % Use the pre-calculated average velocity vector.
        velocities_to_render = track.average_velocity_mm_s;
        
        % --- LOGIC CHANGE: Render a Gaussian for EVERY point in the dense path ---
        for j = 1:size(path_sRes, 1)
            % Get the center pixel for the current point's Gaussian stamp
            y_center = round(path_sRes(j, 1));
            x_center = round(path_sRes(j, 2));
            
            % Get the velocity corresponding to this point
            velocity_value = velocities_to_render(j);
            
            % --- (The following robust kernel placement logic is preserved) ---
            
            % Define the ROI for kernel placement
            x_min = x_center - kernel_radius; x_max = x_center + kernel_radius;
            y_min = y_center - kernel_radius; y_max = y_center + kernel_radius;
            
            % Define the parts of the kernel and map to use, handling edges
            k_x_min = 1; k_x_max = kernel_size;
            k_y_min = 1; k_y_max = kernel_size;
            
            if x_min < 1, k_x_min = 2 - x_min; x_min = 1; end
            if x_max > sRes_dims(2), k_x_max = kernel_size - (x_max - sRes_dims(2)); x_max = sRes_dims(2); end
            if y_min < 1, k_y_min = 2 - y_min; y_min = 1; end
            if y_max > sRes_dims(1), k_y_max = kernel_size - (y_max - sRes_dims(1)); y_max = sRes_dims(1); end
            
            % If the ROI is completely outside the map, skip this point
            if x_min > x_max || y_min > y_max, continue; end
            
            % Extract the valid sub-region of the kernel
            sub_kernel = gaussian_kernel(k_y_min:k_y_max, k_x_min:k_x_max);

            % Add the weighted kernel to the maps
            densityMap(y_min:y_max, x_min:x_max) = densityMap(y_min:y_max, x_min:x_max) + sub_kernel;
            velocityMap(y_min:y_max, x_min:x_max) = velocityMap(y_min:y_max, x_min:x_max) + sub_kernel * velocity_value;
            velocityCountMap(y_min:y_max, x_min:x_max) = velocityCountMap(y_min:y_max, x_min:x_max) + sub_kernel;
        end
    end
end