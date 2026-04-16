function localizations = localizeRadialSymmetry(filteredData, candidateBubbles, locParams, indent_prefix)
% localizeRadialSymmetry - Sub-pixel localization using gradient-based radial symmetry.
%
% DESCRIPTION:
%   This function refines the integer-pixel coordinates of candidate microbubbles 
%   to sub-pixel accuracy. It uses a Radial Symmetry algorithm, which is highly 
%   efficient and robust for perfectly circular or symmetric objects. The algorithm 
%   includes integrated pre-processing and multiple Quality Control (QC) checks 
%   to ensure only high-confidence localizations are retained.
%
% METHODOLOGY:
%   1. Principle: For a radially symmetric object (like a microbubble Point 
%      Spread Function), the image gradient vectors at any point on its boundary 
%      should point directly towards or away from the center.
%   2. ROI Extraction: A Region of Interest (ROI) is extracted around each candidate.
%   3. Gradient Calculation: 2D image gradients are calculated and smoothed within the ROI.
%   4. Center Estimation: The algorithm analytically solves for a single sub-pixel 
%      point that minimizes the perpendicular distance to all lines defined by 
%      the gradient vectors (Weighted Least-Squares fit).
%   5. Quality Control: Candidates are rejected if they violate spatial boundaries, 
%      lack a clear single peak, or yield mathematically unstable gradient matrices.
%
% INPUTS:
%   filteredData     - (H x W x T double) 3D matrix of clutter-filtered ultrasound data.
%   candidateBubbles - (Table) Initial candidate peaks with columns: 'Frame', 'X', 'Y', 'Intensity'.
%   locParams        - (Struct) Parameters structure containing:
%                      .fwhm (1x2 double): Full-width half-max [lateral_x, axial_z] to define ROI.
%                      .enable_roi_maxima_check (logical): Toggle to check for multiple peaks.
%                      .qc_max_roi_maxima (int): Max allowed peaks in an ROI.
%                      .min_gradient_squared (double): Threshold for minimum gradient signal.
%                      .min_determinant (double): Threshold for matrix stability.
%                      .enable_divergence_check (logical): Toggle shift distance limits.
%                      .qc_max_shift_factor (double): Max allowed shift relative to FWHM.
%   indent_prefix    - (String, Optional) Text prefix for console output formatting.
%
% OUTPUTS:
%   localizations    - (Table) Validated, sub-pixel coordinates.
%                      Columns: ['Frame', 'X', 'Y', 'Intensity']
%
% REJECTION / ERROR CODES (QC Stats):
%   edge        - The candidate is too close to the image border to extract a full ROI.
%   roi_maxima  - The extracted ROI contains multiple local peaks (implies noise or overlapping bubbles).
%   gradient    - The sum of squared gradients in the ROI is too low (lacks defined shape).
%   determinant - The determinant of the least-squares matrix is near zero (mathematically unstable/unsolvable).
%   divergence  - The calculated sub-pixel shift is too far from the initial integer estimate.
%
% AUTHOR: Grigori Shapiro

    % --- 1. Initial Checks and Setup ---
    if nargin < 4, indent_prefix = ''; end % Set default indent if not provided

    % If no candidate localizations are provided, return an empty table
    % immediately to avoid errors.
    if isempty(candidateBubbles)
        localizations = table.empty(0, 4);
        localizations.Properties.VariableNames = {'Frame', 'X', 'Y', 'Intensity'};
        return;
    end
    
    num_candidates = height(candidateBubbles);
    
    % Initialize a structure to keep track of why candidates are rejected.
    % This is useful for debugging and parameter tuning.
    qc_stats.edge = 0;
    qc_stats.roi_maxima = 0;
    qc_stats.gradient = 0;
    qc_stats.determinant = 0;
    qc_stats.divergence = 0;
    
    % Initialize the output matrix with NaNs. Any candidate that fails a QC
    % check will retain its NaN value. This allows for easy filtering at the end.
    subpixel_coords = NaN(num_candidates, 2);
    
    % Extract key parameters from the locParams struct for readability.
    fwhm_x = locParams.fwhm(1);
    fwhm_z = locParams.fwhm(2);
    radius_x = round(fwhm_x / 2);
    radius_z = round(fwhm_z / 2);

    % --- 2. Main Localization Loop ---
    % Iterate through each candidate bubble detected in the previous step.
    for i = 1:num_candidates
        % Get the integer coordinates and frame index for the current candidate.
        frame_idx = candidateBubbles.Frame(i);
        x_int = round(candidateBubbles.X(i));
        y_int = round(candidateBubbles.Y(i));
        
        % Extract the relevant 2D image (frame) from the 3D data block.
        frame = filteredData(:,:,frame_idx);
        [H, W] = size(frame);
        
        % --- QC 1: Edge Check ---
        % Reject any candidate that is too close to the image border. This
        % prevents errors that can occur when the ROI is incomplete.
        if (y_int - radius_z < 1) || (y_int + radius_z > H) || ...
           (x_int - radius_x < 1) || (x_int + radius_x > W)
            qc_stats.edge = qc_stats.edge + 1;
            continue; % Skip to the next candidate.
        end
        
        % Define the Region of Interest (ROI): a small image patch
        % centered on the candidate.
        roi_z_range = y_int - radius_z : y_int + radius_z;
        roi_x_range = x_int - radius_x : x_int + radius_x;
        roi = frame(roi_z_range, roi_x_range);
        
        % --- QC 2: ROI Maxima Check ---
        % This check ensures the ROI is "clean". If there are multiple local
        % peaks within the ROI, it could be noise or overlapping bubbles,
        % which violates the radial symmetry assumption.
        if locParams.enable_roi_maxima_check
            if nnz(imregionalmax(roi)) > locParams.qc_max_roi_maxima
                qc_stats.roi_maxima = qc_stats.roi_maxima + 1;
                continue; % Skip if too many peaks are found.
            end
        end
        
        % --- 2b. Radial Symmetry Calculation ---
        [Nz, Nx] = size(roi);
        
        % Step 2b.1: Calculate derivatives along 45-degree shifted coordinates
        dIdu = roi(1:Nz-1, 2:Nx)   - roi(2:Nz, 1:Nx-1);
        dIdv = roi(1:Nz-1, 1:Nx-1) - roi(2:Nz, 2:Nx);
        
        % Step 2b.2: Smooth the gradients
        h = ones(3)/9;
        fdu = conv2(dIdu, h, 'same');
        fdv = conv2(dIdv, h, 'same');
        g_mag_sq = fdu.*fdu + fdv.*fdv;
        
        % --- QC 3: Gradient Check ---
        if sum(g_mag_sq(:)) < locParams.min_gradient_squared 
            qc_stats.gradient = qc_stats.gradient + 1;
            continue;
        end
        
        % Step 2b.3: Calculate slope 'm' and intercept 'b'
        % Define midpoint coordinates relative to the ROI center
        zm_onerow = (-(Nz-1)/2.0+0.5 : (Nz-1)/2.0-0.5)';
        zm = repmat(zm_onerow, 1, Nx-1);
        xm_onecol = (-(Nx-1)/2.0+0.5 : (Nx-1)/2.0-0.5);
        xm = repmat(xm_onecol, Nz-1, 1);
        
        % Slope of the gradient
        m = -(fdv + fdu) ./ (fdu - fdv);
        
        % Handle NaN/Inf values in slope, a robust step from PALA
        if any(isnan(m(:)))
            unsmooth_m = -(dIdv + dIdu) ./ (dIdu - dIdv);
            m(isnan(m)) = unsmooth_m(isnan(m));
        end
        m(isnan(m)) = 0; % Final fallback
        m(isinf(m)) = 1e6; % Replace Inf with a large number
        
        % z-intercept
        b = zm - m.*xm;
        
        % Step 2b.4: Implement Weighted Least-Squares fit
        sdI2 = sum(g_mag_sq(:));
        if sdI2 == 0 % Avoid division by zero if all gradients are zero
            qc_stats.gradient = qc_stats.gradient + 1;
            continue;
        end
        
        % Initial centroid guess (weighted by gradient magnitude)
        zcentroid = sum(sum(g_mag_sq .* zm)) / sdI2;
        xcentroid = sum(sum(g_mag_sq .* xm)) / sdI2;
        
        % Calculate weights: stronger weight for high-gradient, near-centroid points
        dist_from_centroid = sqrt((zm - zcentroid).^2 + (xm - xcentroid).^2);
        dist_from_centroid(dist_from_centroid < eps) = eps; % Avoid division by zero
        w = g_mag_sq ./ dist_from_centroid;
        
        % Step 2b.5: Solve for the center of symmetry (xc, zc)
        wm2p1 = w ./ (m.*m + 1);
        sw = sum(wm2p1(:));
        smmw = sum(sum(m.*m .* wm2p1));
        smw = sum(sum(m .* wm2p1));
        smbw = sum(sum(m.*b .* wm2p1));
        sbw = sum(sum(b .* wm2p1));
        
        detM = smw*smw - smmw*sw;
        
        % --- QC 4: Determinant Check ---
        if abs(detM) < locParams.min_determinant
            qc_stats.determinant = qc_stats.determinant + 1;
            continue;
        end
        
        % Solve for the sub-pixel SHIFT relative to the ROI center
        xc_shift = (smbw*sw - smw*sbw) / detM;
        zc_shift = (smbw*smw - smmw*sbw) / detM;
        
        % A more intuitive way to get the final coordinates:
        % The algorithm finds the center relative to the ROI's coordinate system,
        % where (0,0) is the center. We add this shift to the ROI's center in the
        % main image's coordinate system.
        center_coords = [(x_int + xc_shift), (y_int + zc_shift)];
        
        % --- QC 5: Divergence Check --- (Your original check, still perfectly valid)
        if locParams.enable_divergence_check
            % The calculated shifts are directly what we need to check
            max_shift_x = locParams.qc_max_shift_factor * (fwhm_x / 2);
            max_shift_y = locParams.qc_max_shift_factor * (fwhm_z / 2);
            if abs(xc_shift) > max_shift_x || abs(zc_shift) > max_shift_y
                qc_stats.divergence = qc_stats.divergence + 1;
                continue; % Skip this candidate if it shifted too far.
            end
        end
        
        % If the candidate passed all checks, store its calculated sub-pixel coordinates.
        subpixel_coords(i, :) = [center_coords(1), center_coords(2)]; % [X, Y]

    end
    
    % --- 3. Finalization ---
    
    % Create a logical index of all candidates that passed the QC checks (i.e., are not NaN).
    valid_indices = ~isnan(subpixel_coords(:, 1));
    
    % Filter the coordinate matrix and the original candidate table to keep only the valid ones.
    final_coords = subpixel_coords(valid_indices, :);
    final_candidates = candidateBubbles(valid_indices, :);
    
    % Create the final output table with the validated localizations.
    localizations = table(final_candidates.Frame, final_coords(:,1), final_coords(:,2), final_candidates.Intensity,...
        'VariableNames', {'Frame', 'X', 'Y', 'Intensity'});
        
    % --- 4. Report QC Statistics ---
    % Provide detailed feedback to the user on how many candidates were kept or
    % rejected at each stage of the quality control process.
    total_rejected = qc_stats.edge + qc_stats.roi_maxima + qc_stats.gradient + qc_stats.determinant + qc_stats.divergence;
    total_kept = height(localizations);
    
    fprintf('%s--- Localization QC Summary ---\n', indent_prefix);
    fprintf('%sTotal candidates processed: %d\n', indent_prefix, num_candidates);
    fprintf('%s---------------------------------\n', indent_prefix);
    fprintf('%sRejected by Edge Check:       %6d (%5.1f%%)\n', indent_prefix, qc_stats.edge, 100*qc_stats.edge/num_candidates);
    fprintf('%sRejected by ROI Maxima Check: %6d (%5.1f%%)\n', indent_prefix, qc_stats.roi_maxima, 100*qc_stats.roi_maxima/num_candidates);
    fprintf('%sRejected by Gradient Check:   %6d (%5.1f%%)\n', indent_prefix, qc_stats.gradient, 100*qc_stats.gradient/num_candidates);
    fprintf('%sRejected by Determinant Check:%6d (%5.1f%%)\n', indent_prefix, qc_stats.determinant, 100*qc_stats.determinant/num_candidates);
    fprintf('%sRejected by Divergence Check: %6d (%5.1f%%)\n', indent_prefix, qc_stats.divergence, 100*qc_stats.divergence/num_candidates);
    fprintf('%s---------------------------------\n', indent_prefix);
    fprintf('%sTotal localizations kept:     %6d (%5.1f%%)\n', indent_prefix, total_kept, 100*total_kept/num_candidates);
    fprintf('%sTotal localizations rejected: %6d (%5.1f%%)\n', indent_prefix, total_rejected, 100*total_rejected/num_candidates);
    fprintf('%s---------------------------------\n', indent_prefix);
end
