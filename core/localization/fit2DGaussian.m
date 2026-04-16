function localizations = fit2DGaussian(filteredData, candidateBubbles, locParams, indent_prefix)
% fit2DGaussian - Sub-pixel localization of microbubbles using robust 2D Gaussian fitting.
%
% DESCRIPTION:
%   This function iterates through a provided list of candidate microbubble 
%   locations (identified as integer-pixel peaks) and refines their coordinates 
%   to sub-pixel accuracy. It achieves this by extracting a localized region of 
%   interest (ROI) around each peak and fitting a 2D Gaussian surface to the 
%   pixel intensity values using a robust Non-Linear Least Squares (NLLS) solver.
%
% METHODOLOGY:
%   1. ROI Extraction: For each candidate peak, a small window (defined by a 
%      specified radius) is extracted from the filtered ultrasound frame.
%   2. Pre-Processing: The data within the ROI is offset so the minimum value 
%      is zero, effectively handling background noise and negative signals.
%   3. Modeling: The extracted data is fit to a 2D Gaussian equation:
%      $Z = A \cdot \exp\left(-\left(\frac{(x-x_0)^2}{2\sigma_x^2} + \frac{(y-y_0)^2}{2\sigma_y^2}\right)\right) + B$
%   4. Quality Control (QC): Rigorous checks are applied during and after fitting. 
%      Fits are rejected based on boundary proximity, invalid data, solver 
%      failures, or physically impossible geometric parameters.
%
% INPUTS:
%   filteredData     - (H x W x T double) 3D matrix representing the filtered 
%                      ultrasound data (Height x Width x Frames).
%   candidateBubbles - (Table) Contains the initial integer locations and data 
%                      for candidates. Must include columns: 'Frame', 'X', 'Y', 'Intensity'.
%   locParams        - (Struct) Configuration parameters for the localization:
%                      .gauss_fit_box_radius (int): Half-width of the ROI window (e.g., 3 or 4).
%                      .min_r_squared (double): Minimum R-squared value for an acceptable fit (default: 0.5).
%                      .qc_max_shift_factor (double): Maximum allowed sub-pixel shift relative to the radius (default: 1.0).
%   indent_prefix    - (String, Optional) Text prefix used for formatting console output.
%
% OUTPUTS:
%   localizations    - (Table) Validated, sub-pixel coordinates for successful fits.
%                      Columns: ['Frame', 'X', 'Y', 'Intensity', 'Confidence']
%
% REJECTION / ERROR CODES:
%   1  - Edge Proximity: The candidate is too close to the image boundary to extract a full ROI.
%   10 - Invalid ROI: The extracted ROI contains NaN values or is completely flat (variance is zero).
%   20 - Solver Error: The Curve Fitting solver failed to converge, threw an exception, or returned an invalid R-squared.
%   31 - Bad Amplitude: The fitted amplitude (A) is near zero or physically impossible.
%   32 - Bad Sigma (Shape): The fitted Gaussian is either too narrow ($< 0.3$ pixels) or too wide (exceeds half the box width).
%   40 - Low R-Squared: The goodness-of-fit falls below the required threshold defined in `locParams`.
%   41 - Spatial Divergence: The fitted center ($x_0, y_0$) shifted too far from the initial integer peak estimate.
%
% DEPENDENCIES:
%   - Curve Fitting Toolbox (`fit`, `fittype`, `fitoptions`)
%
% AUTHOR: Grigori Shapiro

    % --- 1. Initialization & Validation ---
    if nargin < 4, indent_prefix = ''; end
    
    if isempty(candidateBubbles)
        localizations = table([], [], [], [], [], 'VariableNames', {'Frame', 'X', 'Y', 'Intensity', 'Confidence'});
        fprintf('%sWarning: No candidates provided.\n', indent_prefix);
        return;
    end
    
    if ~license('test', 'Curve_Fitting_Toolbox')
        error('Error: Curve Fitting Toolbox is required.');
    end
    
    num_candidates = height(candidateBubbles);
    results_matrix = NaN(num_candidates, 5); % [Frame, X, Y, Amp, R2]
    rejection_codes = zeros(num_candidates, 1);
    
    % --- 2. Parameter Setup ---
    radius = locParams.gauss_fit_box_radius;
    box_width = 2 * radius + 1;
    
    % Default safeguards
    min_r_sq = 0.5;
    if isfield(locParams, 'min_r_squared'), min_r_sq = locParams.min_r_squared; end
    
    max_shift = 1.0 * radius;
    if isfield(locParams, 'qc_max_shift_factor'), max_shift = locParams.qc_max_shift_factor * radius; end
    
    % Define Model (Corrected with explicit coefficient order)
    gauss_eq = 'a * exp(-((x-x0)^2/(2*sx^2) + (y-y0)^2/(2*sy^2))) + b';
    gauss_model = fittype(gauss_eq, ...
        'independent', {'x', 'y'}, 'dependent', 'z', ...
        'coefficients', {'a', 'x0', 'y0', 'sx', 'sy', 'b'}); % Crucial for correct mapping
    
    [H, W, ~] = size(filteredData);
    
    % Fitting Options
    opts = fitoptions(gauss_model);
    opts.Algorithm = 'Trust-Region';
    opts.Display = 'off';
    opts.MaxIter = 400;
    
    % Suppress specific warnings for cleaner loop
    warning('off', 'curvefit:fit:noStartPoint');
    warning('off', 'curvefit:fit:iterationLimitReached');
    
    % --- 3. Main Processing Loop ---
    for i = 1:num_candidates
        
        % Progress Update
        if mod(i, 500) == 0
            fprintf('%sFitting candidate %d / %d...\n', indent_prefix, i, num_candidates);
        end
        
        frm_idx = candidateBubbles.Frame(i);
        x_int = round(candidateBubbles.X(i));
        y_int = round(candidateBubbles.Y(i));
        
        % QC 1: Edge Check
        if x_int - radius < 1 || x_int + radius > W || y_int - radius < 1 || y_int + radius > H
            rejection_codes(i) = 1; continue;
        end
        
        % Extract ROI
        img_slice = filteredData(:,:,frm_idx);
        if ~isreal(img_slice), img_slice = abs(img_slice); end % Handle IQ data
        
        roi = double(img_slice(y_int-radius:y_int+radius, x_int-radius:x_int+radius));
        
        % Normalize ROI (Shift to positive)
        roi_min = min(roi(:));
        if roi_min < 0, roi = roi - roi_min; end
        
        % QC 10: ROI Validity Check (Flat or NaN)
        roi_range = max(roi(:)) - min(roi(:));
        if roi_range < 1e-10 || any(isnan(roi(:)))
            rejection_codes(i) = 10; continue;
        end
        
        % Prepare Data for Fit
        [grid_x, grid_y] = meshgrid(-radius:radius, -radius:radius);
        [x_data, y_data, z_data] = prepareSurfaceData(grid_x, grid_y, roi);
        
        % Initial Guesses & Bounds
        % Order: {'a', 'x0', 'y0', 'sx', 'sy', 'b'}
        amp_guess = roi_range;
        bg_guess = min(roi(:));
        
        opts.StartPoint = [amp_guess, 0, 0, 1.0, 1.0, bg_guess];
        opts.Lower = [0, -max_shift, -max_shift, 0.2, 0.2, 0];
        opts.Upper = [roi_range*3, max_shift, max_shift, box_width/2, box_width/2, max(roi(:))];
        
        try
            [fit_res, gof] = fit([x_data, y_data], z_data, gauss_model, opts);
            
            % QC 20: Solver output validity
            if gof.rsquare < 0 || isnan(gof.rsquare)
                rejection_codes(i) = 20; continue; % Catastrophic fit failure
            end
            
            % QC 30-39: Parameter Sanity Checks
            % QC 31: Amplitude Sanity
            if fit_res.a < 1e-6 
                rejection_codes(i) = 31; continue;
            end
            
            % QC 32: Sigma Shape (Too skinny < 0.3px or Too fat > box/2)
            if fit_res.sx < 0.3 || fit_res.sy < 0.3 || fit_res.sx > box_width/2 || fit_res.sy > box_width/2
                rejection_codes(i) = 32; continue;
            end
            
            % QC 40-49: Quality Checks
            % QC 40: R-Squared Threshold
            if gof.rsquare < min_r_sq
                rejection_codes(i) = 40; continue;
            end
            
            % QC 41: Divergence Check (Center shifted too far)
            if abs(fit_res.x0) > max_shift || abs(fit_res.y0) > max_shift
                rejection_codes(i) = 41; continue;
            end
            
            % --- Success ---
            results_matrix(i, :) = [
                frm_idx, ...
                x_int + fit_res.x0, ...  % Absolute X
                y_int + fit_res.y0, ...  % Absolute Y
                fit_res.a, ...           % Amplitude
                gof.rsquare              % Confidence
            ];
            
        catch
            rejection_codes(i) = 20; % Solver threw an exception
            continue;
        end
    end
    
    warning('on', 'curvefit:fit:noStartPoint');
    warning('on', 'curvefit:fit:iterationLimitReached');
    
    % --- 4. Final Output Generation ---
    valid_mask = ~isnan(results_matrix(:,1));
    localizations = array2table(results_matrix(valid_mask,:), ...
        'VariableNames', {'Frame', 'X', 'Y', 'Intensity', 'Confidence'});
        
    % --- 5. Detailed QC Report ---
    print_qc_report(rejection_codes, num_candidates, indent_prefix);
end

function print_qc_report(codes, total, prefix)
    % Helper to print organized stats
    kept = sum(codes == 0);
    
    fprintf('\n%s=== 2D Gaussian Localization QC Report ===\n', prefix);
    fprintf('%sTotal Candidates:     %d\n', prefix, total);
    fprintf('%sSuccessful Fits:      %d (%.1f%%)\n', prefix, kept, 100*kept/total);
    fprintf('%s------------------------------------------\n', prefix);
    
    % Category 1: Pre-Fit Rejections
    c1 = sum(codes == 1);
    c10 = sum(codes == 10);
    if c1+c10 > 0
        fprintf('%s[Pre-Fit] Edge Proximity:      %d\n', prefix, c1);
        fprintf('%s[Pre-Fit] Invalid/Flat ROI:    %d\n', prefix, c10);
    end
    
    % Category 2: Solver Failures
    c20 = sum(codes == 20);
    if c20 > 0
        fprintf('%s[Solver]  Convergence Failed:    %d\n', prefix, c20);
    end
    
    % Category 3: Parameter Sanity (Physics check)
    c31 = sum(codes == 31);
    c32 = sum(codes == 32);
    if c31+c32 > 0
        fprintf('%s[Param]   Bad Amplitude:         %d\n', prefix, c31);
        fprintf('%s[Param]   Bad Sigma (Shape):     %d\n', prefix, c32);
    end
    
    % Category 4: Quality Metrics
    c40 = sum(codes == 40);
    c41 = sum(codes == 41);
    fprintf('%s[Quality] Low R-Squared:         %d\n', prefix, c40);
    fprintf('%s[Quality] Spatial Divergence:    %d\n', prefix, c41);
    fprintf('%s------------------------------------------\n', prefix);
end