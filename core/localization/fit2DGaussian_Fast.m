function localizations = fit2DGaussian_Fast(filteredData, candidateBubbles, locParams, indent_prefix)
% fit2DGaussian_Fast - High-speed sub-pixel localization of microbubbles using optimized 2D Gaussian fitting.
%
% DESCRIPTION:
%   This function serves as a highly optimized alternative to `fit2DGaussian`. 
%   It determines the sub-pixel coordinates of candidate microbubbles by fitting 
%   a 2D Gaussian surface to the image data. By utilizing `lsqcurvefit` alongside 
%   parallel computing (`parfor`), this method achieves localization significantly 
%   faster (typically 10x-50x) while maintaining strict quality control standards.
%
% METHODOLOGY:
%   1. Grid Pre-calculation: A static coordinate grid is generated once prior 
%      to the loop to minimize redundant memory allocation.
%   2. Parallel Processing: The list of candidate bubbles is distributed across 
%      available parallel workers using `parfor`.
%   3. Optimization: An optimized anonymous function handle evaluates the 2D 
%      Gaussian model. `lsqcurvefit` solves the non-linear curve-fitting problem 
%      with strict upper and lower bounds to constrain the search space.
%   4. Quality Control (QC): Candidates are aggressively filtered pre- and 
%      post-fit to ensure only highly confident, physically realistic localizations 
%      are retained.
%
% INPUTS:
%   filteredData     - (H x W x T double) 3D matrix representing the filtered 
%                      ultrasound data (Height x Width x Frames).
%   candidateBubbles - (Table) Contains the initial integer locations and data 
%                      for candidates. Must include columns: 'Frame', 'X', 'Y'.
%   locParams        - (Struct) Configuration parameters for the localization:
%                      .gauss_fit_box_radius (int): Half-width of the ROI window.
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
%   10 - Invalid ROI: The extracted ROI contains NaN values or lacks meaningful intensity variation.
%   20 - Solver Failed: `lsqcurvefit` failed to converge or encountered an exception during the fitting process.
%   30 - Bad Sigma (Shape): The fitted widths ($\sigma_x, \sigma_y$) hit the imposed boundary limits, indicating an unphysical shape (either too narrow or too wide).
%   40 - Low R-Squared: The calculated R-squared value of the fit is below the required minimum threshold.
%   41 - Spatial Divergence: The fitted center shifted further from the initial integer peak than allowed by the maximum shift parameter.
%
% DEPENDENCIES:
%   - Optimization Toolbox (`lsqcurvefit`, `optimoptions`)
%   - Parallel Computing Toolbox (for `parfor` execution)
%
% AUTHOR: Grigori Shapiro


    % --- 1. Initialization ---
    if nargin < 4, indent_prefix = ''; end
    
    if isempty(candidateBubbles)
        localizations = table([], [], [], [], [], 'VariableNames', {'Frame', 'X', 'Y', 'Intensity', 'Confidence'});
        fprintf('%sWarning: No candidates provided.\n', indent_prefix);
        return;
    end
    
    if ~license('test', 'optimization_toolbox')
        error('The fast localization method requires the Optimization Toolbox (lsqcurvefit).');
    end
    
    num_candidates = height(candidateBubbles);
    
    % Prepare Output Arrays (Sliced for parfor efficiency)
    out_frame = candidateBubbles.Frame;
    out_x = zeros(num_candidates, 1);
    out_y = zeros(num_candidates, 1);
    out_int = zeros(num_candidates, 1);
    out_conf = zeros(num_candidates, 1);
    
    % Exit Codes Vector: We will store the rejection reason here
    % 0=Success, 1=Edge, 10=ROI, 20=Solver, 30=Params, 40=LowR2, 41=Divergence
    exit_codes = zeros(num_candidates, 1); 
    
    % --- 2. Parameters ---
    radius = locParams.gauss_fit_box_radius;
    box_width = 2 * radius + 1;
    
    min_r_sq = 0.5;
    if isfield(locParams, 'min_r_squared'), min_r_sq = locParams.min_r_squared; end
    
    max_shift = 1.0 * radius;
    if isfield(locParams, 'qc_max_shift_factor'), max_shift = locParams.qc_max_shift_factor * radius; end
    
    [H, W, ~] = size(filteredData);
    
    % --- 3. Pre-calculate Grid ---
    [grid_x_mat, grid_y_mat] = meshgrid(-radius:radius, -radius:radius);
    x_data = [grid_x_mat(:), grid_y_mat(:)]; % Constant for all fits
    
    % --- 4. Define Gaussian Model (Optimized Function Handle) ---
    % p(1)=Amp, p(2)=x0, p(3)=y0, p(4)=sx, p(5)=sy, p(6)=Bg
    gauss_func = @(p, xy) p(1) * exp( -((xy(:,1)-p(2)).^2/(2*p(4).^2) + (xy(:,2)-p(3)).^2/(2*p(5).^2)) ) + p(6);

    % Solver Options
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxIterations', 50, ...    % Fast fail
        'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6, ...
        'UseParallel', false);      % Parallel is handled by our loop
        
    % --- 5. Main Parallel Loop ---
    fprintf('%sStarting parallel processing of %d candidates...\n', indent_prefix, num_candidates);
    
    parfor i = 1:num_candidates
        
        frm_idx = out_frame(i);
        x_int = round(candidateBubbles.X(i));
        y_int = round(candidateBubbles.Y(i));
        
        % QC 1: Edge Check
        if x_int - radius < 1 || x_int + radius > W || y_int - radius < 1 || y_int + radius > H
            exit_codes(i) = 1; continue;
        end
        
        % Extract ROI
        img_slice = filteredData(:,:,frm_idx);
        if ~isreal(img_slice), img_slice = abs(img_slice); end 
        
        roi = double(img_slice(y_int-radius:y_int+radius, x_int-radius:x_int+radius));
        
        % Normalize ROI
        roi_min = min(roi(:));
        if roi_min < 0, roi = roi - roi_min; end
        roi_range = max(roi(:)) - min(roi(:));
        
        % QC 10: Invalid/Flat ROI
        if roi_range < 1e-10 || any(isnan(roi(:)))
            exit_codes(i) = 10; continue;
        end
        
        z_data = roi(:);
        
        % Setup Solver
        % p = [Amp, x0, y0, sx, sy, Bg]
        p0 = [roi_range, 0, 0, 1.0, 1.0, roi_min];
        
        % Bounds
        lb = [0,          -max_shift, -max_shift, 0.2,            0.2,            0];
        ub = [roi_range*3, max_shift,  max_shift,  box_width/2,    box_width/2,    max(roi(:))];
        
        try
            [p_fit, ~, residual, exitflag] = lsqcurvefit(gauss_func, p0, x_data, z_data, lb, ub, opts);
            
            % QC 20: Solver Convergence Failure
            if exitflag <= 0
                exit_codes(i) = 20; continue;
            end
            
            % QC 30: Bad Sigma (Hit bounds?)
            % Note: lsqcurvefit respects bounds, but if it hits the edge (e.g. 0.2), it's a bad fit
            if p_fit(4) <= 0.21 || p_fit(5) <= 0.21 || p_fit(4) >= (box_width/2)-0.1 || p_fit(5) >= (box_width/2)-0.1
                 exit_codes(i) = 30; continue;
            end

            % Calculate R-Squared
            sst = sum((z_data - mean(z_data)).^2);
            sse = sum(residual.^2);
            r_sq = 1 - (sse / sst);
            
            % QC 40: Low R-Squared
            if r_sq < min_r_sq
                exit_codes(i) = 40; continue;
            end
            
            % QC 41: Spatial Divergence
            if abs(p_fit(2)) > max_shift || abs(p_fit(3)) > max_shift
                exit_codes(i) = 41; continue;
            end
            
            % Store Success
            out_x(i) = x_int + p_fit(2);
            out_y(i) = y_int + p_fit(3);
            out_int(i) = p_fit(1);
            out_conf(i) = r_sq;
            exit_codes(i) = 0; % Explicit success
            
        catch
            exit_codes(i) = 20; % Treat exception as solver fail
            continue;
        end
    end
    
    % --- 6. Report Generation ---
    valid_mask = exit_codes == 0;
    
    localizations = table(out_frame(valid_mask), out_x(valid_mask), out_y(valid_mask), ...
                          out_int(valid_mask), out_conf(valid_mask), ...
                          'VariableNames', {'Frame', 'X', 'Y', 'Intensity', 'Confidence'});

    % Print the detailed report using the tracked exit codes
    print_fast_qc_report(exit_codes, num_candidates, indent_prefix);
end

function print_fast_qc_report(codes, total, prefix)
    kept = sum(codes == 0);
    fprintf('\n%s=== Fast Gaussian QC Report ===\n', prefix);
    fprintf('%sTotal Candidates:     %d\n', prefix, total);
    fprintf('%sSuccessful Fits:      %d (%.1f%%)\n', prefix, kept, 100*kept/total);
    fprintf('%s-------------------------------\n', prefix);
    
    fprintf('%s[1]  Edge Proximity:    %d\n', prefix, sum(codes == 1));
    fprintf('%s[10] Invalid ROI:       %d\n', prefix, sum(codes == 10));
    fprintf('%s[20] Solver Failed:     %d\n', prefix, sum(codes == 20));
    fprintf('%s[30] Bad Sigma (Shape): %d\n', prefix, sum(codes == 30));
    fprintf('%s[40] Low R-Squared:     %d\n', prefix, sum(codes == 40));
    fprintf('%s[41] Divergence:        %d\n', prefix, sum(codes == 41));
    fprintf('%s-------------------------------\n', prefix);
end