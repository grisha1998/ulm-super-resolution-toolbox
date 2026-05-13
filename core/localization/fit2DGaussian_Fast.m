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
%   The fitting window and sigma bounds are derived from the FWHM parameter,
%   which can specify different PSF widths for the lateral (X) and axial (Z)
%   axes. This ensures that the Gaussian fitting is physically consistent
%   with the expected Point Spread Function geometry.
%
% METHODOLOGY:
%   1. FWHM-Driven Setup: The ROI radius, expected sigma values, and solver
%      bounds are all derived from the user-specified FWHM [X, Z]. This
%      ensures physical consistency between the PSF model and the fitting
%      constraints.
%   2. Grid Pre-calculation: A static coordinate grid is generated once prior 
%      to the loop to minimize redundant memory allocation.
%   3. Parallel Processing: The list of candidate bubbles is distributed across 
%      available parallel workers using `parfor`.
%   4. Data-Driven Initial Guess: The initial sigma estimate is computed from
%      the second moment of the ROI intensity, then clamped to a physically
%      plausible range derived from FWHM.
%   5. Quality Control (QC): Candidates are aggressively filtered pre- and 
%      post-fit to ensure only highly confident, physically realistic localizations 
%      are retained.
%
% INPUTS:
%   filteredData     - (H x W x T double) 3D matrix representing the filtered 
%                      ultrasound data (Height x Width x Frames).
%   candidateBubbles - (Table) Contains the initial integer locations and data 
%                      for candidates. Must include columns: 'Frame', 'X', 'Y'.
%   locParams        - (Struct) Configuration parameters for the localization:
%                      .fwhm (1x2 double): PSF full-width at half-max [lateral_X, axial_Z] in pixels.
%                                          Drives the ROI size and all sigma bounds.
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
%   30 - Bad Sigma (Shape): The fitted widths (sigma_x, sigma_y) deviate too far from the FWHM-expected
%        sigma values, indicating an unphysical shape.
%   40 - Low R-Squared: The calculated R-squared value of the fit is below the required minimum threshold.
%   41 - Spatial Divergence: The fitted center shifted further from the initial integer peak than allowed
%        by the maximum shift parameter.
%
% DEPENDENCIES:
%   - Optimization Toolbox (`lsqcurvefit`, `optimoptions`)
%   - Parallel Computing Toolbox (for `parfor` execution)
%
% AUTHOR: Grigori Shapiro


    % --- 1. Initialization ---
    if nargin < 4, indent_prefix = ''; end
    filteredData = abs(filteredData);
    
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
    
    % Exit Codes Vector
    % 0=Success, 1=Edge, 10=ROI, 20=Solver, 30=BadSigma, 40=LowR2, 41=Divergence
    exit_codes = zeros(num_candidates, 1); 
    
    % --- 2. FWHM-Driven Parameter Setup ---
    % Read FWHM [lateral_X, axial_Z] - the master PSF parameter
    if isfield(locParams, 'fwhm') && ~isempty(locParams.fwhm) && ~any(isnan(locParams.fwhm))
        fwhm = locParams.fwhm(:)';  % Force row vector [fwhm_x, fwhm_z]
        if isscalar(fwhm), fwhm = [fwhm, fwhm]; end
    else
        fwhm = [3, 3]; % Fallback
        fprintf('%sWarning: FWHM not set. Using fallback [3, 3] px.\n', indent_prefix);
    end
    
    % Convert FWHM to expected sigma: FWHM = 2*sqrt(2*ln(2)) * sigma ~ 2.355 * sigma
    sigma_expected_x = fwhm(1) / 2.355;
    sigma_expected_y = fwhm(2) / 2.355;
    
    % Derive ROI radius from FWHM (use the larger axis to ensure full coverage)
    % The ROI should extend ~1.5x the FWHM half-width to capture the Gaussian tails
    fwhm_derived_radius = ceil(1.5 * max(fwhm) / 2);
    
    radius = fwhm_derived_radius;
    box_width = 2 * radius + 1;
    
    % Print derived configuration
    fprintf('%sFWHM = [%.1f, %.1f] px  =>  sigma_expected = [%.2f, %.2f] px,  ROI radius = %d px (%dx%d window)\n', ...
        indent_prefix, fwhm(1), fwhm(2), sigma_expected_x, sigma_expected_y, radius, box_width, box_width);

    % R-squared threshold
    min_r_sq = 0.5;
    if isfield(locParams, 'min_r_squared'), min_r_sq = locParams.min_r_squared; end
    
    % Maximum allowed sub-pixel shift from integer peak
    max_shift = 1.0 * radius;
    if isfield(locParams, 'qc_max_shift_factor'), max_shift = locParams.qc_max_shift_factor * radius; end
    
    [H, W, ~] = size(filteredData);
    
    % --- 3. FWHM-Derived Sigma Bounds ---
    % The solver bounds and QC thresholds are now anchored to the expected sigma.
    % This prevents the solver from converging to physically impossible shapes.
    %
    % Lower bound: 30% of expected sigma (sharper than expected but possible)
    % Upper bound: 300% of expected sigma (broader, e.g. overlapping PSFs)
    sigma_lb_x = max(0.15, 0.3 * sigma_expected_x);
    sigma_lb_y = max(0.15, 0.3 * sigma_expected_y);
    sigma_ub_x = min(box_width / 2, 3.0 * sigma_expected_x);
    sigma_ub_y = min(box_width / 2, 3.0 * sigma_expected_y);
    
    % QC 30 thresholds: reject if sigma is within 5% of the solver bounds
    % (means the solver was "fighting" the constraint)
    sigma_qc_lb_x = sigma_lb_x * 1.05;
    sigma_qc_lb_y = sigma_lb_y * 1.05;
    sigma_qc_ub_x = sigma_ub_x * 0.95;
    sigma_qc_ub_y = sigma_ub_y * 0.95;
    
    % --- 4. Pre-calculate Grid ---
    [grid_x_mat, grid_y_mat] = meshgrid(-radius:radius, -radius:radius);
    x_data = [grid_x_mat(:), grid_y_mat(:)]; % Constant for all fits
    
    % --- 5. Define Gaussian Model (Optimized Function Handle) ---
    % p(1)=Amp, p(2)=x0, p(3)=y0, p(4)=sx, p(5)=sy, p(6)=Bg
    gauss_func = @(p, xy) p(1) * exp( -((xy(:,1)-p(2)).^2/(2*p(4).^2) + (xy(:,2)-p(3)).^2/(2*p(5).^2)) ) + p(6);

    % Solver Options
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxIterations', 120, ...
        'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6, ...
        'UseParallel', false);      % Parallel is handled by our parfor loop
        
    % --- 6. Main Parallel Loop ---
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
        
        % Normalize ROI (shift so minimum is >= 0)
        roi_min = min(roi(:));
        if roi_min < 0, roi = roi - roi_min; end
        roi_range = max(roi(:)) - min(roi(:));
        
        % QC 10: Invalid/Flat ROI
        if roi_range < 1e-10 || any(isnan(roi(:)))
            exit_codes(i) = 10; continue;
        end
        
        z_data = roi(:);
        
        % --- Data-Driven Initial Sigma Estimate ---
        % Compute the second moment of the intensity distribution.
        % This gives a better starting point than a fixed value, while
        % the FWHM-derived bounds keep it physically plausible.
        weights = z_data - min(z_data);
        weights = weights / (sum(weights) + eps);
        sigma_est_x = sqrt(sum(weights .* x_data(:,1).^2));
        sigma_est_y = sqrt(sum(weights .* x_data(:,2).^2));
        
        % Clamp the data-driven estimate to the FWHM-derived range
        sigma_est_x = max(sigma_lb_x, min(sigma_est_x, sigma_ub_x));
        sigma_est_y = max(sigma_lb_y, min(sigma_est_y, sigma_ub_y));
        
        % Initial Parameters: p = [Amp, x0, y0, sx, sy, Bg]
        p0 = [roi_range, 0, 0, sigma_est_x, sigma_est_y, roi_min];
        
        % Solver Bounds (FWHM-derived)
        lb = [0,          -max_shift, -max_shift, sigma_lb_x, sigma_lb_y, 0];
        ub = [roi_range*3, max_shift,  max_shift, sigma_ub_x, sigma_ub_y, max(roi(:))];
        
        try
            [p_fit, ~, residual, exitflag] = lsqcurvefit(gauss_func, p0, x_data, z_data, lb, ub, opts);
            
            % QC 20: Solver Convergence Failure
            if exitflag <= 0
                exit_codes(i) = 20; continue;
            end
            
            % QC 30: Bad Sigma - fitted sigma is pinned against bounds
            if p_fit(4) <= sigma_qc_lb_x || p_fit(5) <= sigma_qc_lb_y || ...
               p_fit(4) >= sigma_qc_ub_x || p_fit(5) >= sigma_qc_ub_y
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
    
    % --- 7. Report Generation ---
    valid_mask = exit_codes == 0;
    
    localizations = table(out_frame(valid_mask), out_x(valid_mask), out_y(valid_mask), ...
                          out_int(valid_mask), out_conf(valid_mask), ...
                          'VariableNames', {'Frame', 'X', 'Y', 'Intensity', 'Confidence'});

    % Print the detailed report
    print_fast_qc_report(exit_codes, num_candidates, indent_prefix, ...
        fwhm, sigma_expected_x, sigma_expected_y, radius, ...
        [sigma_lb_x, sigma_lb_y], [sigma_ub_x, sigma_ub_y]);
end

function print_fast_qc_report(codes, total, prefix, fwhm, sig_ex, sig_ey, radius, sig_lb, sig_ub)
    kept = sum(codes == 0);
    fprintf('\n%s=== Fast Gaussian QC Report ===\n', prefix);
    fprintf('%sConfiguration:\n', prefix);
    fprintf('%s  FWHM [X, Z]:        [%.1f, %.1f] px\n', prefix, fwhm(1), fwhm(2));
    fprintf('%s  Expected sigma:      [%.2f, %.2f] px\n', prefix, sig_ex, sig_ey);
    fprintf('%s  Sigma bounds X:      [%.2f, %.2f]\n', prefix, sig_lb(1), sig_ub(1));
    fprintf('%s  Sigma bounds Z:      [%.2f, %.2f]\n', prefix, sig_lb(2), sig_ub(2));
    fprintf('%s  ROI radius:          %d px (%dx%d)\n', prefix, radius, 2*radius+1, 2*radius+1);
    fprintf('%s-------------------------------\n', prefix);
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