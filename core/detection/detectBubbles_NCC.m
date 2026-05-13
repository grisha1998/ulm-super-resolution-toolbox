function candidateBubbles = detectBubbles_NCC(filteredData, locParams, roiMask)
% detectBubbles_NCC  Microbubble candidate detection using Normalised
%                    Cross-Correlation (NCC) between each ultrasound frame
%                    and a PSF template, as described in:
%
%   Corazza et al., "Microbubble Identification Based on Decision Theory for
%   Ultrasound Localization Microscopy," IEEE Open J. UFFC, vol. 3, 2023.
%   DOI: 10.1109/OJUFFC.2023.3274512
%   (NCC method originally from Song et al. 2018 / Bourquin et al. 2022)
%
% SYNTAX:
%   candidateBubbles = detectBubbles_NCC(filteredData, locParams)
%   candidateBubbles = detectBubbles_NCC(filteredData, locParams, roiMask)
%
% INPUTS:
%   filteredData : [H x W x T double]  SVD-filtered (clutter-suppressed)
%                  ultrasound image sequence.
%
%   locParams    : struct with the following required fields:
%       .fwhm                  [1x2 double]  PSF FWHM: [lateral_x, axial_z]
%                                            in pixels.  Defines boundary-
%                                            exclusion margin.
%       .MB_image              [Kz x Kx double]  PSF template (kernel K)
%                                            used for cross-correlation
%                                            (Eq. 2 in paper).  BOTH
%                                            dimensions MUST be odd so that
%                                            the template has a well-defined
%                                            centre pixel and normxcorr2
%                                            output can be cropped
%                                            symmetrically.
%                                            Template options per paper:
%                                              - Gaussian kernel
%                                              - Simulated PSF (e.g. SIMUS)
%                                              - Experimentally measured PSF
%       .crosscor_threshold    [double]      Minimum NCC value tau_NCC for a
%                                            pixel to be considered a
%                                            candidate (Eq. 3).  Typical
%                                            values: 0.5-0.7.  Higher values
%                                            = more restrictive (fewer FP
%                                            but more FN).
%       .max_bubbles_per_frame [integer]     Maximum number of MB candidates
%                                            retained per frame.
%
%   locParams optional fields:
%       .ncc_peak_contrast_h   [double]      H-maxima contrast threshold for
%                                            resolving closely spaced bubbles.
%                                            Controls the minimum NCC dip
%                                            required between two peaks for
%                                            them to be reported as separate
%                                            detections.
%                                            - h = 0 (default): falls back to
%                                              imregionalmax (strict local max,
%                                              original behavior).  Two adjacent
%                                              pixels can never both be detected.
%                                            - h > 0: uses imextendedmax, which
%                                              reports any peak that rises at
%                                              least h above its surrounding
%                                              valley.  This can resolve two
%                                              bubbles whose NCC peaks are only
%                                              1-2 pixels apart, as long as the
%                                              correlation dips by >= h between
%                                              them.
%                                            Typical range: 0.02–0.10.
%                                            Lower h = more detections (higher
%                                            sensitivity, risk of splitting
%                                            noise into false positives).
%                                            Higher h = fewer detections (higher
%                                            specificity, safer but may merge
%                                            close pairs).
%
%   roiMask      : [H x W logical]  OPTIONAL.  Binary spatial mask that
%                  confines detection to a user-defined region of interest.
%                  Pixels outside the mask are zeroed in the NCC map before
%                  peak detection, so they can never produce candidates.
%                  Pass [] or omit to process the entire image.
%
%                  --- Is ROI masking compatible with the NCC algorithm? ---
%                  Yes.  normxcorr2 computes a global correlation map; the
%                  mask is applied *after* the full cross-correlation is
%                  computed, restricting which peaks are considered.  This
%                  is legitimate because the NCC coefficient at pixel (i,j)
%                  depends only on the local image patch centred at (i,j)
%                  (Eq. 2); masking the final map does not alter individual
%                  coefficient values, it only discards positions outside
%                  the ROI.  A typical workflow to build the mask:
%                    avgImg  = mean(filteredData, 3);
%                    roiMask = avgImg > prctile(avgImg(:), 20);
%                  or interactively with roipoly().
%
% OUTPUTS:
%   candidateBubbles : table with columns:
%       Frame     - frame index (1-based)
%       X         - lateral  pixel coordinate (column index)
%       Y         - axial    pixel coordinate (row    index)
%       Intensity - ORIGINAL (pre-threshold) signal amplitude at that
%                   pixel/frame, for sub-pixel localisation (Algorithm 1,
%                   line 18).
%
% ALGORITHM OVERVIEW (paper Section II.B & Eq. 2-3):
%   1.  For each frame, compute the 2-D NCC between the image and the PSF
%       template (Eq. 2).
%   2.  Apply the correlation threshold tau_NCC (Eq. 3) to produce the NCC
%       detection map CNCC(t).
%   3.  Apply boundary exclusion and the optional ROI mask.
%   4.  For each frame independently, find peaks and retain the top
%       max_bubbles_per_frame candidates.
%
% PEAK DETECTION MODES:
%   The peak-finding step (step 4) uses one of two strategies depending on
%   the ncc_peak_contrast_h parameter:
%
%   MODE A — Standard (h = 0, default):
%     Uses imregionalmax(cc_thresh), which requires a pixel to be STRICTLY
%     GREATER than all 8 of its neighbours to be considered a peak.
%     Consequence: two bubbles whose NCC peaks are on adjacent pixels
%     (distance <= 1 px) can never both be detected.  The dimmer one is
%     always suppressed.  This is safe for sparse bubble fields.
%
%   MODE B — H-maxima (h > 0):
%     Uses imextendedmax(cc_thresh, h), which finds peaks that rise at
%     least h units above their surrounding valley floor.
%     Example: if two bubbles produce an NCC profile like
%       [..., 0.60, 0.82, 0.75, 0.84, 0.55, ...]
%     imregionalmax sees only the 0.84 peak (the 0.82 pixel has a higher
%     neighbour at 0.84 to its right, so it fails the strict-maximum test).
%     imextendedmax with h = 0.05 detects BOTH 0.82 and 0.84, because the
%     valley between them (0.75) drops by 0.07 from 0.82 and 0.09 from
%     0.84, both exceeding h.
%     This mode is recommended for high-density bubble fields where PSFs
%     overlap and adjacent detections are physically meaningful.
%
% TEMPLATE SIZE NOTE:
%   The PSF kernel size SK presents a trade-off (paper Section IV.B):
%   too wide => risk of capturing signal from neighbouring MBs;
%   too narrow => insufficient PSF shape information.
%   The paper recommends SK = 5*lambda for in vivo studies.
%
% CHANGES FROM ORIGINAL VERSION:
%   - Added h-maxima peak detection mode (ncc_peak_contrast_h parameter)
%     using imextendedmax for resolving closely spaced bubbles in dense
%     fields.  Falls back to imregionalmax when h = 0 (full backward
%     compatibility).
%   - Added roiMask (3rd argument) support.
%   - Added strict validation of MB_image dimensions (must be odd in both
%     axes to guarantee symmetric normxcorr2 cropping).
%   - Renamed sx_template/sz_template to half_rows/half_cols matching the
%     physical axes (rows=axial=z, cols=lateral=x) so the semantics are
%     explicit and correct for non-square templates.
%   - Moved normxcorr2 loop into parfor for parallel frame processing.
%   - Fixed per-frame imregionalmax: now runs on each frame independently
%     (Algorithm 1, lines 11-12) rather than the reshaped 2-D trick that
%     introduced spurious cross-frame neighbours.
%   - Per-frame localmax + selection combined into a single parfor sweep.
%   - Boundary-exclusion indices made consistent with original formula.
%   - Output table column order: Frame, X, Y, Intensity.
%
% AUTHOR: Grigori Shapiro, updated from Corazza et al. (2023) reference implementation.
% DATE:   March 2026 (h-maxima extension: April 2026)

    % ------------------------------------------------------------------
    % 0. INPUT VALIDATION
    % ------------------------------------------------------------------
    narginchk(2, 3);

    validateattributes(filteredData, {'numeric'}, {'3d','finite'}, ...
        'detectBubbles_NCC', 'filteredData');

    template = locParams.MB_image;
    validateattributes(template, {'numeric'}, {'2d','nonempty'}, ...
        'detectBubbles_NCC', 'locParams.MB_image');

    filteredData = abs(filteredData);

    % Both template dimensions must be odd for symmetric normxcorr2 cropping.
    % normxcorr2 pads the output by (template_size - 1) on each side;
    % with an even template the centre is ambiguous and symmetric cropping
    % would be off by 0.5 pixels, introducing a systematic localisation bias.
    assert( mod(size(template, 1), 2) == 1 && ...
            mod(size(template, 2), 2) == 1, ...
        ['detectBubbles_NCC: locParams.MB_image has dimensions [%d x %d]. ' ...
         'Both dimensions must be odd for symmetric normxcorr2 cropping.  ' ...
         'Crop or pad your template to the nearest odd size.'], ...
        size(template, 1), size(template, 2));

    [height, width, numberOfFrames] = size(filteredData);

    % Parse and validate ROI mask
    if nargin < 3 || isempty(roiMask)
        roiMask = true(height, width);
    else
        roiMask = logical(roiMask);
        assert(isequal(size(roiMask), [height, width]), ...
            ['detectBubbles_NCC: roiMask size [%d x %d] does not match ' ...
             'filteredData spatial dimensions [%d x %d].'], ...
            size(roiMask,1), size(roiMask,2), height, width);
    end

    fwhmz = locParams.fwhm(2);   % axial   (row)    FWHM in pixels
    fwhmx = locParams.fwhm(1);   % lateral (column) FWHM in pixels
    tau   = locParams.crosscor_threshold;
    max_bpf = locParams.max_bubbles_per_frame;

    % ------------------------------------------------------------------
    % H-MAXIMA PARAMETER: controls whether to use imextendedmax (h > 0)
    % or the classic imregionalmax (h = 0) for peak detection.
    %
    % When h = 0, behaviour is identical to the original implementation.
    % When h > 0, two NCC peaks are reported as separate detections if
    % the valley between them drops by at least h.  This resolves closely
    % spaced bubbles that imregionalmax would merge into a single detection.
    % ------------------------------------------------------------------
    if isfield(locParams, 'ncc_peak_contrast_h') && ...
            ~isempty(locParams.ncc_peak_contrast_h)
        h_contrast = max(0, double(locParams.ncc_peak_contrast_h));
    else
        h_contrast = 0;   % default: classic imregionalmax (backward compat)
    end

    use_hmaxima = (h_contrast > 0);

    % Log which peak-detection mode is active
    if use_hmaxima
        fprintf('  NCC detector: h-maxima mode (h = %.4f) — resolves closely spaced peaks\n', h_contrast);
    else
        fprintf('  NCC detector: strict regional-max mode (imregionalmax)\n');
    end

    if isfield(locParams, 'detection_threshold') && ~isempty(locParams.detection_threshold)
        int_thresh = max(0, double(locParams.detection_threshold));
    else
        int_thresh = 0;
    end
    use_int_thresh = (int_thresh > 0);

    % Template half-extents used to crop normxcorr2 output back to [H x W].
    % normxcorr2(T, I) returns a matrix of size
    %   [size(I,1)+size(T,1)-1 x size(I,2)+size(T,2)-1]
    % Removing half_rows from top & bottom (and half_cols left & right)
    % realigns the correlation peak with the original pixel grid.
    half_rows = (size(template, 1) - 1) / 2;   % axial   (z) half-size
    half_cols = (size(template, 2) - 1) / 2;   % lateral (x) half-size

    MatIn_origin = abs(filteredData);

    % ------------------------------------------------------------------
    % 1. BOUNDARY EXCLUSION + ROI MASK (precomputed once)
    % ------------------------------------------------------------------
    % Boundary-exclusion margins ensure sub-pixel localisation kernels
    % (centred on detected maxima) never extend beyond the image border.
    z_start = 2 + round(fwhmz / 2);
    z_end   = height - 1 - round(fwhmz / 2);
    x_start = 2 + round(fwhmx / 2);
    x_end   = width  - 1 - round(fwhmx / 2);

    z_start = max(z_start, 1);  z_end = min(z_end, height);
    x_start = max(x_start, 1);  x_end = min(x_end, width);

    % Precompute a combined spatial validity mask:
    %   boundary_mask  = 1 inside the exclusion margin, 0 on edges
    %   combined_mask  = boundary_mask AND roiMask
    % Applying it once per frame inside parfor avoids repeated index
    % arithmetic and is safe for independent-frame parallelism.
    combined_mask = false(height, width);
    combined_mask(z_start:z_end, x_start:x_end) = true;
    combined_mask = combined_mask & roiMask;

    % ------------------------------------------------------------------
    % 2. PER-FRAME NCC + PEAK DETECTION + CANDIDATE SELECTION (parfor)
    % ------------------------------------------------------------------
    % We merge all steps into a single parfor loop over frames.
    % This is more efficient than separate loops because:
    %   (a) normxcorr2 is the dominant cost and parallelises trivially;
    %   (b) peak-finding is then done on the already-computed map
    %       without requiring a second pass over the data.
    %
    % parfor sliced variables:
    %   MatIn_origin(:,:,t)  — sliced along dimension 3 (valid parfor pattern)
    %   frame_results{t}     — classified as output cell array (valid)
    %
    % The template, tau, max_bpf, combined_mask, half_rows, half_cols,
    % h_contrast, and use_hmaxima are broadcast variables (read-only
    % scalars/arrays, copied to each worker).

    frame_results = cell(numberOfFrames, 1);

    parfor t = 1:numberOfFrames

        orig_frame = MatIn_origin(:, :, t);

        % --- NCC (Eq. 2) ---
        % normxcorr2 normalises for local mean and std, making the result
        % robust to intensity variations across the image (Section II.B).
        cc = normxcorr2(template, orig_frame);   % [(H+2*half_rows) x (W+2*half_cols)]

        % Crop padding added by normxcorr2 to restore original [H x W] size.
        % Rows: remove half_rows from top and bottom.
        % Cols: remove half_cols from left and right.
        cc_cropped = cc(half_rows + 1 : end - half_rows, ...
                        half_cols + 1 : end - half_cols);  % [H x W]

        % --- Threshold (Eq. 3): retain only pixels with NCC >= tau_NCC ---
        cc_thresh = cc_cropped .* (cc_cropped >= tau);

        % --- Apply combined spatial mask (boundary exclusion + ROI) ---
        cc_thresh = cc_thresh .* combined_mask;

        % ==============================================================
        % PEAK DETECTION — two modes depending on h_contrast
        % ==============================================================
        if use_hmaxima
            % ---- MODE B: H-MAXIMA (imextendedmax) ----
            %
            % imextendedmax(I, h) finds all regional maxima of the H-maxima
            % transform of I.  A peak is retained if it rises at least h
            % units above the highest contour line that does NOT encircle
            % any other peak.  Equivalently, two peaks are merged into one
            % only if the valley between them is shallower than h.
            %
            % This allows detection of two bubbles whose NCC peaks are on
            % adjacent or near-adjacent pixels, as long as a valley of
            % depth >= h exists between them.  For NCC maps (values in
            % [0, 1]), h in the range 0.02–0.10 is typically effective.
            %
            % Example with h = 0.05:
            %   NCC profile: [..., 0.60, 0.82, 0.75, 0.84, 0.55, ...]
            %   imregionalmax:  detects only 0.84 (0.82 < 0.84 neighbor)
            %   imextendedmax:  detects BOTH 0.82 and 0.84
            %     because valley 0.75 is 0.07 below 0.82 (> h)
            %     and                   0.09 below 0.84 (> h)
            %
            % Note: imextendedmax requires non-negative input.  cc_thresh
            % is already >= 0 due to the tau threshold above.
            local_max_mask = imextendedmax(cc_thresh, h_contrast);
        else
            % ---- MODE A: STRICT REGIONAL MAXIMA (imregionalmax) ----
            %
            % Classic behaviour: a pixel is a peak only if it is strictly
            % greater than ALL 8 of its neighbours.  Two adjacent pixels
            % can never both be peaks.  This is safe for sparse fields
            % but suppresses closely spaced bubbles.
            local_max_mask = imregionalmax(cc_thresh);
        end

        candidate_idx  = find(local_max_mask);

        if isempty(candidate_idx)
            frame_results{t} = zeros(0, 4, 'double');
            continue;
        end
        
        % Enforce normalized intensity threshold to reject weak candidates
        if use_int_thresh
            frame_max = max(orig_frame(:));
            if frame_max > 0
                frame_norm = orig_frame / frame_max;
            else
                frame_norm = orig_frame;
            end
            
            valid_int_mask = frame_norm(candidate_idx) >= int_thresh;
            candidate_idx = candidate_idx(valid_int_mask);
            
            if isempty(candidate_idx)
                frame_results{t} = zeros(0, 4, 'double');
                continue;
            end
        end
        
        % NCC values at candidate positions (used for ranking)
        ncc_vals = cc_thresh(candidate_idx);

        % Retain top max_bpf candidates by NCC value
        if numel(ncc_vals) > max_bpf
            [~, rank_idx] = sort(ncc_vals, 'descend');
            candidate_idx = candidate_idx(rank_idx(1:max_bpf));
        end

        % Retrieve ORIGINAL (pre-threshold) amplitudes for sub-pixel
        % localisation (Algorithm 1, line 18: localisation uses X(t),
        % the original filtered image, not the detection map).
        orig_vals = orig_frame(candidate_idx);

        % Convert linear indices to (row=z, col=x) subscripts
        [idx_z, idx_x] = ind2sub([height, width], candidate_idx);

        % Accumulate: [Frame, X_lateral, Y_axial, Intensity]
        n = numel(idx_z);
        frame_results{t} = [repmat(t, n, 1), idx_x, idx_z, orig_vals];

    end  % parfor

    % ------------------------------------------------------------------
    % 3. ASSEMBLE OUTPUT TABLE
    % ------------------------------------------------------------------
    all_candidates = cat(1, frame_results{:});

    if isempty(all_candidates)
        candidateBubbles = table( ...
            'Size', [0, 4], ...
            'VariableTypes', {'double','double','double','double'}, ...
            'VariableNames', {'Frame','X','Y','Intensity'});
    else
        candidateBubbles = array2table(all_candidates, ...
            'VariableNames', {'Frame','X','Y','Intensity'});
    end

end