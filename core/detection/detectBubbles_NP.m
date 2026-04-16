function candidateBubbles = detectBubbles_NP(filteredData, locParams, roiMask)
% detectBubbles_NP  Microbubble candidate detection using the Neyman-Pearson
%                   (NP) decision criterion, as described in:
%
%   Corazza et al., "Microbubble Identification Based on Decision Theory for
%   Ultrasound Localization Microscopy," IEEE Open J. UFFC, vol. 3, 2023.
%   DOI: 10.1109/OJUFFC.2023.3274512
%
% SYNTAX:
%   candidateBubbles = detectBubbles_NP(filteredData, locParams)
%   candidateBubbles = detectBubbles_NP(filteredData, locParams, roiMask)
%
% INPUTS:
%   filteredData : [H x W x T double]  SVD-filtered (clutter-suppressed)
%                  ultrasound image sequence.
%
%   locParams    : struct with the following required fields:
%       .fwhm                  [1x2 double]  Full-width at half-maximum of
%                                            the PSF: [lateral_x, axial_z]
%                                            in pixels.  Used to define the
%                                            boundary-exclusion margin so
%                                            that the sub-pixel localisation
%                                            kernel never extends outside the
%                                            image (Algorithm 1, line 13).
%       .NP_alpha0             [double]      False-alarm rate, i.e. the
%                                            probability of declaring H1
%                                            (MB present) when H0 (noise
%                                            only) is true.  Must satisfy
%                                            0 < NP_alpha0 < 0.5  (paper
%                                            Section III.C).  Typical range:
%                                            1e-7 to 0.01.
%       .max_bubbles_per_frame [integer]     Maximum number of MB candidates
%                                            retained per frame after local-
%                                            maxima selection.
%
%   roiMask      : [H x W logical]  OPTIONAL.  Binary spatial mask that
%                  confines detection to a user-defined region of interest
%                  (e.g. a vascular structure segmented on the mean-intensity
%                  image).  Pixels outside the mask are forced to zero in
%                  the detection map, so they can never produce candidates.
%                  Statistical moments (median, MAD) are still estimated
%                  from the full temporal signal of each pixel so that the
%                  NP threshold is unbiased; only the *output* is masked.
%                  Pass [] or omit to process the entire image.
%
%                  --- Is ROI masking compatible with the NP algorithm? ---
%                  Yes.  The NP criterion is applied independently at each
%                  spatio-temporal pixel (xi,j(t)); masking simply discards
%                  the result at unwanted locations without touching the
%                  per-pixel threshold computation.  It is analogous to
%                  restricting the hypothesis test to a spatial subset of
%                  interest.  A typical workflow to generate the mask:
%                    avgImg   = mean(filteredData, 3);
%                    roiMask  = avgImg > prctile(avgImg(:), 20);
%                  or draw it interactively with roipoly().
%
% OUTPUTS:
%   candidateBubbles : table  with columns:
%       Frame     - frame index (1-based)
%       X         - lateral  pixel coordinate (column index)
%       Y         - axial    pixel coordinate (row    index)
%       Intensity - ORIGINAL (pre-threshold) signal amplitude at that
%                   pixel/frame, forwarded to the sub-pixel localisation
%                   step (Algorithm 1, line 18).
%
% ALGORITHM OVERVIEW (Algorithm 1 in the paper):
%   1.  For each spatial pixel (i,j), estimate the noise mean m_ij via the
%       temporal median and the noise std sigma_ij via the scaled MAD
%       (Eq. 16-17).  Both are robust to outliers caused by MB passages.
%   2.  Compute the pixel-wise NP threshold eta_ij (Eq. 15):
%           eta_ij = sqrt(2)*sigma_ij*erfinv(1 - 2*alpha0) + m_ij
%   3.  Build the detection map C (Eq. 18):
%           C_ij(t) = (x_ij(t) - eta_ij) * [x_ij(t) >= eta_ij]
%   4.  Apply boundary exclusion and the optional ROI mask.
%   5.  For each frame independently, find regional maxima on C(t) and
%       retain the top max_bubbles_per_frame candidates.
%
% CHANGES FROM ORIGINAL VERSION:
%   - Added roiMask (3rd argument) support.
%   - Added strict validation of NP_alpha0 range (must be in (0, 0.5)).
%   - Guard against sigma = 0 (dead/constant pixels) by clamping to eps.
%   - Fixed per-frame imregionalmax: now runs on each frame independently
%     (Algorithm 1, lines 11-12) instead of the reshaped 2D trick that
%     caused spurious cross-frame neighbour comparisons.
%   - Replaced manual sort/threshold selection loop with vectorised parfor
%     over frames (steps 3-4 combined per frame).
%   - Boundary-exclusion indices made consistent with original formula.
%   - Output table column order: Frame, X, Y, Intensity (matches
%     detectBubbles.m and detectBubbles_NCC.m).
%
% AUTHOR: Grigori Shapiro, updated from Corazza et al. (2023) reference implementation.
% DATE: March 2026

    % ------------------------------------------------------------------
    % 0. INPUT VALIDATION
    % ------------------------------------------------------------------
    narginchk(2, 3);

    validateattributes(filteredData, {'numeric'}, {'3d','finite'}, ...
        'detectBubbles_NP', 'filteredData');
    validateattributes(locParams.NP_alpha0, {'numeric'}, ...
        {'scalar','real'}, 'detectBubbles_NP', 'locParams.NP_alpha0');

    % Paper Section III.C: alpha0 must lie in the open interval (0, 0.5).
    % Outside this range the detector is either trivially off (>=0.5) or
    % produces near-zero detections that are impractical for ULM (<=0).
    assert(locParams.NP_alpha0 > 0 && locParams.NP_alpha0 < 0.5, ...
        ['detectBubbles_NP: NP_alpha0 = %.2e is outside the valid ' ...
         'interval (0, 0.5) defined in the paper (Section III.C). ' ...
         'Typical values are 1e-7 to 0.01.'], locParams.NP_alpha0);

    [height, width, numberOfFrames] = size(filteredData);

    % Parse and validate ROI mask
    if nargin < 3 || isempty(roiMask)
        roiMask = true(height, width);   % process entire image
    else
        roiMask = logical(roiMask);
        assert(isequal(size(roiMask), [height, width]), ...
            ['detectBubbles_NP: roiMask size [%d x %d] does not match ' ...
             'filteredData spatial dimensions [%d x %d].'], ...
            size(roiMask,1), size(roiMask,2), height, width);
    end

    fwhmz = locParams.fwhm(2);   % axial   (row)    FWHM in pixels
    fwhmx = locParams.fwhm(1);   % lateral (column) FWHM in pixels
    max_bpf = locParams.max_bubbles_per_frame;

    % Work on |filteredData|  (consistent with original implementation)
    MatIn_origin = abs(filteredData);

    % ------------------------------------------------------------------
    % 1. CASORATI MATRIX & ROBUST NOISE STATISTICS (Eqs. 16-17)
    % ------------------------------------------------------------------
    % Reshape to [N_pixels x N_frames] (Casorati matrix).
    % The median and MAD are computed along the time dimension (dim 2),
    % giving one robust noise estimate per spatial pixel.
    % This is the fully vectorised equivalent of the per-pixel loop in
    % Algorithm 1, lines 6-7.
    Casorati = reshape(MatIn_origin, [height * width, numberOfFrames]);

    % Robust noise mean: temporal median (Algorithm 1, line 6)
    moy = median(Casorati, 2);   % [N_pixels x 1]

    % Robust noise std via scaled MAD (Algorithm 1, line 7 & Eq. 17)
    %   sigma = MAD / (sqrt(2) * erfinv(0.5))
    % erfinv(0.5) ≈ 0.4769 => scaling factor ≈ 1.4826
    scale = 1 / (erfinv(0.5) * sqrt(2));
    stdev = scale * mad(Casorati, 1, 2);   % [N_pixels x 1]

    % Guard: constant/dead pixels (sigma = 0) would cause division issues
    % downstream.  Clamping to eps ensures thresh >> 0 and those pixels
    % never produce detections.
    stdev(stdev == 0) = eps;

    % ------------------------------------------------------------------
    % 2. NP THRESHOLD AND DETECTION MAP (Eqs. 15 & 18)
    % ------------------------------------------------------------------
    % Per-pixel NP threshold (Eq. 15):
    %   eta_ij = sqrt(2) * sigma_ij * erfinv(1 - 2*alpha0) + m_ij
    thresh = sqrt(2) .* stdev .* erfinv(1 - 2 * locParams.NP_alpha0) + moy;
    % thresh is [N_pixels x 1]; broadcast across time via bsxfun-style subtraction

    % Detection map C (Eq. 18):  C = (x - eta) * [x >= eta]
    % Result is zero for noise-only samples and positive for putative MBs.
    diff_map = Casorati - thresh;          % [N_pixels x N_frames]
    diff_map(diff_map < 0) = 0;
    diff_map(isnan(diff_map))  = 0;

    % Reshape detection map back to [H x W x T]
    MatIn = reshape(diff_map, [height, width, numberOfFrames]);

    % ------------------------------------------------------------------
    % 3. APPLY ROI MASK TO DETECTION MAP
    % ------------------------------------------------------------------
    % Zero out detection values outside the user-defined region of
    % interest.  This is applied *after* threshold computation so that
    % the per-pixel noise statistics remain unbiased by the masking.
    if ~all(roiMask(:))
        % Expand mask to 3-D and zero out non-ROI detections
        MatIn = MatIn .* repmat(roiMask, [1, 1, numberOfFrames]);
    end

    % ------------------------------------------------------------------
    % 4. BOUNDARY EXCLUSION (matrix cropping)
    % ------------------------------------------------------------------
    % Exclude a margin of ceil(FWHM/2)+1 pixels on each side so that
    % sub-pixel localisation kernels (centred on detected maxima) never
    % extend beyond the image border.
    MatInReduced = zeros(height, width, numberOfFrames, 'like', MatIn_origin);
    z_start = 2 + round(fwhmz / 2);
    z_end   = height - 1 - round(fwhmz / 2);
    x_start = 2 + round(fwhmx / 2);
    x_end   = width  - 1 - round(fwhmx / 2);

    % Safety: ensure indices are valid (can be violated with very large FWHM)
    z_start = max(z_start, 1);  z_end = min(z_end, height);
    x_start = max(x_start, 1);  x_end = min(x_end, width);

    MatInReduced(z_start:z_end, x_start:x_end, :) = ...
        MatIn(z_start:z_end, x_start:x_end, :);

    % ------------------------------------------------------------------
    % 5. PER-FRAME LOCAL MAXIMA DETECTION & CANDIDATE SELECTION
    %    (Algorithm 1, lines 11-13) — parallelised with parfor
    % ------------------------------------------------------------------
    % KEY DESIGN DECISION:
    %   The original code collapsed all frames into a single 2-D matrix
    %   before calling imregionalmax.  This creates spurious spatial
    %   neighbours across frame boundaries (last row of frame t is
    %   adjacent to first row of frame t+1 in reshaped memory).  The
    %   paper's Algorithm 1 explicitly loops over frames (line 11:
    %   "for t = 1 to Nt"); we honour that here by using parfor so
    %   each frame is processed independently and in parallel.
    %
    %   MATLAB parfor requirement: MatIn_origin and MatInReduced are
    %   sliced along their 3rd dimension (t), which is a valid sliced-
    %   variable pattern in MATLAB's parfor.

    frame_results = cell(numberOfFrames, 1);

    parfor t = 1:numberOfFrames

        % Slice current frame from detection map and original data
        det_frame  = MatInReduced(:, :, t);   % detection-map frame
        orig_frame = MatIn_origin(:, :, t);   % original amplitude frame

        % Find regional maxima on the detection map for this frame only
        local_max_mask = imregionalmax(det_frame);

        % Linear indices of all local maxima in this frame
        candidate_idx = find(local_max_mask);

        if isempty(candidate_idx)
            frame_results{t} = zeros(0, 4, 'double');
            continue;
        end

        % Detection-map value at each candidate (used for ranking/limiting)
        det_vals = det_frame(candidate_idx);

        % Keep only the top max_bpf candidates by detection-map value
        % (equivalent to the N_bubbles limit in the intensity-based method)
        if numel(det_vals) > max_bpf
            [~, rank_idx] = sort(det_vals, 'descend');
            candidate_idx = candidate_idx(rank_idx(1:max_bpf));
        end

        % Retrieve ORIGINAL (pre-threshold) intensities at candidate
        % locations — these are passed to sub-pixel localisation which
        % operates on the unmodified image (Algorithm 1, line 18).
        orig_vals = orig_frame(candidate_idx);

        % Convert linear indices to (row=z, col=x) subscripts
        [idx_z, idx_x] = ind2sub([height, width], candidate_idx);

        % Accumulate: [Frame, X_lateral, Y_axial, Intensity]
        n = numel(idx_z);
        frame_results{t} = [repmat(t, n, 1), idx_x, idx_z, orig_vals];

    end  % parfor

    % ------------------------------------------------------------------
    % 6. ASSEMBLE OUTPUT TABLE
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
