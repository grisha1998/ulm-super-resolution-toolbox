function [filtered_Data, diagnostics] = SVD_blockwise(rawData, params_or_framerate, varargin)
% =========================================================================
% FUNCTION: SVD_blockwise
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% ACADEMIC REFERENCE:
%   Main Architecture:
%   Song et al., "Ultrasound Small Vessel Imaging With Block-Wise Adaptive
%   Local Clutter Filtering," IEEE Trans. Med. Imag., vol. 36, no. 1, 2017.
%   DOI: 10.1109/TMI.2016.2605819
%
%   SSM Thresholding Option:
%   J. Baranger, et al. "Adaptive Spatiotemporal SVD Clutter Filtering
%   for Ultrafast Doppler Imaging Using Similarity of Spatial Singular
%   Vectors." IEEE T-MI, 2018.
%
% PURPOSE & ADVANTAGES:
%   Implements a spatially adaptive, block-wise SVD clutter filter. 
%   - Advantages: Global SVD assumes clutter motion (tissue pulsation) is 
%     uniform across the entire field of view, which is false for in-vivo 
%     imaging (e.g., a pulsating artery near static skull tissue). This 
%     function solves this by dividing the image into local overlapping blocks, 
%     finding the optimal SVD threshold for the specific local dynamics, and 
%     seamlessly blending them together. This preserves small vessels near 
%     highly pulsatile regions that a global filter would destroy.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Grid Generation: Calculates block dimensions (e.g., 4x4 mm) and 
%      step sizes based on the requested overlap percentage (typically 75%+).
%   2. Parallel Local Processing: The image is segmented into hundreds of 
%      sub-blocks. Utilizing `parfor`, an independent Casorati matrix and 
%      SVD decomposition is computed for each block simultaneously.
%   3. Local Thresholding: For each block, low/high cutoffs are estimated 
%      using the chosen mathematical method ('DopplerGradient', 'SSM', etc.).
%   4. Local Reconstruction: The block is reconstructed using only its specific 
%      blood components.
%   5. Global Accumulation (Blending): The local blocks are added back into a 
%      global 3D accumulator matrix. A 2D "hit-count" map tracks how many 
%      blocks overlapped at each pixel.
%   6. Normalization: The accumulator is divided by the hit-count map to 
%      produce a smooth, artifact-free final image.
%
% SYNTAX OPTIONS:
%   filt_data = SVD_blockwise(rawData, params)
%   filt_data = SVD_blockwise(rawData, framerate_Hz, 'ThresholdMethod', 'SSM')
%   [filt_data, diag] = SVD_blockwise(rawData, params, 'BlockSizeMm', [4, 4], 'OverlapPct', 80)
%
% EXAMPLES:
%   % Example 1: Standard usage with ULM params struct
%   filtered_Data = SVD_blockwise(IQ_data, params);
%
%   % Example 2: High overlap, specific block size in mm, using Hybrid method
%   [filtered_Data, diagnostics] = SVD_blockwise(IQ_data, params, ...
%       'ThresholdMethod', 'Hybrid', ...
%       'BlockSizeMm', [4.5, 4.5], ...
%       'OverlapPct', 93.75, ...
%       'PlotThresholdMaps', true);
%
% INPUTS:
%   rawData             - (Type: 3D Numeric Array [H x W x T]) Input sequence.
%   params_or_framerate - (Type: Struct OR Double) Either the main ULM parameters 
%                         struct (containing .acq.framerate and pixel sizes) OR 
%                         just a scalar framerate (if pixel sizes are unknown).
%   varargin            - (Type: Name-Value Pairs)
%       * 'ThresholdMethod'       : (String) 'DopplerGradient' (default), 'SSM', 
%                                   'Hybrid', or 'Manual'.
%       * 'BlockSizeMm'           : (1x2 Double) Physical size of blocks [Z, X]. 
%                                   Example: [4.0, 4.0]. Overrides px size.
%       * 'BlockSizePx'           : (1x2 Int) Pixel size of blocks. Default: Auto.
%       * 'OverlapPct'            : (Double, 0-100) Percentage of block overlap. 
%                                   Example: 75.0 (default), 93.75 (paper optimal).
%       * 'ManualCutoff'          : (1x2 Int) [Low, High] for 'Manual' method only.
%       * 'TissueFreqThreshHz'    : (Double) Cutoff for Doppler freq. Default: auto.
%       * 'MPDeviationSigma'      : (Double) Marchenko-Pastur strictness. Default: 2.0.
%       * 'PlotThresholdMaps'     : (Logical) Shows heatmaps of cutoffs. Default: false.
%
% OUTPUTS:
%   filtered_Data  [H x W x T]  Clutter-filtered blood signal.
%
%   diagnostics    struct  with fields:
%     .low_cutoff_map    [H x W]  Interpolated map of low-order cutoffs
%     .high_cutoff_map   [H x W]  Interpolated map of high-order cutoffs
%     .n_blocks          int      Total blocks processed
%     .block_size_px     [1x2]    Block size used [nz, nx]
%     .step_px           [1x2]    Step size used [nz, nx]
%     .block_size_mm     [1x2]    Block size in mm (NaN if pixel size unknown)
%     .overlap_pct       [1x2]    Achieved overlap [z%, x%]
%     .elapsed_sec       double   Total wall-clock time
%     .threshold_method  string   Method used
%
% PHYSICS CONSTRAINT (SVD reliability):
%   nz * nx > T must hold for each block.
%   Rationale: the Casorati matrix must have more rows (spatial pixels) than
%   columns (frames) so that SVD gives K=T full-rank components.
%   If violated, singular values are truncated and thresholding is unreliable.
%   Auto-sizing and warnings enforce this constraint.
%
% PARALLELISM:
%   parfor is used for per-block SVD computations. Results are stored in a
%   cell array and accumulated sequentially — this is the correct pattern
%   for overlapping-write parallelism in MATLAB (avoids race conditions).
%
% =========================================================================

    % =========================================================================
    % 1. INPUT PARSING
    % =========================================================================
    narginchk(2, Inf);

    % Parse options as Name-Value pairs
    p = inputParser();
    addParameter(p, 'ThresholdMethod',        'DopplerGradient');
    addParameter(p, 'BlockSizeMm',            []);
    addParameter(p, 'BlockSizePx',            []);
    addParameter(p, 'OverlapPct',              75.0);
    addParameter(p, 'ManualCutoff',           [10, 200]);
    addParameter(p, 'TissueFreqThreshHz',      -1);     % -1 = auto
    addParameter(p, 'MPDeviationSigma',         2.0);
    addParameter(p, 'GradientInflectionPct',    0.10);
    addParameter(p, 'MinBloodComponents',       3);
    addParameter(p, 'MaxTissueFraction',        0.60);
    addParameter(p, 'PlotThresholdMaps',        false);
    addParameter(p, 'IndentPrefix',            '');
    parse(p, varargin{:});
    opt = p.Results;

    % =========================================================================
    % 2. EXTRACT PHYSICAL PARAMETERS
    % =========================================================================
    if isstruct(params_or_framerate)
        ps          = params_or_framerate;
        framerate   = get_nested(ps, 'acq.framerate',             200);
        pixel_z_mm  = get_nested(ps, 'expParams.pixel_Z_size',    NaN);
        pixel_x_mm  = get_nested(ps, 'expParams.pixel_X_size',    NaN);
    elseif isnumeric(params_or_framerate) && isscalar(params_or_framerate)
        framerate   = double(params_or_framerate);
        pixel_z_mm  = NaN;
        pixel_x_mm  = NaN;
    else
        error('SVD_blockwise:BadInput', ...
            'Second argument must be a params struct or scalar frame rate [Hz].');
    end

    has_pixel_size = ~isnan(pixel_z_mm) && ~isnan(pixel_x_mm) && ...
                      pixel_z_mm > 0   && pixel_x_mm > 0;

    [H, W, T]  = size(rawData);
    is_complex = ~isreal(rawData);
    pfx        = opt.IndentPrefix;
    method     = opt.ThresholdMethod;

    % Auto tissue frequency threshold
    if opt.TissueFreqThreshHz < 0
        freq_thr = max(5.0, min(20.0, framerate / 50));
    else
        freq_thr = opt.TissueFreqThreshHz;
    end

    % =========================================================================
    % 3. VALIDATE THRESHOLD METHOD
    % =========================================================================
    valid_methods = {'DopplerGradient','SSM','Hybrid','Manual'};
    assert(any(strcmpi(method, valid_methods)), ...
        'SVD_blockwise:BadMethod', ...
        'ThresholdMethod must be one of: %s', strjoin(valid_methods, ', '));

    % SSM requires helper function
    use_ssm = any(strcmpi(method, {'SSM','Hybrid'}));
    if use_ssm && ~exist('SVD_SSM_AutoTh', 'file')
        warning('SVD_blockwise:NoSSM', ...
            '%sSVD_SSM_AutoTh.m not found. Falling back to DopplerGradient.', pfx);
        if strcmpi(method,'SSM'),    method = 'DopplerGradient'; end
        if strcmpi(method,'Hybrid'), method = 'DopplerGradient'; end
    end

    % Doppler-based cutoffs require complex IQ data
    uses_doppler = any(strcmpi(method, {'DopplerGradient','Hybrid'}));
    if uses_doppler && ~is_complex
        warning('SVD_blockwise:RealData', ...
            ['%sData is real (not IQ complex). Doppler frequency estimates ' ...
             'will be unreliable. Cutoff 1B disabled; using Cutoff 1A + MP only.'], pfx);
        uses_doppler = false;
    end

    % =========================================================================
    % 4. BLOCK SIZE DETERMINATION (physics-aware)
    % =========================================================================
    if ~isempty(opt.BlockSizeMm) && has_pixel_size
        % Physical mm specification — convert to pixels
        sz_mm = opt.BlockSizeMm;
        if isscalar(sz_mm), sz_mm = [sz_mm, sz_mm]; end
        blk_z = max(1, round(sz_mm(1) / pixel_z_mm));
        blk_x = max(1, round(sz_mm(2) / pixel_x_mm));
        fprintf('%sBlock size from mm spec: [%.1f x %.1f mm] -> [%d x %d px]\n', ...
            pfx, sz_mm(1), sz_mm(2), blk_z, blk_x);

    elseif ~isempty(opt.BlockSizePx)
        sz = opt.BlockSizePx;
        if isscalar(sz), sz = [sz, sz]; end
        blk_z = sz(1);  blk_x = sz(2);

    else
        % Auto-size: smallest multiple of 10 satisfying nz*nx > T
        % The paper used blk=22 as the absolute minimum for T=450 (22^2=484>450)
        min_side = max(22, ceil(sqrt(T)) + 1);
        blk_auto = ceil(min_side / 10) * 10;
        blk_z = blk_auto;  blk_x = blk_auto;
        fprintf('%sBlock size auto-set: [%d x %d px] (min for T=%d: %d px)\n', ...
            pfx, blk_z, blk_x, T, min_side);
    end

    % Clamp to image dimensions
    blk_z = min(blk_z, H);
    blk_x = min(blk_x, W);

    % --- Enforce SVD reliability constraint: nz * nx > T ---
    if blk_z * blk_x <= T
        new_side = ceil(sqrt(T + 1));
        new_side = ceil(new_side / 5) * 5;      % round to nearest 5
        new_side = max(new_side, 22);
        warning('SVD_blockwise:BlockTooSmall', ...
            ['%sBlock [%dx%d] = %d spatial samples <= T=%d frames.\n' ...
             '%sSVD would be rank-deficient. Auto-increasing to [%dx%d].'], ...
            pfx, blk_z, blk_x, blk_z*blk_x, T, pfx, new_side, new_side);
        blk_z = min(new_side, H);
        blk_x = min(new_side, W);
    end

    % =========================================================================
    % 5. STEP SIZE FROM OVERLAP PERCENTAGE
    % =========================================================================
    step_z = max(1, round(blk_z * (1 - opt.OverlapPct / 100)));
    step_x = max(1, round(blk_x * (1 - opt.OverlapPct / 100)));
    ach_ov_z = (1 - step_z / blk_z) * 100;
    ach_ov_x = (1 - step_x / blk_x) * 100;

    % =========================================================================
    % 6. PHYSICS-AWARE CONSOLE REPORT
    % =========================================================================
    bar = repmat('-', 1, 60);
    fprintf('%s%s\n', pfx, bar);
    fprintf('%s SVD_blockwise — configuration summary\n', pfx);
    fprintf('%s%s\n', pfx, bar);
    fprintf('%s  Data:           [%d x %d] px,  T=%d frames,  FR=%.0f Hz\n', ...
        pfx, H, W, T, framerate);
    if has_pixel_size
        fprintf('%s  Pixel size:     z=%.3f mm,  x=%.3f mm\n', ...
            pfx, pixel_z_mm, pixel_x_mm);
        fprintf('%s  Block size:     [%d x %d px]  =  [%.2f x %.2f mm]\n', ...
            pfx, blk_z, blk_x, blk_z*pixel_z_mm, blk_x*pixel_x_mm);
        fprintf('%s  Step size:      [%d x %d px]  =  [%.2f x %.2f mm]\n', ...
            pfx, step_z, step_x, step_z*pixel_z_mm, step_x*pixel_x_mm);
    else
        fprintf('%s  Block size:     [%d x %d px]  (pixel size not available)\n', ...
            pfx, blk_z, blk_x);
        fprintf('%s  Step size:      [%d x %d px]\n', pfx, step_z, step_x);
    end
    fprintf('%s  Overlap:        z=%.1f%%  x=%.1f%%\n', pfx, ach_ov_z, ach_ov_x);
    fprintf('%s  SVD constraint: blk_z*blk_x=%d > T=%d  [%s]\n', ...
        pfx, blk_z*blk_x, T, yesno(blk_z*blk_x > T));
    fprintf('%s  Method:         %s\n', pfx, method);
    if strcmpi(method,'Manual')
        fprintf('%s  Manual cutoff:  [%d, %d]\n', pfx, ...
            opt.ManualCutoff(1), opt.ManualCutoff(2));
    else
        fprintf('%s  Tissue f_thr:   %.1f Hz\n', pfx, freq_thr);
        fprintf('%s  MP sigma:       %.1f\n', pfx, opt.MPDeviationSigma);
    end
    mem_mb = blk_z * blk_x * T * 8 * 2 / 1e6;   % 2x for complex
    fprintf('%s  Est. mem/block: ~%.1f MB\n', pfx, mem_mb);
    fprintf('%s%s\n', pfx, bar);

    % =========================================================================
    % 7. GENERATE BLOCK GRID
    % =========================================================================
    i_starts = unique([1 : step_z : (H - blk_z + 1), H - blk_z + 1]);
    j_starts = unique([1 : step_x : (W - blk_x + 1), W - blk_x + 1]);
    i_starts = i_starts(i_starts >= 1 & i_starts + blk_z - 1 <= H);
    j_starts = j_starts(j_starts >= 1 & j_starts + blk_x - 1 <= W);

    n_row  = numel(i_starts);
    n_col  = numel(j_starts);
    n_blks = n_row * n_col;
    fprintf('%s  Block grid:     %d rows x %d cols = %d blocks\n', ...
        pfx, n_row, n_col, n_blks);

    % Flatten block list for parfor indexing
    block_list = zeros(n_blks, 2, 'int32');
    b = 0;
    for ii = 1:n_row
        for jj = 1:n_col
            b = b + 1;
            block_list(b,:) = [i_starts(ii), j_starts(jj)];
        end
    end

    % =========================================================================
    % 8. PER-BLOCK SVD + THRESHOLD ESTIMATION  (parfor)
    %
    % PARALLELISM STRATEGY:
    %   - parfor processes each block independently (no shared state).
    %   - Results stored in a cell array (parfor-safe output variable).
    %   - Sequential accumulation loop follows — this is the only correct
    %     way to handle overlapping writes in MATLAB parfor.
    % =========================================================================
    block_results = cell(n_blks, 1);

    % Cache options as scalars for parfor broadcast (structs are broadcast vars)
    par_method      = method;
    par_freq_thr    = freq_thr;
    par_mp_sigma    = opt.MPDeviationSigma;
    par_grad_pct    = opt.GradientInflectionPct;
    par_min_blood   = opt.MinBloodComponents;
    par_max_tissue  = opt.MaxTissueFraction;
    par_manual      = opt.ManualCutoff;
    par_use_doppler = uses_doppler;
    par_fr          = framerate;

    t_wall = tic;
    fprintf('%s  Processing blocks...\n', pfx);

    parfor b_idx = 1:n_blks
        r1 = block_list(b_idx, 1);
        c1 = block_list(b_idx, 2);
        r2 = r1 + blk_z - 1;
        c2 = c1 + blk_x - 1;

        % Build Casorati matrix for this block
        casorati = reshape(rawData(r1:r2, c1:c2, :), blk_z * blk_x, T);

        % SVD — 'econ' returns min(spatial, T) components
        [U, S_mat, V] = svd(casorati, 'econ');
        S_vec = diag(S_mat);
        K     = numel(S_vec);

        % Determine cutoff indices
        switch upper(par_method)

            case 'DOPPLERGRADIENT'
                [cut_lo, cut_hi] = blk_DopplerGradient( ...
                    S_vec, V, par_fr, par_freq_thr, par_mp_sigma, ...
                    par_grad_pct, par_use_doppler, par_max_tissue, K);

            case 'SSM'
                [cut_lo, cut_hi] = blk_SSM( ...
                    U, S_vec, V, par_fr, par_mp_sigma, par_max_tissue, K);

            case 'HYBRID'
                [cut_lo, ~]      = blk_SSM( ...
                    U, S_vec, V, par_fr, par_mp_sigma, par_max_tissue, K);
                [~,      cut_hi] = blk_DopplerGradient( ...
                    S_vec, V, par_fr, par_freq_thr, par_mp_sigma, ...
                    par_grad_pct, par_use_doppler, par_max_tissue, K);

            case 'MANUAL'
                cut_lo = par_manual(1);
                cut_hi = par_manual(2);

            otherwise
                cut_lo = max(1, round(0.10 * K));
                cut_hi = round(0.80 * K);
        end

        % Clamp and enforce minimum blood components
        cut_lo = max(1,  min(cut_lo, K - par_min_blood));
        cut_hi = min(K,  max(cut_hi, cut_lo + par_min_blood - 1));

        % Reconstruct blood signal for this block
        idx          = cut_lo : cut_hi;
        blood_matrix = U(:, idx) * S_mat(idx, idx) * V(:, idx)';
        blood_block  = reshape(blood_matrix, blk_z, blk_x, T);

        % Pack result for sequential accumulation
        block_results{b_idx} = {blood_block, int32(cut_lo), int32(cut_hi), ...
                                 int32(r1), int32(r2), int32(c1), int32(c2)};
    end  % parfor

    fprintf('%s  parfor done in %.1f s. Accumulating...\n', pfx, toc(t_wall));

    % =========================================================================
    % 9. SEQUENTIAL ACCUMULATION  (Eq. 1 in the paper)
    % =========================================================================
    accumulator    = zeros(H, W, T, 'double');
    hit_count      = zeros(H, W,    'double');
    low_cut_sum    = zeros(H, W,    'double');   % for threshold maps
    high_cut_sum   = zeros(H, W,    'double');
    cut_hit        = zeros(H, W,    'double');   % weight for map averaging

    for b_idx = 1:n_blks
        res = block_results{b_idx};
        blk_data = res{1};
        cut_lo   = double(res{2});   cut_hi = double(res{3});
        r1 = double(res{4});  r2 = double(res{5});
        c1 = double(res{6});  c2 = double(res{7});

        accumulator(r1:r2, c1:c2, :) = accumulator(r1:r2, c1:c2, :) + blk_data;
        hit_count(r1:r2, c1:c2)      = hit_count(r1:r2, c1:c2) + 1;

        % Store threshold value at block centre
        r_mid = round((r1 + r2) / 2);
        c_mid = round((c1 + c2) / 2);
        low_cut_sum(r_mid,  c_mid) = low_cut_sum(r_mid,  c_mid)  + cut_lo;
        high_cut_sum(r_mid, c_mid) = high_cut_sum(r_mid, c_mid)  + cut_hi;
        cut_hit(r_mid, c_mid)      = cut_hit(r_mid, c_mid) + 1;
    end

    % Pixel-wise average (keep only real part — blood reconstructions are real)
    filtered_Data = real(accumulator ./ repmat(max(hit_count,1), 1, 1, T));

    elapsed = toc(t_wall);
    fprintf('%s  Total time: %.1f s  (%.1f blocks/s)\n', ...
        pfx, elapsed, n_blks / elapsed);

    % =========================================================================
    % 10. DIAGNOSTICS
    % =========================================================================
    % Build sparse threshold centre maps and interpolate to full resolution
    low_sparse  = low_cut_sum  ./ max(cut_hit, 1);
    high_sparse = high_cut_sum ./ max(cut_hit, 1);
    low_sparse(cut_hit == 0)  = NaN;
    high_sparse(cut_hit == 0) = NaN;

    [zz, xx] = ndgrid(1:H, 1:W);
    has_data  = ~isnan(low_sparse);
    if any(has_data(:))
        low_full  = griddata(xx(has_data), zz(has_data), ...
                        low_sparse(has_data),  xx, zz, 'linear');
        high_full = griddata(xx(has_data), zz(has_data), ...
                        high_sparse(has_data), xx, zz, 'linear');
        fill_lo   = median(low_sparse(:), 'omitnan');
        fill_hi   = median(high_sparse(:), 'omitnan');
        low_full(isnan(low_full))   = fill_lo;
        high_full(isnan(high_full)) = fill_hi;
    else
        low_full  = zeros(H, W);
        high_full = zeros(H, W);
    end

    diagnostics.low_cutoff_map    = low_full;
    diagnostics.high_cutoff_map   = high_full;
    diagnostics.n_blocks          = n_blks;
    diagnostics.block_size_px     = [blk_z, blk_x];
    diagnostics.step_px           = [step_z, step_x];
    diagnostics.overlap_pct       = [ach_ov_z, ach_ov_x];
    diagnostics.elapsed_sec       = elapsed;
    diagnostics.threshold_method  = method;
    if has_pixel_size
        diagnostics.block_size_mm = [blk_z*pixel_z_mm, blk_x*pixel_x_mm];
    else
        diagnostics.block_size_mm = [NaN, NaN];
    end

    % =========================================================================
    % 11. OPTIONAL THRESHOLD MAP VISUALISATION
    % =========================================================================
    if opt.PlotThresholdMaps
        figure('Name', 'SVD_blockwise: Adaptive Threshold Maps');
        subplot(1,2,1);
        imagesc(low_full); colorbar; colormap(hot); axis image;
        title({'Low-order cutoff (tissue/blood)', sprintf('Median: %.0f', fill_lo)});
        xlabel('X (px)'); ylabel('Z (px)');
        subplot(1,2,2);
        imagesc(high_full); colorbar; colormap(hot); axis image;
        title({'High-order cutoff (blood/noise)', sprintf('Median: %.0f', fill_hi)});
        xlabel('X (px)'); ylabel('Z (px)');
        sgtitle(sprintf('Method: %s  |  Block: [%dx%d px]  |  Overlap: %.0f%%', ...
            method, blk_z, blk_x, mean([ach_ov_z, ach_ov_x])));
    end

end  % ====== END MAIN FUNCTION ======


% =============================================================================
%  LOCAL: DopplerGradient thresholds (Song et al. 2017)
% =============================================================================
function [cut_lo, cut_hi] = blk_DopplerGradient( ...
        S_vec, V, framerate, freq_thr, mp_sigma, grad_pct, ...
        use_doppler, max_tissue_frac, K)

    % --- Cutoff 1A: inflection of dB-scale singular value curve ---
    sv_db    = 20 * log10(S_vec / (S_vec(1) + eps));
    grad2    = diff(diff(sv_db));
    n_search = max(3, round(0.60 * K));
    g2_abs   = abs(grad2(1:min(n_search, end)));
    flat_idx = find(g2_abs < grad_pct * max(g2_abs), 1, 'first');
    if isempty(flat_idx)
        cutoff_1A = max(2, round(0.15 * K));
    else
        cutoff_1A = flat_idx + 1;   % +1: grad2 shifts index by 2
    end

    % --- Cutoff 1B: Doppler frequency crosses tissue threshold ---
    if use_doppler && size(V, 1) > 2
        dfreqs    = doppler_freqs(V, framerate);
        idx_1B    = find(dfreqs > freq_thr, 1, 'first');
        cutoff_1B = nonempty_or(idx_1B, max(1, round(0.10 * K)));
    else
        cutoff_1B = 1;
    end

    % Low cutoff = most conservative (higher) of the two tissue indicators
    cut_lo = max(cutoff_1A, cutoff_1B);
    cut_lo = min(cut_lo, round(max_tissue_frac * K));
    cut_lo = max(cut_lo, 1);

    % --- High cutoff: Marchenko-Pastur noise floor ---
    cut_hi = mp_high_cutoff(S_vec, V, framerate, cut_lo, mp_sigma, K);
end


% =============================================================================
%  LOCAL: SSM-based thresholds (Baranger 2018/2023, per-block)
% =============================================================================
function [cut_lo, cut_hi] = blk_SSM(U, S_vec, V, framerate, ...
        mp_sigma, max_tissue_frac, K)
    try
        SSM = corr(abs(U));                          % [K x K] Pearson corr
        [blo, bhi] = SVD_SSM_AutoTh(SSM);
        cut_lo = min(blo, bhi);
        cut_hi = max(blo, bhi);
        cut_lo = max(1, min(cut_lo, round(max_tissue_frac * K)));
    catch
        % Fallback if SVD_SSM_AutoTh fails on this block
        cut_lo = max(2, round(0.10 * K));
        cut_hi = round(0.80 * K);
    end

    % Refine high cutoff with Marchenko-Pastur — take the more conservative
    mp_hi  = mp_high_cutoff(S_vec, V, framerate, cut_lo, mp_sigma, K);
    cut_hi = min(cut_hi, mp_hi);
end


% =============================================================================
%  LOCAL: Marchenko-Pastur high-order cutoff (blood/noise boundary)
% =============================================================================
function cut_hi = mp_high_cutoff(S_vec, V, framerate, cut_lo, mp_sigma, K)
% Linear fit to log(SV) tail identifies noise floor; deviation above it = blood.

    % --- Step 1: pre-cutoff from Doppler noise floor ---
    noise_tail = max(cut_lo + 1, round(0.80 * K));
    pre_cut    = noise_tail;

    if noise_tail < K
        try
            dfreqs     = doppler_freqs(V, framerate);
            noise_freq = mean(dfreqs(noise_tail:end));
            for k = K:-1:(cut_lo+1)
                if dfreqs(k) > noise_freq * 1.5
                    pre_cut = k;
                    break;
                end
            end
        catch
            pre_cut = noise_tail;
        end
    end

    pre_cut = min(pre_cut, K);

    % --- Step 2: linear fit to log(SV) from pre_cut to end ---
    fit_idx = (pre_cut:K)';
    if numel(fit_idx) < 4
        cut_hi = max(cut_lo + 1, round(cut_lo + 0.7*(K - cut_lo)));
        return;
    end

    log_sv  = log(S_vec(fit_idx) + eps);
    poly_c  = polyfit(fit_idx, log_sv, 1);
    fitted  = polyval(poly_c, fit_idx);
    nstd    = std(log_sv - fitted) + eps;

    % --- Step 3: last component above the noise linear fit by mp_sigma ---
    sr   = cut_lo : pre_cut;
    res  = log(S_vec(sr) + eps) - polyval(poly_c, sr);
    dev  = find(res > mp_sigma * nstd);

    if isempty(dev)
        cut_hi = max(pre_cut - 1, cut_lo + 1);
    else
        cut_hi = sr(dev(end));
    end
    cut_hi = max(cut_lo + 1, min(cut_hi, K));
end


% =============================================================================
%  LOCAL: Lag-one IQ autocorrelation Doppler frequency estimator
% =============================================================================
function freqs = doppler_freqs(V, framerate)
% Mean Doppler frequency per right singular vector via lag-one autocorrelation.
% Evans & McDicken (2000). V: [T x K], returns [K x 1] in Hz.

    [~, K] = size(V);
    freqs  = zeros(K, 1);
    for k = 1:K
        v = V(:, k);
        R1 = sum(v(2:end) .* conj(v(1:end-1)));
        freqs(k) = abs(framerate / (2*pi) * angle(R1));
    end
end


% =============================================================================
%  UTILITY HELPERS
% =============================================================================
function val = get_nested(s, fpath, default)
% Safe nested struct access using dot-path string.
    parts = strsplit(fpath, '.');
    val   = default;
    cur   = s;
    for i = 1:numel(parts)
        if isstruct(cur) && isfield(cur, parts{i})
            cur = cur.(parts{i});
        else
            return;
        end
    end
    if ~isempty(cur)
        val = cur;
    end
end

function v = nonempty_or(x, fallback)
    if isempty(x), v = fallback; else, v = x; end
end

function s = yesno(tf)
    if tf, s = 'OK'; else, s = 'VIOLATION — block too small'; end
end
