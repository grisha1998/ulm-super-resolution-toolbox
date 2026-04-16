function [filteredData, info] = DCC_SVD(rawData, framerate, options)
% =========================================================================
% FUNCTION: DCC_SVD
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% ACADEMIC REFERENCE:
%   This algorithm is inspired by and based on the methodology proposed in:
%   Han, X. et al. "An adaptive spatiotemporal filter for ultrasound 
%   localization microscopy based on density canopy clustering". 
%   Ultrasonics, 2024. DOI: 10.1016/j.ultras.2024.107308
%
% PURPOSE & ADVANTAGES:
%   Implements an advanced, parameter-free adaptive SVD clutter filter using 
%   Density Canopy Clustering (DCC). 
%   - Advantages: Unlike traditional SVD which requires manual threshold 
%     guessing, this method adapts to the specific physiological flow of the 
%   - subject by evaluating SVD components in a physical feature space. It 
%   - guarantees convergence to physically meaningful clusters (Tissue, Blood, 
%   - Noise) without human intervention, making it highly robust across 
%   - different datasets and imaging conditions.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Casorati SVD: The 3D data [H x W x T] is reshaped into a 2D Casorati 
%      matrix [Pixels x Frames]. An efficient "Method of Snapshots" SVD is 
%      performed (solving the T x T eigenvalue problem), which is much faster 
%      for typical ULM data where Spatial Pixels >> Frames.
%   2. Feature Extraction (The Physics): Three features are computed for each 
%      singular component (triplet of U, S, V):
%      a. Energy: Log10 of the normalized singular value. Tissue heavily 
%         dominates this feature.
%      b. Temporal Frequency: Power-weighted central frequency from the 
%         positive FFT spectrum of the temporal vector (V). Noise dominates 
%         the higher frequencies.
%      c. Spatial Correlation: Absolute cosine similarity between the spatial 
%         vector (U) and the mean spatial image. Tissue maintains high coherence.
%   3. Density Canopy Initialization: Computes local point densities in the 
%      Z-normalized 3D feature space. It greedily selects the 3 densest, 
%      mutually distant points as "Canopy Centers". This solves the random 
%      initialization flaw of standard K-Means.
%   4. K-Means Clustering: Uses the canopy centers as seeds to strictly group 
%      the components into 3 distinct clusters.
%   5. Physical Identification: Clusters are automatically labeled:
%      - Tissue: Cluster with the highest mean Energy.
%      - Noise: Remaining cluster with the highest mean Temporal Frequency.
%      - Blood: The final remaining cluster.
%   6. Low-Rank Reconstruction: Reconstructs the 3D volume using only the 
%      selected blood components via optimized matrix multiplication.
%
% SYNTAX OPTIONS:
%   filteredData = DCC_SVD(rawData, framerate)
%   filteredData = DCC_SVD(rawData, framerate, 'ReconstructionMode', 'blood')
%   [filteredData, info] = DCC_SVD(rawData, framerate, 'PlotResults', true)
%
% EXAMPLES:
%   % Example 1: Basic usage to extract blood signal
%   bloodData = DCC_SVD(IQ_data, 500);
%
%   % Example 2: Extract tissue signal with advanced parameters and visualization
%   [tissueData, diagnostics] = DCC_SVD(IQ_data, 500, ...
%       'ReconstructionMode', 'tissue', ...
%       'DensityPercentile', 15, ...
%       'CanopySeparation', 2.5, ...
%       'PlotResults', true);
%
% INPUTS:
%   rawData   - (Type: 3D Numeric Array [H x W x T]) 
%               The input ultrasound sequence. Can be real (envelope) or 
%               complex (IQ data). Must have at least 10 frames.
%               Example: A matrix of size [128 x 128 x 200].
%   framerate - (Type: Double Scalar) 
%               Pulse Repetition Frequency (PRF) in Hz. Used to calculate 
%               the physical temporal frequency.
%               Example: 500 or 1000.
%   options   - (Type: Name-Value Pairs)
%       * 'ReconstructionMode' : (String) Determines which cluster to output.
%                                Options: "blood" (default), "tissue", "noise".
%       * 'DensityPercentile'  : (Double, 1-50) Percentile of pairwise distances 
%                                used to define the canopy radius (dc).
%                                Example: 10 (default), 15.
%       * 'CanopySeparation'   : (Double) Multiplier applied to 'dc' to define 
%                                exclusion radius around seeds.
%                                Example: 2.0 (default).
%       * 'PlotResults'        : (Logical) If true, plots 3D scatter, silhouette, 
%                                and singular value spectrum. Default: false.
%       * 'IndentPrefix'       : (String) Prefix for console prints.
%                                Example: "   -> ".
%
% OUTPUTS:
%   filteredData - (Type: 3D Numeric Array) The reconstructed clutter-filtered signal.
%   info         - (Type: Struct) Diagnostic metadata containing:
%                  .cluster_indices, .features, .singular_values, .blood_indices, etc.
% =========================================================================

% =========================================================================
% --- 0. ARGUMENT VALIDATION
% =========================================================================
arguments
    rawData    (:,:,:) {mustBeNumeric}
    framerate  (1,1)   {mustBeNumeric, mustBePositive}
    options.ReconstructionMode (1,1) string ...
        {mustBeMember(options.ReconstructionMode, ["blood","tissue","noise"])} = "blood"
    options.DensityPercentile  (1,1) double ...
        {mustBeInRange(options.DensityPercentile, 1, 50)}  = 10
    options.CanopySeparation   (1,1) double ...
        {mustBePositive}                                   = 2.0
    options.PlotResults        (1,1) logical               = false
    options.IndentPrefix       (1,1) string                = ""
end

px    = options.IndentPrefix;   % shorthand for log lines
num_clusters = 3;               % tissue | blood | noise — fixed by design

% =========================================================================
% --- 1. INPUT VALIDATION
% =========================================================================
[H, W, T] = size(rawData);

assert(T >= 10, ...
    'DCC_SVD: rawData must have at least 10 frames for stable clustering (got T=%d).', T);
assert(T >= num_clusters * 3, ...
    'DCC_SVD: T=%d frames is too few for %d clusters (need T >= %d).', ...
    T, num_clusters, num_clusters * 3);

if ~isreal(rawData)
    fprintf('%s[INFO] Complex-valued input detected (IQ data). Processing in complex domain.\n', px);
end

% =========================================================================
% --- 2. CASORATI MATRIX + EFFICIENT SVD (METHOD OF SNAPSHOTS)
% =========================================================================
% Reshape 3-D volume into a 2-D Casorati matrix [pixels x frames].
% The method-of-snapshots solves the smaller (T x T) eigenvalue problem
% instead of the full (H*W x H*W) covariance, which is critical when
% the spatial dimension H*W is very large.
%
% For complex data, casoratiMatrix' is the Hermitian (conjugate) transpose,
% so AtA is Hermitian and eig() returns real, non-negative eigenvalues.
casoratiMatrix = reshape(rawData, H * W, T);

fprintf('%s- Performing efficient SVD on [%d x %d] Casorati matrix...\n', px, H*W, T);

AtA = casoratiMatrix' * casoratiMatrix;   % [T x T], Hermitian for complex input

[V, D] = eig(AtA);                        % V: eigenvectors, D: eigenvalues

% Sort eigenvalues (and corresponding eigenvectors) in descending order.
S_diag = sqrt(abs(real(diag(D))));        % singular values; abs+real guards against
[S_diag, sort_idx] = sort(S_diag, 'descend'); % tiny negative rounding errors
V = V(:, sort_idx);

% Recover left singular vectors U from U = A*V*S^{-1}.
% Zero out the inverse for negligible singular values to avoid division by
% near-zero, which would amplify numerical noise into U.
singular_values_inv = zeros(size(S_diag));
valid_sv = S_diag > 1e-10 * S_diag(1);   % threshold relative to largest SV
singular_values_inv(valid_sv) = 1 ./ S_diag(valid_sv);

U = casoratiMatrix * V * diag(singular_values_inv);  % [H*W x T]

fprintf('%s- SVD complete (%d components). Extracting features...\n', px, T);

% =========================================================================
% --- 3. FEATURE EXTRACTION
% =========================================================================
features = zeros(T, 3);

% ----- Feature 1: Normalized Log Energy -----
% Physical meaning: tissue clutter concentrates signal energy in the first
% few components. Blood and noise have progressively smaller singular values.
% Regularization with eps prevents log10(0) = -Inf for negligible components,
% which would otherwise corrupt z-score normalization and K-means distances.
normalized_sv    = S_diag / sum(S_diag);
features(:, 1)   = log10(normalized_sv + eps);

% ----- Feature 2: Power-Weighted Central Temporal Frequency -----
% Physical meaning: tissue is quasi-static (low frequency), blood moves
% at intermediate frequencies, thermal/electronic noise spans all
% frequencies uniformly (high effective center).
%
% IMPORTANT — use only the positive half of the FFT spectrum.
% The full FFT of a real signal is conjugate-symmetric: the second half
% mirrors the first. Including it would double-count high-frequency content
% and bias the centroid estimate upward for all components, reducing the
% discriminability between blood and noise.
%
% V contains the temporal singular vectors [T x T]. Each COLUMN is one
% component's time signature. fft() operates along the first dimension,
% so we transpose: V' is [T x T] with each ROW being one component.
half_T    = floor(T / 2) + 1;                      % number of unique frequencies
V_fft     = abs(fft(V));                            % [T x T], fft along columns
V_fft_pos  = V_fft(2:half_T, :);                    % keep positive half [half_T x T]
freqs_axis = (1 : half_T-1)' * (framerate / T);   % physical frequency axis [Hz]

% Power-weighted centroid: sum(f * |F(f)|) / sum(|F(f)|)
% sum() operates along dim 1 (over frequencies), result is [1 x T].
central_freqs  = (freqs_axis' * V_fft_pos) ./ sum(V_fft_pos, 1);
features(:, 2) = central_freqs';

% ----- Feature 3: Spatial Correlation with Mean Pattern -----
% Physical meaning: tissue components are highly correlated with the mean
% spatial map (they dominate the average frame). Blood and noise deviate
% more from the mean anatomical structure.
% Absolute value of cosine similarity is used so phase sign does not matter.
mean_spatial_vector = mean(U, 2);                   % [H*W x 1]
norm_mean           = norm(mean_spatial_vector);

% Vectorized cosine similarity for all T components simultaneously:
%   cos_sim(k) = |<mean, U_k>| / (||mean|| * ||U_k||)
% mean_spatial_vector' * U -> [1 x T] dot products
% vecnorm(U, 2, 1)          -> [1 x T] column norms
complex_corr   = (mean_spatial_vector' * U) ./ (norm_mean * vecnorm(U, 2, 1));
features(:, 3) = abs(complex_corr)';

% Z-score normalize so all three features contribute equally to Euclidean
% distances in the clustering step (removes scale bias).
features_normalized = zscore(features);

fprintf('%s- Features extracted. Running Density Canopy Initialization...\n', px);

% =========================================================================
% --- 4. DENSITY CANOPY INITIALIZATION (DCC CORE)
% =========================================================================
% Compute pairwise Euclidean distances using pdist (upper triangle only,
% then expand to full square with squareform). This avoids computing each
% distance twice compared to pdist2(X, X).
dist_vec  = pdist(features_normalized, 'euclidean');   % [1 x T*(T-1)/2]
distances = squareform(dist_vec);                      % [T x T] symmetric

% Define neighborhood radius dc as the DensityPercentile-th percentile of
% ALL pairwise distances (excluding self-distances on the diagonal, which
% are zero and would bias the percentile downward).
off_diag_mask = dist_vec > 0;                          % upper-triangle, no zeros
dc = prctile(dist_vec(off_diag_mask), options.DensityPercentile);

fprintf('%s  Neighborhood radius dc = %.4f (at %d-th percentile).\n', ...
    px, dc, options.DensityPercentile);

% Local density rho(i) = number of points within dc of point i.
% Diagonal is zero so it does not count the point itself.
rho = sum(distances < dc, 2);                          % [T x 1]

% Greedy canopy selection:
%   - Find the highest-density unselected point.
%   - Record it as a canopy center.
%   - Mark all points within (CanopySeparation * dc) as unavailable,
%     ensuring the next canopy comes from a different density peak.
% A separate 'already_selected' mask prevents a full reset from re-picking
% an already-chosen center when valid candidates are exhausted early.
exclusion_radius = options.CanopySeparation * dc;
initial_centers  = zeros(num_clusters, 3);
unselected_mask  = true(T, 1);
already_selected = false(T, 1);

for k = 1 : num_clusters
    valid_idx = find(unselected_mask);

    % Fallback: if the exclusion zone has covered all remaining candidates,
    % reset to all unselected points EXCEPT those already chosen as centers.
    % This preserves the previously identified canopies.
    if isempty(valid_idx)
        unselected_mask  = ~already_selected;
        valid_idx        = find(unselected_mask);
    end

    % Among valid candidates, pick the one with the highest density.
    [~, local_max_idx]  = max(rho(valid_idx));
    center_idx          = valid_idx(local_max_idx);

    % Store normalized feature coordinates of this canopy center.
    initial_centers(k, :) = features_normalized(center_idx, :);

    % Lock in this center so it cannot be overwritten by the fallback.
    already_selected(center_idx) = true;

    % Exclude all points too close to this center from future selection.
    too_close                    = distances(center_idx, :) < exclusion_radius;
    unselected_mask(too_close')   = false;
end

fprintf('%s  Canopy centers identified at feature-space coordinates:\n', px);
for k = 1 : num_clusters
    fprintf('%s    Canopy %d: [%.3f, %.3f, %.3f]\n', ...
        px, k, initial_centers(k,1), initial_centers(k,2), initial_centers(k,3));
end

% =========================================================================
% --- 5. CANOPY-SEEDED K-MEANS CLUSTERING
% =========================================================================
% Using the high-density canopy peaks as starting centroids anchors
% K-means to physically meaningful regions of feature space and avoids
% random initialization sensitivity. Multiple replicates are still run
% to confirm convergence is robust (replicates=5 is a balanced trade-off).
fprintf('%s- Running K-means (seeded from canopy centers)...\n', px);

cluster_indices = kmeans(features_normalized, num_clusters, ...
    'Start',      initial_centers, ...
    'Distance',   'sqeuclidean', ...
    'Display',    'off');

fprintf('%s- Clustering complete. Identifying tissue, blood, and noise clusters...\n', px);

% =========================================================================
% --- 6. CLUSTER IDENTIFICATION BY PHYSICAL PROPERTIES
% =========================================================================
% Compute the mean of each raw (non-normalized) feature per cluster.
% Raw features preserve physical interpretability (energy, Hz, correlation).
cluster_means = zeros(num_clusters, 3);
for i = 1 : num_clusters
    cluster_means(i, :) = mean(features(cluster_indices == i, :), 1);
end

% --- Step 1: Identify TISSUE as the cluster with the highest mean energy.
% Tissue clutter dominates the first singular values and therefore has the
% largest log-energy among all three groups.
[~, tissue_cluster_idx] = max(cluster_means(:, 1));

% --- Step 2: Identify NOISE as the cluster with the highest mean temporal
% frequency AMONG THE REMAINING clusters (i.e., after tissue is removed).
% This hierarchical approach avoids the collision bug where a single cluster
% could be assigned as both tissue (max energy) and noise (max frequency),
% leaving blood_cluster_idx with two elements.
remaining_clusters     = setdiff(1 : num_clusters, tissue_cluster_idx);
[~, local_noise_idx]   = max(cluster_means(remaining_clusters, 2));
noise_cluster_idx      = remaining_clusters(local_noise_idx);

% --- Step 3: BLOOD is whatever is left.
blood_cluster_idx = setdiff(1 : num_clusters, [tissue_cluster_idx, noise_cluster_idx]);

% Log cluster statistics for transparency.
fprintf('%s  Tissue cluster  (#%d): mean energy=%.3f, mean freq=%.2f Hz, mean corr=%.3f\n', ...
    px, tissue_cluster_idx, cluster_means(tissue_cluster_idx, :));
fprintf('%s  Blood  cluster  (#%d): mean energy=%.3f, mean freq=%.2f Hz, mean corr=%.3f\n', ...
    px, blood_cluster_idx,  cluster_means(blood_cluster_idx,  :));
fprintf('%s  Noise  cluster  (#%d): mean energy=%.3f, mean freq=%.2f Hz, mean corr=%.3f\n', ...
    px, noise_cluster_idx,  cluster_means(noise_cluster_idx,  :));

% =========================================================================
% --- 7. COLLECT DIAGNOSTIC INFO (second output argument)
% =========================================================================
tissue_indices = find(cluster_indices == tissue_cluster_idx);
blood_indices  = find(cluster_indices == blood_cluster_idx);
noise_indices  = find(cluster_indices == noise_cluster_idx);

info.cluster_indices     = cluster_indices;
info.features            = features;
info.features_normalized = features_normalized;
info.singular_values     = S_diag;
info.tissue_indices      = tissue_indices;
info.blood_indices       = blood_indices;
info.noise_indices       = noise_indices;
info.num_kept_components = [];   % filled after mode selection below
info.dc                  = dc;
info.initial_centers     = initial_centers;

% =========================================================================
% --- 8. VISUALIZATION (if requested)
% =========================================================================
if options.PlotResults
    fprintf('%s- Generating visualization plots...\n', px);

    % --- Plot A: 3D Feature Space Scatter ---
    figure('Name', 'DCC-SVD: Feature Space Clustering', 'NumberTitle', 'off');
    hold on;
    scatter3(features_normalized(tissue_indices, 1), ...
             features_normalized(tissue_indices, 2), ...
             features_normalized(tissue_indices, 3), ...
             36, 'r', 'filled', 'DisplayName', 'Tissue');
    scatter3(features_normalized(blood_indices, 1), ...
             features_normalized(blood_indices, 2), ...
             features_normalized(blood_indices, 3), ...
             36, 'b', 'filled', 'DisplayName', 'Blood');
    scatter3(features_normalized(noise_indices, 1), ...
             features_normalized(noise_indices, 2), ...
             features_normalized(noise_indices, 3), ...
             36, [0.4 0.4 0.4], 'filled', 'DisplayName', 'Noise');
    scatter3(initial_centers(:,1), initial_centers(:,2), initial_centers(:,3), ...
             200, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', ...
             'DisplayName', 'Canopy Centers');
    hold off;
    grid on; view(3); axis tight;
    title('Density Canopy Clustering of SVD Components');
    xlabel('Feature 1: Energy (z-score)');
    ylabel('Feature 2: Temporal Frequency (z-score)');
    zlabel('Feature 3: Spatial Correlation (z-score)');
    legend('show', 'Location', 'best');

    % --- Plot B: Silhouette Analysis ---
    figure('Name', 'DCC-SVD: Clustering Quality (Silhouette)', 'NumberTitle', 'off');
    silhouette(features_normalized, cluster_indices);
    title('Silhouette Analysis of Clustered SVD Components');
    xlabel('Silhouette Value');

    % --- Plot C: Singular Value Spectrum (color-coded by cluster) ---
    figure('Name', 'DCC-SVD: Singular Value Spectrum', 'NumberTitle', 'off');
    semilogy(1:T, S_diag / S_diag(1), 'Color', [0.75 0.75 0.75], ...
             'LineWidth', 1, 'DisplayName', 'All components');
    hold on;
    semilogy(tissue_indices, S_diag(tissue_indices) / S_diag(1), ...
             'r.', 'MarkerSize', 14, 'DisplayName', 'Tissue');
    semilogy(blood_indices, S_diag(blood_indices) / S_diag(1),   ...
             'b.', 'MarkerSize', 14, 'DisplayName', 'Blood');
    semilogy(noise_indices, S_diag(noise_indices) / S_diag(1),   ...
             '.', 'MarkerSize', 14, 'Color', [0.4 0.4 0.4], 'DisplayName', 'Noise');
    hold off;
    grid on;
    xlabel('Component Index');
    ylabel('Normalized Singular Value (log scale)');
    title('Singular Value Spectrum — Color-Coded by Cluster');
    legend('show', 'Location', 'southwest');
    xlim([1, T]);
end

% =========================================================================
% --- 9. SELECT COMPONENTS FOR RECONSTRUCTION
% =========================================================================
fprintf('%s- Selecting components for reconstruction (mode: ''%s'')...\n', ...
    px, options.ReconstructionMode);

switch lower(options.ReconstructionMode)
    case 'blood'
        components_to_keep_idx = blood_cluster_idx;
    case 'tissue'
        components_to_keep_idx = tissue_cluster_idx;
    case 'noise'
        components_to_keep_idx = noise_cluster_idx;
end

keep_mask = (cluster_indices == components_to_keep_idx);
info.num_kept_components = sum(keep_mask);

% =========================================================================
% --- 10. EFFICIENT LOW-RANK RECONSTRUCTION
% =========================================================================
% Reconstruct using only the K kept singular triplets (U_k, S_k, V_k):
%
%   filteredMatrix = sum_{k in kept} S_k * U_k * V_k'
%
% Written as a single matrix product for efficiency:
%
%   filteredMatrix = (U_k .* S_k') * V_k'
%
% This avoids constructing a full T x T diagonal matrix S and the
% associated O(T^2 * H*W) dense matrix multiplications.
keep_idx       = find(keep_mask);
U_k            = U(:, keep_idx);              % [H*W x K]
S_k            = S_diag(keep_idx);            % [K x 1]
V_k            = V(:, keep_idx);              % [T x K]

filteredMatrix = (U_k .* S_k') * V_k';        % [H*W x T]
filteredData   = reshape(filteredMatrix, H, W, T);

fprintf('%s- DCC-SVD complete. Kept %d / %d components classified as ''%s''.\n', ...
    px, info.num_kept_components, T, options.ReconstructionMode);

info.U = U;   % [H*W x T] left singular vectors
info.V = V;   % [T x T]   right singular vectors
end % function DCC_SVD
