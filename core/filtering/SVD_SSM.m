function filtered_Data = SVD_SSM(rawData, options)
% =========================================================================
% FUNCTION: SVD_SSM
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% ACADEMIC REFERENCE:
%   [1] J. Baranger, et al. "Fast Thresholding of SVD Clutter Filter Using
%       the Spatial Similarity Matrix and a Sum-Table Algorithm." IEEE
%       T-UFFC, 2023.
%   [2] J. Baranger, et al. "Adaptive Spatiotemporal SVD Clutter Filtering
%       for Ultrafast Doppler Imaging Using Similarity of Spatial Singular
%       Vectors." IEEE T-MI, 2018.
%
% PURPOSE & ADVANTAGES:
%   Implements an adaptive global SVD clutter filter using a Spatial 
%   Similarity Matrix (SSM).
%   - Advantages: Standard global SVD requires the user to manually guess 
%     an arbitrary index (e.g., 15) to separate tissue from blood. This 
%     method calculates the physical, structural coherence of the image 
%     to automatically find the exact cutoff point. It adapts dynamically 
%     if the amount of clutter changes between datasets.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. SVD Decomposition: The sequence is reshaped and decomposed (U, S, V).
%   2. SSM Calculation: Computes the Pearson correlation matrix between the 
%      absolute values of all spatial singular vectors (U). 
%      SSM(i,j) = corr(|U_i|, |U_j|).
%   3. Thresholding Logic (`SVD_SSM_AutoTh`): Tissue components represent the 
%      static, dominant anatomy and are therefore highly spatially correlated 
%      with one another. Blood and noise represent transient/random signals 
%      and lack this correlation. The algorithm analyzes the SSM to find the 
%      sharp drop in similarity, defining the exact low/high cutoff indices.
%   4. Classification: Components are split into tissue (1 to LowTh), blood 
%      (LowTh to HighTh), and noise (HighTh to End).
%   5. Signal Reconstruction: The data is reconstructed using only the target 
%      components.
%
% SYNTAX OPTIONS:
%   filt_data = SVD_SSM(rawData)
%   filt_data = SVD_SSM(rawData, 'ReconstructionMode', 'tissue')
%
% EXAMPLES:
%   % Example 1: Auto-filter blood
%   filteredData = SVD_SSM(IQ_data);
%
%   % Example 2: Extract noise profile
%   noiseData = SVD_SSM(IQ_data, 'ReconstructionMode', 'noise', 'IndentPrefix', '>> ');
%
% INPUTS:
%   rawData - (Type: 3D Numeric Array [H x W x T]) The raw/IQ ultrasound sequence.
%   options - (Type: Name-Value Pairs) 
%       * 'ReconstructionMode' : (String) "blood" (default), "tissue", or "noise".
%       * 'IndentPrefix'       : (String) Prefix for formatting console output.
%
% OUTPUTS:
%   filtered_Data - (Type: 3D Numeric Array) The reconstructed data.
% =========================================================================

    arguments
        rawData (:,:,:) {mustBeNumeric}
        options.ReconstructionMode (1,1) string {mustBeMember(options.ReconstructionMode, ["blood", "tissue", "noise"])} = "blood"
        options.PlotResults (1,1) logical = false
        options.IndentPrefix (1,1) string = ""
    end

    % --- Step 1: Perform a full SVD on the raw data to get U, S, and V ---
    [H, W, T] = size(rawData);
    casoratiMatrix = reshape(rawData, H * W, T);
    [U,S,V] = svd(casoratiMatrix,0);
    S_diag = diag(S);

    % --- Step 2: Calculate SSM ---
    SSM = corr(abs(U)); %

    % --- Step 3: Call the core thresholding function with the calculated SSM ---
    [bloodLowTh, bloodHighTh] = SVD_SSM_AutoTh(SSM);

    % --- Step 4: Classify all components based on the two thresholds ---
    % Ensure thresholds are sorted correctly for comparison
    if bloodLowTh > bloodHighTh
        [bloodLowTh, bloodHighTh] = deal(bloodHighTh, bloodLowTh); % Swap them
    end
    
    tissue_indices = 1:bloodLowTh;
    blood_indices  = bloodLowTh:bloodHighTh;
    noise_indices  = bloodHighTh:length(S_diag);
    
    switch lower(options.ReconstructionMode)
        case 'blood'
            components_to_keep_idx = blood_indices;
        case 'tissue'
            components_to_keep_idx = tissue_indices;
        case 'noise'
            components_to_keep_idx = noise_indices;
    end
    
    % --- 8. Reconstruct Data Using ONLY the Selected Components ---
    if isempty(components_to_keep_idx)
        filtered_Data = zeros(H, W, T, 'like', rawData);
        return;
    end
    
    U_filtered = U(:, components_to_keep_idx);
    S_filtered_small = S(components_to_keep_idx, components_to_keep_idx);
    V_filtered = V(:, components_to_keep_idx);
    
    filteredMatrix = U_filtered * S_filtered_small * V_filtered';
    filtered_Data = reshape(filteredMatrix, H, W, T);
    
end