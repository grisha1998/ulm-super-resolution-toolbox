function filtered_Data = SVD_filter(rawData, cutoff, options)
% =========================================================================
% FUNCTION: SVD_filter
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Applies a standard, global Singular Value Decomposition (SVD) clutter 
%   filter using hard-coded, user-defined cutoffs.
%   - Advantages: Extremely fast, highly predictable, and mathematically 
%     transparent. It is the gold standard for testing and establishing 
%     a baseline before moving to adaptive or block-wise methods. 
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Reshapes the 3D ultrasound data [H x W x T] into a 2D Casorati 
%      matrix [H*W x T].
%   2. Performs an economical SVD (`svd(..., 'econ')`), separating the 
%      data into Spatial components (U), Energy weights (S), and Temporal 
%      components (V).
%   3. Applies the manually provided `cutoff` array (e.g., [10, 200]).
%      - Indices 1 to (cutoff(1)-1) = Tissue clutter (Discarded).
%      - Indices cutoff(1) to cutoff(2) = Blood signal (Retained).
%      - Indices above cutoff(2) = Thermal noise (Discarded).
%   4. Efficiently reconstructs the 3D volume using reduced-rank matrix 
%      multiplication of only the retained components.
%
% SYNTAX OPTIONS:
%   filt_data = SVD_filter(rawData, [low_cut, high_cut])
%   filt_data = SVD_filter(rawData, [15, 150], 'ReconstructionMode', 'tissue')
%
% EXAMPLES:
%   % Example 1: Standard blood filtering dropping first 15 components
%   bloodData = SVD_filter(IQ_data, [15, 200]);
%
%   % Example 2: Extracting only the tissue background
%   tissueData = SVD_filter(IQ_data, [15, 200], 'ReconstructionMode', 'tissue');
%
% INPUTS:
%   rawData - (Type: 3D Numeric Array) The input ultrasound sequence.
%   cutoff  - (Type: 1x2 Double Array) [low_index, high_index]. 
%             Example: [15, 200] or [8, 180].
%   options - (Type: Name-Value Pairs) 
%       * 'ReconstructionMode' : (String) "blood", "tissue", "noise".
%
% OUTPUTS:
%   filtered_Data - (Type: 3D Numeric Array) The manually filtered sequence.
% =========================================================================

    arguments
        rawData (:,:,:) {mustBeNumeric}
        cutoff (1,2) double
        options.ReconstructionMode (1,1) string {mustBeMember(options.ReconstructionMode, ["blood", "tissue", "noise"])} = "blood"
        options.PlotResults (1,1) logical = false
        options.IndentPrefix (1,1) string = ""
    end
    [H, W, T] = size(rawData);
    
    % --- Store original cutoff values for index definition ---
    cutoff_start = cutoff(1);
    cutoff_end = cutoff(2);
    
    if cutoff_end > T
        cutoff_end = T;
    end
        
    if or(isequal(cutoff_start:cutoff_end, 1:T), cutoff_start < 2)
        filtered_Data = rawData;
        return
    end
    
    casoratiMatrix = reshape(rawData, H * W, T);
    [U,S,V] = svd(casoratiMatrix,'econ'); % 'econ' is more efficient than '0'
    S_diag = diag(S);
    
    % --- Define component indices ---
    tissue_indices = 1:(cutoff_start - 1);
    blood_indices  = cutoff_start:cutoff_end;
    noise_indices  = (cutoff_end + 1):length(S_diag);
    
    switch lower(options.ReconstructionMode)
        case 'blood'
            components_to_keep_idx = blood_indices;
        case 'tissue'
            components_to_keep_idx = tissue_indices;
        case 'noise'
            components_to_keep_idx = noise_indices;
    end

    % --- If no components are selected, return zero matrix ---
    if isempty(components_to_keep_idx)
        filtered_Data = zeros(H, W, T, 'like', rawData);
        return;
    end
    
    % --- Reconstruct Data Using ONLY the Selected Components (Efficient Method) ---
    U_filtered = U(:, components_to_keep_idx);
    S_filtered_small = S(components_to_keep_idx, components_to_keep_idx);
    V_filtered = V(:, components_to_keep_idx);
    
    filteredMatrix = U_filtered * S_filtered_small * V_filtered';
    filtered_Data = reshape(filteredMatrix, H, W, T);
end
