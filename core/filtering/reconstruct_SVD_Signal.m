function filteredData = reconstruct_SVD_Signal(U, S_diag, V, dataDims, cutoff, mode)
% =========================================================================
% FUNCTION: reconstruct_SVD_Signal
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Rapidly reconstructs a 3D spatiotemporal signal from pre-computed SVD 
%   matrices based on a specific cutoff range.
%   - Advantages: Highly optimized for speed. It acts as the "renderer" 
%     for the matrices produced by `run_SVD_Decomposition`.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   To maximize speed and minimize RAM footprint, it does NOT reconstruct 
%   the full T x T diagonal matrix `S`. Instead, it extracts only the 
%   requested columns from `U` and `V`, and performs a reduced-rank matrix 
%   multiplication: X_filt = U(:, idx) * S_reduced * V(:, idx)'.
%   Finally, it reshapes the 2D result back into the 3D video space.
%
% SYNTAX OPTIONS:
%   filt_data = reconstruct_SVD_Signal(U, S_diag, V, dims, [low, high])
%   filt_data = reconstruct_SVD_Signal(U, S_diag, V, dims, [15, 100], 'tissue')
%
% EXAMPLES:
%   bloodData = reconstruct_SVD_Signal(U, S_diag, V, dims, [10, 150], 'blood');
%
% INPUTS:
%   U, S_diag, V - Outputs from run_SVD_Decomposition.
%   dataDims     - (Type: 1x3 Double) Original dimensions [H, W, T].
%   cutoff       - (Type: 1x2 Double) Desired [low, high] cutoff indices.
%   mode         - (Type: String) 'blood' (default), 'tissue', or 'noise'.
%
% OUTPUTS:
%   filteredData - (Type: 3D Numeric Array) The instantly reconstructed sequence.
% =========================================================================

    if nargin < 6, mode = 'blood'; end

    T = dataDims(3);
    cutoff_start = cutoff(1);
    cutoff_end = cutoff(2);
    
    % 1. Define Indices based on Cutoff
    tissue_indices = 1:(cutoff_start - 1);
    blood_indices  = cutoff_start:cutoff_end;
    noise_indices  = (cutoff_end + 1):length(S_diag);
    
    % 2. Select Components based on Mode
    switch lower(mode)
        case 'blood'
            idx_to_keep = blood_indices;
        case 'tissue'
            idx_to_keep = tissue_indices;
        case 'noise'
            idx_to_keep = noise_indices;
        otherwise
            error('Unknown reconstruction mode');
    end
    
    % Safety check for indices bounds
    idx_to_keep = idx_to_keep(idx_to_keep <= length(S_diag) & idx_to_keep >= 1);

    if isempty(idx_to_keep)
        filteredData = zeros(dataDims);
        return;
    end

    % 3. Fast Reconstruction
    %    We reconstruct only using the relevant columns.
    %    Formula: X_filt = U(:, idx) * S(idx, idx) * V(:, idx)'
    
    S_reduced = diag(S_diag(idx_to_keep));
    
    % Perform matrix multiplication
    filteredMatrix = U(:, idx_to_keep) * S_reduced * V(:, idx_to_keep)';
    
    % 4. Reshape back to 3D image
    filteredData = reshape(filteredMatrix, dataDims(1), dataDims(2), dataDims(3));
end