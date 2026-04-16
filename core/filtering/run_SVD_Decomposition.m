function [U, S_diag, V, dataDims] = run_SVD_Decomposition(rawData)
% =========================================================================
% FUNCTION: run_SVD_Decomposition
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Performs the heavy computational lifting of SVD independently of the 
%   reconstruction phase.
%   - Advantages: When hyperparameter tuning (trying to find the perfect 
%     manual cutoff), running the full SVD takes 95% of the time. This 
%     function calculates the decomposition once. You can then pass the 
%     outputs to `reconstruct_SVD_Signal` 100 times per second to visualize 
%     different cutoffs instantly.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Extracts the data dimensions and reshapes to a Casorati matrix.
%   2. Runs `svd(..., 'econ')` to get U, S, V.
%   3. Extracts just the diagonal of S to save memory.
%
% SYNTAX OPTIONS:
%   [U, S_diag, V, dims] = run_SVD_Decomposition(rawData)
%
% EXAMPLES:
%   [U, S_diag, V, dims] = run_SVD_Decomposition(IQ_data);
%   % Then use reconstruct_SVD_Signal...
%
% INPUTS:
%   rawData - (Type: 3D Numeric Array) The ultrasound sequence.
%
% OUTPUTS:
%   U        - (Type: 2D Double) Spatial singular vectors.
%   S_diag   - (Type: 1D Double) Vector containing the singular values.
%   V        - (Type: 2D Double) Temporal singular vectors.
%   dataDims - (Type: 1x3 Double) Original dimensions [H, W, T].
% =========================================================================

    % 1. Get dimensions
    [H, W, T] = size(rawData);
    dataDims = [H, W, T];

    % 2. Reshape to Casorati Matrix (2D)
    %    Rows = Spatial pixels, Cols = Time frames
    casoratiMatrix = reshape(rawData, H * W, T);

    % 3. Perform Efficient SVD ('econ')
    %    This is the computationally expensive step.
    [U, S, V] = svd(casoratiMatrix, 'econ');
    
    % 4. Extract diagonal values as a vector (saves memory)
    S_diag = diag(S);
    
end