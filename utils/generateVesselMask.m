function mask = generateVesselMask(filteredData, params)
% =========================================================================
% FUNCTION: generateVesselMask
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Automatically creates a binary Region of Interest (ROI) mask by highlighting 
%   vascular structures from the clutter-filtered ultrasound data.
%   - Advantages: Enhances detection performance and tracking safety by 
%     restricting operations strictly to areas with physical blood flow. It 
%     incorporates advanced morphological and contrast enhancement techniques 
%     to isolate vessels even in low-SNR regions.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Base Map Generation: Takes the absolute envelope of the filtered data, 
%      averages it over time (Temporal Mean), and applies a square-root 
%      compression to boost weaker vessels.
%   2. Image Enhancement: Applies the user-selected algorithm:
%      - CLAHE: Adaptive histogram equalization to enhance local contrast.
%      - Top-Hat: Morphological filtering using a disk structuring element 
%        to highlight vessel-like tubular structures and suppress flat backgrounds.
%      - Sharpen: Edge enhancement.
%   3. Gamma Correction: Non-linear scaling (`map .^ gamma`) to manipulate 
%      overall contrast and suppress residual noise floors.
%   4. Thresholding: Normalizes the final enhanced map to [0, 1] and generates 
%      a boolean binary mask based on `params.threshold`.
%
% SYNTAX OPTIONS:
%   mask = generateVesselMask(filteredData, params)
%
% EXAMPLES:
%   % Generate mask using CLAHE enhancement:
%   vesselParams.method = 'CLAHE';
%   vesselParams.strength = 0.5;
%   vesselParams.gamma = 1.0;
%   vesselParams.threshold = 0.15;
%   roiMask = generateVesselMask(filteredData, vesselParams);
%
% INPUTS:
%   filteredData - (Type: 3D Numeric Matrix) Clutter-suppressed ultrasound sequence.
%   params       - (Type: Struct) Contains enhancement parameters:
%                  .method, .strength, .gamma, .threshold.
%
% OUTPUTS:
%   mask         - (Type: 2D Logical Matrix [H x W]) Binary map defining valid vessel regions.
%
% AUTHOR: Grigori Shapiro
% =========================================================================

    % 1. Calculate Base Vessel Map (Temporal Average + Sqrt)
    if ~isreal(filteredData)
        absData = abs(filteredData);
    else
        absData = filteredData;
    end
    
    % Temporal Average
    vMap = mean(absData, 3);
    
    % Square Root for initial dynamic range compression
    vMap = vMap .^ 0.5;
    
    % Normalize to 0-1
    mx = max(vMap(:));
    if mx > 0, vMap = vMap / mx; end
    
    % 2. Apply Selected Enhancement Algorithm
    val = params.strength;
    
    switch params.method
        case 'None'
            procImg = vMap;
            
        case {'CLAHE', 'CLAHE (Local Contrast)'}
            % ClipLimit: 0.001 to 0.041
            clipLim = 0.001 + (val * 0.04); 
            procImg = adapthisteq(vMap, 'ClipLimit', clipLim, 'Distribution', 'rayleigh');
            
        case {'Top-Hat', 'Top-Hat (Vesselness)'}
            % Structural Element Radius: 1 to 11 pixels
            radius = 1 + round(val * 10); 
            se = strel('disk', radius);
            procImg = imtophat(vMap, se);
            % Re-normalize after top-hat
            mx = max(procImg(:)); if mx>0, procImg=procImg/mx; end
            
        case 'Sharpen'
            % Amount: 0 to 2
            amount = val * 2;
            procImg = imsharpen(vMap, 'Radius', 1, 'Amount', amount);
            
        otherwise
            warning('Unknown enhancement method: %s. Using raw map.', params.method);
            procImg = vMap;
    end
    
    % 3. Apply Gamma Correction
    procImg = procImg .^ params.gamma;
    
    % 4. Normalize Final Image (Ensure 0-1 range for thresholding)
    mn = min(procImg(:));
    mx = max(procImg(:));
    if mx > mn
        procImg = (procImg - mn) / (mx - mn);
    end
    
    % 5. Generate Binary Mask
    mask = procImg >= params.threshold;

end