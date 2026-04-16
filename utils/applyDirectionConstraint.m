function validTracks = applyDirectionConstraint(tracks, params)
% =========================================================================
% FUNCTION: applyDirectionConstraint
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   A Quality Control (QC) filter that rejects tracks exhibiting erratic, 
%   physically impossible angular turns (e.g., sharp 180-degree reversals 
%   between consecutive frames).
%   - Advantages: Highly effective at eliminating false tracks formed by 
%     noise spikes or "identity swapping" (where the tracker jumps between 
%     two passing bubbles).
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Validation: Aborts if the constraint is disabled in params.
%   2. Vectorization: For each track, step-by-step 2D spatial vectors are 
%      computed using `diff(path)`.
%   3. Angular Calculation: The angle theta between consecutive vectors (v1, v2) 
%      is calculated using the normalized dot product formula: 
%      theta = acos(dot(v1, v2) / (norm(v1) * norm(v2))).
%   4. Numerical Stability: The dot product ratio is explicitly clamped between 
%      [-1, 1] to prevent complex-number outputs from `acos` due to floating-point 
%      inaccuracies.
%   5. Rejection: If the angle exceeds `params.qc.max_angle_change_deg`, the 
%      track is instantly flagged for removal and evaluation breaks early.
%
% SYNTAX OPTIONS:
%   validTracks = applyDirectionConstraint(tracks, params)
%
% EXAMPLES:
%   % Remove zigzagging tracks caused by false positive noise linking:
%   smooth_tracks = applyDirectionConstraint(rawTracks, params);
%
% INPUTS:
%   tracks - (Type: Struct Array) Array containing the trajectory '.path' (Nx2).
%   params - (Type: Struct) QC configuration. Requires:
%            .qc.max_angle_change_deg (double) Maximum turn angle in degrees.
%
% OUTPUTS:
%   validTracks - (Type: Struct Array) Subset of tracks meeting directional criteria.
%
% AUTHOR: Grigori Shapiro
% =========================================================================

    % Filters tracks based on a maximum allowed change in direction.
    if isempty(tracks) || ~params.qc.enable_direction_constraint
        validTracks = tracks;
        return;
    end
    
    max_angle_rad = deg2rad(params.qc.max_angle_change_deg);
    keep_mask = true(1, length(tracks));
    
    for i = 1:length(tracks)
        if tracks(i).length < 3
            continue; % Need at least 3 points to define two vectors
        end
        
        path = tracks(i).path;
        vectors = diff(path, 1, 1); % Get step-by-step velocity vectors
        
        for j = 1:(size(vectors, 1) - 1)
            v1 = vectors(j, :);
            v2 = vectors(j+1, :);
            
            % Calculate angle using dot product formula
            cos_theta = dot(v1, v2) / (norm(v1) * norm(v2));
            angle = acos(max(-1, min(1, cos_theta))); % Clamp for numerical stability
            
            if angle > max_angle_rad
                keep_mask(i) = false;
                break; % No need to check the rest of this track
            end
        end
    end
    
    validTracks = tracks(keep_mask);
end