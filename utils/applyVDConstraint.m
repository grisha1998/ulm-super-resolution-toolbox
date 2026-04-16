function validTracks = applyVDConstraint(tracks, params)
% =========================================================================
% FUNCTION: applyVDConstraint
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   A Quality Control (QC) filter based on the Velocity Dispersion (VD) ratio.
%   - Advantages: Highly effective at differentiating true microbubble flow 
%     from stationary, vibrating noise. True bubbles travel a net distance 
%     through the vessel (VD near 1.0). Noise points vibrate randomly in place, 
%     accumulating high path length but near-zero net displacement (VD >> 1.0).
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Validation: Requires a minimum of 2 points to define a path.
%   2. Total Path Length: Calculates the sum of all step-by-step Euclidean 
%      distances across the entire trajectory.
%   3. Net Displacement: Calculates the straight-line Euclidean distance 
%      between the absolute first and absolute last point of the track.
%   4. Ratio Calculation: VD = Total Path Length / Net Displacement. 
%      Safe-guards against division by zero (if net displacement < 1e-6, VD = inf).
%   5. Rejection: Tracks exceeding `params.qc.max_vd_ratio` are removed.
%
% SYNTAX OPTIONS:
%   validTracks = applyVDConstraint(tracks, params)
%
% EXAMPLES:
%   % Remove static jitter noise from the data:
%   flowing_tracks = applyVDConstraint(tracks, params);
%
% INPUTS:
%   tracks - (Type: Struct Array) Contains the trajectory '.path' (Nx2).
%   params - (Type: Struct) QC configuration. Requires:
%            .qc.max_vd_ratio (double) Maximum allowed tortuosity/dispersion.
%
% OUTPUTS:
%   validTracks - (Type: Struct Array) Tracks exhibiting true directional flow.
%
% AUTHOR: Grigori Shapiro
% =========================================================================

    % Filters tracks based on the Velocity Dispersion (VD) ratio.
    if isempty(tracks) || ~params.qc.enable_vd_constraint
        validTracks = tracks;
        return;
    end
    
    max_vd = params.qc.max_vd_ratio;
    keep_mask = true(1, length(tracks));

    for i = 1:length(tracks)
        if tracks(i).length < 2
            continue;
        end
        
        path = tracks(i).path;
        
        % Sum of step-by-step distances
        step_lengths = vecnorm(diff(path, 1, 1), 2, 2);
        pathLength = sum(step_lengths);
        
        % Net displacement from start to end
        netDisplacement = norm(path(end,:) - path(1,:));
        
        if netDisplacement < 1e-6 % Avoid division by zero
            vd_ratio = inf;
        else
            vd_ratio = pathLength / netDisplacement;
        end
        
        if vd_ratio > max_vd
            keep_mask(i) = false;
        end
    end

    validTracks = tracks(keep_mask);
end