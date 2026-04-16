function validTracks = applyAccelerationConstraint(tracks, params)
% =========================================================================
% FUNCTION: applyAccelerationConstraint
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   A Quality Control (QC) filter that removes physically impossible tracks 
%   based on their acceleration. 
%   - Advantages: Uses an *adaptive* threshold. Instead of a hard-coded 
%     acceleration limit, the maximum allowed acceleration scales with the 
%     track's mean velocity. This prevents fast-moving bubbles (which naturally 
%     exhibit higher absolute acceleration due to vessel tortuosity) from 
%     being unfairly penalized compared to slow-moving bubbles.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Validation: Skips processing if tracks are empty or if the constraint 
%      is disabled in `params.qc.enable_acceleration_constraint`.
%   2. Iteration: Evaluates each track individually (requires minimum 3 points 
%      to calculate at least 2 velocities and 1 acceleration).
%   3. Kinematics: Derives step-by-step acceleration (diff(velocity) / dt).
%   4. Adaptive Thresholding: Calculates the limit as: 
%      threshold = C_factor * Mean_Velocity / dt.
%   5. Rejection: If *any* single acceleration step within the track exceeds 
%      the threshold, the entire track is marked for deletion.
%
% SYNTAX OPTIONS:
%   validTracks = applyAccelerationConstraint(tracks, params)
%
% EXAMPLES:
%   % Clean tracking output by removing high-acceleration outliers:
%   filtered_tracks = applyAccelerationConstraint(rawTracks, params);
%
% INPUTS:
%   tracks - (Type: Struct Array) Contains microbubble trajectories with 
%            pre-calculated '.velocities_mm_s' and '.average_velocity_mm_s'.
%   params - (Type: Struct) QC configuration. Requires:
%            .dt (double) Time delta between frames.
%            .qc.acceleration_C_factor (double) Multiplier for the adaptive limit.
%
% OUTPUTS:
%   validTracks - (Type: Struct Array) Subset of the input tracks that passed 
%                 the acceleration test.
%
% AUTHOR: Grigori Shapiro
% =========================================================================

    % Filters tracks based on an adaptive maximum acceleration.
    if isempty(tracks) || ~params.qc.enable_acceleration_constraint
        validTracks = tracks;
        return;
    end
    
    C = params.qc.acceleration_C_factor;
    dt = params.dt;
    keep_mask = true(1, length(tracks));
    
    for i = 1:length(tracks)
        if tracks(i).length < 3
            continue;
        end
        
        velocities = tracks(i).velocities_mm_s;
        mean_velocity = tracks(i).average_velocity_mm_s(1);
        
        % Adaptive threshold based on the track's mean speed
        accel_threshold = C * mean_velocity / dt;
        
        accelerations = diff(velocities) / dt;
        
        if any(abs(accelerations) > accel_threshold)
            keep_mask(i) = false;
        end
    end
    
    validTracks = tracks(keep_mask);
end