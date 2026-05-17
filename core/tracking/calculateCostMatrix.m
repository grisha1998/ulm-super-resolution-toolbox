function costMatrix = calculateCostMatrix(activeTracks, currentLocs, params)
% =========================================================================
% FUNCTION: calculateCostMatrix
% AUTHOR:   Grigori Shapiro
% =========================================================================
%
% PURPOSE:
%   Centralized engine for computing the assignment cost matrix between
%   active tracks and new localizations. Supports Kalman-based and
%   Hungarian tracking methods with optional advanced penalties.
%   
%   - Uses the Kalman-smoothed correctedPath for direction estimation
%     instead of raw noisy localizations, yielding more stable direction
%     vectors and more accurate angle gating.
%   - Regularizes the Mahalanobis distance (S + epsilon*I) to prevent
%     numerical failures when measurement noise is near-zero.
%   - Gracefully handles missing 'Intensity' columns in either the track
%     history or current localizations (brightness penalty is skipped).
%   - Uses the incrementally cached meanBrightness for performance.
%
% COST CALCULATION (3 Stages):
%
%   STAGE 1: Base Cost (Method Dependent)
%     Kalman / Kalman_v2 / Kalman_Advanced:
%       Mahalanobis distance — evaluates the residual (predicted vs actual)
%       scaled by the filter's innovation covariance. Dynamically forgives
%       larger jumps when the filter's predictive uncertainty is high.
%     Hungarian:
%       Squared Euclidean distance from the track's last smoothed position.
%
%   STAGE 2: Gating (Hard Constraints → Inf)
%     Distance Gate: links exceeding max_linking_distance are rejected.
%     Angle Gate: direction vector is estimated via linear regression
%       (polyfit) over the last N points of the Kalman-smoothed path.
%       Links exceeding gating_max_angle_change_deg are rejected.
%
%   STAGE 3: Advanced Penalties (Soft Multiplicative Factors >= 1.0)
%     Active only when use_advanced_cost_matrix is enabled.
%     Direction Penalty:
%       Ramped above max_angle_change_deg to discourage zigzagging:
%       P = 1 + W_dir * W_slope * (angle - threshold)
%     Brightness Penalty:
%       Penalizes intensity mismatch to prevent identity swaps:
%       P = 1 + W_int * |I_current - I_avg| / (I_avg + eps)
%
%   Final: Cost(i,j) = Base_Cost * Direction_Penalty * Brightness_Penalty
%
% INPUTS:
%   activeTracks - (Cell Array) Track structs with .kalmanFilter,
%                  .correctedPath, .localizations, .meanBrightness.
%   currentLocs  - (Table) Columns: 'X','Y' (required), 'Intensity' (optional).
%   params       - (Struct) Master parameters. Required fields:
%                  params.track.method, .max_linking_distance,
%                  .use_advanced_cost_matrix, .kalman.motion_model,
%                  .kalman.gating_max_angle_change_deg,
%                  .kalman.max_angle_change_deg, .kalman.angle_penalty_slope,
%                  .kalman.direction_penalty_weight,
%                  .kalman.brightness_penalty_weight.
%
% OUTPUTS:
%   costMatrix - (M x N double) Assignment costs; Inf for gated links.
% =========================================================================

    % --- Initialization ---
    numActiveTracks = length(activeTracks);
    numCurrentLocs  = height(currentLocs);
    costMatrix      = zeros(numActiveTracks, numCurrentLocs);

    if isempty(activeTracks) || isempty(currentLocs)
        return;
    end

    currentPositions = [currentLocs.X, currentLocs.Y];
    trackParams      = params.track;

    % Check whether the Intensity column is present in current detections
    has_intensity_current = ismember('Intensity', currentLocs.Properties.VariableNames);

    % --- Main Logic Switch ---
    switch trackParams.method
        case {'Kalman', 'Kalman_v2', 'Kalman_Advanced'}

            % Extract predicted positions from Kalman filter states
            all_states   = cellfun(@(track) track.kalmanFilter.State, activeTracks, 'UniformOutput', false);
            state_matrix = [all_states{:}];
            if strcmpi(trackParams.kalman.motion_model, 'ConstantAcceleration')
                predictedPositions = state_matrix([1, 4], :)';
            else % ConstantVelocity
                predictedPositions = state_matrix([1, 3], :)';
            end

            if trackParams.use_advanced_cost_matrix
                % --- Advanced Kalman Cost ---
                % Mahalanobis base + direction penalty + brightness penalty
                for i = 1:numActiveTracks
                    track = activeTracks{i};
                    kf    = track.kalmanFilter;

                    [valid_mask, angles_deg] = performGating(track, currentPositions, trackParams);
                    combined_cost = Inf(1, numCurrentLocs);
                    if ~any(valid_mask)
                        costMatrix(i, :) = combined_cost;
                        continue;
                    end

                    % Factor 1: Regularized Mahalanobis distance
                    mahalanobis_dist_sq = zeros(1, numCurrentLocs);
                    for k = find(valid_mask)
                        try
                            [res, S] = residual(kf, currentPositions(k, :)');
                            S_reg = S + eye(size(S)) * 1e-6;
                            mahalanobis_dist_sq(k) = res' * (S_reg \ res);
                        catch
                            % Fallback to squared Euclidean if residual fails
                            delta = currentPositions(k, :) - predictedPositions(i, :);
                            mahalanobis_dist_sq(k) = sum(delta.^2);
                        end
                    end

                    % Factor 2: Directional penalty
                    directional_penalty = calculateDirectionalPenalty(angles_deg, valid_mask, trackParams);

                    % Factor 3: Brightness penalty
                    brightness_penalty = calculateBrightnessPenalty(track, currentLocs, valid_mask, trackParams, has_intensity_current);

                    % Combine all factors
                    valid_costs = mahalanobis_dist_sq(valid_mask) .* ...
                                  directional_penalty(valid_mask) .* ...
                                  brightness_penalty(valid_mask);
                    combined_cost(valid_mask) = valid_costs;
                    costMatrix(i, :) = combined_cost;
                end
            else
                % --- Simple Kalman Cost ---
                % Squared Euclidean from predicted position, hard distance gate
                costMatrix = pdist2(predictedPositions, currentPositions, 'squaredeuclidean');
                costMatrix(costMatrix > trackParams.max_linking_distance^2) = Inf;
            end

        case 'Hungarian'
            % Use last smoothed position from correctedPath if available
            lastPositions = cell2mat(cellfun(@(c) getLastPosition(c), activeTracks, 'UniformOutput', false)');

            if trackParams.use_advanced_cost_matrix
                % --- Advanced Hungarian Cost ---
                euc_dist_sq_base = pdist2(lastPositions, currentPositions, 'squaredeuclidean');

                for i = 1:numActiveTracks
                    track = activeTracks{i};

                    [valid_mask, angles_deg] = performGating(track, currentPositions, trackParams);
                    combined_cost = Inf(1, numCurrentLocs);
                    if ~any(valid_mask)
                        costMatrix(i, :) = combined_cost;
                        continue;
                    end

                    base_cost           = euc_dist_sq_base(i, :);
                    directional_penalty = calculateDirectionalPenalty(angles_deg, valid_mask, trackParams);
                    brightness_penalty  = calculateBrightnessPenalty(track, currentLocs, valid_mask, trackParams, has_intensity_current);

                    valid_costs = base_cost(valid_mask) .* ...
                                  directional_penalty(valid_mask) .* ...
                                  brightness_penalty(valid_mask);
                    combined_cost(valid_mask) = valid_costs;
                    costMatrix(i, :) = combined_cost;
                end
            else
                % --- Simple Hungarian Cost ---
                costMatrix = pdist2(lastPositions, currentPositions, 'squaredeuclidean');
                costMatrix(costMatrix > trackParams.max_linking_distance^2) = Inf;
            end
    end
end

%% ========================================================================
%  Local Helper Functions
% =========================================================================

function pos = getLastPosition(track)
% GETLASTPOSITION  Returns the track's most recent position, preferring
%   the Kalman-smoothed correctedPath over raw localizations.

    if isfield(track, 'correctedPath') && ~isempty(track.correctedPath)
        pos = track.correctedPath(end, :);
    else
        pos = [track.localizations.X(end), track.localizations.Y(end)];
    end
end

function [valid_mask, angles_deg] = performGating(track, currentPositions, trackParams)
% PERFORMGATING  Hard gating based on distance and angle constraints.
%   Uses the Kalman-smoothed correctedPath for direction estimation,
%   yielding more stable direction vectors than raw localizations.
%   Direction is estimated via linear regression (polyfit) over the last
%   N points of the smoothed path to suppress sub-pixel jitter.

    numCurrentLocs = size(currentPositions, 1);

    % Use smoothed path for position and direction if available
    if isfield(track, 'correctedPath') && ~isempty(track.correctedPath)
        last_known_position = track.correctedPath(end, :);
        direction_history   = track.correctedPath;
    else
        last_known_position = [track.localizations.X(end), track.localizations.Y(end)];
        direction_history   = [track.localizations.X, track.localizations.Y];
    end

    % Distance gate
    euc_dist_sq   = sum((currentPositions - last_known_position).^2, 2)';
    distance_mask = euc_dist_sq <= trackParams.max_linking_distance^2;

    % Angle gate: estimate direction via polyfit on last N smoothed points
    angle_mask = true(1, numCurrentLocs);
    angles_deg = zeros(1, numCurrentLocs);

    num_pts = size(direction_history, 1);
    last_motion_vector = [];

    if num_pts >= 5
        hist_points = direction_history(end-4:end, :);
        k = 5;
    elseif num_pts >= 2
        hist_points = direction_history;
        k = num_pts;
    else
        k = 0;
    end

    if k >= 2
        t   = (1:k)';
        p_x = polyfit(t, hist_points(:,1), 1);
        p_y = polyfit(t, hist_points(:,2), 1);
        last_motion_vector = [p_x(1), p_y(1)];
    end

    if ~isempty(last_motion_vector) && norm(last_motion_vector) > eps
        candidate_vectors = currentPositions - last_known_position;
        dot_products  = candidate_vectors * last_motion_vector';
        norms_product = vecnorm(candidate_vectors, 2, 2) * norm(last_motion_vector);
        cos_theta     = dot_products ./ (norms_product + eps);
        cos_theta     = max(-1, min(1, cos_theta));
        angles_deg    = acosd(cos_theta');

        angle_mask = angles_deg <= trackParams.kalman.gating_max_angle_change_deg;
    end

    valid_mask = distance_mask & angle_mask;
end

function penalty = calculateDirectionalPenalty(angles_deg, valid_mask, trackParams)
% CALCULATEDIRECTIONALPENALTY  Ramped soft penalty for turning angles that
%   exceed the soft threshold (max_angle_change_deg). Below the threshold,
%   no penalty is applied. Above it, penalty ramps linearly:
%   P = 1 + W_dir * W_slope * (angle - threshold)

    numCurrentLocs = length(angles_deg);
    penalty = ones(1, numCurrentLocs);

    if ~any(valid_mask)
        return;
    end

    angle_threshold_deg = trackParams.kalman.max_angle_change_deg;
    angle_penalty_slope = trackParams.kalman.angle_penalty_slope;
    direction_weight    = trackParams.kalman.direction_penalty_weight;

    valid_angles = angles_deg(valid_mask);
    exceeds_mask = valid_angles > angle_threshold_deg;

    if any(exceeds_mask)
        excess = valid_angles(exceeds_mask) - angle_threshold_deg;
        ramped = 1 + direction_weight * angle_penalty_slope * excess;

        temp_penalty = ones(1, sum(valid_mask));
        temp_penalty(exceeds_mask) = ramped;
        penalty(valid_mask) = temp_penalty;
    end
end

function penalty = calculateBrightnessPenalty(track, currentLocs, valid_mask, trackParams, has_intensity_current)
% CALCULATEBRIGHTNESSPENALTY  Soft penalty for acoustic intensity mismatch.
%   Prevents identity swaps between spatially close bubbles of different
%   brightness. Uses the incrementally cached meanBrightness if available.
%   Gracefully returns no penalty (all ones) when the Intensity column is
%   absent from either the track history or current detections.
%
%   P = 1 + W_int * |I_current - I_avg| / (I_avg + eps)

    numCurrentLocs = height(currentLocs);
    penalty = ones(1, numCurrentLocs);

    % Skip if Intensity data is unavailable
    has_intensity_track = ismember('Intensity', track.localizations.Properties.VariableNames);
    if ~has_intensity_current || ~has_intensity_track
        return;
    end

    % Use cached mean brightness for performance (maintained incrementally
    % during tracking), with fallback to full-history computation
    if isfield(track, 'meanBrightness') && ~isempty(track.meanBrightness) && ...
            isfinite(track.meanBrightness) && track.meanBrightness > 0
        avg_brightness = track.meanBrightness;
    else
        avg_brightness = mean(track.localizations.Intensity, 'omitnan');
    end

    if avg_brightness <= 0 || ~isfinite(avg_brightness)
        return;
    end

    currentBrightnesses = currentLocs.Intensity';
    brightness_diff     = abs(currentBrightnesses(valid_mask) - avg_brightness);
    normalized_diff     = brightness_diff / (avg_brightness + eps);

    brightness_weight   = trackParams.kalman.brightness_penalty_weight;
    penalty(valid_mask)  = 1 + brightness_weight * normalized_diff;
end