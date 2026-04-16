function costMatrix = calculateCostMatrix(activeTracks, currentLocs, params)
% =========================================================================
% FUNCTION: calculateCostMatrix
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Acts as the centralized, highly optimized engine for calculating the 
%   assignment cost matrix between currently active tracks and new localizations.
%   - Advantages: By consolidating cost logic, it ensures uniformity across 
%     all tracking algorithms (Kalman, Hungarian, NN). It introduces an 
%     "Advanced Cost Matrix" that goes beyond simple spatial distance, 
%     incorporating motion predictability, path smoothness, and intensity 
%     conservation to resolve tracking ambiguities in dense microbubble fields.
%
% DETAILED METHODOLOGY & COST CALCULATION:
%   The final cost matrix (M x N) is computed in three main stages:
%
%   STAGE 1: Base Cost Calculation (Method Dependent)
%   -------------------------------------------------
%   The core distance metric adapts based on the active tracking algorithm:
%   * Kalman / Kalman_Advanced:
%     Uses the Mahalanobis distance. Instead of just measuring physical 
%     Euclidean distance, it evaluates the Kalman filter's `residual` 
%     (the difference between the predicted position and the actual candidate 
%     position) scaled by the filter's uncertainty covariance matrix 
%     (`MeasurementNoise`). This means the algorithm dynamically forgives 
%     larger physical jumps if the filter's predictive uncertainty is high.
%   * Hungarian / Nearest Neighbor:
%     Uses Squared Euclidean Distance between the track's last known 
%     position and the new candidate's position.
%
%   STAGE 2: Gating (Hard Constraints)
%   -------------------------------------------------
%   Invalid or physically impossible links are immediately set to Infinity 
%   (Inf) to save computation time and prevent catastrophic assignments:
%   * Distance Gate: Any link exceeding `params.track.max_linking_distance`.
%   * Angle Gate: Sub-pixel localization jitter can cause wildly inaccurate 
%     angle calculations between just two points. To solve this, `performGating` 
%     uses linear regression (`polyfit`) on up to the last 5 points of a track 
%     to find a stable "Direction Vector". Links exceeding `max_angle_change_deg` 
%     relative to this stable vector are gated (set to Inf).
%
%   STAGE 3: Advanced Penalties (Soft Constraints as Multipliers)
%   -------------------------------------------------
%   If `use_advanced_cost_matrix` is enabled, the Base Cost is multiplied 
%   by specific penalty factors (Values >= 1.0). This artificially inflates 
%   the cost of undesirable (but theoretically possible) links:
%
%   * Direction Penalty: 
%     - Why? Forces the tracker to prefer smooth, continuous vascular flow 
%       and strongly discourages erratic zigzagging.
%     - How? If the candidate's angle exceeds a defined soft threshold, a 
%       ramped scalar penalty is applied based on the slope:
%       Penalty = 1 + (Weight * Slope * (Angle - Threshold))
%       Variables used: `trackParams.kalman.direction_penalty_weight`.
%
%   * Brightness (Intensity) Penalty:
%     - Why? Microbubbles generally maintain consistent acoustic intensity 
%       over short timeframes. This penalty prevents "identity swapping" 
%       (e.g., a faint track suddenly jumping to a bright bubble nearby just 
%       because it is spatially close).
%     - How? Evaluates the fractional difference between the candidate's 
%       intensity and the track's historical mean intensity.
%       Penalty = 1 + (Weight * abs(I_current - I_avg) / I_avg)
%       Variables used: `trackParams.kalman.brightness_penalty_weight`.
%
%   FINAL EQUATION:
%   Cost(i,j) = Base_Cost * Direction_Penalty * Brightness_Penalty
%
% SYNTAX OPTIONS:
%   C = calculateCostMatrix(activeTracks, currentLocs, params)
%
% EXAMPLES:
%   % Generate the cost matrix during a tracking loop iteration:
%   cost_matrix = calculateCostMatrix(activeTracks, currentDetections, pipelineParams);
%
% INPUTS:
%   activeTracks - (Type: Cell Array) Array of active track structures. 
%                  Must contain `.kalmanFilter` (if using Kalman modes) or 
%                  `.localizations` / `.lastPos` (for Hungarian/NN).
%   currentLocs  - (Type: Table) Localizations detected in the current frame.
%   params       - (Type: Struct) Master parameters containing `.track.method` 
%                  and `.track.kalman...` weights for the advanced penalties.
%
% OUTPUTS:
%   costMatrix   - (Type: M x N Double Matrix) Where M = active tracks, 
%                  N = new localizations. Values represent the calculated 
%                  assignment cost; Inf represents gated/invalid links.
% =========================================================================
    
    % --- 1. Initialization ---
    numActiveTracks = length(activeTracks);
    numCurrentLocs = height(currentLocs);
    costMatrix = zeros(numActiveTracks, numCurrentLocs);

    if isempty(activeTracks) || isempty(currentLocs)
        return;
    end
    
    currentPositions = [currentLocs.X, currentLocs.Y];
    trackParams = params.track;

    % --- 2. Main Logic Switch ---
    switch trackParams.method
        case {'Kalman', 'Kalman_Advanced'}
            
            % --- Get Predicted Positions (common for both simple & advanced Kalman) ---
            all_states = cellfun(@(track) track.kalmanFilter.State, activeTracks, 'UniformOutput', false);
            state_matrix = [all_states{:}];
            if strcmpi(trackParams.kalman.motion_model, 'ConstantAcceleration')
                predictedPositions = state_matrix([1, 4], :)';
            else % ConstantVelocity
                predictedPositions = state_matrix([1, 3], :)';
            end
            
            if trackParams.use_advanced_cost_matrix
                % --- Advanced Kalman Cost (Mahalanobis + Direction + Brightness) ---
                for i = 1:numActiveTracks
                    track = activeTracks{i};
                    kf = track.kalmanFilter;
                    
                    [valid_links_mask, angles_deg] = performGating(track, currentPositions, trackParams);
                    combined_cost = Inf(1, numCurrentLocs);
                    if ~any(valid_links_mask), costMatrix(i, :) = combined_cost; continue; end

                    % Factor 1: Mahalanobis Distance
                    mahalanobis_dist_sq = zeros(1, numCurrentLocs);
                    for k = find(valid_links_mask)
                        [res, S] = residual(kf, currentPositions(k, :)');
                        mahalanobis_dist_sq(k) = res' * A \ res;
                    end

                    % Factor 2: Directional Penalty
                    directional_penalty = calculateDirectionalPenalty(angles_deg, valid_links_mask, trackParams);
                    
                    % Factor 3: Brightness Penalty
                    brightness_penalty = calculateBrightnessPenalty(track, currentLocs, valid_links_mask, trackParams);

                    % Final Cost Combination
                    valid_costs = mahalanobis_dist_sq(valid_links_mask) .* ...
                                  directional_penalty(valid_links_mask) .* ...
                                  brightness_penalty(valid_links_mask);
                    combined_cost(valid_links_mask) = valid_costs;
                    costMatrix(i, :) = combined_cost;
                end
            else
                % --- Simple Kalman Cost (Squared Euclidean from Prediction) ---
                costMatrix = pdist2(predictedPositions, currentPositions, 'squaredeuclidean');
                costMatrix(costMatrix > trackParams.max_linking_distance^2) = Inf;
            end

        
        case 'Hungarian'
            % --- Hungarian Cost (Euclidean + Direction + Brightness) ---
            lastPositions = cell2mat(cellfun(@(c) [c.localizations.X(end), c.localizations.Y(end)], activeTracks, 'UniformOutput', false)');
            
            if trackParams.use_advanced_cost_matrix
                % --- Advanced Hungarian Cost ---
                euc_dist_sq_base = pdist2(lastPositions, currentPositions, 'squaredeuclidean');
                
                for i = 1:numActiveTracks
                    track = activeTracks{i};
                    
                    % --- Gating (Distance & Angle) ---
                    [valid_links_mask, angles_deg] = performGating(track, currentPositions, trackParams);
                    combined_cost = Inf(1, numCurrentLocs);
                    if ~any(valid_links_mask), costMatrix(i, :) = combined_cost; continue; end

                    % --- Factor 1: Euclidean Distance (Base Cost) ---
                    base_cost = euc_dist_sq_base(i, :);

                    % --- Factor 2: Directional Penalty ---
                    directional_penalty = calculateDirectionalPenalty(angles_deg, valid_links_mask, trackParams);
                    
                    % --- Factor 3: Brightness Penalty ---
                    brightness_penalty = calculateBrightnessPenalty(track, currentLocs, valid_links_mask, trackParams);

                    % --- Final Cost Combination ---
                    valid_costs = base_cost(valid_links_mask) .* ...
                                  directional_penalty(valid_links_mask) .* ...
                                  brightness_penalty(valid_links_mask);
                    combined_cost(valid_links_mask) = valid_costs;
                    costMatrix(i, :) = combined_cost;
                end
            else
                % --- Simple Hungarian Cost (Squared Euclidean) ---
                costMatrix = pdist2(lastPositions, currentPositions, 'squaredeuclidean');
                costMatrix(costMatrix > trackParams.max_linking_distance^2) = Inf;
            end
    end
end

%% ========================================================================
%  Local Helper Functions
% =========================================================================

function [valid_links_mask, angles_deg] = performGating(track, currentPositions, trackParams)
    % Performs hard gating based on distance and angle.
    % Uses robust polyfit-based direction vector.
    
    numCurrentLocs = size(currentPositions, 1);
    last_known_position = [track.localizations.X(end), track.localizations.Y(end)];
    
    % 1a: Distance Gate
    euc_dist_sq_from_last = sum((currentPositions - last_known_position).^2, 2)';
    distance_mask = euc_dist_sq_from_last <= trackParams.max_linking_distance^2;
    
    % 1b: Angle Gate
    angle_mask = true(1, numCurrentLocs);
    angles_deg = zeros(1, numCurrentLocs);
    
    num_locs = height(track.localizations);
    last_motion_vector = []; % Initialize
    k = 0;
    
    if num_locs >= 5
        % --- Smart method: Use last 5 points ---
        k = 5;
        hist_points = track.localizations((end-k+1):end, :);
        
    elseif num_locs >= 2 % User said "less than 5", so this catches 2, 3, 4
        % --- Smart method: Use ALL available points (k = num_locs) ---
        k = num_locs;
        hist_points = track.localizations; % Use all of them
    end
    
    if k >= 2
        % Perform linear regression
        x_hist = hist_points.X;
        y_hist = hist_points.Y;
        t = (1:k)'; % Time vector
        
        % Fit line: y = p(1)*t + p(2)
        % The slope, p(1), is the velocity (direction vector component)
        p_x = polyfit(t, x_hist, 1); 
        p_y = polyfit(t, y_hist, 1);
        
        last_motion_vector = [p_x(1), p_y(1)]; % Vector is [vx, vy]
    end
    
    % --- Calculate angles (same as before, but using the new robust vector) ---
    if ~isempty(last_motion_vector) && norm(last_motion_vector) > eps
        candidate_vectors = currentPositions - last_known_position;
        dot_products = (candidate_vectors * last_motion_vector');
        norms_product = vecnorm(candidate_vectors, 2, 2) * norm(last_motion_vector);

        cos_theta = dot_products ./ (norms_product + eps);
        cos_theta = max(-1, min(1, cos_theta)); % Protect against numerical errors
        angles_deg = acosd(cos_theta');
        
        angle_mask = angles_deg <= trackParams.kalman.gating_max_angle_change_deg;
    end
    
    valid_links_mask = distance_mask & angle_mask;
end

function penalty = calculateDirectionalPenalty(angles_deg, valid_mask, trackParams)
    % Calculates a soft, ramped penalty for turning angles that exceed a threshold.
    numCurrentLocs = length(angles_deg);
    penalty = ones(1, numCurrentLocs);
    
    % Check if any links are valid to begin with
    if ~any(valid_mask)
        return;
    end
    
    % Get parameters from the params struct for configurability
    angle_threshold_deg = trackParams.kalman.max_angle_change_deg;      % e.g., 120
    angle_penalty_slope = trackParams.kalman.angle_penalty_slope;       % e.g., 0.5

    % Isolate the angles for valid potential links
    valid_angles = angles_deg(valid_mask);
    
    % Find which of these valid angles exceed the penalty threshold
    exceeds_threshold_mask = valid_angles > angle_threshold_deg;
    
    if any(exceeds_threshold_mask)
        % Isolate only the angles that require a penalty
        angles_that_exceed = valid_angles(exceeds_threshold_mask);
        
        % Calculate the ramped penalty for these specific angles
        ramped_penalty = 1 + trackParams.kalman.direction_penalty_weight * angle_penalty_slope * (angles_that_exceed - angle_threshold_deg);
        
        % Create a temporary penalty vector for all valid links (defaulting to 1)
        temp_penalty_values = ones(1, sum(valid_mask));
        % Place the calculated penalties into the correct positions
        temp_penalty_values(exceeds_threshold_mask) = ramped_penalty;
        
        % Assign the calculated penalties back to the main penalty vector
        penalty(valid_mask) = temp_penalty_values;
    end
end

function penalty = calculateBrightnessPenalty(track, currentLocs, valid_mask, trackParams)
    % Calculates a soft penalty for changes in brightness.
    numCurrentLocs = height(currentLocs);
    penalty = ones(1, numCurrentLocs);
    
    avg_track_brightness = mean(track.localizations.Intensity);
    currentBrightnesses = currentLocs.Intensity';
    
    brightness_diff = abs(currentBrightnesses(valid_mask) - avg_track_brightness);
    normalized_brightness_diff = brightness_diff / (avg_track_brightness + eps);
    
    penalty_values = 1 + trackParams.kalman.brightness_penalty_weight * normalized_brightness_diff;
    penalty(valid_mask) = penalty_values;
end