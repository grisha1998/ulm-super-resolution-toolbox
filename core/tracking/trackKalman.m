function tracks = trackKalman(localizations, params, indent_prefix)
% =========================================================================
% FUNCTION: trackKalman
% =========================================================================
%
% PURPOSE:
%   Implements predictive microbubble tracking utilizing a 2D Kalman Filter. 
%   By maintaining a kinematic state (position and velocity/acceleration) for 
%   each microbubble, this tracker can accurately predict where a bubble will 
%   appear in the next frame. This drastically reduces identity-swapping in dense 
%   regions and smooths out sub-pixel localization jitter.
%
% DETAILED ALGORITHM WORKFLOW:
%   1. Mask Initialization: Checks for `params.proc.vesselMask`. If present, 
%      it enables physical boundary constraints.
%   2. Predict Phase: Every active track's `kalmanFilter` is stepped forward 
%      (`predict()`) to estimate the [X, Y] position in the current frame.
%   3. Cost & Assignment: Builds a cost matrix comparing *predicted* track 
%      locations to *measured* localizations. Solves via Hungarian or NN.
%   4. Update & Constrain (The Core Logic):
%      - For assigned links, the Kalman filter is updated (`correct()`) with the 
%        new measurement.
%      - Constraint Check: The corrected [X, Y] state is checked against the 
%        `vesselMask`. If the Kalman filter "drifts" outside the defined vessels, 
%        the state is overridden and clamped strictly back to the raw measurement 
%        coordinate to prevent track divergence.
%   5. Track Path Storage: The newly corrected (smoothed) state is saved to 
%      `.correctedPath`. This path is used for all downstream velocity mapping.
%   6. Two-Step Track Initiation: To prevent noise spikes from starting fake tracks, 
%      a new track is ONLY created if an unassigned localization in frame `T` can 
%      be successfully linked to an unassigned localization in frame `T-1`.
%   7. Track Management & Finalization: Gap closing logic is applied. Dead tracks 
%      are evaluated, velocities are calculated in mm/s using physical parameters, 
%      and mask-safety checks are performed globally.
%
% INPUT PARAMETERS & VARIABLES:
%   localizations - (Table) 'Frame', 'X', 'Y' columns.
%   params        - (Struct) Tracking configuration:
%                   .track.kalman.motion_model (String): 'ConstantVelocity' (4-state) 
%                         or 'ConstantAcceleration' (6-state).
%                   .track.kalman.process_noise (double): Covariance Q; defines how 
%                         much the velocity/acceleration is allowed to change.
%                   .track.kalman.assignment_method (String): 'hungarian' or 'nn'.
%                   .loc.fwhm (double array): Used to dynamically set the measurement 
%                         noise (R) of the Kalman filter based on bubble size.
%                   .proc.vesselMask.vesselMask (Logical Matrix): Spatial constraint.
%
% INTERNAL TRACK STATE VARIABLES:
%   .kalmanFilter              - (trackingKF object) Stores state, covariances, etc.
%   .correctedPath             - [Nx2 double] The smoothed trajectory output from KF.
%   .age                       - (Integer) Total lifespan frames.
%   .consecutiveInvisibleCount - (Integer) Current gap counter.
%   unassignedLocs_prevFrame   - (Table) Cross-loop variable enabling two-step initiation.
%
% FINAL OUTPUT STRUCTURE (tracks array):
%   Matches trackHungarian, but critically, `.path` is populated using the 
%   smoothed `.correctedPath` to ensure high-fidelity super-resolution rendering.
%
% DEPENDENCIES:
%   - Sensor Fusion and Tracking Toolbox (`trackingKF`, `predict`, `correct`)
%
% AUTHOR: Grigori Shapiro
% =========================================================================
    
    % --- 1. Initialization and Pre-checks ---
    if nargin < 3, indent_prefix = ''; end 
    
    % Setup Mask for Boundary Constraints
    use_vessel_mask = false;
    if isfield(params, 'proc') && isfield(params.proc, 'vesselMask') && ...
            isfield(params.proc.vesselMask, 'enable') && params.proc.vesselMask.enable && ...
            isfield(params.proc.vesselMask, 'vesselMask') && ~isempty(params.proc.vesselMask.vesselMask)
        vesselMask = params.proc.vesselMask.vesselMask;
        use_vessel_mask = true;
        [maskH, maskW] = size(vesselMask);
    end

    if isempty(localizations)
        tracks = [];
        return;
    end
    
    if ~license('test', 'Sensor_Fusion_and_Tracking')
        error('The "kalman" tracking method requires the Sensor Fusion and Tracking Toolbox.');
    end

    locsByFrame = splitapply(@(varargin) {table(varargin{:}, 'VariableNames', localizations.Properties.VariableNames)},...
        localizations, findgroups(localizations.Frame));
    numFrames = length(locsByFrame);

    activeTracks = {};
    terminatedTracks = {};
    nextTrackID = 1;
    
    % Get the assignment method from params ---
    if isfield(params.track.kalman, 'assignment_method')
        assignment_method = lower(params.track.kalman.assignment_method);
    else
        assignment_method = 'hungarian'; % Default to Hungarian
    end

    fprintf('%sStarting standard Kalman tracking...\n', indent_prefix);
    
    % === Initialize a table to store unassigned localizations from the previous frame ===
    unassignedLocs_prevFrame = table();

    % --- 2. Main Tracking Loop ---
    for t = 1:numFrames
        currentLocs = locsByFrame{t};
        
        % --- PREDICTION STEP ---
        numActiveTracks = length(activeTracks);
        if numActiveTracks > 0
            % Optimized state extraction and prediction
            cellfun(@(track) predict(track.kalmanFilter), activeTracks);
        end
        
        % --- ASSIGNMENT STEP ---
        unassignedDetections = true(height(currentLocs), 1);
        if ~isempty(activeTracks) && ~isempty(currentLocs)
            costMatrix = calculateCostMatrix(activeTracks, currentLocs, params);

            switch assignment_method
                case 'hungarian'
                    if exist('munkres', 'file')
                        [assignments, ~] = munkres(costMatrix);
                    else
                        error('Hungarian algorithm implementation (e.g., ''munkres'') not found on MATLAB path.');
                    end
                    
                case 'nn'
                    % Use the greedy nearest neighbor assignment
                    max_link_dist_sq = params.track.max_linking_distance^2;
                    assignments = greedyNN_assign(costMatrix, max_link_dist_sq);
                    
                otherwise
                    error('Unknown assignment method: %s. Use ''hungarian'' or ''nn''.', params.track.kalman.assignment_method);
            end

            % --- UPDATE STEP ---
            assignedTracks = false(numActiveTracks, 1);
            currentPositions = [currentLocs.X, currentLocs.Y]; % Get positions for correction
            for i = 1:numActiveTracks
                if assignments(i) > 0
                    detectionIdx = assignments(i);
                    
                    % 1. Correct the filter (This is good)
                    correct(activeTracks{i}.kalmanFilter, currentPositions(detectionIdx, :));
                    
                    % 2. Extract the corrected state (This is the solution)
                    corrected_state = activeTracks{i}.kalmanFilter.State;
                    
                    % 3. Extract the [x, y] position and apply Mask Constraint
                    if strcmpi(params.track.kalman.motion_model, 'ConstantAcceleration')
                         state_idx_x = 1; state_idx_y = 4;
                    else % ConstantVelocity
                         state_idx_x = 1; state_idx_y = 3;
                    end
                    
                    corrected_pos = [corrected_state(state_idx_x), corrected_state(state_idx_y)];

                    % --- CONSTRAINT CHECK ---
                    if use_vessel_mask
                        % Convert continuous position to grid indices
                        cx = round(corrected_pos(1));
                        cy = round(corrected_pos(2));
                        
                        % Check boundaries
                        is_outside = (cx < 1 || cx > maskW || cy < 1 || cy > maskH);
                        
                        % If inside bounds, check the specific mask pixel
                        if ~is_outside && ~vesselMask(cy, cx)
                            is_outside = true;
                        end
                        
                        if is_outside
                            % VIOLATION DETECTED: Force position back to the raw localization
                            raw_loc = [currentLocs.X(detectionIdx), currentLocs.Y(detectionIdx)];
                            corrected_pos = raw_loc;
                            
                            % CRITICAL: Update the Kalman Filter state to stop it from "drifting"
                            % If we don't do this, the filter will keep predicting further outside next time.
                            activeTracks{i}.kalmanFilter.State(state_idx_x) = corrected_pos(1);
                            activeTracks{i}.kalmanFilter.State(state_idx_y) = corrected_pos(2);
                        end
                    end
                    
                    % 4. Store the corrected position in a *new* field
                    activeTracks{i}.correctedPath = [activeTracks{i}.correctedPath; corrected_pos];
                    
                    % 5. Store the raw localization (still good practice)
                    activeTracks{i}.localizations = [activeTracks{i}.localizations; currentLocs(detectionIdx,:)];
                    
                    activeTracks{i}.age = activeTracks{i}.age + 1;
                    activeTracks{i}.lastFrame = t;
                    activeTracks{i}.consecutiveInvisibleCount = 0;
                    
                    unassignedDetections(detectionIdx) = false;
                    assignedTracks(i) = true;
                end
            end
        else
            assignedTracks = false(numActiveTracks, 1);
        end

        % --- 4. Track Management ---
        % Handle tracks that were not assigned (gap closing or termination)
        for i = numActiveTracks:-1:1
            if ~assignedTracks(i)
                activeTracks{i}.consecutiveInvisibleCount = activeTracks{i}.consecutiveInvisibleCount + 1;
                if activeTracks{i}.consecutiveInvisibleCount > params.track.max_gap_closing_frames
                    terminatedTracks{end+1} = activeTracks{i};
                    activeTracks(i) = [];
                end
            end
        end
        
        % === Two-Step Track Initiation Logic ===
        
        % 1. Get all localizations from the current frame that were not assigned to any active track.
        unassignedLocs_currentFrame = currentLocs(unassignedDetections, :);
        
        % 2. Try to link current unassigned locs with unassigned locs from the PREVIOUS frame.
        if ~isempty(unassignedLocs_currentFrame) && ~isempty(unassignedLocs_prevFrame)
            costMatrix_init = pdist2([unassignedLocs_prevFrame.X, unassignedLocs_prevFrame.Y], ...
                                      [unassignedLocs_currentFrame.X, unassignedLocs_currentFrame.Y], 'squaredeuclidean');
            
            init_dist_gate_sq = (params.track.max_linking_distance * 1.5)^2;
            costMatrix_init(costMatrix_init > init_dist_gate_sq) = Inf;
            
            [assignments_init, ~] = munkres(costMatrix_init);
            
            used_current_locs_mask = false(height(unassignedLocs_currentFrame), 1);
            
            % 3. Create new tracks ONLY from valid two-frame links.
            for i = 1:length(assignments_init)
                if assignments_init(i) > 0
                    detectionIdx = assignments_init(i);
                    
                    loc1 = unassignedLocs_prevFrame(i, :);
                    loc2 = unassignedLocs_currentFrame(detectionIdx, :);
                    
                    % Initialize filter with the first point
                    kf = initializeKalmanFilter(loc1, params);
                    initial_state = kf.State; % Get state for point 1
                    
                    % Immediately correct it with the second point
                    correct(kf, [loc2.X, loc2.Y]);
                    corrected_state_2 = kf.State; % Get corrected state for point 2
                    
                    % Extract positions
                    if strcmpi(params.track.kalman.motion_model, 'ConstantAcceleration')
                         state_idx_x = 1; state_idx_y = 4;
                    else % ConstantVelocity
                         state_idx_x = 1; state_idx_y = 3;
                    end
                    
                    pos1 = [initial_state(state_idx_x), initial_state(state_idx_y)];
                    pos2 = [corrected_state_2(state_idx_x), corrected_state_2(state_idx_y)];

                    % Apply Mask Constraint to the Second Point (pos2)
                    if use_vessel_mask
                        cx = round(pos2(1));
                        cy = round(pos2(2));
                        
                        is_outside = (cx < 1 || cx > maskW || cy < 1 || cy > maskH);
                        
                        if ~is_outside && ~vesselMask(cy, cx)
                            is_outside = true;
                        end
                        
                        if is_outside
                            % Force pos2 back to the raw localization loc2
                            pos2 = [loc2.X, loc2.Y];
                            
                            % Update the filter state so it doesn't start with a bad state
                            kf.State(state_idx_x) = pos2(1);
                            kf.State(state_idx_y) = pos2(2);
                        end
                    end
                    
                    newTrack.id = nextTrackID; nextTrackID = nextTrackID + 1;
                    newTrack.kalmanFilter = kf;
                    newTrack.localizations = [loc1; loc2]; % Store raw data
                    newTrack.correctedPath = [pos1; pos2];
                    newTrack.age = 2;
                    newTrack.lastFrame = t;
                    newTrack.consecutiveInvisibleCount = 0;
                    activeTracks{end+1} = newTrack;
                    
                    used_current_locs_mask(detectionIdx) = true;
                end
            end
            
            % 4. Update the list of unassigned locs for the next frame.
            unassignedLocs_prevFrame = unassignedLocs_currentFrame(~used_current_locs_mask, :);
            
        else
            unassignedLocs_prevFrame = unassignedLocs_currentFrame;
        end

    end

    % --- 5. Finalization ---
    terminatedTracks = [terminatedTracks, activeTracks];
    
    terminatedTracks = validate_Tracks(terminatedTracks, params);
    tracks = formatTracks(terminatedTracks, params, indent_prefix); 
    
    fprintf('%sStandard Kalman tracking complete. Found %d valid tracks after QC.\n', indent_prefix, length(tracks));

    total_initial_locs = height(localizations);
    if ~isempty(tracks)
        total_locs_in_tracks = sum([tracks.length]);
    else
        total_locs_in_tracks = 0;
    end
    untracked_locs = total_initial_locs - total_locs_in_tracks;
    
    fprintf('%s  - Total localizations in tracks: %d\n', indent_prefix, total_locs_in_tracks);
    if total_initial_locs > 0
        untracked_percent = 100 * untracked_locs / total_initial_locs;
        fprintf('%s  - Un-tracked localizations: %d (%.1f%% of total)\n', indent_prefix, untracked_locs, untracked_percent);
    else
        fprintf('%s  - Un-tracked localizations: 0\n', indent_prefix);
    end
end

%% ========================================================================
%  Local Helper Functions
% =========================================================================

function kf = initializeKalmanFilter(loc, params)
    % This function configures and initializes a new Kalman filter for a new track.
    
    kalmanParams = params.track.kalman;
    motionModel = kalmanParams.motion_model;
    
    % Set up the filter based on the chosen motion model.
    if strcmpi(motionModel, 'ConstantAcceleration')
        % State is [x; vx; ax; y; vy; ay]
        initialState = [loc.X; 0; 0; loc.Y; 0; 0];
        kf = trackingKF('MotionModel', '2D Constant Acceleration', 'State', initialState);
        kf.ProcessNoise = kalmanParams.process_noise;
    else % Default to 'ConstantVelocity'
        % State is [x; vx; y; vy]
        initialState = [loc.X; 0; loc.Y; 0];
        kf = trackingKF('MotionModel', '2D Constant Velocity', 'State', initialState);
        kf.ProcessNoise = kalmanParams.process_noise;
    end
    
    localization_uncertainty_pixels = 0.5 * mean(params.loc.fwhm);
    measurement_noise_variance = localization_uncertainty_pixels^2;
    kf.MeasurementNoise = eye(2) * measurement_noise_variance;
end


function finalTracks = formatTracks(terminatedTracks, params, indent_prefix)
    % This helper function converts raw track objects into the final output structure,
    % calculates velocities in mm/s, and filters by length.
    
    % Extract the tracking-specific parameters.
    trackParams = params.track;
    min_track_length = trackParams.min_track_length;
    
    % Pre-allocate a cell array for temporary storage for performance.
    tempTracksCell = cell(1, length(terminatedTracks));
    track_count = 0;
    
    for i = 1:length(terminatedTracks)
        track = terminatedTracks{i};
        if track.age >= min_track_length
            track_count = track_count + 1;
            
            % Store basic track info
            outTrack = struct();
            outTrack.id = track.id;
            outTrack.localizations = track.localizations; % Store raw data
            
            % Build the path from the CORRECTED states, not the raw localizations
            if isfield(track, 'correctedPath') && ~isempty(track.correctedPath)
                outTrack.path = track.correctedPath;
                % Ensure path length matches age, in case of logic error
                if size(outTrack.path, 1) ~= track.age
                     fprintf('Warning: Track %d path/age mismatch. Using localizations as fallback.\n', track.id);
                     outTrack.path = [track.localizations.X, track.localizations.Y];
                end
            else
                 fprintf('Warning: Track %d missing correctedPath. Using localizations as fallback.\n', track.id);
                 outTrack.path = [track.localizations.X, track.localizations.Y];
            end
            
            outTrack.frames = track.localizations.Frame;
            outTrack.length = track.age;

            % --- Velocity Calculation (in mm/s) ---
            % This calculation will now use the SMOOTHED path
            if track.age > 1
                dx_pixels = diff(outTrack.path(:,1));
                dy_pixels = diff(outTrack.path(:,2));
                dt_frames = diff(outTrack.frames);

                displacement_mm = sqrt((dx_pixels * trackParams.pixel_X_size).^2 + (dy_pixels * trackParams.pixel_Z_size).^2);
                dt_sec = dt_frames * trackParams.dt;

                velocities_mm_s = zeros(size(dt_sec));
                valid_dt = dt_sec > 0;
                velocities_mm_s(valid_dt) = displacement_mm(valid_dt) ./ dt_sec(valid_dt);
                outTrack.velocities_mm_s = [velocities_mm_s; velocities_mm_s(end)];
                
                mean_vel = mean(velocities_mm_s(isfinite(velocities_mm_s)), 'omitnan');
                if isempty(mean_vel) || isnan(mean_vel), mean_vel = 0; end
                outTrack.average_velocity_mm_s = repmat(mean_vel, track.age, 1);
            else
                outTrack.velocities_mm_s = 0;
                outTrack.average_velocity_mm_s = 0;
            end
            tempTracksCell{track_count} = outTrack;
        end
    end

    % Convert the cell array to a structure array, removing empty cells.
    initial_filtered_tracks = [tempTracksCell{1:track_count}];

    finalTracks = applyQualityControl(initial_filtered_tracks, trackParams, indent_prefix);
end

function assignments = greedyNN_assign(costMatrix, max_cost)
    % Solves the assignment problem using a "best-first" greedy algorithm.
    [numActiveTracks, numDetections] = size(costMatrix);
    assignments = zeros(numActiveTracks, 1);
    if numDetections == 0, return; end
    
    availableDetections = true(1, numDetections);
    
    % Find all costs below the gate
    valid_links = costMatrix <= max_cost;
    [trackIndices_all, detIndices_all] = find(valid_links);
    
    if isempty(trackIndices_all), return; end % No valid links found
    
    % Get the linear indices of these valid links to find their costs
    linear_indices = sub2ind(size(costMatrix), trackIndices_all, detIndices_all);
    costs_all = costMatrix(linear_indices);
    
    % Sort them from best (lowest cost) to worst
    [~, sortedOrder] = sort(costs_all, 'ascend');
    
    trackIndices = trackIndices_all(sortedOrder);
    detIndices = detIndices_all(sortedOrder);
    
    % Iterate through the sorted list
    for i = 1:length(trackIndices)
        trackIdx = trackIndices(i);
        detIdx = detIndices(i);
        
        % If this track is not yet assigned AND this detection is not yet taken
        if assignments(trackIdx) == 0 && availableDetections(detIdx)
            assignments(trackIdx) = detIdx;
            availableDetections(detIdx) = false; % This detection is now taken
        end
    end
end

function validTracks = validate_Tracks(tracks, params)
    % Removes ENTIRE tracks if they contain points outside the active masks.
    
    if isempty(tracks)
        validTracks = tracks;
        return;
    end

    % --- 1. Combine Masks (Updated to use 'params' not 'obj') ---
    H = params.expParams.size(1);
    W = params.expParams.size(2);
    combinedMask = true(H, W); 
    maskActive = false;

    % A. Check Vessel Mask
    if isfield(params, 'proc') && isfield(params.proc, 'vesselMask') && ...
       isfield(params.proc.vesselMask, 'enable') && params.proc.vesselMask.enable && ...
       isfield(params.proc.vesselMask, 'vesselMask') && ~isempty(params.proc.vesselMask.vesselMask)
   
        combinedMask = combinedMask & logical(params.proc.vesselMask.vesselMask);
        maskActive = true;
    end

    % B. Check Interactive ROI Mask
    % Note: In trackKalman, metadata might not be fully available in params, 
    % but if you passed it, this works. Usually, for tracking, we rely on vesselMask.
    if isfield(params, 'metadata') && isfield(params.metadata, 'roiMask') && ~isempty(params.metadata.roiMask)
        combinedMask = combinedMask & logical(params.metadata.roiMask);
        maskActive = true;
    end
    
    if ~maskActive
        validTracks = tracks;
        return;
    end

    % --- 2. Iterate Backwards to Filter Tracks ---
    % We loop backwards because we are deleting elements from the array
    tracks_to_remove = false(1, length(tracks));
    
    for i = 1:length(tracks)
        currTrack = tracks{i};
        
        % Determine which path to check (Corrected or Raw)
        if isfield(currTrack, 'correctedPath') && ~isempty(currTrack.correctedPath)
            pts = currTrack.correctedPath;
        else
            pts = [currTrack.localizations.X, currTrack.localizations.Y];
        end
        
        x_idx = round(pts(:, 1));
        y_idx = round(pts(:, 2));
        
        % Check bounds
        in_bounds = (x_idx >= 1 & x_idx <= W & y_idx >= 1 & y_idx <= H);
        
        % If any point is out of image bounds, mark for removal
        if any(~in_bounds)
            tracks_to_remove(i) = true;
            continue;
        end
        
        % Check mask
        indices = sub2ind([H, W], y_idx, x_idx);
        if any(~combinedMask(indices))
            tracks_to_remove(i) = true;
        end
    end

    % --- 3. Remove Invalid Tracks ---
    validTracks = tracks(~tracks_to_remove);
    
    num_removed = sum(tracks_to_remove);
    if num_removed > 0
        fprintf('      -> Mask Safety Check: Removed %d tracks that drifted outside the mask.\n', num_removed);
    end
end