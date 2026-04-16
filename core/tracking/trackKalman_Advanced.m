function finalTracks = trackKalman_Advanced(localizations, params, indent_prefix)
% =========================================================================
% FUNCTION: trackKalman_Advanced
% =========================================================================
%
% PURPOSE:
%   Implements a highly robust, multi-pass Hierarchical Kalman (HK) tracking 
%   architecture. This algorithm is specifically designed to overcome the 
%   limitations of standard tracking in ultra-dense, high-concentration 
%   microbubble datasets by decomposing the tracking problem into distinct 
%   velocity-dependent bands.
%
% ACADEMIC REFERENCE:
%   The core methodology—including the subtractive hierarchical approach, 
%   dynamic parameter scaling, and forward-backward tracking—is directly 
%   adapted from the research of Dr. Iman Taghavi:
%   
%   1. Primary Paper: 
%      Taghavi, I., et al. (2022). "Ultrasound super-resolution imaging with 
%      a hierarchical Kalman tracker." Ultrasonics, 124, 106742.
%   2. Doctoral Thesis:
%      Taghavi, Iman. Specific architectural implementations are based on the 
%      thesis, including Velocity-Range Filtering (Section 3.3.1.3), Forward-
%      Backward Temporal Tracking (Section 3.3.6), and Dynamic Noise Scaling 
%      Formulas (Figure 3.17).
%
% DETAILED ALGORITHM WORKFLOW:
%   1. Pre-Processing & Subtractive Framework: The problem is divided into 
%      iterative passes over defined `velocity_levels` (e.g., slow to fast). 
%      Every localization is assigned a unique `LocalizationID` to track 
%      its usage across these passes.
%   2. Dynamic Parameter Scaling: For each velocity level (n), the Kalman 
%      filter's noise matrices adapt strictly to the expected kinematics:
%      - Process Noise (sigma_e): Scales linearly with the level's max 
%        velocity (sigma_e = alpha * v_max).
%      - Measurement Noise (sigma_nu): Decays exponentially across levels 
%        (sigma_nu = beta / 2^(n-1)). The filter inherently trusts physical 
%        measurements more at higher speeds to prevent divergence.
%   3. Track Initiation & Assignment: Employs a strict "Two-Step Initiation", 
%      requiring an unassigned detection to appear in two consecutive frames 
%      within a distance gate to spawn a track (suppresses noise). Global 
%      optimal assignment (Hungarian Algorithm / Munkres) links the tracks.
%   4. Velocity-Range Filtering (Critical Check): After tracking a level, 
%      the physical mean velocity of every newly formed track is computed. 
%      If a track's speed falls strictly outside the current level's 
%      [v_min, v_max] range, it is rejected, and its localizations are 
%      returned to the available pool.
%   5. Subtractive Updates: The `LocalizationID`s of validated tracks are 
%      permanently deleted from the global pool. Subsequent (faster) levels 
%      only attempt to track the remaining, unassigned bubbles, preventing 
%      fast bubbles from being falsely linked to slow trajectories.
%   6. Forward-Backward Dual Pass & Deduplication: If enabled, tracking 
%      executes temporally forward, then backward using `parfor`. To prevent 
%      artificial density inflation, backward tracks are only kept if they 
%      trace unique paths (i.e., less than 50% overlap with forward tracks 
%      based on consumed LocalizationIDs).
%   7. Quality Control (QC): Finally, tracks shorter than `min_track_length` 
%      or those drifting outside the user-defined `combinedMask` are purged.
%
% INPUT PARAMETERS & VARIABLES:
%   localizations - (Table) Must contain: 'Frame', 'X', 'Y', 'Intensity'.
%                   (A 'LocalizationID' is generated internally if missing).
%   params        - (Struct) Tracking configuration:
%                   .track.dt, .pixel_X_size, .pixel_Z_size (double): Physical 
%                         scaling parameters for mm/s conversions.
%                   .track.min_track_length, .max_gap_closing_frames (integer).
%                   .track.kalman.velocity_levels (Cell Array): Defines the 
%                         hierarchical bands, e.g., {[1,5], [5,15]}.
%                   .track.kalman.hk_alpha (double): Process noise scale factor.
%                   .track.kalman.hk_beta (double): Measurement noise base factor.
%                   .track.kalman.hk_forward_backward (Logical): Enables dual-pass.
%                   .track.kalman.motion_model (String): e.g., 'ConstantVelocity'.
%                   .mask.combinedMask (Logical): Global spatial constraint map.
%
% FINAL OUTPUT STRUCTURE (finalTracks array):
%   Produces a standardized output struct compatible with downstream rendering:
%   .id                    - Unique Track ID.
%   .localizations         - The complete table of raw localizations.
%   .path                  - [Nx2 double] The smoothed trajectory output from KF.
%   .correctedPath         - Maintained identically to .path for compatibility.
%   .frames                - [Nx1 double] Frame indices.
%   .length                - (Integer) Total number of points in the track.
%   .velocities_mm_s       - [Nx1 double] Instantaneous physical velocities.
%   .average_velocity_mm_s - [Nx1 double] Robust mean velocity of the track.
%
% DEPENDENCIES:
%   - Sensor Fusion and Tracking Toolbox (`trackingKF`)
%   - Parallel Computing Toolbox (`parfor`)
%   - munkres.m (Hungarian algorithm implementation)
%   - calculateCostMatrix.m
%
% AUTHOR: Grigori Shapiro
% DATE: September 2025 (Revised April 2026)
% =========================================================================

    % --- 1. Initialization and Pre-checks ---
    if nargin < 3, indent_prefix = ''; end 

    % --- Validate Toolbox ---
    if ~license('test', 'Sensor_Fusion_and_Tracking')
        error('trackKalman_Advanced requires the Sensor Fusion and Tracking Toolbox.');
    end

    % --- Add a unique LocalizationID if not present (crucial for subtractive phase) ---
    if ~ismember('LocalizationID', localizations.Properties.VariableNames)
        localizations.LocalizationID = (1:height(localizations))';
    end

    fprintf('\n=======================================================\n');
    fprintf('   STARTING ADVANCED HIERARCHICAL KALMAN TRACKING\n');
    fprintf('=======================================================\n');

    % --- PREPARE DATA FOR PARALLEL EXECUTION ---
    numPasses = 1;
    datasets = {localizations};
    directions = {'forward'};

    if isfield(params.track.kalman, 'hk_forward_backward') && params.track.kalman.hk_forward_backward
        numPasses = 2;
        fprintf('   [INFO] Dual-pass enabled. Preparing Backward data...\n');

        maxFrame = max(localizations.Frame);
        locsBackward = localizations;
        locsBackward.Frame = maxFrame - locsBackward.Frame + 1;
        locsBackward = sortrows(locsBackward, 'Frame', 'ascend');

        datasets{2} = locsBackward;
        directions{2} = 'backward';
    end

    passResults = cell(1, numPasses);

    % --- PARALLEL TRACKING EXECUTION ---
    fprintf('\n---> Executing Tracking Passes (utilizing parfor if available)...\n');

    parfor passIdx = 1:numPasses
        passResults{passIdx} = runHierarchicalPass(datasets{passIdx}, params, directions{passIdx});
    end

    % --- POST-PROCESSING & MERGING ---
    tracksForward = passResults{1};

    if numPasses == 2
        tracksBackward = passResults{2};
        maxFrame = max(localizations.Frame);

        % Re-align backward tracks to match original timeline
        for i = 1:length(tracksBackward)
            tracksBackward(i).path = flipud(tracksBackward(i).path);
            if isfield(tracksBackward(i), 'correctedPath') && ~isempty(tracksBackward(i).correctedPath)
                tracksBackward(i).correctedPath = flipud(tracksBackward(i).correctedPath);
            end
            tracksBackward(i).frames = maxFrame - tracksBackward(i).frames + 1;

            [tracksBackward(i).frames, sortIdx] = sort(tracksBackward(i).frames, 'ascend');
            tracksBackward(i).path = tracksBackward(i).path(sortIdx, :);
            if isfield(tracksBackward(i), 'correctedPath') && ~isempty(tracksBackward(i).correctedPath)
                tracksBackward(i).correctedPath = tracksBackward(i).correctedPath(sortIdx, :);
            end
        end

        % --- [FIX-5] DEDUPLICATION: Keep backward tracks only if they ---
        %     contribute localizations NOT already used by forward tracks. ---
        if ~isempty(tracksForward) && ~isempty(tracksBackward)
            % Collect all LocalizationIDs consumed by forward tracks
            fwdIDs = arrayfun(@(t) t.localizations.LocalizationID, tracksForward, 'UniformOutput', false);
            fwdIDs = vertcat(fwdIDs{:});
            fwdIDSet = unique(fwdIDs);

            keepBackward = true(length(tracksBackward), 1);
            for i = 1:length(tracksBackward)
                bwdIDs_i = tracksBackward(i).localizations.LocalizationID;
                overlapFraction = sum(ismember(bwdIDs_i, fwdIDSet)) / length(bwdIDs_i);
                % If more than 50% of this backward track's localizations
                % were already used by forward tracks, discard it.
                if overlapFraction > 0.5
                    keepBackward(i) = false;
                end
            end
            tracksBackward = tracksBackward(keepBackward);

            fprintf('\n[DEDUP] Kept %d of %d backward tracks (%.1f%% overlap filtered).\n', ...
                sum(keepBackward), length(keepBackward), 100*(1 - sum(keepBackward)/length(keepBackward)));
        end

        % Merge
        if isempty(tracksForward)
            finalTracks = tracksBackward;
        elseif isempty(tracksBackward)
            finalTracks = tracksForward;
        else
            finalTracks = [tracksForward; tracksBackward];
        end
        fprintf('[RESULT] Combined %d Forward + %d Backward = %d total tracks.\n', ...
            length(tracksForward), length(tracksBackward), length(finalTracks));
    else
        finalTracks = tracksForward;
        fprintf('\n[RESULT] Generated %d Forward tracks.\n', length(finalTracks));
    end
    fprintf('=======================================================\n\n');
end

% =========================================================================
% HELPER FUNCTION: runHierarchicalPass
% =========================================================================
% Runs the full hierarchical tracking pipeline for one temporal direction.
% For each velocity level:
%   1. Configure noise parameters scaled to the level's v_max.
%   2. Run the Kalman tracker on remaining localizations.
%   3. Estimate each track's mean velocity; reject tracks outside [v_min, v_max].
%   4. Subtract valid tracks' localizations from the remaining pool.
% =========================================================================
function allTracks = runHierarchicalPass(localizations, params, directionStr)

    velocityLevels = params.track.kalman.velocity_levels;
    numLevels = length(velocityLevels);
    allTracksCells = cell(numLevels, 1);  % Accumulate per-level results
    remainingLocs = localizations;

    avg_pixel_size_mm = mean([params.track.pixel_X_size, params.track.pixel_Z_size]);
    dt = params.track.dt;
    pixel_X = params.track.pixel_X_size;
    pixel_Z = params.track.pixel_Z_size;

    for level = 1:numLevels
        v_min_mm_s = velocityLevels{level}(1);
        v_max_mm_s = velocityLevels{level}(2);

        levelParams = params;
        
        %   sigma_epsilon = alpha * v_max
        %   sigma_nu      = beta / 2^(level-1)
        levelParams.track.kalman.process_noise  = params.track.kalman.hk_alpha * v_max_mm_s;
        levelParams.track.kalman.measurement_noise = params.track.kalman.hk_beta / (2^(level - 1));

        % Maximum linking distance for this level
        max_displacement_mm = v_max_mm_s * dt;
        levelParams.track.max_linking_distance = ceil(max_displacement_mm / avg_pixel_size_mm);
        levelParams.track.kalman.assignment_method = 'hungarian';

        fprintf('   [%s] Level %d/%d (v: %.1f-%.1f mm/s) | sigma_e=%.5f, sigma_v=%.5f | max_link=%d px | Remaining MBs: %d\n', ...
            upper(directionStr), level, numLevels, v_min_mm_s, v_max_mm_s, ...
            levelParams.track.kalman.process_noise, levelParams.track.kalman.measurement_noise, ...
            levelParams.track.max_linking_distance, height(remainingLocs));

        if height(remainingLocs) == 0
            break;
        end

        levelTracks = trackSingleLevelKalman(remainingLocs, levelParams);

        % --- [FIX-2] Velocity-Range Filtering ---
        % Estimate each track's mean velocity and reject tracks outside
        % this level's [v_min, v_max] range.
        if ~isempty(levelTracks)
            trackVelocities = zeros(length(levelTracks), 1);
            for i = 1:length(levelTracks)
                t = levelTracks(i);
                if size(t.path, 1) >= 2
                    dx_px = diff(t.path(:, 1));
                    dy_px = diff(t.path(:, 2));
                    dt_frames = diff(t.frames);
                    displacement_mm = sqrt((dx_px * pixel_X).^2 + (dy_px * pixel_Z).^2);
                    dt_sec = dt_frames * dt;
                    valid_dt = dt_sec > 0;
                    if any(valid_dt)
                        trackVelocities(i) = mean(displacement_mm(valid_dt) ./ dt_sec(valid_dt));
                    end
                end
            end

            % Keep only tracks within this level's velocity range
            inRange = (trackVelocities >= v_min_mm_s*0.8) & (trackVelocities <= v_max_mm_s*1.2);
            rejectedTracks = levelTracks(~inRange);
            levelTracks = levelTracks(inRange);

            fprintf('   [%s]   -> %d tracks formed, %d passed velocity filter [%.1f-%.1f mm/s], %d rejected.\n', ...
                upper(directionStr), sum(inRange) + sum(~inRange), sum(inRange), ...
                v_min_mm_s, v_max_mm_s, sum(~inRange));
        end

        % Subtractive Step: Remove localizations used by VALID tracks only
        if ~isempty(levelTracks)
            usedIDs = arrayfun(@(t) t.localizations.LocalizationID, levelTracks, 'UniformOutput', false);
            usedIDs = vertcat(usedIDs{:});
            remainingLocs(ismember(remainingLocs.LocalizationID, usedIDs), :) = [];
            allTracksCells{level} = levelTracks;
        end
    end

    % Concatenate all levels
    nonEmpty = ~cellfun(@isempty, allTracksCells);
    if any(nonEmpty)
        allTracks = vertcat(allTracksCells{nonEmpty});
    else
        allTracks = [];
    end
end

% =========================================================================
% HELPER FUNCTION: trackSingleLevelKalman
% =========================================================================
% Runs the Kalman tracking loop for a single velocity level.
%
% Key features:
%   - Uses trackingKF (Sensor Fusion Toolbox) for full Kalman state access,
%     including residual() for Mahalanobis distance in the advanced cost matrix.
%   - Two-step track initiation: detections must appear unassigned in two
%     consecutive frames to spawn a track.
%   - Pre-groups localizations by frame for efficient iteration.
%   - Uses cell-based track storage to avoid repeated struct array growth.
% =========================================================================
function validTracks = trackSingleLevelKalman(localizations, params)

    frames = unique(localizations.Frame);
    if isempty(frames), validTracks = []; return; end

    % --- [FIX-7] Pre-group localizations by frame for performance ---
    frameGroups = findgroups(localizations.Frame);
    uniqueFrames = unique(localizations.Frame);
    locsByFrame = splitapply(@(varargin) {table(varargin{:}, ...
        'VariableNames', localizations.Properties.VariableNames)}, ...
        localizations, frameGroups);

    % Build a frame-number-to-index map
    frameMap = containers.Map(double(uniqueFrames), num2cell(1:length(uniqueFrames)));

    min_frame = min(uniqueFrames);
    max_frame = max(uniqueFrames);

    % --- Track storage as cell array to avoid struct array growth ---
    activeTracks = {};
    terminatedTracks = {};
    nextId = 1;
    max_gap = params.track.max_gap_closing_frames;

    % Kalman model configuration
    motionModel = params.track.kalman.motion_model;
    processNoise = params.track.kalman.process_noise;
    measurementNoise = params.track.kalman.measurement_noise;

    % Determine state indices based on motion model
    if strcmpi(motionModel, 'ConstantAcceleration')
        stateIdxX = 1; stateIdxY = 4;
        kfMotionStr = '2D Constant Acceleration';
    else
        stateIdxX = 1; stateIdxY = 3;
        kfMotionStr = '2D Constant Velocity';
    end

    % --- [FIX-4] Two-step initiation: store previous frame's unassigned locs ---
    unassignedLocs_prevFrame = table();

    for f = min_frame:max_frame
        % Retrieve current frame's localizations efficiently
        if frameMap.isKey(f)
            current_locs = locsByFrame{frameMap(f)};
        else
            current_locs = localizations([], :);  % Empty table with correct columns
        end
        numDetections = height(current_locs);
        numTracks = length(activeTracks);

        % 1. PREDICT
        for i = 1:numTracks
            predict(activeTracks{i}.kalmanFilter);
        end

        unassignedTracksList = (1:numTracks)';
        unassignedDetections = (1:numDetections)';

        % 2. ASSIGNMENT using external calculateCostMatrix
        if numTracks > 0 && numDetections > 0

            % Execute the unified cost matrix function
            costMatrix = calculateCostMatrix(activeTracks, current_locs, params);

            % Edge Case: Protect Munkres from all-Inf matrix
            if all(costMatrix(:) == Inf)
                assignments = zeros(numTracks, 1);
            else
                assignments = munkres(costMatrix);
            end

            assignedTrackIdxs = [];
            assignedDetIdxs = [];
            for i = 1:length(assignments)
                j = assignments(i);
                if j > 0 && costMatrix(i, j) ~= Inf
                    assignedTrackIdxs(end+1) = i; %#ok<AGROW>
                    assignedDetIdxs(end+1) = j;   %#ok<AGROW>

                    detLoc = [current_locs.X(j), current_locs.Y(j)];

                    % --- [FIX-3] correct() works with trackingKF ---
                    correct(activeTracks{i}.kalmanFilter, detLoc);
                    corrState = activeTracks{i}.kalmanFilter.State;
                    correctedPos = [corrState(stateIdxX), corrState(stateIdxY)];

                    activeTracks{i}.age = activeTracks{i}.age + 1;
                    activeTracks{i}.totalVisibleCount = activeTracks{i}.totalVisibleCount + 1;
                    activeTracks{i}.consecutiveInvisibleCount = 0;
                    activeTracks{i}.lastFrame = f;

                    activeTracks{i}.localizations = [activeTracks{i}.localizations; current_locs(j, :)];
                    activeTracks{i}.path = [activeTracks{i}.path; detLoc];
                    activeTracks{i}.correctedPath = [activeTracks{i}.correctedPath; correctedPos];
                    activeTracks{i}.frames = [activeTracks{i}.frames; f];
                end
            end
            unassignedTracksList = setdiff(unassignedTracksList, assignedTrackIdxs);
            unassignedDetections = setdiff(unassignedDetections, assignedDetIdxs);
        end

        % 3. MANAGE LOST TRACKS
        for i = 1:length(unassignedTracksList)
            idx = unassignedTracksList(i);
            activeTracks{idx}.age = activeTracks{idx}.age + 1;
            activeTracks{idx}.consecutiveInvisibleCount = ...
                activeTracks{idx}.consecutiveInvisibleCount + 1;
        end

        % 4. DELETE DEAD TRACKS (iterate backward to avoid index shifting)
        for i = length(activeTracks):-1:1
            if activeTracks{i}.consecutiveInvisibleCount > max_gap
                terminatedTracks{end+1} = activeTracks{i}; %#ok<AGROW>
                activeTracks(i) = [];
            end
        end

        % 5. TWO-STEP TRACK INITIATION [FIX-4]
        % Get unassigned detections for this frame
        if ~isempty(unassignedDetections)
            unassignedLocs_currentFrame = current_locs(unassignedDetections, :);
        else
            unassignedLocs_currentFrame = current_locs([], :);
        end

        if ~isempty(unassignedLocs_currentFrame) && ~isempty(unassignedLocs_prevFrame)
            % Try linking current unassigned locs with previous frame's unassigned locs
            costMatrix_init = pdist2( ...
                [unassignedLocs_prevFrame.X, unassignedLocs_prevFrame.Y], ...
                [unassignedLocs_currentFrame.X, unassignedLocs_currentFrame.Y], ...
                'squaredeuclidean');

            init_dist_gate_sq = (params.track.max_linking_distance * 1.5)^2;
            costMatrix_init(costMatrix_init > init_dist_gate_sq) = Inf;

            if ~all(costMatrix_init(:) == Inf)
                assignments_init = munkres(costMatrix_init);
            else
                assignments_init = zeros(height(unassignedLocs_prevFrame), 1);
            end

            used_current_locs_mask = false(height(unassignedLocs_currentFrame), 1);

            for i = 1:length(assignments_init)
                if assignments_init(i) > 0
                    detIdx = assignments_init(i);
                    loc1 = unassignedLocs_prevFrame(i, :);
                    loc2 = unassignedLocs_currentFrame(detIdx, :);

                    % --- [FIX-3] Initialize with trackingKF ---
                    kf = initializeTrackingKF(loc1, motionModel, kfMotionStr, ...
                        processNoise, measurementNoise);

                    initialState = kf.State;

                    % Immediately correct with second point
                    correct(kf, [loc2.X, loc2.Y]);
                    corrState2 = kf.State;

                    pos1 = [initialState(stateIdxX), initialState(stateIdxY)];
                    pos2 = [corrState2(stateIdxX), corrState2(stateIdxY)];

                    newTrack = struct();
                    newTrack.id = nextId; nextId = nextId + 1;
                    newTrack.kalmanFilter = kf;
                    newTrack.localizations = [loc1; loc2];
                    newTrack.path = [loc1.X, loc1.Y; loc2.X, loc2.Y];
                    newTrack.correctedPath = [pos1; pos2];
                    newTrack.frames = [f-1; f];
                    newTrack.age = 2;
                    newTrack.totalVisibleCount = 2;
                    newTrack.consecutiveInvisibleCount = 0;
                    newTrack.lastFrame = f;

                    activeTracks{end+1} = newTrack; %#ok<AGROW>
                    used_current_locs_mask(detIdx) = true;
                end
            end

            % Remaining unassigned current locs become candidates for next frame
            unassignedLocs_prevFrame = unassignedLocs_currentFrame(~used_current_locs_mask, :);
        else
            unassignedLocs_prevFrame = unassignedLocs_currentFrame;
        end

    end  % End frame loop

    % Move remaining active tracks to terminated
    terminatedTracks = [terminatedTracks, activeTracks];

    % --- PHASE 3: QUALITY CONTROL (QC) FILTERING ---
    if isempty(terminatedTracks)
        validTracks = [];
        return;
    end

    % Filter by minimum track length
    keepMask = cellfun(@(t) t.totalVisibleCount >= params.track.min_track_length, terminatedTracks);
    terminatedTracks = terminatedTracks(keepMask);

    % --- [FIX-8] Mask filtering using params.mask.combinedMask ---
    maskActive = isfield(params, 'mask') && isfield(params.mask, 'combinedMask') ...
                 && ~isempty(params.mask.combinedMask);
    if maskActive && ~isempty(terminatedTracks)
        combinedMask = params.mask.combinedMask;
        [H, W] = size(combinedMask);
        tracks_to_remove = false(1, length(terminatedTracks));

        for i = 1:length(terminatedTracks)
            currTrack = terminatedTracks{i};
            if isfield(currTrack, 'correctedPath') && ~isempty(currTrack.correctedPath)
                pts = currTrack.correctedPath;
            else
                pts = currTrack.path;
            end

            x_idx = round(pts(:, 1));
            y_idx = round(pts(:, 2));

            in_bounds = (x_idx >= 1 & x_idx <= W & y_idx >= 1 & y_idx <= H);

            if any(~in_bounds)
                tracks_to_remove(i) = true;
                continue;
            end

            indices = sub2ind([H, W], y_idx, x_idx);
            if any(~combinedMask(indices))
                tracks_to_remove(i) = true;
            end
        end
        terminatedTracks(tracks_to_remove) = [];

        num_removed = sum(tracks_to_remove);
        if num_removed > 0
            fprintf('      -> Mask QC: Removed %d tracks outside combinedMask.\n', num_removed);
        end
    end

    % --- [FIX-6] FORMAT OUTPUT (compatible with standard trackKalman.m) ---
    validTracks = formatTracksAdvanced(terminatedTracks, params);
end

% =========================================================================
% HELPER FUNCTION: initializeTrackingKF
% =========================================================================
% Creates a trackingKF object (Sensor Fusion Toolbox) instead of
% configureKalmanFilter (Computer Vision Toolbox).
% This ensures compatibility with residual() for Mahalanobis distance
% in the advanced cost matrix.
% =========================================================================
function kf = initializeTrackingKF(loc, motionModel, kfMotionStr, processNoise, measurementNoise)
    if strcmpi(motionModel, 'ConstantAcceleration')
        initialState = [loc.X; 0; 0; loc.Y; 0; 0];
    else  % ConstantVelocity
        initialState = [loc.X; 0; loc.Y; 0];
    end

    kf = trackingKF('MotionModel', kfMotionStr, 'State', initialState);
    kf.ProcessNoise = processNoise;
    kf.MeasurementNoise = eye(2) * measurementNoise;
end

% =========================================================================
% HELPER FUNCTION: formatTracksAdvanced
% =========================================================================
% Converts raw cell-array tracks into a standardized struct array with
% velocity calculations in mm/s, matching the output format of trackKalman.m.
% =========================================================================
function outTracks = formatTracksAdvanced(terminatedTracks, params)

    if isempty(terminatedTracks)
        outTracks = [];
        return;
    end

    pixel_X = params.track.pixel_X_size;
    pixel_Z = params.track.pixel_Z_size;
    dt = params.track.dt;

    numTracks = length(terminatedTracks);
    outCell = cell(numTracks, 1);

    for i = 1:numTracks
        t = terminatedTracks{i};

        out = struct();
        out.id = t.id;
        out.localizations = t.localizations;

        % Use correctedPath if available, fall back to raw path
        if isfield(t, 'correctedPath') && ~isempty(t.correctedPath)
            out.path = t.correctedPath;
            out.correctedPath = t.correctedPath;
        else
            out.path = t.path;
            out.correctedPath = t.path;
        end

        out.frames = t.frames;
        out.length = t.totalVisibleCount;

        % --- Velocity Calculation (mm/s) ---
        if out.length > 1
            dx_px = diff(out.path(:, 1));
            dy_px = diff(out.path(:, 2));
            dt_frames = diff(out.frames);

            displacement_mm = sqrt((dx_px * pixel_X).^2 + (dy_px * pixel_Z).^2);
            dt_sec = dt_frames * dt;

            velocities = zeros(size(dt_sec));
            valid_dt = dt_sec > 0;
            velocities(valid_dt) = displacement_mm(valid_dt) ./ dt_sec(valid_dt);
            % Replicate last velocity for equal-length vector
            out.velocities_mm_s = [velocities; velocities(end)];

            mean_vel = mean(velocities(isfinite(velocities)), 'omitnan');
            if isempty(mean_vel) || isnan(mean_vel), mean_vel = 0; end
            out.average_velocity_mm_s = repmat(mean_vel, out.length, 1);
        else
            out.velocities_mm_s = 0;
            out.average_velocity_mm_s = 0;
        end

        outCell{i} = out;
    end

    outTracks = [outCell{:}]';
end