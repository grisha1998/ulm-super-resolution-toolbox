function tracks = trackKalman(localizations, params, indent_prefix)
% =========================================================================
% FUNCTION: trackKalman
% AUTHOR:   Grigori Shapiro
% =========================================================================
%
% PURPOSE:
%   Implements predictive microbubble tracking using a 2D Kalman Filter,
%   designed for demanding vascular beds such as the rat kidney. By
%   maintaining kinematic state (position + velocity, optionally
%   acceleration) for each microbubble, the tracker predicts where a bubble
%   will appear in the next frame(s). This reduces identity swaps in dense
%   regions, smooths sub-pixel localization jitter, and allows bubbles to
%   be recovered after temporary disappearance at branching points or due
%   to momentary SNR drops.
%
% KEY DESIGN PRINCIPLES:
%   - Velocity-bootstrapped initiation: New tracks require two consecutive
%     detections (two-step initiation). The initial velocity is computed
%     from the loc1->loc2 displacement, giving the filter a physically
%     meaningful first prediction instead of defaulting to zero velocity.
%   - Frame-synchronous prediction: The main loop iterates over real frame
%     numbers (min_frame:max_frame) rather than frame indices. This ensures
%     predict() is called exactly once per physical frame, correctly
%     handling empty frames without desynchronizing the filter's time.
%   - Anisotropic measurement noise: Separate variances for lateral (X)
%     and axial (Z) dimensions, reflecting the inherent resolution
%     asymmetry of ultrasound imaging.
%   - Adaptive process noise: Q is dynamically scaled by recent track
%     speed, accommodating the wide velocity range of kidney vasculature
%     (arcuate arteries ~20 mm/s down to peritubular capillaries <1 mm/s).
%   - Mask-tolerant validation: Tracks are rejected only when a configurable
%     fraction of their path falls outside the vessel mask, rather than
%     discarding entire tracks for single boundary-edge points.
%   - Gap interpolation: Predicted positions during gap frames are stored
%     separately and can be merged into a dense output path for density-map
%     rendering without contaminating velocity calculations.
%
% ALGORITHM WORKFLOW:
%   1. Mask Initialization — checks params.proc.vesselMask; enables
%      physical boundary constraints if present.
%   2. Core Loop (iterates frame numbers min_frame:max_frame):
%      a. Predict: every active track's filter is stepped forward.
%      b. Cost & Assignment: cost matrix (from calculateCostMatrix_v2) is
%         solved via Hungarian or greedy nearest-neighbor.
%      c. Update: for assigned links, correct() the filter, clamp to mask
%         if needed, update adaptive process noise.
%      d. Unassigned track management: store predicted position for gap
%         interpolation, increment gap counter, terminate if exceeded.
%      e. Two-step initiation: unassigned detections from consecutive
%         verified frames spawn new velocity-bootstrapped tracks.
%   3. Finalization — mask-tolerance validation, path formatting, velocity
%      computation in mm/s, optional QC.
%
% INPUT PARAMETERS:
%   localizations - (Table) Required columns: 'Frame', 'X', 'Y'.
%                   Optional column: 'Intensity'.
%   params        - (Struct) Tracking configuration. Relevant fields:
%     params.track.method                       'Kalman' or 'Kalman_v2'
%     params.track.max_linking_distance         Base distance gate (pixels)
%     params.track.max_gap_closing_frames       Max gap frames before termination
%     params.track.min_track_length             Min corrected frames for valid track
%     params.track.pixel_X_size, pixel_Z_size   mm/pixel in X and Z
%     params.track.dt                           Seconds between frames
%     params.track.use_advanced_cost_matrix     Enable cost matrix penalties
%     params.track.kalman.motion_model          'ConstantVelocity' | 'ConstantAcceleration'
%     params.track.kalman.process_noise         Base process noise Q (scalar or matrix)
%     params.track.kalman.assignment_method     'hungarian' | 'nn'
%     params.track.kalman.process_noise_velocity_scale Adaptive Q scaling (0=off).
%                                                      Q_eff = Q * (1 + scale * speed).
%     params.track.kalman.mask_tolerance_fraction      Max outside-mask path fraction
%     params.track.kalman.interpolate_gaps             Include gap predictions in output
%     params.track.kalman.min_track_length_velocity    Min length for reliable velocity
%     params.track.kalman.gating_max_angle_change_deg  Hard angle gate (degrees)
%     params.track.kalman.max_angle_change_deg         Soft angle threshold (degrees)
%     params.track.kalman.angle_penalty_slope          Angle penalty ramp slope
%     params.track.kalman.direction_penalty_weight     Direction penalty weight
%     params.track.kalman.brightness_penalty_weight    Brightness penalty weight
%     params.loc.fwhm                           PSF FWHM [lateral, axial] (pixels)
%     params.proc.vesselMask                    (optional) vessel boundary mask
%
% INTERNAL TRACK STATE:
%   .kalmanFilter              - (trackingKF) Kalman filter handle
%   .correctedPath             - [Nx2] Kalman-smoothed trajectory
%   .predictedPath_gap         - [Mx2] predicted positions during gap frames
%   .predictedFrames_gap       - [Mx1] frame numbers for gap predictions
%   .localizations             - (table) raw measurements assigned to track
%   .age                       - count of corrected (measured) frames
%   .lastFrame                 - last frame with a corrected measurement
%   .consecutiveInvisibleCount - current gap frame counter
%   .meanBrightness            - incrementally maintained mean intensity
%
% OUTPUT:
%   tracks - Struct array. Each track contains:
%     .id                       Integer ID
%     .path                     [Nx2] Kalman-smoothed trajectory (age rows)
%     .path_interpolated        [Kx2] dense path including gap predictions
%     .frames                   Frame numbers of corrected positions
%     .frames_interpolated      Frame numbers including gap frames
%     .length                   Number of corrected frames (= age)
%     .localizations            Raw localization table
%     .velocities_mm_s          Per-step velocity in mm/s
%     .average_velocity_mm_s    Mean velocity (scalar broadcast to length)
%     .has_reliable_velocity    true if track meets min length for velocity
%
% DEPENDENCIES:
%   - Sensor Fusion and Tracking Toolbox (trackingKF, predict, correct)
%   - munkres.m on path (for Hungarian assignment)
%   - calculateCostMatrix_v2.m (companion cost matrix function)
%   - Optional: applyQualityControl.m
% =========================================================================

    % --- Argument defaults ---
    if nargin < 3, indent_prefix = ''; end

    if isempty(localizations)
        tracks = [];
        return;
    end

    if ~license('test', 'Sensor_Fusion_and_Tracking')
        error('The Kalman tracker requires the Sensor Fusion and Tracking Toolbox.');
    end

    % --- Run the core tracking pass ---
    terminatedTracks = runTrackingPass(localizations, params, indent_prefix);

    % --- Validation & Formatting ---
    terminatedTracks = validateTracks(terminatedTracks, params);
    tracks           = formatTracks(terminatedTracks, params, indent_prefix);

    % --- Summary ---
    printTrackingSummary(tracks, localizations, indent_prefix);
end

%% ========================================================================
%  CORE TRACKING LOOP
% =========================================================================

function terminatedTracks = runTrackingPass(localizations, params, indent_prefix)
% RUNTRACKINGPASS  Executes the main temporal tracking loop.
%
%   Iterates over every physical frame in the dataset (including empty
%   frames). For each frame: predicts all active tracks, computes the
%   assignment cost matrix, solves the assignment, updates matched tracks,
%   manages unmatched tracks (gap counting, termination), and initiates
%   new tracks from consecutive unassigned detections.

    % --- Mask setup ---
    use_vessel_mask = false;
    vesselMask = [];
    maskH = 0; maskW = 0;
    if isfield(params, 'proc') && isfield(params.proc, 'vesselMask') && ...
            isfield(params.proc.vesselMask, 'enable') && params.proc.vesselMask.enable && ...
            isfield(params.proc.vesselMask, 'vesselMask') && ~isempty(params.proc.vesselMask.vesselMask)
        vesselMask = params.proc.vesselMask.vesselMask;
        use_vessel_mask = true;
        [maskH, maskW] = size(vesselMask);
    end

    % --- State-vector indices for the active motion model ---
    if strcmpi(params.track.kalman.motion_model, 'ConstantAcceleration')
        idx_x = 1; idx_vx = 2;   % State: [x; vx; ax; y; vy; ay]
        idx_y = 4; idx_vy = 5;
    else
        idx_x = 1; idx_vx = 2;   % State: [x; vx; y; vy]
        idx_y = 3; idx_vy = 4;
    end

    activeTracks       = {};
    terminatedTracks   = {};
    nextTrackID        = 1;
    assignment_method  = lower(params.track.kalman.assignment_method);

    fprintf('%s  Kalman v2 tracking pass...\n', indent_prefix);

    % --- Build a frame-number -> localization-table lookup ---
    % This allows iterating over the full frame range (including empty
    % frames) so that predict() is called exactly once per physical frame.
    all_frames = unique(localizations.Frame);
    min_frame  = min(all_frames);
    max_frame  = max(all_frames);

    [G, frame_ids] = findgroups(localizations.Frame);
    locsByFrameMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for g = 1:max(G)
        locsByFrameMap(frame_ids(g)) = localizations(G == g, :);
    end

    % --- Two-step initiation state ---
    % Tracks the previous frame's unassigned detections and frame number
    % to enforce frame continuity (consecutive frames only).
    unassignedLocs_prevFrame = table();
    prevFrameNumber          = -Inf;

    % --- Main temporal loop ---
    for currentFrame = min_frame:max_frame

        % Retrieve current frame's localizations (may be empty)
        if isKey(locsByFrameMap, currentFrame)
            currentLocs = locsByFrameMap(currentFrame);
        else
            currentLocs = [];
        end

        % === PREDICT =========================================================
        % Step every active track's Kalman filter forward by one time step.
        numActiveTracks = length(activeTracks);
        if numActiveTracks > 0
            cellfun(@(track) predict(track.kalmanFilter), activeTracks, 'UniformOutput', false);
        end

        % === ASSIGN ==========================================================
        % Build cost matrix and solve the assignment problem.
        hasLocs        = ~isempty(currentLocs) && height(currentLocs) > 0;
        if hasLocs
            unassignedDet = true(height(currentLocs), 1);
        else
            unassignedDet = false(0, 1);
        end
        assignedTracks = false(numActiveTracks, 1);

        if numActiveTracks > 0 && hasLocs
            costMatrix = calculateCostMatrix(activeTracks, currentLocs, params);

            switch assignment_method
                case 'hungarian'
                    if ~exist('munkres', 'file')
                        error('Hungarian algorithm (munkres) not found on MATLAB path.');
                    end
                    [assignments, ~] = munkres(costMatrix);
                case 'nn'
                    max_link_dist_sq = params.track.max_linking_distance^2;
                    assignments = greedyNN_assign(costMatrix, max_link_dist_sq);
                otherwise
                    error('Unknown assignment method: %s', assignment_method);
            end

            % === UPDATE ======================================================
            % For each assigned track-detection pair: correct the Kalman
            % filter, enforce mask constraints, update track state.
            currentPositions = [currentLocs.X, currentLocs.Y];
            for i = 1:numActiveTracks
                if assignments(i) > 0
                    detectionIdx = assignments(i);

                    % Correct the Kalman filter with the new measurement
                    correct(activeTracks{i}.kalmanFilter, currentPositions(detectionIdx, :));

                    % Extract the corrected (smoothed) position
                    s = activeTracks{i}.kalmanFilter.State;
                    corrected_pos = [s(idx_x), s(idx_y)];

                    % Enforce vessel mask boundary constraint: if the
                    % corrected position drifts outside the mask, snap
                    % it back to the raw measurement position.
                    if use_vessel_mask
                        [corrected_pos, was_clamped] = applyMaskConstraint( ...
                            corrected_pos, currentPositions(detectionIdx, :), vesselMask, maskH, maskW);
                        if was_clamped
                            activeTracks{i}.kalmanFilter.State(idx_x) = corrected_pos(1);
                            activeTracks{i}.kalmanFilter.State(idx_y) = corrected_pos(2);
                        end
                    end

                    % Append to the Kalman-smoothed trajectory
                    activeTracks{i}.correctedPath = [activeTracks{i}.correctedPath; corrected_pos];
                    activeTracks{i}.localizations = [activeTracks{i}.localizations; currentLocs(detectionIdx, :)];
                    activeTracks{i}.age                       = activeTracks{i}.age + 1;
                    activeTracks{i}.lastFrame                 = currentFrame;
                    activeTracks{i}.consecutiveInvisibleCount = 0;

                    % Incrementally update mean brightness (avoids
                    % recomputing over the full history each frame)
                    if ismember('Intensity', activeTracks{i}.localizations.Properties.VariableNames)
                        n = activeTracks{i}.age;
                        new_I = currentLocs.Intensity(detectionIdx);
                        activeTracks{i}.meanBrightness = ...
                            ((n-1) * activeTracks{i}.meanBrightness + new_I) / n;
                    end

                    % Adaptive process noise: scale Q by recent track speed
                    % so fast tracks can turn more freely while slow tracks
                    % remain tightly smoothed. Disabled when scale = 0.
                    if params.track.kalman.process_noise_velocity_scale > 0
                        activeTracks{i}.kalmanFilter = updateAdaptiveProcessNoise(...
                            activeTracks{i}.kalmanFilter, activeTracks{i}.correctedPath, params);
                    end

                    unassignedDet(detectionIdx) = false;
                    assignedTracks(i)           = true;
                end
            end
        end

        % === HANDLE UNASSIGNED TRACKS (gap management) =======================
        % For unassigned tracks: increment their gap counter and store the
        % predicted position for optional gap interpolation in the output.
        for i = 1:numActiveTracks
            if ~assignedTracks(i)
                activeTracks{i}.consecutiveInvisibleCount = activeTracks{i}.consecutiveInvisibleCount + 1;

                % Store predicted position during this gap frame
                s = activeTracks{i}.kalmanFilter.State;
                pred_pos = [s(idx_x), s(idx_y)];
                activeTracks{i}.predictedPath_gap   = [activeTracks{i}.predictedPath_gap;   pred_pos];
                activeTracks{i}.predictedFrames_gap = [activeTracks{i}.predictedFrames_gap; currentFrame];
            end
        end

        % === TERMINATION =====================================================
        % Tracks exceeding max_gap_closing_frames are terminated.
        to_remove = false(1, numActiveTracks);
        max_gap   = params.track.max_gap_closing_frames;

        for i = 1:numActiveTracks
            if activeTracks{i}.consecutiveInvisibleCount > max_gap
                terminatedTracks{end+1} = activeTracks{i}; %#ok<AGROW>
                to_remove(i) = true;
            end
        end
        activeTracks(to_remove) = [];

        % === TWO-STEP INITIATION =============================================
        % New tracks are only created when an unassigned detection in the
        % current frame can be linked to an unassigned detection in the
        % immediately preceding frame (frame continuity enforced). The
        % initial velocity is bootstrapped from the displacement between
        % the two detections, giving the Kalman filter a physically
        % meaningful first prediction.
        if hasLocs
            unassignedLocs_currentFrame = currentLocs(unassignedDet, :);
        else
            unassignedLocs_currentFrame = table();
        end

        frame_delta  = currentFrame - prevFrameNumber;
        can_initiate = (frame_delta == 1) && ...
                       ~isempty(unassignedLocs_currentFrame) && ...
                       ~isempty(unassignedLocs_prevFrame);

        if can_initiate
            costMatrix_init = pdist2([unassignedLocs_prevFrame.X, unassignedLocs_prevFrame.Y], ...
                                     [unassignedLocs_currentFrame.X, unassignedLocs_currentFrame.Y], ...
                                     'squaredeuclidean');
            init_gate_sq = (params.track.max_linking_distance * 1.5)^2;
            costMatrix_init(costMatrix_init > init_gate_sq) = Inf;

            % Use the same assignment method as the main tracking loop
            switch assignment_method
                case 'nn'
                    assignments_init = greedyNN_assign(costMatrix_init, init_gate_sq);
                otherwise
                    [assignments_init, ~] = munkres(costMatrix_init);
            end

            used_current = false(height(unassignedLocs_currentFrame), 1);

            for i = 1:length(assignments_init)
                if assignments_init(i) > 0
                    detectionIdx = assignments_init(i);
                    loc1 = unassignedLocs_prevFrame(i, :);
                    loc2 = unassignedLocs_currentFrame(detectionIdx, :);

                    % Initialize a Kalman filter at loc1's position with
                    % anisotropic measurement noise
                    kf = initializeKalmanFilter(loc1, params);

                    % Bootstrap velocity from the loc1->loc2 displacement.
                    % Without this, the filter starts at vx=vy=0 and
                    % correct() alone cannot infer velocity from a single
                    % position measurement.
                    vx = loc2.X - loc1.X;
                    vy = loc2.Y - loc1.Y;
                    kf.State(idx_vx) = vx;
                    kf.State(idx_vy) = vy;

                    % Advance filter one step and correct with loc2
                    predict(kf);
                    correct(kf, [loc2.X, loc2.Y]);
                    s = kf.State;

                    pos1 = [loc1.X, loc1.Y];
                    pos2 = [s(idx_x), s(idx_y)];

                    % Apply mask constraint to the corrected position
                    if use_vessel_mask
                        [pos2, was_clamped] = applyMaskConstraint( ...
                            pos2, [loc2.X, loc2.Y], vesselMask, maskH, maskW);
                        if was_clamped
                            kf.State(idx_x) = pos2(1);
                            kf.State(idx_y) = pos2(2);
                        end
                    end

                    % Create the new track structure
                    newTrack = createNewTrack(nextTrackID, kf, loc1, loc2, pos1, pos2, currentFrame);
                    nextTrackID = nextTrackID + 1;
                    activeTracks{end+1} = newTrack; %#ok<AGROW>

                    used_current(detectionIdx) = true;
                end
            end

            unassignedLocs_prevFrame = unassignedLocs_currentFrame(~used_current, :);
        else
            % Previous frame was not immediately prior — reset to current
            unassignedLocs_prevFrame = unassignedLocs_currentFrame;
        end

        prevFrameNumber = currentFrame;
    end

    % Append any tracks still active at the end of the sequence
    terminatedTracks = [terminatedTracks, activeTracks];
end

%% ========================================================================
%  KALMAN FILTER INITIALIZATION & HELPERS
% =========================================================================

function kf = initializeKalmanFilter(loc, params)
% INITIALIZEKALMANFILTER  Creates a new trackingKF with anisotropic
%   measurement noise reflecting the different lateral (X) vs axial (Z)
%   resolution of ultrasound imaging.
%
%   The measurement noise R is derived from the PSF FWHM: sigma = 0.5*FWHM.
%   When FWHM is a 2-element vector [lateral, axial], R is diagonal with
%   separate variances per axis. When scalar, R is isotropic.

    kalmanParams = params.track.kalman;
    motionModel  = kalmanParams.motion_model;

    if strcmpi(motionModel, 'ConstantAcceleration')
        initialState = [loc.X; 0; 0; loc.Y; 0; 0];
        kf = trackingKF('MotionModel', '2D Constant Acceleration', 'State', initialState);
    else
        initialState = [loc.X; 0; loc.Y; 0];
        kf = trackingKF('MotionModel', '2D Constant Velocity', 'State', initialState);
    end

    kf.ProcessNoise = kalmanParams.process_noise;

    % Anisotropic measurement noise from PSF FWHM
    fwhm = params.loc.fwhm;
    if numel(fwhm) == 1
        sigma_x = 0.5 * fwhm;
        sigma_z = 0.5 * fwhm;
    else
        sigma_x = 0.5 * fwhm(1);
        sigma_z = 0.5 * fwhm(2);
    end
    kf.MeasurementNoise = diag([sigma_x^2, sigma_z^2]);
end

function kf = updateAdaptiveProcessNoise(kf, correctedPath, params)
% UPDATEADAPTIVEPROCESSNOISE  Scales ProcessNoise by recent inter-frame
%   displacement so that fast tracks can accelerate/turn more freely while
%   slow tracks remain tightly smoothed. This accommodates the wide dynamic
%   velocity range encountered in kidney vasculature.
%
%   Formula: Q_effective = Q_base * (1 + scale * recent_speed)

    if size(correctedPath, 1) < 2
        return;
    end

    recent_velocity = norm(correctedPath(end, :) - correctedPath(end-1, :));
    base_Q = params.track.kalman.process_noise;
    scale  = params.track.kalman.process_noise_velocity_scale;

    try
        kf.ProcessNoise = base_Q * (1 + scale * recent_velocity);
    catch
        % Graceful fallback if base_Q dimensionality is incompatible
    end
end

function [corrected_pos, was_clamped] = applyMaskConstraint(corrected_pos, raw_pos, vesselMask, maskH, maskW)
% APPLYMASKCONSTRAINT  Enforces vessel boundary constraints.
%   If the Kalman-corrected position falls outside the vessel mask (either
%   out of image bounds or on a non-vessel pixel), it is snapped back to
%   the raw measurement position to prevent track divergence at vessel
%   boundaries.

    cx = round(corrected_pos(1));
    cy = round(corrected_pos(2));

    is_outside = (cx < 1 || cx > maskW || cy < 1 || cy > maskH);
    if ~is_outside && ~vesselMask(cy, cx)
        is_outside = true;
    end

    was_clamped = false;
    if is_outside
        corrected_pos = raw_pos;
        was_clamped   = true;
    end
end

function newTrack = createNewTrack(id, kf, loc1, loc2, pos1, pos2, currentFrame)
% CREATENEWTRACK  Assembles a new track structure from two initial
%   detections and a velocity-bootstrapped Kalman filter.

    newTrack = struct();
    newTrack.id                        = id;
    newTrack.kalmanFilter              = kf;
    newTrack.localizations             = [loc1; loc2];
    newTrack.correctedPath             = [pos1; pos2];
    newTrack.predictedPath_gap         = zeros(0, 2);
    newTrack.predictedFrames_gap       = zeros(0, 1);
    newTrack.age                       = 2;
    newTrack.lastFrame                 = currentFrame;
    newTrack.consecutiveInvisibleCount = 0;
    if ismember('Intensity', newTrack.localizations.Properties.VariableNames)
        newTrack.meanBrightness = mean(newTrack.localizations.Intensity, 'omitnan');
    else
        newTrack.meanBrightness = NaN;
    end
end

%% ========================================================================
%  ASSIGNMENT HELPERS
% =========================================================================

function assignments = greedyNN_assign(costMatrix, max_cost)
% GREEDYNN_ASSIGN  Greedy nearest-neighbor assignment.
%   Finds the globally best (lowest-cost) valid link, assigns it, then
%   removes both the track and detection from the pool. Repeats until no
%   valid links remain. Faster than Hungarian for sparse cost matrices.

    [numTracks, numDet] = size(costMatrix);
    assignments = zeros(numTracks, 1);
    if numDet == 0, return; end

    availableDet = true(1, numDet);

    valid_links = costMatrix <= max_cost & isfinite(costMatrix);
    [trackIdx_all, detIdx_all] = find(valid_links);
    if isempty(trackIdx_all), return; end

    lin = sub2ind(size(costMatrix), trackIdx_all, detIdx_all);
    costs_all = costMatrix(lin);
    [~, order] = sort(costs_all, 'ascend');

    trackIdx = trackIdx_all(order);
    detIdx   = detIdx_all(order);

    for i = 1:length(trackIdx)
        tI = trackIdx(i);
        dI = detIdx(i);
        if assignments(tI) == 0 && availableDet(dI)
            assignments(tI) = dI;
            availableDet(dI) = false;
        end
    end
end

%% ========================================================================
%  TRACK VALIDATION & FORMATTING
% =========================================================================

function validTracks = validateTracks(tracks, params)
% VALIDATETRACKS  Removes tracks where more than mask_tolerance_fraction of
%   the smoothed path falls outside the combined vessel/ROI mask.
%
%   This tolerance-based approach prevents losing valid tracks that have
%   a few Kalman-corrected points drifting slightly past vessel boundaries,
%   which is common at sharp mask edges.

    if isempty(tracks)
        validTracks = tracks;
        return;
    end

    H = params.expParams.size(1);
    W = params.expParams.size(2);
    combinedMask = true(H, W);
    maskActive   = false;

    if isfield(params, 'proc') && isfield(params.proc, 'vesselMask') && ...
       isfield(params.proc.vesselMask, 'enable') && params.proc.vesselMask.enable && ...
       isfield(params.proc.vesselMask, 'vesselMask') && ~isempty(params.proc.vesselMask.vesselMask)
        combinedMask = combinedMask & logical(params.proc.vesselMask.vesselMask);
        maskActive   = true;
    end

    if isfield(params, 'metadata') && isfield(params.metadata, 'roiMask') && ~isempty(params.metadata.roiMask)
        combinedMask = combinedMask & logical(params.metadata.roiMask);
        maskActive   = true;
    end

    if ~maskActive
        validTracks = tracks;
        return;
    end

    tolerance = params.track.kalman.mask_tolerance_fraction;
    tracks_to_remove = false(1, length(tracks));

    for i = 1:length(tracks)
        currTrack = tracks{i};

        if isfield(currTrack, 'correctedPath') && ~isempty(currTrack.correctedPath)
            pts = currTrack.correctedPath;
        else
            pts = [currTrack.localizations.X, currTrack.localizations.Y];
        end

        x_idx = round(pts(:, 1));
        y_idx = round(pts(:, 2));

        in_bounds   = (x_idx >= 1 & x_idx <= W & y_idx >= 1 & y_idx <= H);
        outside_any = ~in_bounds;

        if any(in_bounds)
            ib = find(in_bounds);
            lin_idx = sub2ind([H, W], y_idx(ib), x_idx(ib));
            outside_any(ib(~combinedMask(lin_idx))) = true;
        end

        outside_frac = sum(outside_any) / length(outside_any);
        if outside_frac > tolerance
            tracks_to_remove(i) = true;
        end
    end

    validTracks = tracks(~tracks_to_remove);

    n_removed = sum(tracks_to_remove);
    if n_removed > 0
        fprintf('      -> Mask validation: Removed %d tracks (>%.0f%% outside mask).\n', ...
            n_removed, 100 * tolerance);
    end
end

function finalTracks = formatTracks(terminatedTracks, params, indent_prefix)
% FORMATTRACKS  Converts raw track structs to the final output format.
%   Computes velocities in mm/s from the Kalman-smoothed path, builds
%   interpolated (gap-filled) paths for density rendering, and optionally
%   delegates to applyQualityControl for additional filtering.
%
%   Tracks shorter than min_track_length_velocity receive zeroed velocity
%   fields and are flagged with has_reliable_velocity = false. They still
%   contribute to density maps but are excluded from velocity maps.

    trackParams              = params.track;
    min_track_length         = trackParams.min_track_length;
    min_track_length_velocity = trackParams.kalman.min_track_length_velocity;

    tempTracksCell = cell(1, length(terminatedTracks));
    track_count    = 0;

    for i = 1:length(terminatedTracks)
        track = terminatedTracks{i};
        if track.age < min_track_length
            continue;
        end

        track_count = track_count + 1;
        outTrack = struct();
        outTrack.id            = track.id;
        outTrack.localizations = track.localizations;

        % Enforce strict consistency between correctedPath and age
        if isfield(track, 'correctedPath') && ~isempty(track.correctedPath)
            assert(size(track.correctedPath, 1) == track.age, ...
                ['Track %d: correctedPath length (%d) != age (%d). ', ...
                 'This indicates a logic error in the tracking loop.'], ...
                track.id, size(track.correctedPath, 1), track.age);
            outTrack.path = track.correctedPath;
        else
            outTrack.path = [track.localizations.X, track.localizations.Y];
        end

        outTrack.frames = track.localizations.Frame;
        outTrack.length = track.age;

        % Build the dense, gap-filled path for density rendering.
        % Gap positions come from Kalman predictions during invisible
        % frames and are stored separately from the corrected path.
        if trackParams.kalman.interpolate_gaps && ...
           isfield(track, 'predictedPath_gap') && ~isempty(track.predictedPath_gap)
            [outTrack.path_interpolated, outTrack.frames_interpolated] = ...
                buildInterpolatedPath(outTrack.path, outTrack.frames, ...
                                      track.predictedPath_gap, track.predictedFrames_gap);
        else
            outTrack.path_interpolated   = outTrack.path;
            outTrack.frames_interpolated = outTrack.frames;
        end

        % --- Velocity computation (mm/s) on the corrected path ---
        if track.age >= min_track_length_velocity && track.age > 1
            dx = diff(outTrack.path(:, 1));
            dy = diff(outTrack.path(:, 2));
            dt_frames = diff(outTrack.frames);

            disp_mm = sqrt((dx * trackParams.pixel_X_size).^2 + (dy * trackParams.pixel_Z_size).^2);
            dt_sec  = dt_frames * trackParams.dt;

            v = zeros(size(dt_sec));
            valid = dt_sec > 0;
            v(valid) = disp_mm(valid) ./ dt_sec(valid);
            outTrack.velocities_mm_s = [v; v(end)];

            mean_v = mean(v(isfinite(v)), 'omitnan');
            if isempty(mean_v) || isnan(mean_v), mean_v = 0; end
            outTrack.average_velocity_mm_s = repmat(mean_v, track.age, 1);
            outTrack.has_reliable_velocity = true;
        else
            % Track is too short for a reliable velocity estimate
            outTrack.velocities_mm_s       = zeros(track.age, 1);
            outTrack.average_velocity_mm_s = zeros(track.age, 1);
            outTrack.has_reliable_velocity = false;
        end

        tempTracksCell{track_count} = outTrack;
    end

    initial_filtered_tracks = [tempTracksCell{1:track_count}];

    if exist('applyQualityControl', 'file')
        finalTracks = applyQualityControl(initial_filtered_tracks, trackParams, indent_prefix);
    else
        finalTracks = initial_filtered_tracks;
    end
end

function [interpPath, interpFrames] = buildInterpolatedPath(corrPath, corrFrames, gapPath, gapFrames)
% BUILDINTERPOLATEDPATH  Combines corrected and gap-predicted positions
%   into a single chronologically sorted dense path. Useful for density-map
%   rendering where gaps would otherwise create holes.

    all_frames = [corrFrames(:); gapFrames(:)];
    all_pts    = [corrPath;      gapPath];
    [interpFrames, order] = sort(all_frames);
    interpPath            = all_pts(order, :);
end

%% ========================================================================
%  SUMMARY OUTPUT
% =========================================================================

function printTrackingSummary(tracks, localizations, indent_prefix)
% PRINTTRACKINGSUMMARY  Prints post-tracking statistics to the console.

    fprintf('%sKalman v2 tracking complete. Found %d valid tracks after QC.\n', ...
        indent_prefix, length(tracks));

    total_initial = height(localizations);
    if ~isempty(tracks)
        total_in_tracks = sum([tracks.length]);
    else
        total_in_tracks = 0;
    end
    untracked = total_initial - total_in_tracks;

    fprintf('%s  - Total localizations in tracks: %d\n', indent_prefix, total_in_tracks);
    if total_initial > 0
        pct = 100 * untracked / total_initial;
        fprintf('%s  - Un-tracked localizations: %d (%.1f%% of total)\n', ...
            indent_prefix, untracked, pct);
    else
        fprintf('%s  - Un-tracked localizations: 0\n', indent_prefix);
    end
end