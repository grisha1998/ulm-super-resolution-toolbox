function tracks = trackHungarian(localizations, params, indent_prefix)
% =========================================================================
% FUNCTION: trackHungarian
% =========================================================================
%
% PURPOSE:
%   Performs global optimal microbubble tracking using the Hungarian algorithm 
%   (Munkres assignment). This function evaluates the entire cost landscape 
%   of potential links between existing tracks and new localizations in each 
%   frame, ensuring that the final set of assignments minimizes the total 
%   tracking cost, rather than settling for local, greedy nearest-neighbor matches.
%
% DETAILED ALGORITHM WORKFLOW:
%   1. Data Preparation: Localizations are pre-grouped by 'Frame' using 
%      `splitapply` and `findgroups` for high-performance sequential access.
%   2. Cost Matrix Construction (`calculateCostMatrix`): 
%      For every frame, an [M x N] cost matrix is built where M = active tracks 
%      and N = new localizations.
%      - Standard Mode: Cost is purely squared spatial distance.
%      - Advanced Mode: Cost integrates distance, directional penalty (angle 
%        deviation from historical velocity), and intensity similarity.
%   3. Global Assignment (`munkres`): The Hungarian algorithm solves the linear 
%      assignment problem on the cost matrix.
%   4. Track Update: 
%      - Assigned tracks append the new localization and reset their gap counter.
%      - Unassigned tracks increment their gap counter (`gapFrames`). If 
%        `gapFrames > max_gap_closing_frames`, the track is moved to `terminatedTracks`.
%   5. Track Initiation: Any new localization not assigned to an active track 
%      immediately spawns a brand new track (Single-step initiation).
%   6. Finalization & Kinematics: All tracks are finalized. The function calculates 
%      frame-to-frame displacement using `pixel_X_size` and `pixel_Z_size`, divides 
%      by `dt` to yield `velocities_mm_s`, and filters out any track shorter than 
%      `min_track_length`.
%
% INPUT PARAMETERS & VARIABLES:
%   localizations - (Table) Must contain:
%                   'Frame': Integer frame index.
%                   'X', 'Y': Sub-pixel spatial coordinates.
%                   'Intensity': Amplitude/Brightness of the localization.
%   params        - (Struct) Core tracking configuration:
%                   .track.max_linking_distance (double): Maximum allowed spatial 
%                         jump between frames in pixels. Forms the gating threshold.
%                   .track.max_gap_closing_frames (integer): Number of consecutive 
%                         frames a track can miss a detection before dying.
%                   .track.min_track_length (integer): Tracks with fewer points 
%                         than this are deleted at the end.
%                   .track.pixel_X_size / .pixel_Z_size (double): Physical pixel 
%                         dimensions in mm, used for velocity calculation.
%                   .track.dt (double): Time delta between frames in seconds.
%                   .track.kalman... : Sub-struct containing weights for the 
%                         advanced cost matrix (direction/brightness penalties).
%
% INTERNAL TRACK STATE VARIABLES (Inside activeTracks cell array):
%   .id            - (Integer) Unique monotonically increasing track identifier.
%   .localizations - (Table) Accumulated raw localizations for this track.
%   .lastFrame     - (Integer) The index of the last frame this track was seen.
%   .gapFrames     - (Integer) Current counter of consecutive missing frames.
%
% FINAL OUTPUT STRUCTURE (tracks array):
%   .id                    - Unique Track ID.
%   .localizations         - The complete table of raw localizations.
%   .path                  - [Nx2 double] array of the [X, Y] coordinates.
%   .frames                - [Nx1 double] array of the frame indices.
%   .length                - (Integer) Total number of points in the track.
%   .velocities_mm_s       - [Nx1 double] Instantaneous velocity vector.
%   .average_velocity_mm_s - [Nx1 double] Robust mean velocity (replicated for 
%                            matrix compatibility during rendering).
%
% DEPENDENCIES:
%   - munkres.m (Hungarian algorithm)
%   - calculateCostMatrix.m
%   - applyQualityControl.m
%
% AUTHOR: Grigori Shapiro
% =========================================================================

    % --- 1. Initialization ---
    if nargin < 3, indent_prefix = ''; end
    if isempty(localizations)
        tracks = [];
        fprintf('%sNo localizations found to track. Returning empty.\n', indent_prefix);
        return;
    end
    
    trackParams = params.track;
    % Map user-defined parameter names for clarity
    trackParams.direction_penalty_weight = params.track.kalman.direction_penalty_weight;
    trackParams.brightness_penalty_weight = params.track.kalman.brightness_penalty_weight;
    trackParams.max_angle_gate_deg = params.track.kalman.max_angle_change_deg;

    locsByFrame = splitapply(@(varargin) {table(varargin{:}, 'VariableNames', localizations.Properties.VariableNames)},...
        localizations, findgroups(localizations.Frame));
    
    numFrames = length(locsByFrame);
    activeTracks = {};
    terminatedTracks = {};
    nextTrackID = 1;
    
    fprintf('%sStarting Hungarian tracking...\n', indent_prefix);
    
    % --- 2. Main Tracking Loop ---
    for t = 1:numFrames
        currentLocs = locsByFrame{t};
        
        if isempty(activeTracks) || isempty(currentLocs)
            if ~isempty(activeTracks)
                for i = length(activeTracks):-1:1
                    if (t - activeTracks{i}.lastFrame) > trackParams.max_gap_closing_frames
                        terminatedTracks{end+1} = activeTracks{i};
                        activeTracks(i) = [];
                    end
                end
            end
            for i = 1:height(currentLocs)
                newTrack.id = nextTrackID; nextTrackID = nextTrackID + 1;
                newTrack.localizations = currentLocs(i,:);
                newTrack.lastFrame = t;
                newTrack.gapFrames = 0;
                activeTracks{end+1} = newTrack;
            end
            continue;
        end
        
        % --- 3. Build Cost Matrix ---
        % The cost matrix calculation is now centralized.
        % The function will use the 'Hungarian' method logic internally.
        costMatrix = calculateCostMatrix(activeTracks, currentLocs, params);
        
        if exist('munkres', 'file')
            [assignments, ~] = munkres(costMatrix);
        else
            error('Hungarian algorithm implementation (e.g., ''munkres'') not found on MATLAB path.'); 
        end
        
        % --- 4 & 5. Update Tracks Based on Assignments ---
        unassignedDetections = true(height(currentLocs), 1);
        assignedTracks = false(length(activeTracks), 1);
        
        for i = 1:length(activeTracks)
            if assignments(i) > 0
                detectionIdx = assignments(i);
                activeTracks{i}.localizations = [activeTracks{i}.localizations; currentLocs(detectionIdx,:)];
                activeTracks{i}.lastFrame = t;
                activeTracks{i}.gapFrames = 0;
                unassignedDetections(detectionIdx) = false;
                assignedTracks(i) = true;
            end
        end
        
        for i = length(activeTracks):-1:1
            if ~assignedTracks(i)
                activeTracks{i}.gapFrames = activeTracks{i}.gapFrames + 1;
                if activeTracks{i}.gapFrames > trackParams.max_gap_closing_frames
                    terminatedTracks{end+1} = activeTracks{i};
                    activeTracks(i) = [];
                end
            end
        end
        
        newDetectionIndices = find(unassignedDetections);
        for i = 1:length(newDetectionIndices)
            idx = newDetectionIndices(i);
            newTrack.id = nextTrackID; nextTrackID = nextTrackID + 1;
            newTrack.localizations = currentLocs(idx,:);
            newTrack.lastFrame = t;
            newTrack.gapFrames = 0;
            activeTracks{end+1} = newTrack;
        end
    end
    
    % --- 6. Finalization ---
    terminatedTracks = [terminatedTracks, activeTracks];
    tracks = formatTracks(terminatedTracks, params, indent_prefix);

    fprintf('%sHungarian tracking complete. Found %d valid tracks after QC.\n', indent_prefix, length(tracks));
    
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

% Local helper function to format the final tracks
function finalTracks = formatTracks(terminatedTracks, params, indent_prefix)
    if isempty(terminatedTracks)
        finalTracks = struct('id', {}, 'localizations', {}, 'path', {}, 'frames', {}, ...
                             'velocities_mm_s', {}, 'average_velocity_mm_s', {}, 'length', {});
        return;
    end
    trackParams = params.track;
    min_track_length = trackParams.min_track_length;
    
    % Use a cell array for efficient pre-allocation
    tempTracksCell = cell(1, length(terminatedTracks));
    track_count = 0;
    
    for i = 1:length(terminatedTracks)
        track = terminatedTracks{i};
        track_len = height(track.localizations); % Get length
        
        if track_len >= min_track_length
            track_count = track_count + 1;
            
            % Create a new struct for this track
            outTrack = struct();
            outTrack.id = track.id;
            outTrack.localizations = track.localizations;
            outTrack.path = [track.localizations.X, track.localizations.Y];
            outTrack.frames = track.localizations.Frame;
            outTrack.length = track_len; % Use the calculated length

            if track_len > 1
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
                
                outTrack.average_velocity_mm_s = repmat(mean_vel, outTrack.length, 1);
            else
                outTrack.velocities_mm_s = 0;
                outTrack.average_velocity_mm_s = 0;
            end
            
            % Store the struct in the cell array
            tempTracksCell{track_count} = outTrack;
        end
    end
    
    % Convert the cell array to a struct array in one go (very fast)
    initial_filtered_tracks = [tempTracksCell{1:track_count}];

    % Call the final QC function (assuming applyQualityControl.m is in your path)
    finalTracks = applyQualityControl(initial_filtered_tracks, trackParams, indent_prefix);
end

function D_sq = calculate_dist_sq(A, B)
    % Calculates the squared Euclidean distance between each row in A and each row in B.
    % A is Mx2, B is Nx2. Output D_sq is MxN.
    % This is a toolbox-independent replacement for pdist2(A, B, 'squaredeuclidean').
    
    m = size(A, 1);
    n = size(B, 1);
    D_sq = zeros(m, n);
    
    for i = 1:m
        % Subtract the i-th row of A from all rows of B
        diff_coords = B - repmat(A(i,:), n, 1);
        % Calculate squared distance and store as the i-th row of the output
        D_sq(i,:) = sum(diff_coords.^2, 2)';
    end
end