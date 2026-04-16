function tracks = trackNearestNeighbor(localizations, params, indent_prefix)
% =========================================================================
% FUNCTION: trackNearestNeighbor
% =========================================================================
%
% PURPOSE:
%   Executes tracking using a greedy "Nearest Neighbor" (NN) heuristic. 
%   This is the fastest and computationally lightest tracking method. It relies 
%   entirely on spatial proximity and makes local, best-first decisions without 
%   looking at the global cost matrix or utilizing predictive kinematics.
%
% DETAILED ALGORITHM WORKFLOW:
%   1. Setup: Localizations are grouped by frame. Maximum linking distance is 
%      squared (`max_link_dist_sq`) to avoid calculating costly square roots 
%      during the distance matrix computation.
%   2. Cost Matrix (Distance Map): For the current frame, an [M x N] matrix is 
%      computed using `calculate_dist_sq`, representing the Euclidean squared 
%      distance between the `lastPos` of all M active tracks and all N new localizations.
%   3. Greedy Sorting (`greedyNN_assign`): 
%      - All valid links (where cost <= `max_link_dist_sq`) are identified.
%      - These links are extracted into a linear list and strictly sorted from 
%        lowest distance (best) to highest distance (worst).
%   4. Assignment: The algorithm iterates down the sorted list. If the track 
%      in the pairing is currently unassigned AND the localization in the pairing 
%      is currently unassigned, the link is finalized. 
%      *Note:* Because it is greedy, a slightly sub-optimal early link might prevent 
%      a globally better link later down the list.
%   5. Track Update & Management: 
%      - Assigned tracks update their `.lastPos` and append data.
%      - Unassigned tracks increment `.gapFrames`.
%      - Unassigned localizations instantly spawn new tracks (no two-step logic).
%   6. Finalization: Tracks are formatted and velocities calculated similarly to 
%      the other algorithms.
%
% INPUT PARAMETERS & VARIABLES:
%   localizations - (Table) 'Frame', 'X', 'Y' columns.
%   params        - (Struct) Basic tracking configuration:
%                   .track.max_linking_distance (double): Maximum radius.
%                   .track.max_gap_closing_frames (int): Missing frame tolerance.
%                   .track.min_track_length (int): Final QC length threshold.
%
% INTERNAL TRACK STATE VARIABLES:
%   .lastPos   - [1x2 double] The [X,Y] coordinate of the track in its most 
%                recently detected frame. Used purely as the anchor point for 
%                the next frame's distance search.
%   .gapFrames - (Integer) Counts how long the track has been "blind".
%
% FINAL OUTPUT STRUCTURE (tracks array):
%   Standard format. Note that because there is no smoothing filter (Kalman), 
%   `.path` is a direct copy of the raw `.localizations` coordinates, which 
%   preserves all original sub-pixel localization jitter.
%
% DEPENDENCIES:
%   - None (Uses custom `calculate_dist_sq` to avoid dependency on Stats toolbox).
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
    max_link_dist_sq = trackParams.max_linking_distance^2;
    
    % Group localizations by frame
    locsByFrame = splitapply(@(varargin) {table(varargin{:}, 'VariableNames', localizations.Properties.VariableNames)},...
        localizations, findgroups(localizations.Frame));
    
    numFrames = length(locsByFrame);
    activeTracks = {};
    terminatedTracks = {};
    nextTrackID = 1;
    
    fprintf('%sStarting Nearest Neighbor (NN) tracking...\n', indent_prefix);
    
    % --- 2. Main Tracking Loop ---
    for t = 1:numFrames
        currentLocs = locsByFrame{t};
        
        if isempty(activeTracks) || isempty(currentLocs)
            % --- No active tracks OR no new detections ---
            % Handle gaps for existing tracks
            if ~isempty(activeTracks)
                for i = length(activeTracks):-1:1
                    activeTracks{i}.gapFrames = activeTracks{i}.gapFrames + 1;
                    if activeTracks{i}.gapFrames > trackParams.max_gap_closing_frames
                        terminatedTracks{end+1} = activeTracks{i};
                        activeTracks(i) = [];
                    end
                end
            end
            % Start new tracks for any new detections
            for i = 1:height(currentLocs)
                newTrack.id = nextTrackID; nextTrackID = nextTrackID + 1;
                newTrack.localizations = currentLocs(i,:);
                newTrack.lastFrame = t;
                newTrack.gapFrames = 0;
                newTrack.lastPos = [currentLocs.X(i), currentLocs.Y(i)];
                activeTracks{end+1} = newTrack;
            end
            continue; % Move to the next frame
        end
        
        % --- 3. Build Cost Matrix (Simple Euclidean Distance) ---
        numActiveTracks = length(activeTracks);
        numDetections = height(currentLocs);
        
        % Get last positions of all active tracks
        allPos = cellfun(@(t) t.lastPos, activeTracks, 'UniformOutput', false);
        trackLastPos = vertcat(allPos{:});
        % Get positions of all new detections
        detectionPos = [currentLocs.X, currentLocs.Y];
        
        % Calculate all-to-all squared distances
        costMatrix = calculate_dist_sq(trackLastPos, detectionPos);
        
        % --- 4. Assignment Step (Greedy Nearest Neighbor) ---
        assignments = greedyNN_assign(costMatrix, max_link_dist_sq);
        
        % --- 5. Update Tracks Based on Assignments ---
        unassignedDetections = true(numDetections, 1);
        assignedTracks = false(numActiveTracks, 1);
        
        for i = 1:numActiveTracks
            if assignments(i) > 0
                detectionIdx = assignments(i);
                
                % Update track
                activeTracks{i}.localizations = [activeTracks{i}.localizations; currentLocs(detectionIdx,:)];
                activeTracks{i}.lastFrame = t;
                activeTracks{i}.gapFrames = 0;
                activeTracks{i}.lastPos = [currentLocs.X(detectionIdx), currentLocs.Y(detectionIdx)];
                
                % Mark as assigned
                unassignedDetections(detectionIdx) = false;
                assignedTracks(i) = true;
            end
        end
        
        % --- 6. Track Management ---
        % Handle unassigned tracks (gap or terminate)
        for i = numActiveTracks:-1:1
            if ~assignedTracks(i)
                activeTracks{i}.gapFrames = activeTracks{i}.gapFrames + 1;
                if activeTracks{i}.gapFrames > trackParams.max_gap_closing_frames
                    terminatedTracks{end+1} = activeTracks{i};
                    activeTracks(i) = [];
                end
            end
        end
        
        % Start new tracks for unassigned detections
        newDetectionIndices = find(unassignedDetections);
        for i = 1:length(newDetectionIndices)
            idx = newDetectionIndices(i);
            newTrack.id = nextTrackID; nextTrackID = nextTrackID + 1;
            newTrack.localizations = currentLocs(idx,:);
            newTrack.lastFrame = t;
            newTrack.gapFrames = 0;
            newTrack.lastPos = [currentLocs.X(idx), currentLocs.Y(idx)];
            activeTracks{end+1} = newTrack;
        end
    end
    
    % --- 7. Finalization ---
    terminatedTracks = [terminatedTracks, activeTracks];
    % Format tracks (same helper function as other trackers)
    tracks = formatTracks(terminatedTracks, params, indent_prefix);

    fprintf('%sNN tracking complete. Found %d valid tracks after QC.\n', indent_prefix, length(tracks));
end

%% ========================================================================
%  Local Helper Functions
% =========================================================================

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

function D_sq = calculate_dist_sq(A, B)
    % Calculates the squared Euclidean distance between each row in A and each row in B.
    m = size(A, 1);
    n = size(B, 1);
    D_sq = zeros(m, n);
    for i = 1:m
        diff_coords = B - repmat(A(i,:), n, 1);
        D_sq(i,:) = sum(diff_coords.^2, 2)';
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
