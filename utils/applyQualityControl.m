function finalTracks = applyQualityControl(initialTracks, trackParams, indent_prefix)
% =========================================================================
% FUNCTION: applyQualityControl
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Acts as the centralized orchestrator for the post-tracking Quality Control 
%   (QC) pipeline. It applies a standardized sequence of physical constraints 
%   to raw tracks.
%   - Advantages: Enforces a strict, reproducible filtering pipeline. It 
%     provides extensive, formatted console logging and text-based histograms, 
%     giving the user immediate visual feedback on how aggressively each filter 
%     (Direction, Acceleration, VD) is cutting the data.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Initial Log: Outputs the baseline statistics for tracks that passed 
%      the minimum length constraint during tracking.
%   2. Sequential Filtering:
%      a. Calls `applyDirectionConstraint` to remove sharp turns.
%      b. Calls `applyAccelerationConstraint` to remove impossible speed jumps.
%      c. Calls `applyVDConstraint` to remove highly tortuous, jittery tracks.
%   3. Logging: After each specific filter, it calls the local helper 
%      `displayTrackHistograms` to print the updated Mean/Median/Max of length 
%      and velocity, alongside an ASCII art histogram to visualize distributions.
%
% SYNTAX OPTIONS:
%   finalTracks = applyQualityControl(initialTracks, trackParams, indent_prefix)
%
% EXAMPLES:
%   % Finalize tracks inside a tracking loop with nested log indentation:
%   clean_tracks = applyQualityControl(rawTracks, params.track, '      -> ');
%
% INPUTS:
%   initialTracks - (Type: Struct Array) The unrefined output from a tracker.
%   trackParams   - (Type: Struct) The sub-struct `params.track` containing 
%                   all specific `.qc` settings.
%   indent_prefix - (Type: String) Formatting prefix (e.g., '  ') for clean logs.
%
% OUTPUTS:
%   finalTracks   - (Type: Struct Array) The finalized, physically valid tracks.
%
% AUTHOR: Grigori Shapiro
% =========================================================================
    
    % If no tracks were passed in, or if QC is disabled, return immediately.
    if isempty(initialTracks) || ~isfield(trackParams, 'qc')
        if isempty(initialTracks) && isfield(trackParams, 'min_track_length')
             fprintf('%s       - Initial Filtering: Kept 0 tracks (length >= %d).\n', indent_prefix, trackParams.track.min_track_length);
        end
        finalTracks = initialTracks;
        return;
    end
    
    fprintf('%s   - Applying post-tracking Quality Control filters...\n', indent_prefix);
    
    % --- DETAILED QC LOGGING BLOCK STARTS HERE ---
    min_track_length = trackParams.min_track_length;
    n_before_qc = length(initialTracks);
    fprintf('%s       - Initial Filtering: Kept %d tracks (length >= %d).\n', indent_prefix, n_before_qc, min_track_length);
    displayTrackHistograms(initialTracks, 'Initial', indent_prefix);
    
    % Apply Direction Constraint
    tracks_after_dir = applyDirectionConstraint(initialTracks, trackParams);
    n_after_dir = length(tracks_after_dir);
    fprintf('%s       - QC Direction: Kept %d of %d tracks (%.1f%% pass).\n', indent_prefix, n_after_dir, n_before_qc, 100*n_after_dir/max(1, n_before_qc));
    displayTrackHistograms(tracks_after_dir, 'After Dir. QC', indent_prefix);
    
    % Apply Acceleration Constraint
    tracks_after_accel = applyAccelerationConstraint(tracks_after_dir, trackParams);
    n_after_accel = length(tracks_after_accel);
    fprintf('%s       - QC Acceleration: Kept %d of %d tracks (%.1f%% pass).\n', indent_prefix, n_after_accel, n_after_dir, 100*n_after_accel/max(1, n_after_dir));
    displayTrackHistograms(tracks_after_accel, 'After Accel. QC', indent_prefix);
    
    % Apply Velocity Dispersion Constraint
    tracks_after_vd = applyVDConstraint(tracks_after_accel, trackParams);
    n_after_vd = length(tracks_after_vd);
    fprintf('%s       - QC Velocity Dispersion: Kept %d of %d tracks (%.1f%% pass).\n', indent_prefix, n_after_vd, n_after_accel, 100*n_after_vd/max(1, n_after_accel));
    displayTrackHistograms(tracks_after_vd, 'Final', indent_prefix);
    
    finalTracks = tracks_after_vd;
    % --- DETAILED QC LOGGING BLOCK ENDS HERE ---
end


%% ========================================================================
%  Local Helper Function: Display Track Statistics and Histograms
% =========================================================================
function displayTrackHistograms(tracks, stageTitle, indent_prefix)
    if isempty(tracks)
        fprintf('%s         > No tracks to analyze for stage: %s\n', indent_prefix, stageTitle);
        return;
    end
    lengths = [tracks.length];
    avg_velocities = arrayfun(@(t) t.average_velocity_mm_s(1), tracks);
    
    mean_len = mean(lengths); median_len = median(lengths); max_len = max(lengths);
    mean_vel = mean(avg_velocities); median_vel = median(avg_velocities); max_vel = max(avg_velocities);

    fprintf('%s         > Stats (%s): Len(μ=%.1f, M=%d, max=%d), Vel(μ=%.1f, M=%.1f, max=%.1f mm/s)\n', ...
            indent_prefix, stageTitle, mean_len, median_len, max_len, mean_vel, median_vel, max_vel);
            
    len_bins = [3, 6, 11, 21, inf]; len_counts = histcounts(lengths, len_bins);
    len_labels = {'3-5', '6-10', '11-20', '21+'};
    print_text_histogram('Len Hist', len_labels, len_counts, indent_prefix);
    
    vel_bins = [0, 5, 10, 20, 40, inf]; vel_counts = histcounts(avg_velocities, vel_bins);
    vel_labels = {'0-5', '5-10', '10-20', '20-40', '40+'};
    print_text_histogram('Vel Hist', vel_labels, vel_counts, indent_prefix);
end

function print_text_histogram(title, labels, counts, indent_prefix)
    max_bar_width = 40; total_counts = sum(counts);
    if total_counts == 0, return; end
    
    str = sprintf('%s         > %s: ', indent_prefix, title);
    
    for i = 1:length(labels)
        if counts(i) > 0
            bar_width = round((counts(i) / max(counts)) * max_bar_width);
            bar = repmat('#', 1, bar_width);
            str = [str, sprintf('[%s]: %s (%d)  ', labels{i}, bar, counts(i))];
        end
    end
    fprintf('%s\n', str);
end