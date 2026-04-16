function printTrackMetrics(tracks, available_locs_count, indent_prefix)
% =========================================================================
% FUNCTION: printTrackMetrics
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   A lightweight diagnostic utility that prints immediate statistical feedback 
%   to the console regarding the quality and composition of tracked trajectories.
%   - Advantages: Provides developers and users with instant visibility into 
%     how efficiently localizations are being grouped into tracks (Localization 
%     Usage Percentage) and the distribution of track lengths, without needing 
%     to wait for the heavy GUI analysis stage.
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Safety Check: Silently returns if the track array is empty.
%   2. Metric Calculation: Extracts all lengths into an array. Calculates 
%      mean, median, and max length.
%   3. Efficiency Metric: Computes the percentage of total raw localizations 
%      that were successfully bundled into tracks vs discarded as noise.
%   4. Text Histogram: Uses `histcounts` with predefined bins (e.g., 3-4, 5-9) 
%      to generate and format a clean, inline text representation of the distribution.
%
% SYNTAX OPTIONS:
%   printTrackMetrics(tracks, available_locs_count, indent_prefix)
%
% EXAMPLES:
%   % Print metrics after a tracking pass:
%   printTrackMetrics(finalTracks, height(rawLocalizationsTable), '    -> ');
%
% INPUTS:
%   tracks               - (Type: Struct Array) Output from the tracking algorithms.
%   available_locs_count - (Type: Integer) The total number of points fed into 
%                          the tracker (used for efficiency percentage).
%   indent_prefix        - (Type: String) Spaces for console alignment.
%
% OUTPUTS:
%   None. (Outputs directly to the MATLAB console).
%
% AUTHOR: Grigori Shapiro
% =========================================================================

    % Prints a summary of track quality metrics for a given set of tracks.
    if isempty(tracks)
        fprintf('%s   -> Metrics: No valid tracks found at this level.\n', indent_prefix);
        return;
    end
    
    lengths = [tracks.length];
    total_locs_in_tracks = sum(lengths);
    
    % --- 1. Basic Stats ---
    mean_len = mean(lengths);
    median_len = median(lengths);
    max_len = max(lengths);
    loc_usage_percent = 100 * total_locs_in_tracks / available_locs_count;
    
    fprintf('%s   -> Metrics: Mean Length: %.1f | Median: %d | Max: %d | Locs Used: %.1f%%\n', ...
            indent_prefix, mean_len, median_len, max_len, loc_usage_percent);
            
    % --- 2. Text-based Histogram ---
    bins = [3, 5, 10, 20, 50, Inf]; % Define the edges of the length bins
    counts = histcounts(lengths, bins);
    
    fprintf('%s      Length Dist: ', indent_prefix);
    for i = 1:length(counts)
        if i < length(bins) - 1
            label = sprintf('[%d-%d]: %d', bins(i), bins(i+1)-1, counts(i));
        else
            label = sprintf('[%d+]: %d', bins(i), counts(i));
        end
        
        if i < length(counts)
            fprintf('%s, ', label);
        else
            fprintf('%s\n', label); % Newline at the end
        end
    end
end