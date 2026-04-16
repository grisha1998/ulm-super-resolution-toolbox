classdef DisplayManager < handle
    % DISPLAYMANAGER - Handles all visualization and display logic
    %
    % This class separates display concerns from the main GUI logic,
    % making the code more maintainable and testable.
    
    properties (Access = private)
        lastDisplayTime = 0;
        debounceTimer = [];
        trackColormap   = [];   % cached lines(256) — avoids recomputation every frame
    end
    
    methods
        function obj = DisplayManager()
            % Constructor
            obj.lastDisplayTime = 0;
            obj.trackColormap   = lines(256);
        end
        
        function displayFrame(obj, app)
            % Main display function - routes to appropriate display method
            %
            % Args:
            %   app: Application structure
            
            if ~isvalid(app.fig) || isempty(app)
                return;
            end
            
            % Clear axes
            cla(app.ui.ax);
            hold(app.ui.ax, 'on');
            axis(app.ui.ax, 'image');
            
            % Check for ROI preview mode
            if app.state.isROIPreview && ~isempty(app.data.vesselMap)
                obj.displayROIPreview(app);
                return;
            end
            
            % Display based on current state
            switch app.state.currentState
                case 0
                    obj.displayRawData(app);
                case 1
                    obj.displayFilteredData(app);
                case 2
                    obj.displayDetectionResults(app);
                case 3
                    obj.displayLocalizations(app);
                case 4
                    obj.displayRawTracks(app);
                case 5
                    obj.displayProcessedTracks(app);
                otherwise
                    title(app.ui.ax, 'Please load data to begin...');
            end
            
            % Add mask contour if exists and not in ROI preview
            if ~isempty(app.data.mask) && ~app.state.isROIPreview
                obj.addMaskContour(app);
            end
            
            set(app.ui.ax, 'YDir', 'reverse');
            hold(app.ui.ax, 'off');
            drawnow limitrate;
        end
        
        function displayRawData(obj, app)
            % Display raw ultrasound data
            t = app.state.currentFrame;
            
            if ~isempty(app.data.rawData)
                % Process the image using the GUI parameters
                img = obj.processImageForDisplay(app, app.data.rawData(:,:,t));
                imagesc(app.ui.ax, img);
                obj.applyAxesDisplaySettings(app); % Apply CLim & Colormap
            end
            
            if isfield(app.data, 'isIQData') && app.data.isIQData
                dataType = 'IQ Complex';
            else
                dataType = 'Real';
            end
            title(app.ui.ax, sprintf('Raw Data [%s] | Frame %d / %d', ...
                dataType, t, app.state.maxFrame));
        end
        
        function displayFilteredData(obj, app)
            % Display filtered data
            t = app.state.currentFrame;
            
            if ~isempty(app.data.filteredData)
                % Process the image using the GUI parameters
                img = obj.processImageForDisplay(app, app.data.filteredData(:,:,t));
                imagesc(app.ui.ax, img);
                obj.applyAxesDisplaySettings(app); % Apply CLim & Colormap
            end
            
            title(app.ui.ax, sprintf('Filtered Data | Frame %d / %d', ...
                t, app.state.maxFrame));
        end
        
        function displayDetectionResults(obj, app)
            % Display detection results with overlays
            t = app.state.currentFrame;
            count = 0;
            
            % Background image
            if ~isempty(app.data.filteredData)
                img = obj.processImageForDisplay(app, app.data.filteredData(:,:,t));
                imagesc(app.ui.ax, img);
                obj.applyAxesDisplaySettings(app);
            end
            
            % Overlay detections
            if ~isempty(app.data.candidateBubbles) && istable(app.data.candidateBubbles)
                frameData = app.data.candidateBubbles(app.data.candidateBubbles.Frame == t, :);
                if ~isempty(frameData)
                    plot(app.ui.ax, frameData.X, frameData.Y, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
                    count = height(frameData);
                end
            end
            
            title(app.ui.ax, sprintf('Detection: %d Candidates | Frame %d / %d', ...
                count, t, app.state.maxFrame));
        end
        
        function displayLocalizations(obj, app)
            % Display localization results
            t = app.state.currentFrame;
            count = 0;
            
            % Background image
            if ~isempty(app.data.filteredData)
                img = obj.processImageForDisplay(app, app.data.filteredData(:,:,t));
                imagesc(app.ui.ax, img);
                obj.applyAxesDisplaySettings(app);
            end
            
            % Overlay localizations
            if ~isempty(app.data.localizations) && istable(app.data.localizations)
                frameData = app.data.localizations(app.data.localizations.Frame == t, :);
                if ~isempty(frameData)
                    plot(app.ui.ax, frameData.X, frameData.Y, 'r.', 'MarkerSize', 12);
                    count = height(frameData);
                end
            end
            
            title(app.ui.ax, sprintf('Localization: %d Particles | Frame %d / %d', ...
                count, t, app.state.maxFrame));
        end
        
        function displayRawTracks(obj, app)
            % Display raw tracking results
            t = app.state.currentFrame;
            trackCount = 0;
            
            % Background image (mean)
            if ~isempty(app.data.filteredMeanBG)
                img = obj.processImageForDisplay(app, app.data.filteredMeanBG);
                imagesc(app.ui.ax, img);
                obj.applyAxesDisplaySettings(app);
            end
            
            % Overlay tracks
            if ~isempty(app.data.tracks_raw)
                cmap = obj.trackColormap;
                
                for i = 1:length(app.data.tracks_raw)
                    track = app.data.tracks_raw(i);
                    idx = track.frames <= t;
                    
                    if any(idx)
                        color = cmap(mod(track.id, 256) + 1, :);
                        plot(app.ui.ax, track.path(idx,1), track.path(idx,2), ...
                            'Color', color, 'LineWidth', 1.5);
                        trackCount = trackCount + 1;
                    end
                end
            end
            
            title(app.ui.ax, sprintf('Tracking (Raw): %d Active Tracks | Frame %d / %d', ...
                trackCount, t, app.state.maxFrame));
        end
        
        function displayProcessedTracks(obj, app)
            % Display processed/smoothed tracks
            
            % Background image
            if ~isempty(app.data.filteredMeanBG)
                img = obj.processImageForDisplay(app, app.data.filteredMeanBG);
                imagesc(app.ui.ax, img);
                obj.applyAxesDisplaySettings(app);
            end
            
            % Overlay tracks (filtered by length)
            trackCount = 0;
            if ~isempty(app.data.tracks_final)
                min_len = app.ui.DisplayMinLengthField.Value;
                tracks = app.data.tracks_final([app.data.tracks_final.original_length] >= min_len);
                trackCount = length(tracks);
                
                cmap = obj.trackColormap;
                for i = 1:length(tracks)
                    color = cmap(mod(tracks(i).id, 256) + 1, :);
                    plot(app.ui.ax, tracks(i).path(:,1), tracks(i).path(:,2), ...
                        'Color', color, 'LineWidth', 1.5);
                end
            end
            
            title(app.ui.ax, sprintf('Final Processing: %d Total Tracks (Len>=%d)', ...
                trackCount, app.ui.DisplayMinLengthField.Value));
        end
        
        function displayROIPreview(obj, app)
            % Display ROI selection preview with vessel map and mask overlay
            
            % 1. Display vessel map
            imagesc(app.ui.ax, app.data.vesselMap);
            obj.applyAxesDisplaySettings(app); % Apply CLim & Colormap only
            
            % 2. Red semi-transparent overlay for mask
            if ~isempty(app.data.mask)
                % Pre-compute red overlay template if needed
                if isempty(app.data.redOverlayTemplate) || ...
                   ~isequal(size(app.data.redOverlayTemplate(:,:,1)), size(app.data.mask))
                    app.data.redOverlayTemplate = cat(3, ...
                        ones(size(app.data.mask)), ...
                        zeros(size(app.data.mask)), ...
                        zeros(size(app.data.mask)));
                end
                
                hOv = image(app.ui.ax, app.data.redOverlayTemplate);
                set(hOv, 'AlphaData', double(app.data.mask) * 0.3);
            end
            
            title(app.ui.ax, 'ROI Selection: Average Intensity (Vessel Map)');
        end
        
        function addMaskContour(obj, app)
            % Add mask contour to current display
            if ~isempty(app.data.mask)
                contour(app.ui.ax, double(app.data.mask), [0.5 0.5], ...
                    'r--', 'LineWidth', 1);
            end
        end
        
        function debouncedDisplay(obj, app, delay)
            % Debounced display update to prevent excessive redraws
            %
            % Args:
            %   app: Application structure
            %   delay: Debounce delay in seconds (default: from constants)
            
            if nargin < 3
                delay = ULM_Constants.DEBOUNCE_DELAY;
            end
            
            % Stop existing timer
            if ~isempty(obj.debounceTimer) && isvalid(obj.debounceTimer)
                stop(obj.debounceTimer);
                delete(obj.debounceTimer);
            end
            
            % Create new timer
            obj.debounceTimer = timer('StartDelay', delay, ...
                'TimerFcn', @(~,~) obj.displayFrame(app));
            start(obj.debounceTimer);
        end

        function img = processImageForDisplay(obj, app, rawImg)
            % Applies Gamma, Log, and Mat2Gray based on UI settings
            img = abs(rawImg);
            
            % 1. Log compression
            if app.ui.disp_log.Value
                img = 20 * log10(img + max(img(:))*1e-6); 
            end
            
            % 2. Mat2Gray (Normalization)
            if app.ui.disp_mat2gray.Value
                img = mat2gray(img);
            end
            
            % 3. Gamma stretch
            gamma = app.ui.disp_gamma.Value;
            if gamma ~= 1.0
                if ~app.ui.disp_mat2gray.Value
                    % Safely apply gamma if not normalized 0-1
                    img_min = min(img(:)); img_max = max(img(:));
                    if img_max > img_min
                        img_norm = (img - img_min) / (img_max - img_min);
                        img = (img_norm .^ gamma) * (img_max - img_min) + img_min;
                    end
                else
                    img = img .^ gamma;
                end
            end
        end
        
        function applyAxesDisplaySettings(obj, app)
            % Applies colormap and CLim limits
            colormap(app.ui.ax, app.ui.disp_cmap.Value);
            
            if ~app.ui.disp_clim_auto.Value
                cmin = app.ui.disp_clim_min.Value;
                cmax = app.ui.disp_clim_max.Value;
                if cmax > cmin
                    clim(app.ui.ax, [cmin, cmax]);
                end
            else
                clim(app.ui.ax, 'auto');
            end
        end
        
        function cleanup(obj)
            % Clean up resources
            if ~isempty(obj.debounceTimer) && isvalid(obj.debounceTimer)
                stop(obj.debounceTimer);
                delete(obj.debounceTimer);
            end
        end
    end
    
    methods (Static)
        function [cLow, cHigh] = calculateColorLimits(data, pLow, pHigh)
            % Calculate robust color limits using percentiles
            %
            % Args:
            %   data: Image data
            %   pLow: Lower percentile (default: 1)
            %   pHigh: Upper percentile (default: 99.9)
            %
            % Returns:
            %   cLow, cHigh: Color limit values
            
            if nargin < 2
                pLow = ULM_Constants.DEFAULT_PERCENTILE_LOW;
            end
            if nargin < 3
                pHigh = ULM_Constants.DEFAULT_PERCENTILE_HIGH;
            end
            
            absData = abs(data(:));
            cLow = prctile(absData, pLow);
            cHigh = prctile(absData, pHigh);
            
            % Prevent identical limits
            if cLow == cHigh
                cHigh = cLow + 1;
            end
        end
        
        function overlayImage = createMaskOverlay(mask, color, alpha)
            % Create colored semi-transparent overlay for mask
            %
            % Args:
            %   mask: Binary mask
            %   color: RGB color vector [R G B] (default: red)
            %   alpha: Transparency (default: 0.3)
            %
            % Returns:
            %   overlayImage: RGB image with alpha channel
            
            if nargin < 2
                color = [1 0 0]; % Red
            end
            if nargin < 3
                alpha = 0.3;
            end
            
            overlayImage = cat(3, ...
                ones(size(mask)) * color(1), ...
                ones(size(mask)) * color(2), ...
                ones(size(mask)) * color(3));
        end
    end
end
