classdef SessionManager < handle
    % SESSIONMANAGER - Handles saving and loading complete GUI sessions
    %
    % This class manages serialization of the entire application state
    % including data, parameters, and UI state for later restoration.
    
    properties
        version = '3.0';
        lastSaveFile = '';
    end
    
    methods
        function sessionData = createSessionData(obj, app)
            % Create a complete snapshot of current session
            
            fprintf('Creating session snapshot...\n');
            
            sessionData = struct();
            sessionData.version = obj.version;
            sessionData.timestamp = datetime('now');
            
            % === PARAMETERS ===
            sessionData.params = app.data.params;
            
            % === DATA (selectively save based on size) ===
            sessionData.data = struct();
            
            % Always save these if they exist
            if ~isempty(app.data.rawData)
                sessionData.data.rawDataSize = size(app.data.rawData);
                sessionData.data.rawDataHash = app.data.rawDataHash;
                % Save raw data (this can be large)
                sessionData.data.rawData = app.data.rawData;
            end
            
            if ~isempty(app.data.mask)
                sessionData.data.mask = app.data.mask;
            end
            
            % Save processing results (tables are efficient)
            if ~isempty(app.data.candidateBubbles)
                sessionData.data.candidateBubbles = app.data.candidateBubbles;
            end
            
            if ~isempty(app.data.localizations)
                sessionData.data.localizations = app.data.localizations;
            end
            
            if ~isempty(app.data.tracks_raw)
                sessionData.data.tracks_raw = app.data.tracks_raw;
            end
            
            if ~isempty(app.data.tracks_final)
                sessionData.data.tracks_final = app.data.tracks_final;
            end
            
            % Optionally save filtered data (can be regenerated)
            if ~isempty(app.data.filteredData)
                % Store dimensions and hash, but allow option to not save full data
                sessionData.data.filteredDataSize = size(app.data.filteredData);
                % Uncomment to save full filtered data:
                % sessionData.data.filteredData = app.data.filteredData;
            end
            
            % SVD components (can be regenerated but save for speed)
            if ~isempty(app.data.U)
                sessionData.data.hasSVD = true;
                sessionData.data.svdDims = app.data.svdDims;
                % Save SVD components (these can be large)
                % sessionData.data.U = app.data.U;
                % sessionData.data.S_diag = app.data.S_diag;
                % sessionData.data.V = app.data.V;
            end
            
            % DCC indices
            if ~isempty(app.data.tissue_indices)
                sessionData.data.tissue_indices = app.data.tissue_indices;
                sessionData.data.blood_indices = app.data.blood_indices;
                sessionData.data.noise_indices = app.data.noise_indices;
            end
            
            % === STATE ===
            sessionData.state = struct();
            sessionData.state.currentState = app.state.currentState;
            sessionData.state.currentFrame = app.state.currentFrame;
            sessionData.state.maxFrame = app.state.maxFrame;
            
            % === UI STATE ===
            sessionData.ui = struct();
            sessionData.ui.selectedTab = app.ui.tabGroup.SelectedTab.Title;
            sessionData.ui.frameSliderLimits = app.ui.FrameSlider.Limits;
            
            fprintf('Session snapshot created successfully.\n');
            fprintf('  - State level: %d\n', sessionData.state.currentState);
            fprintf('  - Total frames: %d\n', sessionData.state.maxFrame);
            if isfield(sessionData.data, 'tracks_final')
                fprintf('  - Final tracks: %d\n', length(sessionData.data.tracks_final));
            end
        end
        
        function app = restoreSessionData(obj, app, sessionData)
            % Restore application state from session data
            
            fprintf('Restoring session...\n');
            
            % Check version compatibility
            if ~isfield(sessionData, 'version')
                warning('Session version unknown. May have compatibility issues.');
            elseif ~strcmp(sessionData.version, obj.version)
                warning('Session version mismatch: %s vs %s', sessionData.version, obj.version);
            end
            
            % === RESTORE PARAMETERS ===
            if isfield(sessionData, 'params')
                app.data.params = sessionData.params;
                fprintf('  - Parameters restored\n');
            end
            
            % === RESTORE DATA ===
            if isfield(sessionData, 'data')
                sd = sessionData.data;
                
                % Raw data
                if isfield(sd, 'rawData')
                    app.data.rawData = sd.rawData;
                    app.data.rawDataHash = sd.rawDataHash;
                    [~, ~, T] = size(app.data.rawData);
                    app.state.maxFrame = T;
                    
                    % Update slider limits
                    app.ui.FrameSlider.Limits = [1 T];
                    app.ui.FrameField.Limits = [1 T];
                    
                    % Calculate color limits
                    abs_data = abs(app.data.rawData(:));
                    app.data.rawClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
                    if app.data.rawClim(1) == app.data.rawClim(2)
                        app.data.rawClim(2) = app.data.rawClim(1) + 1;
                    end
                    
                    fprintf('  - Raw data restored (%dx%dx%d)\n', size(app.data.rawData));
                end
                
                % Mask
                if isfield(sd, 'mask')
                    app.data.mask = sd.mask;
                    app.ui.maskStatusLabel.Text = 'Status: Loaded';
                    fprintf('  - Mask restored\n');
                end
                
                % Processing results
                if isfield(sd, 'candidateBubbles')
                    app.data.candidateBubbles = sd.candidateBubbles;
                    fprintf('  - Candidates restored (%d)\n', height(sd.candidateBubbles));
                end
                
                if isfield(sd, 'localizations')
                    app.data.localizations = sd.localizations;
                    fprintf('  - Localizations restored (%d)\n', height(sd.localizations));
                end
                
                if isfield(sd, 'tracks_raw')
                    app.data.tracks_raw = sd.tracks_raw;
                    fprintf('  - Raw tracks restored (%d)\n', length(sd.tracks_raw));
                end
                
                if isfield(sd, 'tracks_final')
                    app.data.tracks_final = sd.tracks_final;
                    fprintf('  - Final tracks restored (%d)\n', length(sd.tracks_final));
                end
                
                % Filtered data (if saved)
                if isfield(sd, 'filteredData')
                    app.data.filteredData = sd.filteredData;
                    
                    % Recalculate color limits
                    abs_data = app.data.filteredData(:);
                    app.data.filteredClim = [prctile(abs_data, 1), prctile(abs_data, 99.9)];
                    if app.data.filteredClim(1) == app.data.filteredClim(2)
                        app.data.filteredClim(2) = app.data.filteredClim(1) + 1;
                    end
                    
                    fprintf('  - Filtered data restored\n');
                end
                
                % SVD components (if saved)
                if isfield(sd, 'U')
                    app.data.U = sd.U;
                    app.data.S_diag = sd.S_diag;
                    app.data.V = sd.V;
                    app.data.svdDims = sd.svdDims;
                    fprintf('  - SVD components restored\n');
                end
                
                % DCC indices
                if isfield(sd, 'tissue_indices')
                    app.data.tissue_indices = sd.tissue_indices;
                    app.data.blood_indices = sd.blood_indices;
                    app.data.noise_indices = sd.noise_indices;
                    fprintf('  - DCC indices restored\n');
                end
            end
            
            % === RESTORE STATE ===
            if isfield(sessionData, 'state')
                app.state.currentState = sessionData.state.currentState;
                app.state.currentFrame = sessionData.state.currentFrame;
                
                % Update frame controls
                app.ui.FrameSlider.Value = app.state.currentFrame;
                app.ui.FrameField.Value = app.state.currentFrame;
                
                fprintf('  - Application state restored (level %d)\n', app.state.currentState);
            end
            
            % === RESTORE UI STATE ===
            if isfield(sessionData, 'ui')
                % Find and select the saved tab
                if isfield(sessionData.ui, 'selectedTab')
                    tabTitle = sessionData.ui.selectedTab;
                    tabs = app.ui.tabGroup.Children;
                    for i = 1:length(tabs)
                        if strcmp(tabs(i).Title, tabTitle)
                            app.ui.tabGroup.SelectedTab = tabs(i);
                            break;
                        end
                    end
                end
            end
            
            fprintf('Session restored successfully.\n');
        end
        
        function success = quickSave(obj, app, filename)
            % Quick save without dialog
            if nargin < 3
                if isempty(obj.lastSaveFile)
                    success = false;
                    warning('No previous save file. Use regular save.');
                    return;
                end
                filename = obj.lastSaveFile;
            end
            
            try
                sessionData = obj.createSessionData(app);
                save(filename, 'sessionData', '-v7.3'); % Use -v7.3 for large files
                obj.lastSaveFile = filename;
                success = true;
                fprintf('Quick save to: %s\n', filename);
            catch ME
                success = false;
                warning('Quick save failed: %s', ME.message);
            end
        end
        
        function [sessionData, filename] = autoSave(obj, app, autosaveDir)
            % Automatic save with timestamp
            if nargin < 3
                autosaveDir = tempdir;
            end
            
            % Create autosave filename with timestamp
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            filename = fullfile(autosaveDir, sprintf('ULM_autosave_%s.mat', timestamp));
            
            try
                sessionData = obj.createSessionData(app);
                save(filename, 'sessionData', '-v7.3');
                fprintf('Autosave created: %s\n', filename);
            catch ME
                warning('Autosave failed: %s', ME.message);
                sessionData = [];
                filename = '';
            end
        end
        
        function info = getSessionInfo(obj, filename)
            % Get information about a session file without loading all data
            try
                vars = whos('-file', filename);
                if ~ismember('sessionData', {vars.name})
                    error('Not a valid session file');
                end
                
                % Load only metadata
                temp = load(filename, 'sessionData');
                sd = temp.sessionData;
                
                info = struct();
                info.version = sd.version;
                info.timestamp = sd.timestamp;
                info.currentState = sd.state.currentState;
                info.maxFrame = sd.state.maxFrame;
                
                % Data summary
                info.hasRawData = isfield(sd.data, 'rawData');
                info.hasMask = isfield(sd.data, 'mask');
                info.hasLocalizations = isfield(sd.data, 'localizations');
                info.hasTracks = isfield(sd.data, 'tracks_final');
                
                if info.hasTracks
                    info.numTracks = length(sd.data.tracks_final);
                end
                
            catch ME
                warning('Could not read session info: %s', ME.message);
                info = struct('error', ME.message);
            end
        end

        function session = createSession(obj, data, params)
            % V3 wrapper: packs data + params into a pseudo-app struct,
            % delegates to the existing createSessionData, and returns
            % a flat session struct (no nested 'sessionData' key).
            fakeApp = struct();
            fakeApp.data   = data;
            fakeApp.data.params = params;
 
            % Build a minimal ui/state stub so createSessionData won't crash
            % on fields it optionally reads.
            fakeApp.ui    = struct();
            fakeApp.state = struct();
            if isfield(data, 'params'), fakeApp.data.params = data.params; end
 
            % Copy state fields that createSessionData expects
            stateFields = {'currentState','currentFrame','maxFrame'};
            for i = 1:numel(stateFields)
                if isfield(data, stateFields{i})
                    fakeApp.state.(stateFields{i}) = data.(stateFields{i});
                end
            end
 
            % We need a figure handle check — supply a dummy
            fakeApp.fig = gobjects(0);
 
            % --- Build the session directly (simplified) ---
            session = struct();
            session.version   = obj.version;
            session.timestamp = datetime('now');
            session.params    = params;
 
            % Serialize data fields
            session.data = struct();
            dataFields = {'rawData','rawDataHash','mask','candidateBubbles', ...
                          'localizations','tracks_raw','tracks_final', ...
                          'filteredData','U','S_diag','V','svdDims', ...
                          'tissue_indices','blood_indices','noise_indices', ...
                          'blockwiseDiag','vesselMap','baseVesselMap', ...
                          'rawClim','filteredClim','filteredMeanBG','MeanBG'};
            for i = 1:numel(dataFields)
                f = dataFields{i};
                if isfield(data, f) && ~isempty(data.(f))
                    session.data.(f) = data.(f);
                end
            end
 
            % State
            session.state = struct();
            if isfield(data, 'currentState'), session.state.currentState = data.currentState;
            else,                             session.state.currentState = -1; end
            if isfield(data, 'currentFrame'), session.state.currentFrame = data.currentFrame;
            else,                             session.state.currentFrame = 1; end
            if isfield(data, 'maxFrame'),     session.state.maxFrame = data.maxFrame;
            else
                if isfield(data, 'rawData') && ~isempty(data.rawData)
                    session.state.maxFrame = size(data.rawData, 3);
                else
                    session.state.maxFrame = 1;
                end
            end
 
            fprintf('Session snapshot created (v3 wrapper).\n');
        end
 
        function [data, params] = restoreSession(obj, session)
            % V3 wrapper: unpacks a session struct into separate data and
            % params outputs, matching the v3 calling convention:
            %   [app.data, app.data.params] = sessionManager.restoreSession(session)
 
            fprintf('Restoring session (v3 wrapper)...\n');
 
            % Version check
            if isfield(session, 'version') && ~strcmp(session.version, obj.version)
                warning('Session version mismatch: file=%s, current=%s', ...
                    session.version, obj.version);
            end
 
            % Params
            if isfield(session, 'params')
                params = session.params;
            else
                params = struct();
            end
 
            % Data
            data = struct();
            if isfield(session, 'data')
                flds = fieldnames(session.data);
                for i = 1:numel(flds)
                    data.(flds{i}) = session.data.(flds{i});
                end
            end
 
            % Restore state fields onto data so the GUI can pick them up
            if isfield(session, 'state')
                sf = fieldnames(session.state);
                for i = 1:numel(sf)
                    data.(sf{i}) = session.state.(sf{i});
                end
            end
 
            % Attach params into data (v3 expects data.params)
            data.params = params;
 
            fprintf('Session restored (v3 wrapper).\n');
        end
    end
end
