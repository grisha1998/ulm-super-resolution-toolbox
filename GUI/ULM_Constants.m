classdef ULM_Constants
    % ULM_CONSTANTS - Centralized constants for ULM processing
    %
    % This class contains all constant values used throughout the ULM
    % system to ensure consistency and easy maintenance.
    
    properties (Constant)
        % === DEFAULT VALUES ===
        DEFAULT_FRAMERATE = 1000;           % Hz
        DEFAULT_PIXEL_SIZE_X = 0.1;         % mm
        DEFAULT_PIXEL_SIZE_Z = 0.1;         % mm
        
        % === FILTER PARAMETERS ===
        DEFAULT_SVD_CUTOFF = [1, 100];
        DEFAULT_BUTTER_CUTOFF = [10, 100]; % Hz
        DEFAULT_BUTTER_ORDER = 4;
        
        % === DETECTION PARAMETERS ===
        MIN_DETECTION_THRESHOLD = 0.01;
        MAX_DETECTION_THRESHOLD = 1.0;
        DEFAULT_DETECTION_THRESHOLD = 0.3;
        DEFAULT_MAX_BUBBLES_PER_FRAME = 2000;
        MAX_BUBBLES_LIMIT = 5000;
        
        % === LOCALIZATION PARAMETERS ===
        DEFAULT_FWHM = [1.5, 1.5];          % pixels
        DEFAULT_GAUSS_BOX_RADIUS = 3;       % pixels
        MIN_GAUSS_BOX_RADIUS = 2;
        MAX_GAUSS_BOX_RADIUS = 10;
        DEFAULT_QC_MAX_SHIFT_FACTOR = 2.0;
        
        % === TRACKING PARAMETERS ===
        MIN_LINKING_DISTANCE = 0.1;         % pixels
        MAX_LINKING_DISTANCE = 10.0;        % pixels
        DEFAULT_LINKING_DISTANCE = 2.0;     % pixels
        
        MIN_GAP_CLOSING = 0;                % frames
        MAX_GAP_CLOSING = 10;               % frames
        DEFAULT_GAP_CLOSING = 2;            % frames
        
        MIN_TRACK_LENGTH = 2;               % localizations
        MAX_TRACK_LENGTH = 20;              % localizations
        DEFAULT_TRACK_LENGTH = 3;           % localizations
        
        DEFAULT_KALMAN_PROCESS_NOISE = 10;
        DEFAULT_KALMAN_MODEL = 'ConstantVelocity';
        DEFAULT_ASSIGNMENT_METHOD = 'hungarian';
        
        % === QC PARAMETERS ===
        DEFAULT_MAX_ANGLE_CHANGE = 90;      % degrees
        DEFAULT_ACCELERATION_C_FACTOR = 3.0;
        DEFAULT_VD_RATIO = 2.0;
        
        % === POST-PROCESSING ===
        MIN_SMOOTHING_WINDOW = 3;
        MAX_SMOOTHING_WINDOW = 21;
        DEFAULT_SMOOTHING_WINDOW = 5;
        DEFAULT_INTERPOLATION_STEP = 0.5;
        
        % === RENDERING ===
        MIN_UPSAMPLING_FACTOR = 1;
        MAX_UPSAMPLING_FACTOR = 10;
        DEFAULT_UPSAMPLING_FACTOR = 3;
        DEFAULT_RENDER_METHOD = 'histogram';
        
        % === ROI PARAMETERS ===
        DEFAULT_ROI_GAMMA = 1.0;
        MIN_ROI_GAMMA = 0.1;
        MAX_ROI_GAMMA = 3.0;
        DEFAULT_ROI_THRESHOLD = 0.0;
        
        % === UI PARAMETERS ===
        DEBOUNCE_DELAY = 0.05;              % seconds
        MEMORY_UPDATE_INTERVAL = 2.0;       % seconds
        PLAYBACK_FPS = 20;                  % frames per second
        
        % === UNDO/REDO ===
        MAX_UNDO_STATES = 20;               % Maximum undo history
        
        % === PERFORMANCE ===
        LARGE_DATA_THRESHOLD = 1e8;         % elements (trigger warnings)
        PROGRESS_UPDATE_INTERVAL = 100;     % iterations
        
        % === FILE I/O ===
        DEFAULT_DATA_FOLDER = pwd;
        SUPPORTED_DATA_FORMATS = {'*.mat'};
        SESSION_FILE_EXT = '.mat';
        
        % === COLORMAP RANGES ===
        DEFAULT_PERCENTILE_LOW = 1;
        DEFAULT_PERCENTILE_HIGH = 99.9;
        DENSITY_PERCENTILE_HIGH = 99.5;
        VELOCITY_PERCENTILE_HIGH = 99.5;
        
        % === SCALE BAR ===
        SCALE_BAR_LENGTH_MM = 1.0;
        SCALE_BAR_X_POSITION = 0.05;        % fraction of width
        SCALE_BAR_Y_POSITION = 0.95;        % fraction of height
        
        % === ENHANCEMENT ===
        CLAHE_MIN_CLIP = 0.001;
        CLAHE_MAX_CLIP = 0.05;
        CLAHE_NUM_TILES = 8;
        TOPHAT_MIN_RADIUS = 1;
        TOPHAT_MAX_RADIUS = 11;
        SHARPEN_MAX_AMOUNT = 2.0;

        % === PLAYBACK / TIMERS ===
        PLAYBACK_TIMER_PERIOD = 0.05;       % seconds per frame advance

        % === PROCESSING ===
        GAUSSIAN_FIT_BATCH_SIZE = 500;
        SGOLAY_POLY_ORDER = 3;

        % === DISPLAY ===
        TRACK_COLORMAP_SIZE = 256;
        RUN_BUTTON_COLOR = [0.6 1 0.6];
    end
    
    methods (Static)
        function val = validateRange(val, minVal, maxVal)
            % Clamp value to valid range
            val = max(minVal, min(maxVal, val));
        end
        
        function val = roundToOdd(val)
            % Round to nearest odd number (for window sizes)
            val = round(val);
            if mod(val, 2) == 0
                val = val + 1;
            end
        end
        
        function isValid = isValidSVDCutoff(cutoff)
            % Validate SVD cutoff range
            isValid = length(cutoff) == 2 && cutoff(1) < cutoff(2) && ...
                      cutoff(1) >= 1 && all(cutoff > 0);
        end
        
        function isValid = isValidButterCutoff(cutoff, framerate)
            % Validate Butterworth cutoff (must be below Nyquist)
            nyquist = framerate / 2;
            isValid = length(cutoff) == 2 && cutoff(1) < cutoff(2) && ...
                      all(cutoff > 0) && cutoff(2) < nyquist;
        end
        
        function mem = estimateMemoryUsage(dataSize)
            % Estimate memory usage in MB for given data size
            % Assumes double precision (8 bytes per element)
            bytesPerElement = 8;
            mem = (prod(dataSize) * bytesPerElement) / (1024^2);
        end
        
        function fps = calculatePlaybackFPS(framerate, desiredFPS)
            % Calculate appropriate playback speed
            % Returns actual FPS to use based on data framerate
            if nargin < 2
                desiredFPS = ULM_Constants.PLAYBACK_FPS;
            end
            
            if framerate > desiredFPS
                % Slow down playback
                fps = desiredFPS;
            else
                % Play at native framerate
                fps = framerate;
            end
        end
    end
end
