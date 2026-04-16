classdef UndoRedoManager < handle
    % UNDOREDOMANAGER - Manages undo/redo functionality for parameter changes
    %
    % This class maintains a history of parameter states allowing users
    % to undo and redo changes.
    
    properties (Access = private)
        undoStack       % Cell array of previous states
        redoStack       % Cell array of undone states
        maxStates      % Maximum number of states to keep
        operations     % Description of each operation
    end
    
    methods
        function obj = UndoRedoManager(maxStates)
            % Constructor
            if nargin < 1
                maxStates = ULM_Constants.MAX_UNDO_STATES;
            end
            
            obj.maxStates = maxStates;
            obj.undoStack = {};
            obj.redoStack = {};
            obj.operations = {};
        end
        
        function saveState(obj, params, operation)
            % Save current parameter state
            %
            % Args:
            %   params: Current parameter structure
            %   operation: String describing the operation (e.g., 'filter', 'detection')
            
            if nargin < 3
                operation = 'unknown';
            end
            
            % Create a deep copy of params
            state = struct();
            state.params = obj.deepCopy(params);
            state.timestamp = datetime('now');
            state.operation = operation;
            
            % Add to undo stack
            obj.undoStack{end+1} = state;
            obj.operations{end+1} = operation;
            
            % Limit stack size
            if length(obj.undoStack) > obj.maxStates
                obj.undoStack(1) = [];
                obj.operations(1) = [];
            end
            
            % Clear redo stack (new action invalidates redo history)
            obj.redoStack = {};
            
            fprintf('[Undo] State saved: %s (stack size: %d)\n', operation, length(obj.undoStack));
        end
        
        function params = undo(obj)
            if isempty(obj.undoStack)
                warning('Nothing to undo');
                params = [];
                return;
            end
            
            % Pop current state to redo stack
            currentState = obj.undoStack{end};
            obj.undoStack(end) = [];
            if ~isempty(obj.operations)          % keep in sync
                obj.operations(end) = [];
            end
            obj.redoStack{end+1} = currentState;
            
            % Return the state now at the top of the stack (the "previous" state)
            if ~isempty(obj.undoStack)
                prevState = obj.undoStack{end};
                params = prevState.params;
                fprintf('[Undo] Reversed: %s (back to: %s)\n', ...
                    currentState.operation, prevState.operation);
            else
                % stack is now empty — nothing left to restore.
                % Return empty so caller knows undo reached the beginning.
                fprintf('[Undo] Reached beginning of history.\n');
                params = [];
            end
        end
        
        function params = redo(obj)
            % Redo previously undone operation
            %
            % Returns:
            %   params: Next parameter state
            
            if isempty(obj.redoStack)
                warning('Nothing to redo');
                params = [];
                return;
            end
            
            % Pop from redo stack
            nextState = obj.redoStack{end};
            obj.redoStack(end) = [];
            
            % Push back to undo stack
            obj.undoStack{end+1} = nextState;
            
            params = nextState.params;
            fprintf('[Redo] Applied: %s\n', nextState.operation);
        end
        
        function tf = canUndo(obj)
            % Check if undo is available
            tf = length(obj.undoStack) > 1; % Need at least 2 states to undo
        end
        
        function tf = canRedo(obj)
            % Check if redo is available
            tf = ~isempty(obj.redoStack);
        end
        
        function clear(obj)
            % Clear all undo/redo history
            obj.undoStack = {};
            obj.redoStack = {};
            obj.operations = {};
            fprintf('[Undo] History cleared\n');
        end
        
        function info = getInfo(obj)
            % Get information about undo/redo state
            info = struct();
            info.undoStackSize = length(obj.undoStack);
            info.redoStackSize = length(obj.redoStack);
            info.canUndo = obj.canUndo();
            info.canRedo = obj.canRedo();
            
            if ~isempty(obj.undoStack)
                info.lastOperation = obj.undoStack{end}.operation;
                info.lastTimestamp = obj.undoStack{end}.timestamp;
            end
            
            if ~isempty(obj.redoStack)
                info.nextRedoOperation = obj.redoStack{end}.operation;
            end
        end
        
        function history = getHistory(obj, maxEntries)
            % Get undo history
            %
            % Args:
            %   maxEntries: Maximum number of entries to return (default: all)
            %
            % Returns:
            %   history: Cell array of operation descriptions
            
            if nargin < 2 || maxEntries > length(obj.operations)
                maxEntries = length(obj.operations);
            end
            
            startIdx = max(1, length(obj.operations) - maxEntries + 1);
            history = obj.operations(startIdx:end);
        end
        
        function displayHistory(obj)
            % Print undo history to console
            fprintf('\n=== Undo History ===\n');
            if isempty(obj.undoStack)
                fprintf('  (empty)\n');
            else
                for i = 1:length(obj.undoStack)
                    state = obj.undoStack{i};
                    fprintf('  [%d] %s - %s\n', i, char(state.timestamp), state.operation);
                end
            end
            
            fprintf('\n=== Redo Queue ===\n');
            if isempty(obj.redoStack)
                fprintf('  (empty)\n');
            else
                for i = 1:length(obj.redoStack)
                    state = obj.redoStack{i};
                    fprintf('  [%d] %s - %s\n', i, char(state.timestamp), state.operation);
                end
            end
            fprintf('===================\n\n');
        end
    end
    
    methods (Access = private)
        function copy = deepCopy(obj, original)
            % Create a deep copy of a structure
            %
            % This is needed because MATLAB's default copy for structures
            % is shallow for handle objects.
            
            if isstruct(original)
                copy = struct();
                fields = fieldnames(original);
                for i = 1:length(fields)
                    fieldName = fields{i};
                    fieldValue = original.(fieldName);
                    
                    if isstruct(fieldValue)
                        % Recursively copy nested structures
                        copy.(fieldName) = obj.deepCopy(fieldValue);
                    elseif iscell(fieldValue)
                        % Copy cell arrays
                        copy.(fieldName) = cell(size(fieldValue));
                        for j = 1:numel(fieldValue)
                            if isstruct(fieldValue{j})
                                copy.(fieldName){j} = obj.deepCopy(fieldValue{j});
                            else
                                copy.(fieldName){j} = fieldValue{j};
                            end
                        end
                    else
                        % Direct copy for primitives and arrays
                        copy.(fieldName) = fieldValue;
                    end
                end
            else
                % Not a structure, just copy
                copy = original;
            end
        end
    end
end
