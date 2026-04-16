% FILENAME: getExpParams.m
%
% PURPOSE:
%   Parses a detailed experiment text file to extract key physical and
%   setup parameters for the ULM acquisition.
%
%   This utility function reads a specified text file and uses regular
%   expressions to extract a rich set of metadata, including bubble type,
%   phantom geometry, flow rates, and acquisition settings. This allows for
%   automated, context-aware configuration of the processing pipeline.
%
% SYNTAX:
%   = getExpParams(filePath)
%
% INPUTS:
%   filePath: (char) The path to the experiment's.txt parameter file.
%
% OUTPUTS:
%   A struct containing all extracted experimental parameters.
%
% NOTES:
%   This function was provided by the user to replace the previous
%   filename-based parsing, enabling more sophisticated parameter tuning.
%
% AUTHOR: Grigori Shapiro
% DATE: July 14, 2025
%
% -------------------------------------------------------------------------

function expParams = getExpParams(filePath)
    fileContent = fileread(filePath);
    expParams = struct(...
        'bubbleType', '',...
        'mainChannelDiameter', NaN,...
        'secondaryChannelDiameter', NaN,...
        'angle', 45,... % Default angle if not found
        'shootingRate', NaN,...
        'flowSpeed', NaN,...
        'fovZ', NaN,...
        'fovX', NaN,...
        'frequency', NaN...
    );

    % --- DEFINING REGEX PATTERNS ---
    
    % 1. Bubble Type
    bubbleTypePattern = '(?<bubbleType>NanoDropletsC6|NBs|MBs)';
    
    % 2. Channels (Simpler - just looks for the 300/300 part)
    channelPattern = '(\d+)\/(\d+)\s*\[?um\]?';
    
    % 3. Angle (New - Independent search)
    % Matches a number followed optionally by brackets and 'deg'
    anglePattern = '(?<angle>\d+\.?\d*)\s*\[?deg\]?';

    % 4. Other Patterns
    shootingRatePattern = 'FR: (?<shootingRate>\d+) ?Hz';
    flowSpeedPattern = 'Flow Rate: (?<flowSpeed>\d+\.?\d*) ?ml/min';
    fovZPattern = 'FOV_Z: (\d+\.?\d*)-(\d+\.?\d*)( ?mm)?'; 
    fovXPattern = 'FOV_X: (?<fovX>\d+\.?\d*)( ?mm)?';
    frequencyPattern = 'f1? = (?<frequency>\d+\.?\d*) ?MHz';

    % --- EXTRACTION ---

    % Extract bubble type
    bubbleTypeMatch = regexp(fileContent, bubbleTypePattern, 'names');
    if ~isempty(bubbleTypeMatch), expParams.bubbleType = bubbleTypeMatch.bubbleType; end
    
    % Extract channel diameters
    channelMatch = regexp(fileContent, channelPattern, 'tokens');
    if ~isempty(channelMatch)
        % channelMatch{1} contains {Main, Secondary}
        tokens = channelMatch{1};
        expParams.mainChannelDiameter = str2double(tokens{1});
        expParams.secondaryChannelDiameter = str2double(tokens{2});
    end

    % Extract Angle (Independent Logic)
    angleMatch = regexp(fileContent, anglePattern, 'names');
    if ~isempty(angleMatch)
        expParams.angle = str2double(angleMatch.angle);
    end

    % Extract shooting rate
    shootingRateMatch = regexp(fileContent, shootingRatePattern, 'names');
    if ~isempty(shootingRateMatch), expParams.shootingRate = str2double(shootingRateMatch.shootingRate); end
    
    % Extract flow speed
    flowSpeedMatch = regexp(fileContent, flowSpeedPattern, 'names');
    if ~isempty(flowSpeedMatch), expParams.flowSpeed = str2double(flowSpeedMatch.flowSpeed); end
    
    % Extract FOV_Z
    fovZMatch = regexp(fileContent, fovZPattern, 'tokens');
    if ~isempty(fovZMatch)
        z1 = str2double(fovZMatch{1}{1}); z2 = str2double(fovZMatch{1}{2});
        expParams.fovZ = abs(z2 - z1);
    else
        % Fallback for generic FOV
        genericFovRangePattern = 'FOV: ?(\d+\.?\d*)-(\d+\.?\d*)';
        genericFovRangeMatch = regexp(fileContent, genericFovRangePattern, 'tokens');
        if ~isempty(genericFovRangeMatch)
            z1 = str2double(genericFovRangeMatch{1}{1}); z2 = str2double(genericFovRangeMatch{1}{2});
            expParams.fovZ = abs(z2 - z1);
        end
    end
    
    % Extract FOV_X
    fovXMatch = regexp(fileContent, fovXPattern, 'names');
    if ~isempty(fovXMatch), expParams.fovX = str2double(fovXMatch.fovX); end
    
    % Extract frequency
    frequencyMatch = regexp(fileContent, frequencyPattern, 'names');
    if ~isempty(frequencyMatch), expParams.frequency = str2double(frequencyMatch.frequency); end
end