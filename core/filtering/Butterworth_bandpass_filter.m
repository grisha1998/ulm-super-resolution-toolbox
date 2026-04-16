function filteredData = Butterworth_bandpass_filter(Data, cutoffFreq, samplingFreq, filterOrder)
% =========================================================================
% FUNCTION: Butterworth_bandpass_filter
% AUTHOR: Grigori Shapiro
% =========================================================================
%
% PURPOSE & ADVANTAGES:
%   Applies a temporal frequency Infinite Impulse Response (IIR) Butterworth 
%   bandpass filter along the slow-time (frame) dimension.
%   - Advantages: SVD separates components by mathematical coherence, which 
%     can sometimes leave behind noise that shares similar spatial structures 
%     but moves at the wrong speed. A Butterworth filter rigidly enforces 
%     physical speed limits: filtering out anything moving too slow 
%     (residual tissue) or too fast (electronic noise).
%
% DETAILED METHODOLOGY (LOGIC & WORKFLOW):
%   1. Safety Check: Calculates the Nyquist frequency (framerate / 2). If the 
%      requested upper cutoff exceeds this, it silently clamps it to 99% of 
%      the Nyquist limit to prevent filter instability.
%   2. Coefficient Generation: Uses MATLAB's `butter` function to generate 
%      the (b, a) transfer function coefficients for the specified order.
%   3. Filtering: Applies `filter` along the 3rd dimension (time) of the 
%      3D matrix.
%   4. Clean up: Replaces any non-finite values (NaN/Inf) created during 
%      filtering with zeros.
%
% SYNTAX OPTIONS:
%   filt_data = Butterworth_bandpass_filter(Data, [low_Hz, high_Hz], PRF, order)
%
% EXAMPLES:
%   % Example: Keep frequencies between 30 Hz and 250 Hz (3rd order filter)
%   cleanedData = Butterworth_bandpass_filter(bloodData, [30, 250], 1000, 3);
%
% INPUTS:
%   Data         - (Type: 3D Numeric Array) The input data sequence, typically 
%                  already processed by an SVD filter.
%   cutoffFreq   - (Type: 1x2 Double) [low_Hz, high_Hz] boundary frequencies.
%                  Example: [20, 300].
%   samplingFreq - (Type: Double) The PRF / framerate in Hz. Example: 500.
%   filterOrder  - (Type: Integer) The order of the filter. Higher orders 
%                  create steeper roll-offs but can introduce phase distortion.
%                  Example: 2 or 3.
%
% OUTPUTS:
%   filteredData - (Type: 3D Numeric Array) The temporally filtered sequence.
% =========================================================================

    % Calculate the Nyquist frequency
    nyquistFreq = samplingFreq / 2;
    if cutoffFreq(2) > nyquistFreq
        cutoffFreq(2) = floor(0.99 * nyquistFreq);
    end
    % Get the filter coefficients (b, a) for the bandpass filter
    [b, a] = butter(filterOrder, cutoffFreq / nyquistFreq, 'bandpass');
    
    % Apply the filter along the 3rd dimension (time/frames)
    % The result overwrites the SVD-filtered data
    filteredData = filter(b, a, Data, [], 3);
    filteredData(~isfinite(filteredData)) = 0;

end