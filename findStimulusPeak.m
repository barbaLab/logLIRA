function [peakIdx, varargout] = findStimulusPeak(data, sampleRate, blankingPeriod, varargin)
%FINDSTIMULUSPEAK   Find peak location and detect wether the signal is clipped or not.
%
%   peakIdx = FINDSTIMULUSPEAK(data, sampleRate, blankingPeriod) returns
%   the index where the peak of the stimulus is found. It always belongs to
%   the blanking period, that is expressed in seconds.
%
%   [peakIdx, isClipped] = FINDSTIMULUSPEAK(data, sampleRate, blankingPeriod)
%   returns a boolean flag telling if data are clipped or not.
%
%   [peakIdx, isClipped, clippedSamples] = FINDSTIMULUSPEAK(data, sampleRate, blankingPeriod)
%   returns the samples that are actually clipped due to amplifier saturation.
%
%   [peakIdx, isClipped, clippedSamples, polarity] = FINDSTIMULUSPEAK(data, sampleRate, blankingPeriod)
%   returns the polarity of the artifact. A positive polarity means that
%   the signal decays towards the baseline from larger values, while a
%   negative polarity means that the signal goes back to baseline from
%   smaller values.
%
%   [...] = FINDSTIMULUSPEAK(..., saturationVoltage) specifies the recording
%   system operating range in mV as specified in the datasheet. This is useful
%   to properly detect saturation. If a scalare is provided, then the operating
%   range is assumed to be symmetric with respect to 0, otherwise specify lower
%   and upper boundaries through an array.
%
%   [...] = FINDSTIMULUSPEAK(..., saturationVoltage, minClippedNSamples)
%   specifies the minimum number of consecutive clipped samples to flag the
%   artifact as a clipped one.

    %% 0) Check and parse input arguments
    if nargin < 3
        throw(MException('SAR:NotEnoughParameters', 'The parameters data, sampleRate, and blankingPeriod are required.'));
    end
    
    if nargin < 4
        saturationVoltage = 5;
    end

    if nargin < 5
        minClippedNSamples = 2;
    end

    if isscalar(saturationVoltage)
        saturationVoltage = [-saturationVoltage, saturationVoltage];
    end

    saturationVoltage = [min(saturationVoltage), max(saturationVoltage)] * 1e3;

    %% 1) Detect all clipped intervals
    startClippingIdxs = [];
    endClippingIdxs = [];

    % Lower clipping
    clippedIdxs = find(data < saturationVoltage(1));
    if ~isempty(clippedIdxs)
        clippingIntervals = [true, diff(clippedIdxs) ~= 1];
        startClippingIdxs = [startClippingIdxs, clippedIdxs(clippingIntervals)];
        endClippingIdxs = [endClippingIdxs, clippedIdxs(clippingIntervals([2:end, 1]))];
    end

    % Upper clipping
    clippedIdxs = find(data > saturationVoltage(2));
    if ~isempty(clippedIdxs)
        clippingIntervals = [true, diff(clippedIdxs) ~= 1];
        startClippingIdxs = [startClippingIdxs, clippedIdxs(clippingIntervals)];
        endClippingIdxs = [endClippingIdxs, clippedIdxs(clippingIntervals([2:end, 1]))];
    end
    
    if ~isempty(startClippingIdxs) && ~isempty(endClippingIdxs)
        clippingCheck = abs(endClippingIdxs - startClippingIdxs) + 1 >= minClippedNSamples;
        startClippingIdxs = startClippingIdxs(clippingCheck);
        endClippingIdxs = endClippingIdxs(clippingCheck);
    end

    %% 2) Determine the polarity of the artifact, flipped data are useful
    % to properly identify upper and lower peaks, for polarity
    % and also to account for clipped data.
    blankingSamples = 1:round(blankingPeriod * sampleRate);
    [~, maxIdx] = max(flip(data(blankingSamples)));
    [~, minIdx] = min(flip(data(blankingSamples)));
    maxIdx = length(blankingSamples) + 1 - maxIdx;
    minIdx = length(blankingSamples) + 1 - minIdx;

    peakIdx = [minIdx, maxIdx];
    polarity = [-1, 1];
    peakCheck = islocalmax(data) | islocalmin(data);
    peakCheck = peakCheck(peakIdx);

    for idx = 1:length(peakIdx)    
        if ~peakCheck(idx) && ~isempty(startClippingIdxs) && ~isempty(endClippingIdxs)
            if sum((startClippingIdxs <= peakIdx(idx)) & (endClippingIdxs >= peakIdx(idx))) >= 1
                peakCheck(idx) = true;
            end
        end
    end

    peakIdx = peakIdx(peakCheck);
    polarity = polarity(peakCheck);
    
    if length(peakIdx) > 1    
        polarity = maxIdx > minIdx;
        peakIdx = peakIdx(polarity + 1);

        if polarity == 0
            polarity = -1;
        end
    end

    %% 3) Check if peak is clipped
    clippedSamples = [];
    isClipped = false;

    if ~isempty(startClippingIdxs) && ~isempty(endClippingIdxs)
        peakClippingIdx = find((peakIdx >= startClippingIdxs & peakIdx <= endClippingIdxs) == 1);
        
        if ~isempty(peakClippingIdx) && length(startClippingIdxs(peakClippingIdx):endClippingIdxs(peakClippingIdx)) >= 3
            startClippingIdx = startClippingIdxs(peakClippingIdx);
            endClippingIdx = endClippingIdxs(peakClippingIdx);
    
            clippedSamples = startClippingIdx:endClippingIdx;
            isClipped = true;
        end
    end

    %% 4) Return output values
    varargout{1} = isClipped;
    varargout{2} = clippedSamples;
    varargout{3} = polarity;

end

