function [peakIdx, varargout] = findStimulusPeak(data, blankingPeriod, sampleRate, clippingThreshold, minClippedNSamples, minClippedVoltage)
%FINDSTIMULUSPEAK Summary of this function goes here
%   Detailed explanation goes here

    %% 1) Detect all clipped intervals
    % Lower clipping
    clippedVoltage = min(data);
    if clippedVoltage > -minClippedVoltage
        clippedIdxs = 1;
    else
        clippedIdxs = find(abs(data - clippedVoltage) < clippingThreshold);
    end
    clippingIntervals = [true, diff(clippedIdxs) ~= 1];
    startClippingIdxsL = clippedIdxs(clippingIntervals);
    endClippingIdxsL = clippedIdxs(clippingIntervals([2:end, 1]));

    % Upper clipping
    clippedVoltage = max(data);
    if clippedVoltage < minClippedVoltage
        clippedIdxs = 1;
    else
        clippedIdxs = find(abs(data - clippedVoltage) < clippingThreshold);
    end
    clippingIntervals = [true, diff(clippedIdxs) ~= 1];
    startClippingIdxsU = clippedIdxs(clippingIntervals);
    endClippingIdxsU = clippedIdxs(clippingIntervals([2:end, 1]));

    startClippingIdxs = [startClippingIdxsL, startClippingIdxsU];
    endClippingIdxs = [endClippingIdxsL, endClippingIdxsU];

    clippingCheck = abs(endClippingIdxs - startClippingIdxs) + 1 >= minClippedNSamples;
    startClippingIdxs = startClippingIdxs(clippingCheck);
    endClippingIdxs = endClippingIdxs(clippingCheck);

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
        if ~peakCheck(idx)
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
    peakClippingIdx = find((peakIdx >= startClippingIdxs & peakIdx <= endClippingIdxs) == 1);
    
    if ~isempty(peakClippingIdx) && length(startClippingIdxs(peakClippingIdx):endClippingIdxs(peakClippingIdx)) >= 3
        startClippingIdx = startClippingIdxs(peakClippingIdx);
        endClippingIdx = endClippingIdxs(peakClippingIdx);

        clippedSamples = startClippingIdx:endClippingIdx;
        isClipped = true;
    else
        clippedSamples = [];
        isClipped = false;
    end

    %% 4) Build varargout
    varargout{1} = isClipped;
    varargout{2} = clippedSamples;
    varargout{3} = polarity;

end

