function [peakIdx, varargout] = findArtifactPeak(data, sampleRate, blankingPeriod, varargin)
%FINDARTIFACTPEAK   Find peak location and detect wether the signal is clipped or not.
%
%   peakIdx = FINDARTIFACTPEAK(data, sampleRate, blankingPeriod) returns
%   the index where the peak of the artifact is found. The blanking period
%   should be expressed in seconds.
%
%   [peakIdx, isClipped] = FINDARTIFACTPEAK(data, sampleRate, blankingPeriod)
%   returns a boolean flag telling if data are clipped or not.
%
%   [peakIdx, isClipped, clippedSamples] = FINDARTIFACTPEAK(data, sampleRate, blankingPeriod)
%   returns the samples that are actually clipped due to amplifier saturation.
%
%   [peakIdx, isClipped, clippedSamples, polarity] = FINDARTIFACTPEAK(data, sampleRate, blankingPeriod)
%   returns the polarity of the artifact. A positive polarity means that
%   the signal decays towards the baseline from larger values, while a
%   negative polarity means that the signal goes back to baseline from
%   smaller values.
%
%   [...] = FINDARTIFACTPEAK(..., saturationVoltage) specifies the recording
%   system operating range in mV as specified in the datasheet. This is useful
%   to properly detect saturation. If a scalare is provided, then the operating
%   range is assumed to be symmetric with respect to 0, otherwise specify lower
%   and upper boundaries through an array. By default, the 95% of the maximum
%   absolute value of the data is employed.
%
%   [...] = FINDARTIFACTPEAK(..., saturationVoltage, minClippedNSamples)
%   specifies the minimum number of consecutive clipped samples to flag the
%   artifact as a clipped one.

    %% 0) Check and parse input arguments
    if nargin < 3
        throw(MException('SAR:NotEnoughParameters', 'The parameters data, sampleRate, and blankingPeriod are required.'));
    end
    
    if nargin < 4 || isempty(varargin{1})
        saturationVoltage = 0.95 * max(abs(data)) / 1e3;
    else
        saturationVoltage = varargin{1};
    end

    if nargin < 5 || isempty(varargin{2})
        minClippedNSamples = 2;
    else
        minClippedNSamples = varargin{2};
    end

    if isscalar(saturationVoltage)
        saturationVoltage = [-saturationVoltage, saturationVoltage];
    end

    saturationVoltage = [min(saturationVoltage), max(saturationVoltage)] * 1e3;

    %% 1) Detect all clipped intervals
    startClippingIdxs = [];
    endClippingIdxs = [];

    for i = 1:numel(saturationVoltage)
        clippedIdxs = find((-1)^i * data > (-1)^i * saturationVoltage(i));

        if ~isempty(clippedIdxs)
            clippedIntervals = [true, diff(clippedIdxs) ~= 1];
            startClippingIdxs = [startClippingIdxs, clippedIdxs(clippedIntervals)];
            endClippingIdxs = [endClippingIdxs, clippedIdxs(clippedIntervals([2:end, 1]))];
        end
    end
    
    if ~isempty(startClippingIdxs) && ~isempty(endClippingIdxs)
        samplesCheck = abs(endClippingIdxs - startClippingIdxs) + 1 >= minClippedNSamples;

        startClippingIdxs = unique(startClippingIdxs(samplesCheck));
        endClippingIdxs = unique(endClippingIdxs(samplesCheck));
    end

    %% 2) Determine the polarity of the artifact, flipped data are useful
    % to properly identify upper and lower peaks, for polarity
    % and also to account for clipped data.
    blankingSamples = 1:round(blankingPeriod * sampleRate);

    maxValue = max(data(blankingSamples)) * 0.975;
    minValue = min(data(blankingSamples)) * 0.975;

    [~, maxIdx] = findpeaks([0, flip(data(blankingSamples))], 'NPeaks', 1, 'MinPeakHeight', maxValue);
    [~, minIdx] = findpeaks([0, -flip(data(blankingSamples))], 'NPeaks', 1, 'MinPeakHeight', abs(minValue));
    maxIdx = (length(blankingSamples) + 1) - maxIdx + 1;
    minIdx = (length(blankingSamples) + 1) - minIdx + 1;

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

    if ~isempty(peakIdx) && ~isempty(startClippingIdxs) && ~isempty(endClippingIdxs)
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

    %% 5) Plot
    % fig = figure();
    % hold('on');
    % plot(data);
    % scatter(startClippingIdxs, data(startClippingIdxs), 'green');
    % scatter(endClippingIdxs, data(endClippingIdxs), 'red');
    % uiwait(fig);

end

