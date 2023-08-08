function [peakIdx, varargout] = findArtifactPeak(data, sampleRate, blankingPeriod, varargin)
%FINDARTIFACTPEAK Find peak location and detect wether the signal is clipped or not.
%
%   peakIdx = FINDARTIFACTPEAK(data, sampleRate, blankingPeriod) returns
%   the index where the peak of the artifact is found. The peak is defined as the
%   last point after which it becomes possible to recover data. Blanking period
%   should be expressed in seconds.
%
%   [peakIdx, isClipped] = FINDARTIFACTPEAK(data, sampleRate, blankingPeriod)
%   returns a boolean flag true if data are clipped, false otherwise.
%
%   [peakIdx, isClipped, polarity] = FINDARTIFACTPEAK(data, sampleRate, blankingPeriod)
%   returns the polarity of the artifact. A positive polarity means that
%   the signal decays towards the baseline from larger values, while a
%   negative polarity means that the signal goes back to baseline from
%   smaller values.
%
%   [...] = FINDARTIFACTPEAK(..., saturationVoltage) specifies the recording
%   system operating range in mV as specified in the datasheet. This is useful
%   to properly detect saturation. If a scalare is provided, then the operating
%   range is assumed to be symmetric with respect to 0, otherwise specify lower
%   and upper boundaries through an array. By default, 95% of the input data absolute
%   value maximum is employed.
%
%   [...] = FINDARTIFACTPEAK(..., saturationVoltage, minClippedNSamples)
%   specifies the minimum number of consecutive clipped samples to mark the
%   artifact as a clipped one. It should be a 1x1 positive integer. By default,
%   it is 2.

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

    %% 1) Find peakIdx
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

    peakIdx = peakIdx(peakCheck);
    polarity = polarity(peakCheck);
    
    if length(peakIdx) > 1    
        polarity = maxIdx > minIdx;
        peakIdx = peakIdx(polarity + 1);

        if polarity == 0
            polarity = -1;
        end
    end

    %% 2) Detect clipping
    startClippingIdxs = [];
    endClippingIdxs = [];
    isClipped = false;

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

        if ~isempty(startClippingIdxs) && ~isempty(endClippingIdxs)
            isClipped = true;
            peakIdx = max([peakIdx, max(endClippingIdxs)]);
            polarity = sign(data(peakIdx) - median(data));
        end
    end

    %% 3) Return output values
    varargout{1} = isClipped;
    varargout{2} = polarity;

    %% 4) Plot
    % fig = figure();
    % hold('on');
    % plot(data);
    % scatter(startClippingIdxs, data(startClippingIdxs), 'green');
    % scatter(endClippingIdxs, data(endClippingIdxs), 'red');
    % scatter(peakIdx, data(peakIdx), 'black', 'Marker', '*');
    % uiwait(fig);

end