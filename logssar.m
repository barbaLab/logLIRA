function output = logssar(signal, stimIdxs, sampleRate, varargin)
%LOGSSAR LOGarithmically distributed intervals Spline Stimulus Artifacts Rejection.
%   output = LOGSSAR(signal, stimIdxs, sampleRate) returns the input signal
%   without the artifacts caused by electrical stimulation. The stimIdxs
%   are the indexes of stimulation onsets.
%
%   output = LOGSSAR(..., blankingPeriod) specifies the time after the
%   stimulus onset that is discarded. It must be expressed in seconds. By
%   default it is 1 ms.
%
%   output = LOGSSAR(..., blankingPeriod, searchWindow) specifies the time
%   window used to find when the artifact disappears as it returns to
%   baseline. It must be expressed in seconds. By default it is 3 ms.
%
%   output = LOGSSAR(..., blankingPeriod, searchWindow, correctionWindow)
%   specifies the time window used to correct the signal without artifacts
%   to avoid discontinuities. It must be expressed in seconds. By default
%   it is 0.2 ms.
%
%   output = LOGSSAR(..., blankingPeriod, searchWindow, correctionWindow, correctionMethod)
%   specifies the interpolation method used to shift the signal once the
%   artifacts are rejected to avoid discontinuities. Possible values are
%   'linear', 'cubic' or 'spline'. By default it is 'linear'.

    %% 0) Check and parse input arguments
    addpath(genpath('./src'));


    blankingPeriod = 1e-3;
    searchWindow = 3e-3;
    correctionWindow = 0.2e-3;
    correctionMethod = 'linear';
    sFraction = 0.05;
    paddingDuration = 1e-3;

    validNumPosCheck = @(x) isnumeric(x) && (x >= 0);
    
    parser = inputParser;
    addRequired(parser, 'signal', @isnumeric);
    addRequired(parser, 'stimIdxs', @(x) isnumeric(x) && all(x > 0));
    addRequired(parser, 'sampleRate', validNumPosCheck);
    addOptional(parser, 'blankingPeriod', blankingPeriod, validNumPosCheck);
    addParameter(parser, 'searchWindow', searchWindow, validNumPosCheck);
    addParameter(parser, 'correctionWindow', correctionWindow, validNumPosCheck);
    addParameter(parser, 'correctionMethod', correctionMethod, ...
        @(x) any(validatestring(x, {'linear', 'pchip', 'cubic', 'v5cubic', 'makima', 'spline'})));
    addParameter(parser, 'sFraction', sFraction, @(x) isnumeric(x) && (x > 0) && (x <= 1));
    addParameter(parser, 'paddingDuration', paddingDuration, validNumPosCheck);
    addParameter(parser, 'saturationVoltage', [], @isnumeric);
    addParameter(parser, 'minClippedNSamples', [], validNumPosCheck);

    parse(parser, signal, stimIdxs, sampleRate, varargin{:});

    signal = parser.Results.signal;
    stimIdxs = parser.Results.stimIdxs;
    sampleRate = parser.Results.sampleRate;
    blankingPeriod = parser.Results.blankingPeriod;
    searchWindow = parser.Results.searchWindow;
    correctionWindow = parser.Results.correctionWindow;
    correctionMethod = parser.Results.correctionMethod;
    sFraction = parser.Results.sFraction;
    paddingDuration = parser.Results.paddingDuration;
    saturationVoltage = parser.Results.saturationVoltage;
    minClippedNSamples = parser.Results.minClippedNSamples;

    output = signal;

    blankingNSamples = round(blankingPeriod * sampleRate);
    searchNSamples = round(searchWindow * sampleRate);
    correctionNSamples = round(correctionWindow * sampleRate);

    waitbarFig = waitbar(0, 'Starting...');

    %% 1) Find signal baseline and ISI
    [~, baselinePercentiles] = getBaselineFromStim(signal, stimIdxs, sampleRate, 10e-3, 0.5e-3, [25, 75]);
    lowerBaseline = baselinePercentiles(1);
    upperBaseline = baselinePercentiles(end);

    IAI = [diff(stimIdxs), length(signal) - stimIdxs(end)];
    minArtifactDuration = 0.03;
    minArtifactNSamples = min(min(IAI), round(minArtifactDuration * sampleRate));

    %% 2) Clean each artifact iteratively
    for idx = 1:numel(stimIdxs)
        % Identify the samples to clean
        data = signal((1:IAI(idx)) + stimIdxs(idx) - 1);
        hasReachedBaseline = false;
        searchOffset = blankingNSamples;

        while ~hasReachedBaseline && (searchNSamples + searchOffset) < length(data)
            searchMedian = median(data((1:searchNSamples) + searchOffset));

            if searchMedian > lowerBaseline && searchMedian < upperBaseline
                hasReachedBaseline = true;
            else
                searchOffset = searchOffset + 1;
            end
        end

        if searchNSamples + searchOffset > minArtifactNSamples
            data = data(1:(searchNSamples + searchOffset));
        else
            data = data(1:minArtifactNSamples);
        end

        % Find the artifact shape
        artifact = fitArtifact(data, sampleRate, blankingPeriod, ...
            'sFraction', sFraction, 'paddingDuration', paddingDuration, ...
            'saturationVoltage', saturationVoltage, 'minClippedNSamples', minClippedNSamples);

        % Correct discontinuities
        correctionX = [-correctionNSamples:-1, (1:correctionNSamples) + length(artifact) - 1] + stimIdxs(idx);
        correctionY = output(correctionX);
        correctionArtifact = interp1(correctionX, correctionY, (1:length(artifact)) + stimIdxs(idx) - 1, correctionMethod);

        % Update output signal
        output((1:length(artifact)) + stimIdxs(idx) - 1) = data - artifact + correctionArtifact;

        % Update progress bar
        waitbar(idx / numel(stimIdxs), waitbarFig, 'Removing artifacts...');
    end

    close(waitbarFig);
end