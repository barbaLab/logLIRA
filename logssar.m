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
    if nargin < 3
    throw(MException('SAR:NotEnoughParameters', 'The parameters signal, stimIdxs, and sampleRate are required.'));
    end
    
    if nargin < 4
        blankingPeriod = 1e-3;
    else
        blankingPeriod = varargin{1};
    end

    if nargin < 5
        searchWindow = 3e-3;
    else
        searchWindow = varargin{2};
    end

    if nargin < 6
        correctionWindow = 0.2e-3;
    else
        correctionWindow = varargin{3};
    end

    if nargin < 7
        correctionMethod = 'linear';
    else
        correctionMethod = varargin{4};
    end

    output = signal;

    blankingNSamples = round(blankingPeriod * sampleRate);
    searchNSamples = round(searchWindow * sampleRate);
    correctionNSamples = round(correctionWindow * sampleRate);

    waitbarFig = waitbar(0, 'Starting...');

    %% 1) Find signal baseline and ISI
    [~, baselinePercentiles] = getBaselineFromStim(signal, stimIdxs, sampleRate, 10e-3, 0.5e-3, [25, 75]);
    lowerBaseline = baselinePercentiles(1);
    upperBaseline = baselinePercentiles(end);

    ISI = [diff(stimIdxs), length(signal) - stimIdxs(end)];

    %% 2) Clean each artifact iteratively
    for idx = 1:numel(stimIdxs)
        % Identify the samples to clean
        data = signal((1:ISI(idx)) + stimIdxs(idx) - 1);
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

        data = data(1:(searchNSamples + searchOffset));

        % Find the artifact shape
        artifact = fitArtifact(data, sampleRate, blankingPeriod);

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