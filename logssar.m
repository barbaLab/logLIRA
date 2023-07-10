function output = logssar(signal, stimIdxs, sampleRate, blankingPeriod)
%LOGSSAR LOGarithmically distributed intervals Spline Stimulus Artifacts Rejection.
%   Detailed explanation goes here
    %% 0) Check arguments and init variables
    output = signal;
    blankingNSamples = round(blankingPeriod * sampleRate);

    searchWindow = 3e-3;
    searchNSamples = round(searchWindow * sampleRate);

    correctionWindow = 0.2e-3;
    correctionNSamples = round(correctionWindow * sampleRate);

    waitbarFig = waitbar(0, 'Starting...');

    %% 1) Find signal baseline and ISI
    [~, baselinePercentiles] = getBaselineFromStim(signal, stimIdxs, sampleRate, 10e-3, 0.5e-3, [25, 75]);
    lowerBaseline = baselinePercentiles(1);
    upperBaseline = baselinePercentiles(end);

    ISI = [diff(stimIdxs), length(signal) - stimIdxs(end)];

    %% 2) Clean each artifact iteratively
    for idx = 1:numel(stimIdxs)
        data = signal((1:ISI(idx)) + stimIdxs(idx) - 1);

        % Identify the samples to clean
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
        correctionArtifact = interp1(correctionX, correctionY, (1:length(artifact)) + stimIdxs(idx) - 1, 'linear');

        % Update output signal
        output((1:length(artifact)) + stimIdxs(idx) - 1) = data - artifact + correctionArtifact;

        % Update progress bar
        waitbar(idx / numel(stimIdxs), waitbarFig, 'Removing artifacts...');
    end

    close(waitbarFig);
end