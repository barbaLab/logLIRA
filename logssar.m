function output = logssar(signal, stimIdxs, sampleRate, blankingPeriod)
%LOGSSAR LOGarithmically distributed intervals Spline Stimulus Artifacts Rejection.
%   Detailed explanation goes here
    %% 0) Check arguments and init variables
    output = signal;
    blankingNSamples = round(blankingPeriod * sampleRate);

    searchWindow = 3e-3;
    searchSamples = 1:round(searchWindow * sampleRate);

    correctionWindow = 0.2e-3;
    correctionNSamples = round(correctionWindow * sampleRate);

    waitbarFig = waitbar(0, 'Starting...');

    %% 1) Find signal baseline and ISI
    [~, baselinePercentiles] = getBaselineFromStim(signal, stimIdxs, sampleRate, 10e-3, 0.5e-3, [25, 75]);
    lowerBaseline = baselinePercentiles(1);
    upperBaseline = baselinePercentiles(end);

    ISI = getIEI(stimIdxs);

    %% 2) Clean each artifact iteratively
    for idx = 1:numel(stimIdxs)
        % Identify the samples to clean
        hasReachedBaseline = false;
        searchOffset = blankingNSamples;

        while ~hasReachedBaseline && (stimIdxs(idx) + length(searchSamples) + searchOffset - 1) < length(signal)
            searchMedian = median(signal(searchSamples + stimIdxs(idx) + searchOffset - 1));

            if searchMedian > lowerBaseline && searchMedian < upperBaseline
                hasReachedBaseline = true;
            else
                searchOffset = searchOffset + 1;
            end
        end

        searchOffset = length(searchSamples) + searchOffset;

        if idx ~= numel(stimIdxs)
            if searchOffset < ISI(idx)
                artifactNSamples = searchOffset;
            else
                artifactNSamples = stimIdxs(idx + 1) - stimIdxs(idx);
            end
        else
            if hasReachedBaseline
                artifactNSamples = searchOffset;
            else
                artifactNSamples = length(signal) - stimIdxs(idx);
            end
        end

        artifactSamples = (1:artifactNSamples) + stimIdxs(idx) - 1;

        % Find the artifact shape
        artifact = fitArtifact(signal(artifactSamples), sampleRate, blankingPeriod);

        % Correct discontinuities
        correctionX = [-correctionNSamples:-1, (1:correctionNSamples) + length(artifact) - 1] + stimIdxs(idx);
        correctionY = output(correctionX);
        fullX = correctionX(1):correctionX(end);
        fullY = interp1(correctionX, correctionY, fullX, 'linear');

        % Update output signal
        output(artifactSamples) = output(artifactSamples) - artifact + fullY((correctionNSamples + 1):(length(fullY) - correctionNSamples));

        % Update progress bar
        waitbar(idx / numel(stimIdxs), waitbarFig, 'Removing artifacts...');
    end

    close(waitbarFig);
end