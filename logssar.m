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
    warning('off');

    validNumPosCheck = @(x) isnumeric(x) && (x >= 0);
    
    parser = inputParser();
    addRequired(parser, 'signal', @isnumeric);
    addRequired(parser, 'stimIdxs', @(x) isnumeric(x) && all(x > 0));
    addRequired(parser, 'sampleRate', validNumPosCheck);
    addOptional(parser, 'blankingPeriod', 1e-3, validNumPosCheck);
    addParameter(parser, 'saturationVoltage', [], @isnumeric);
    addParameter(parser, 'minClippedNSamples', [], validNumPosCheck);

    parse(parser, signal, stimIdxs, sampleRate, varargin{:});

    signal = parser.Results.signal;
    stimIdxs = parser.Results.stimIdxs;
    sampleRate = parser.Results.sampleRate;
    blankingPeriod = parser.Results.blankingPeriod;
    saturationVoltage = parser.Results.saturationVoltage;
    minClippedNSamples = parser.Results.minClippedNSamples;

    output = signal;

    waitbarFig = waitbar(0, 'Starting...');

    %% 1) Find signal IAI, minimum artifact length, and check if artifacts requires correction
    IAI = [diff(stimIdxs), length(signal) - stimIdxs(end)];
    minArtifactDuration = 0.04;
    minArtifactNSamples = min([min(IAI), round(minArtifactDuration * sampleRate)]);

    blankingNSamples = round(blankingPeriod * sampleRate);
    artifactSamples = reshape(repmat(stimIdxs, [minArtifactNSamples, 1]), 1, []) + repmat(0:(minArtifactNSamples - 1), [1, numel(stimIdxs)]);
    artifact = signal(artifactSamples);
    artifact = reshape(artifact, minArtifactNSamples, []);
    artifact = mean(artifact, 2)';
    artifact = artifact((blankingNSamples + 1):end);
    isArtifact = abs(mean(artifact(1:blankingNSamples)) - mean(artifact(end - (1:blankingNSamples) + 1))) > 10;

    if isArtifact
        FPRemovalDuration = 0.002;
        FPRemovalNSamples = round(FPRemovalDuration * sampleRate);
        FPRemovalData = zeros(numel(stimIdxs), FPRemovalNSamples);
        FPRemovalSamples = zeros(numel(stimIdxs), FPRemovalNSamples);

        nSkippedTrials = 0;

        %% 2) Clean each artifact iteratively
        for idx = 1:numel(stimIdxs)
            % Identify the samples to clean
            data = signal((1:IAI(idx)) + stimIdxs(idx) - 1);

            minArtifactNSamples = min([IAI(idx), round(minArtifactDuration * sampleRate)]);
            derivative = diff(data(minArtifactNSamples:end), 2);
            derivative = smoothdata(derivative, 'movmedian', round(2 * 1e-3 * sampleRate));
            endIdx = find(derivative >= -0.2 & derivative <= 0.2, 1) + minArtifactNSamples;
            
            data = data(1:min([endIdx, length(data)]));

            % Find the artifact shape
            [artifact, blankingNSamples] = fitArtifact(data, sampleRate, blankingPeriod, ...
                'saturationVoltage', saturationVoltage, 'minClippedNSamples', minClippedNSamples);

            if ~isempty(blankingNSamples) && (length(artifact) - FPRemovalNSamples) > blankingNSamples
                % Get data for false positives removal at artifact beginning
                FPRemovalSamples(idx, :) = (1:FPRemovalNSamples) + blankingNSamples;
                FPRemovalData(idx, :) = data(FPRemovalSamples(idx, :)) - artifact(FPRemovalSamples(idx, :));
            else
                nSkippedTrials = nSkippedTrials + 1;
            end

            % Correct discontinuities
            correctionX = [0, length(artifact) + 1];
            correctionY = [signal(correctionX(1) + stimIdxs(idx) - 1), signal(correctionX(end) + stimIdxs(idx) - 1)];
            correction = interp1(correctionX, correctionY, 1:length(artifact), 'linear');

            % Update output signal
            output((1:length(artifact)) + stimIdxs(idx) - 1) = data - artifact + correction;

            % Update progress bar
            waitbar(idx / numel(stimIdxs), waitbarFig, 'Removing artifacts...');
        end

        if nSkippedTrials > 0
            warning('logssar:logssar:skippedTrials', 'Some trials were skipped and blanked completely: %d/%d.', nSkippedTrials, numel(stimIdxs));
        end

        %% 3) Remove false positives at artifacts beginning
        waitbar(0, waitbarFig, 'Checking signal...');
        minClusterSize = 50;
        FPRemovalDataReduced = pca(FPRemovalData', 'NumComponents', 3);
        FPClustersEvaluation = evalclusters(FPRemovalDataReduced, 'kmeans', 'silhouette', 'KList', 1:floor(numel(stimIdxs) / (2 * minClusterSize)));
        labels = FPClustersEvaluation.OptimalY;
        nClusters = FPClustersEvaluation.OptimalK;
                    
        for clusterIdx = 1:nClusters
            if sum(labels == clusterIdx) >= minClusterSize
                selectedFPRemovalSamples = FPRemovalSamples(labels == clusterIdx, :) + stimIdxs(labels == clusterIdx)' - 1;
                selectedFPRemovalSamples = reshape(selectedFPRemovalSamples', [1, numel(selectedFPRemovalSamples)]);

                % fig = figure();
                % tiledlayout(3, 1);
                % nexttile();
                % plot(reshape(output(selectedFPRemovalSamples), [], sum(labels == clusterIdx)));
                % nexttile();
                % plot(mean(FPRemovalData(labels == clusterIdx, :), 1));
                % nexttile();
                % plot(reshape(output(selectedFPRemovalSamples) - repmat(mean(FPRemovalData(labels == clusterIdx, :), 1), [1, sum(labels == clusterIdx)]), [], sum(labels == clusterIdx)));
                % uiwait(fig);

                output(selectedFPRemovalSamples) = output(selectedFPRemovalSamples) - repmat(mean(FPRemovalData(labels == clusterIdx, :), 1), [1, sum(labels == clusterIdx)]);
            end

            waitbar(clusterIdx / nClusters, waitbarFig, 'Checking signal...');
        end
    else
        artifact = zeros(1, length(signal));
        blankingSamples = reshape(repmat(stimIdxs, [blankingNSamples, 1]), 1, []) + repmat(0:(blankingNSamples - 1), [1, numel(stimIdxs)]);
        artifact(blankingSamples) = signal(blankingSamples);
        
        for idx = 1:numel(stimIdxs)
            correctionX = [0, blankingNSamples + 1];
            correctionY = signal(correctionX  + stimIdxs(idx) - 1);
            correction = interp1(correctionX, correctionY, 1:blankingNSamples, 'linear');
            artifact((1:blankingNSamples) + stimIdxs(idx) - 1) = artifact((1:blankingNSamples) + stimIdxs(idx) - 1) - correction;

            % Update progress bar
            waitbar(idx / numel(stimIdxs), waitbarFig, 'Removing artifacts...');
        end

        output = output - artifact;
    end

    close(waitbarFig);
    warning('on');
end