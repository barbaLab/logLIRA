function [output, varargout] = logssar(signal, stimIdxs, sampleRate, varargin)
%LOGSSAR LOGarithmically distributed intervals Spline Stimulus Artifacts Rejection.
%   output = LOGSSAR(signal, stimIdxs, sampleRate) returns the input signal
%   without the artifacts caused by electrical stimulation. The stimIdxs
%   are the indexes of stimulation onsets. The sampleRate should be expressed
%   in Hz.
%
%   [output, blankingPeriods] = LOGSSAR(signal, stimIdxs, sampleRate) returns a vector
%   containing all the blanking periods determined by the algorithm for each stimulus
%   onset, in samples.
%
%   output = LOGSSAR(..., blankingPeriod) specifies the minimum time after the
%   stimulus onset that is discarded. It must be expressed in seconds. By
%   default it is 1 ms.
%
%   [...] = LOGSSAR(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%   parameter name/value pairs. Parameters are:
%
%       'SaturationVoltage' - It specifies the recording system operating range
%                             in mV as specified in the datasheet. This is useful
%                             to properly detect saturation. Choices are:
%                   default - 95% of the input signal absolute value maximum.
%                1x1 scalar - The operating range is assumed to be symmetric with
%                             respect to 0.
%          1x2 or 2x1 array - The operating range is the specified one.
%
%      'minClippedNSamples' - It is the minimum number of consecutive clipped samples
%                             to mark the artifact as a clipped one. It should be a
%                             1x1 positive integer. By default, it is 2.

    %% 0) Check and parse input arguments
    warning('off', 'signal:findpeaks:largeMinPeakHeight');
    warning('off', 'stats:kmeans:FailedToConvergeRep');

    validNumPosCheck = @(x) isnumeric(x) && (x >= 0);
    
    parser = inputParser();
    addRequired(parser, 'signal', @isnumeric);
    addRequired(parser, 'stimIdxs', @(x) isnumeric(x) && all(x > 0));
    addRequired(parser, 'sampleRate', validNumPosCheck);
    addOptional(parser, 'blankingPeriod', 1e-3, validNumPosCheck);
    addParameter(parser, 'saturationVoltage', 0.95 * max(abs(signal)) / 1e3, @isnumeric);
    addParameter(parser, 'minClippedNSamples', [], validNumPosCheck);
    addParameter(parser, 'randomSeed', 'default', @(x) x >= 0);

    parse(parser, signal, stimIdxs, sampleRate, varargin{:});

    signal = parser.Results.signal;
    stimIdxs = parser.Results.stimIdxs;
    sampleRate = parser.Results.sampleRate;
    blankingPeriod = parser.Results.blankingPeriod;
    saturationVoltage = parser.Results.saturationVoltage;
    minClippedNSamples = parser.Results.minClippedNSamples;
    randomSeed = parser.Results.randomSeed;

    rng(randomSeed);

    output = signal;
    varargout{1} = zeros(size(stimIdxs));
    varargout{2} = false(size(stimIdxs));

    waitbarFig = waitbar(0, 'Starting...');

    %% 1) Find signal IAI and check if artifacts requires correction
    minArtifactDuration = 0.04;
    FPRemovalDuration = 0.002;
    checkDuration = 0.005;
    checkThreshold = 30;

    blankingNSamples = round(blankingPeriod * sampleRate);
    IAI = [diff(stimIdxs), length(signal) - stimIdxs(end)];
    

    checkNSamples = round(checkDuration * sampleRate);
    checkSamples = repmat(0:(checkNSamples - 1), [1, numel(stimIdxs)]);
    artifactSamples = reshape(repmat(stimIdxs, [checkNSamples, 1]), 1, []);
    
    preArtifacts = signal(artifactSamples - flip(checkSamples) - 1);
    preArtifacts = reshape(preArtifacts, checkNSamples, []);
    postArtifacts = signal(artifactSamples + checkSamples + blankingNSamples);
    postArtifacts = reshape(postArtifacts, checkNSamples, []);

    hasArtifact = (abs(preArtifacts(end, :) - postArtifacts(1, :)) > checkThreshold) | ...
                    std(postArtifacts, 0, 1) >= 3 * std(preArtifacts, 0, 1) | ...
                    (blankingNSamples >= IAI);

    FPRemovalNSamples = round(FPRemovalDuration * sampleRate);
    FPRemovalData = zeros(numel(stimIdxs), FPRemovalNSamples);
    FPRemovalSamples = zeros(numel(stimIdxs), FPRemovalNSamples);

    %% 2) Clean each artifact iteratively
    minArtifactNSamples = round(minArtifactDuration * sampleRate) + blankingNSamples;

    for idx = 1:numel(stimIdxs)
        % Identify samples to clean
        data = signal((1:IAI(idx)) + stimIdxs(idx) - 1);

        if hasArtifact(idx)
            endIdx = [];
            if minArtifactNSamples < IAI(idx) 
                smoothData = smoothdata(data(minArtifactNSamples:end), 'movmean', round(5 * 1e-3 * sampleRate));
                endIdx = find(abs(smoothData - median(data(minArtifactNSamples:end))) < 1, 1) + minArtifactNSamples - 1;
            end
                
            endIdx = min([IAI(idx), endIdx]);

            % Find artifact shape
            [artifact, blankingNSamples] = fitArtifact(data(1:endIdx), sampleRate, blankingPeriod, ...
                'saturationVoltage', saturationVoltage, 'minClippedNSamples', minClippedNSamples);
        else
            blankingNSamples = round(blankingPeriod * sampleRate);
            artifact = data(1:blankingNSamples);
        end

        if ~isempty(blankingNSamples) && (length(artifact) - FPRemovalNSamples) > blankingNSamples
            % Get data for false positives removal at artifact beginning
            FPRemovalSamples(idx, :) = (1:FPRemovalNSamples) + blankingNSamples;
            FPRemovalData(idx, :) = data(FPRemovalSamples(idx, :)) - artifact(FPRemovalSamples(idx, :));
            varargout{1}(idx) = blankingNSamples;
        elseif ~hasArtifact(idx)
            varargout{1}(idx) = blankingNSamples;
        else
            artifact = data;
            varargout{1}(idx) = length(artifact);
            varargout{2}(idx) = true;
        end

        % Correct artifact to avoid discontinuities
        if ~hasArtifact(idx) || IAI(idx) > endIdx
            correctionX = [0, length(artifact) + 1];
            correctionY = [output(correctionX(1) + stimIdxs(idx) - 1), output(correctionX(end) + stimIdxs(idx) - 1)];
            correction = interp1(correctionX, correctionY, 1:length(artifact), 'linear');
        else
            correction = output(stimIdxs(idx) - 1) * ones(1, length(artifact));
        end

        % Update output signal
        output((1:length(artifact)) + stimIdxs(idx) - 1) = data(1:length(artifact)) - artifact + correction;

        % Update progress bar
        waitbar(idx / numel(stimIdxs), waitbarFig, 'Removing artifacts...');
    end

    varargout{2} = find(varargout{2} == true);
    if ~isempty(varargout{2})
        warning('logssar:logssar:skippedTrials', 'Some trials were skipped and blanked completely: %d/%d.', numel(varargout{2}), numel(stimIdxs));
    end

    %% 3) Remove false positives at artifacts beginning
    waitbar(0, waitbarFig, 'Checking signal...');
    
    minClusterSize = 50;
    explainedThreshold = 70;
    nRepetitions = 1;
    KList = 2:6;
    cutoffDistance = 0.8;

    warning('off', 'stats:pca:ColRankDefX');
    [FPRemovalDataReduced, ~, ~, ~, explained] = pca(FPRemovalData');
    FPRemovalDataReduced = FPRemovalDataReduced(:, 1:max([2, find(cumsum(explained) >= explainedThreshold, 1)]));
    warning('on', 'stats:pca:ColRankDefX');

    consensusMatrix = zeros([numel(stimIdxs), numel(stimIdxs)]);
    indicatorMatrix = 0;
    idxs = 1:numel(stimIdxs);

    nRuns = nRepetitions * numel(KList);
    runIdx = 0;

    for repetitionIdx = 1:nRepetitions
        for idx = 1:numel(KList)
            rng((randomSeed + KList(idx)) * repetitionIdx);

            % KMeans Clustering
            labels = kmeans(FPRemovalDataReduced, KList(idx), 'Replicates', 5);
            for clusterIdx = 1:numel(unique(labels))
                [rows, cols] = meshgrid(idxs(labels == clusterIdx), idxs(labels == clusterIdx));
                linearIdxs = sub2ind([numel(idxs), numel(idxs)], rows(:), cols(:));
                consensusMatrix(linearIdxs) = consensusMatrix(linearIdxs) + 1;
            end
            indicatorMatrix = indicatorMatrix + 1;

            % Spectral Clustering
            labels = spectralcluster(FPRemovalDataReduced, KList(idx));
            for clusterIdx = 1:numel(unique(labels))
                [rows, cols] = meshgrid(idxs(labels == clusterIdx), idxs(labels == clusterIdx));
                linearIdxs = sub2ind([numel(idxs), numel(idxs)], rows(:), cols(:));
                consensusMatrix(linearIdxs) = consensusMatrix(linearIdxs) + 1;
            end
            indicatorMatrix = indicatorMatrix + 1;

            runIdx = runIdx + 1;
            waitbar(runIdx / nRuns, waitbarFig, 'Checking signal...');
        end
    end

    rng(randomSeed);
    consensusMatrix = consensusMatrix ./ indicatorMatrix;
    labels = cluster(linkage(1 - consensusMatrix, 'average'), 'cutoff', cutoffDistance);
    [clusterSize, ~] = histcounts(labels, unique(labels));
    nClusters = max([1, sum(clusterSize >= minClusterSize)]);

    warning('off', 'stats:gmdistribution:FailedToConvergeReps');
    warning('off', 'stats:gmdistribution:IllCondCov');
    GMModel = [];

    while isempty(GMModel) && nClusters > 0
        try
            GMModel = fitgmdist(FPRemovalDataReduced, nClusters, 'Replicates', 5, 'Options', statset('MaxIter', 1000));
        catch
            nClusters = nClusters - 1;
            GMModel = [];
        end
    end
    warning('on', 'stats:gmdistribution:FailedToConvergeReps');
    warning('on', 'stats:gmdistribution:IllCondCov');

    labels = GMModel.cluster(FPRemovalDataReduced);
    for clusterIdx = 1:numel(unique(labels))
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
    end

    close(waitbarFig);
    warning('on', 'all');
end