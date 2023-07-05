function output = FraAlgorithm(sig, p)
    output = sig;
    nStims = length(p.StimI);
    blankingNSamples = p.fs * 1e-3;

    ISI = getIEI(p.StimI);

    [~, baselinePercentiles] = getBaselineFromStim(sig, p.StimI, p.fs, 10e-3, 0.5e-3, [25, 75]);
    baselineBottom = baselinePercentiles(1);
    baselineTop = baselinePercentiles(end);

    stimNSamples = zeros(1, nStims);

    searchingWindow = 3e-3;
    searchingSamples = 1:round(searchingWindow * p.fs);

    waitbarFig = waitbar(0);

    for idx=1:length(p.StimI)
        % Identify the samples for each artifact
        flag = false;
        searchingOffset = blankingNSamples;

        while ~flag && (p.StimI(idx) + length(searchingSamples) - 1 + searchingOffset) < length(sig)
            searchingVector = sig(p.StimI(idx) + (searchingSamples - 1) + searchingOffset);

            if median(searchingVector) > baselineBottom && median(searchingVector) < baselineTop
                flag = true;
            else
                searchingOffset = searchingOffset + 1;
            end
        end

        searchingOffset = searchingOffset + length(searchingSamples);

        if idx ~= nStims
            if searchingOffset < ISI(idx)
                stimNSamples(idx) = searchingOffset;
            else
                stimNSamples(idx) = p.StimI(idx + 1) - p.StimI(idx);
            end
        else
            if flag
               stimNSamples(idx) = searchingOffset;
            else
                stimNSamples(idx) = length(sig) - p.StimI(idx);
            end
        end

        artifactSamples = (1:stimNSamples(idx)) - 1 + p.StimI(idx);

        % Remove the artifacts
        [data, ~, ~] = fitArtifact(sig(artifactSamples), p.fs, 1e-3, 200);

        % Correct discontinuities for blanked samples
        correctionDuration = 3e-3;
        correctionNSamples = round(correctionDuration * p.fs);
        correctionX = [-correctionNSamples:-1, (1:correctionNSamples) - 1 + length(data)] + p.StimI(idx);
        correctionY = output(correctionX);
        fullX = [-correctionNSamples:-1, (1:length(data)) - 1, (1:correctionNSamples) - 1 + length(data)] + p.StimI(idx);
        fullY = interp1(correctionX, correctionY, fullX, 'linear');

        % Update the output signal
        output(artifactSamples) = output(artifactSamples) - data + fullY((correctionNSamples + 1):(length(fullY) - correctionNSamples));

        % Update progress bar
        waitbar(idx / nStims, waitbarFig, 'Removing artifacts...');
    end

    close(waitbarFig);
end