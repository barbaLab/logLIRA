function [output, peakIdx, blankingSamples] = fitArtifact(data, sampleRate, blankingPeriod)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %% 0) Check input arguments and initialize
    data = double(data);
    blankingNSamples = round(blankingPeriod * sampleRate);

    %% 1) Find peakIdx and adjust it if clipped
    [peakIdx, isClipped, clippedSamples] = findStimulusPeak(data, sampleRate, blankingPeriod);

    if isempty(peakIdx)
        % No peak detect, basically it means that there is no clear
        % artifact.
        output = zeros(size(data));
        peakIdx = 1;
        blankingSamples = 1:blankingNSamples;
        return;
    end

    if isClipped
        peakIdx = clippedSamples(end);
    end

    %% 2) Pad data
    paddingDuration = 2e-3;    
    paddingNSamples = round(paddingDuration * sampleRate);

    paddingAfter = flip(data((end - paddingNSamples + 1):end));
    output = [data, paddingAfter];

    %% 2) Get the smooth data with the lowpass filter applied to the whole
    % raw data with correction of the boundary effects at the end of the
    % data. The boundary effect at the beginning are disregarded as the
    % filtered function is later substituted by the raw data during the
    % blanking period.

    splineX = exp(linspace(log(peakIdx), log(length(output)), length(output) / 8));
    splineX = sort(unique(rmmissing(interp1(peakIdx:length(output), peakIdx:length(output), splineX, 'previous'))));
    splineY = output(splineX);
    output = spline(splineX, splineY, 1:length(output));
   
    %% Remove padding and restore original data in blanking period
    output = output(1:(length(output) - paddingNSamples));

    if  blankingNSamples < peakIdx
        blankingNSamples = peakIdx;
    end

    blankingSamples = 1:blankingNSamples;

    output(blankingSamples) = data(blankingSamples);
    
end