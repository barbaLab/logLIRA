function [baseline, baselinePercentiles] = getBaselineFromStim(data, stim, sampleRate, baselineDuration, baselineOffset, percentiles)
%GETBASELINEFROMSTIM Summary of this function goes here
%   Detailed explanation goes here

baselineSamples = round(baselineDuration * sampleRate):-1:1;
baselineSamples = baselineSamples + round(baselineOffset * sampleRate);

baselineVector = zeros(length(stim), length(baselineSamples));

for idx=1:length(stim)
    baselineVector(idx, :) = data(stim(idx) - baselineSamples);
end

baselineVector = reshape(baselineVector', 1, []);

baseline = median(baselineVector);
baselinePercentiles = prctile(baselineVector, percentiles);
end

