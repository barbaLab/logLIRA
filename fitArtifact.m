function [artifact, varargout] = fitArtifact(data, sampleRate, varargin)
%FITARTIFACT Fit the artifact shape.
%   artifact = FITARTIFACT(data, sampleRate) computes the stimulus artifact
%   shape from the input data, discarding all the high frequency components
%   usually related to spiking activity. By default, a 1 ms blanking period
%   is assumed after the stimulus onset. The input data are expected to
%   start from the stimulus onset.
%
%   [artifact, peakIdx] = FITARTIFACT(data, sampleRate) returns
%   the index where the peak of the artifact is found.
%
%   [artifact, peakIdx, blankingSamples] = FITARTIFACT(data, sampleRate) returns
%   the blanking samples, where the input data were not modified.
%
%   [...] = FITARTIFACT(..., blankingPeriod) specifies the time after the
%   stimulus onset that is not affected by this function. It must be
%   expressed in seconds. By default it is 1 ms.
%
%   [...] = FITARTIFACT(..., blankingPeriod, sFraction) specifies the fraction
%   of data points to use for spline interpolation. It is expected to be a value
%   between 0 and 1. By default it is 0.1.
%
%   [...] = FITARTIFACT(..., blankingPeriod, sFraction, paddingDuration) specifies
%   the padding that must be addedd at the end of the signal to compensate for
%   boundary effects. It must be expressed in seconds. By default it is 1 ms.

    %% 0) Check and parse input arguments
    if nargin < 2
        throw(MException('SAR:NotEnoughParameters', 'The parameters data and sampleRate are required.'));
    end
    
    if nargin < 3
        blankingPeriod = 1e-3;
    else
        blankingPeriod = varargin{1};
    end

    if nargin < 4
        sFraction = 0.1;    % We want to interpolate 10% of the overall points
    else
        sFraction = varargin{2};
    end

    if nargin < 5
        paddingDuration = 1e-3;
    else
        paddingDuration = varargin{3};
    end

    data = double(data);
    blankingNSamples = round(blankingPeriod * sampleRate);

    %% 1) Find peakIdx and adjust it if clipped
    [peakIdx, isClipped, clippedSamples] = findArtifactPeak(data, sampleRate, blankingPeriod);

    if isempty(peakIdx)
        % No peak detect, it implies that there is no artifact.
        artifact = zeros(size(data));
        varargout{1} = 1;
        varargout{2} = 1:blankingNSamples;
        return;
    end

    if isClipped
        peakIdx = clippedSamples(end);
    end

    %% 2) Pad data 
    paddingNSamples = round(paddingDuration * sampleRate);

    paddingAfter = flip(data((end - paddingNSamples + 1):end));
    output = [data, paddingAfter];

    %% 3) Extract the artifact shape
    splineX = exp(linspace(log(peakIdx), log(length(output)), round(length(output) * sFraction)));
    splineX = rmmissing(interp1(peakIdx:length(output), peakIdx:length(output), splineX, 'previous'));
    splineX = sort(unique([splineX, length(output) - paddingNSamples]));
    splineY = output(splineX);
    output = spline(splineX, splineY, 1:length(output));
   
    %% 4) Remove padding and restore original data in the blanking period
    output = output(1:(length(output) - paddingNSamples));

    if  blankingNSamples < peakIdx
        blankingNSamples = peakIdx;
    end

    blankingSamples = 1:blankingNSamples;
    output(blankingSamples) = data(blankingSamples);

    %% 5) Return output values
    artifact = output;
    varargout{1} = peakIdx;
    varargout{2} = blankingSamples;
    
end