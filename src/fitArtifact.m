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
    blankingPeriod = 1e-3;
    sFraction = 0.05;
    paddingDuration = 1e-3;

    validNumPosCheck = @(x) isnumeric(x) && (x >= 0);
    
    parser = inputParser;
    addRequired(parser, 'data', @isnumeric);
    addRequired(parser, 'sampleRate', validNumPosCheck);
    addOptional(parser, 'blankingPeriod', blankingPeriod, validNumPosCheck);
    addParameter(parser, 'sFraction', sFraction, @(x) isnumeric(x) && (x > 0) && (x <= 1));
    addParameter(parser, 'paddingDuration', paddingDuration, validNumPosCheck);
    addParameter(parser, 'saturationVoltage', [], @(x) isempty(x) || isnumeric(x));
    addParameter(parser, 'minClippedNSamples', [], @(x) isempty(x) || (isnumeric(x) && (x >= 0)));

    parse(parser, data, sampleRate, varargin{:});

    data = double(parser.Results.data);
    sampleRate = parser.Results.sampleRate;
    blankingPeriod = parser.Results.blankingPeriod;
    sFraction = parser.Results.sFraction;
    paddingDuration = parser.Results.paddingDuration;
    saturationVoltage = parser.Results.saturationVoltage;
    minClippedNSamples = parser.Results.minClippedNSamples;

    blankingNSamples = round(blankingPeriod * sampleRate);

    %% 1) Find peakIdx and adjust it if clipped
    [peakIdx, isClipped, clippedSamples] = findArtifactPeak(data, sampleRate, blankingPeriod, saturationVoltage, minClippedNSamples);

    if isempty(peakIdx)
        % No peak detect, it implies that there is no artifact.
        varargout{1} = 1;
        varargout{2} = 1:blankingNSamples;
        artifact = zeros(size(data));
        artifact(1:blankingNSamples) = data(1:blankingNSamples);
    else
        if isClipped && blankingNSamples < clippedSamples(end)
            blankingNSamples = clippedSamples(end);
        end
    
        %% 2) Pad data 
        paddingNSamples = round(paddingDuration * sampleRate);
    
        paddingAfter = flip(data((end - paddingNSamples + 1):end));
        output = [data, paddingAfter];
    
        %% 3) Extract the artifact shape
        startInterpX = 1;
        interpX = exp(linspace(log(startInterpX), log(length(output)), round(length(output) * sFraction)));
        interpX = rmmissing(interp1(startInterpX:length(output), startInterpX:length(output), interpX, 'previous'));
        interpX = sort(unique([interpX, blankingNSamples + 1, length(output) - paddingNSamples]));
        interpY = output(interpX);
        output = interp1(interpX, interpY, 1:length(output), 'linear');
        % output = spline(interpX, interpY, 1:length(output));
        
        % [interpX, interpY] = samplePoints(output(blankingNSamples:end));
        % output = interp1(interpX + blankingNSamples - 1, interpY, 1:length(output), 'linear');
       
        %% 4) Remove padding and restore original data in the blanking period
        output = output(1:(length(output) - paddingNSamples));
    
        blankingSamples = 1:blankingNSamples;
        output(blankingSamples) = data(blankingSamples);
    
        %% 5) Return output values
        artifact = output;
        varargout{1} = peakIdx;
        varargout{2} = blankingSamples;
    end

    %% 6) Plot
    % t = 0:1/sampleRate:(length(data)/sampleRate - 1/sampleRate);
    % t = t*1e3;
    % 
    % fig = figure();
    % tiledlayout(2, 1);
    % 
    % nexttile()
    % hold('on');
    % plot(t, data);
    % plot([0, 0], [min(data), max(data)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1)
    % plot([1, 1], [min(data), max(data)], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1)
    % plot(t, output, 'Color', 'magenta')
    % scatter(1e3*(peakIdx/sampleRate - 1/sampleRate), output(peakIdx), 25, 'black', 'Marker', '*');
    % title('Raw Data');
    % xlabel('Time (ms)');
    % ylabel('Voltage (\mu{V})');
    % 
    % residuals = data-artifact;
    % nexttile()
    % hold('on')
    % plot([0, 0], [min(residuals), max(residuals)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1)
    % plot([1, 1], [min(residuals), max(residuals)], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1)
    % plot(t, residuals, 'Color', 'b')
    % scatter(1e3*(peakIdx/sampleRate - 1/sampleRate), residuals(peakIdx), 25, 'black', 'Marker', '*');
    % title('Residuals');
    % xlabel('Time (ms)');
    % ylabel('Voltage (\mu{V})');
    % uiwait(fig);
    
end