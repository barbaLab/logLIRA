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
    magicNPoints = 600;

    validNumPosCheck = @(x) isnumeric(x) && (x >= 0);
    
    parser = inputParser;
    addRequired(parser, 'data', @isnumeric);
    addRequired(parser, 'sampleRate', validNumPosCheck);
    addOptional(parser, 'blankingPeriod', blankingPeriod, validNumPosCheck);
    addParameter(parser, 'sFraction', sFraction, @(x) isnumeric(x) && (x > 0) && (x <= 1));
    addParameter(parser, 'saturationVoltage', [], @(x) isempty(x) || isnumeric(x));
    addParameter(parser, 'minClippedNSamples', [], @(x) isempty(x) || (isnumeric(x) && (x >= 0)));

    parse(parser, data, sampleRate, varargin{:});

    data = double(parser.Results.data);
    sampleRate = parser.Results.sampleRate;
    blankingPeriod = parser.Results.blankingPeriod;
    sFraction = parser.Results.sFraction;
    saturationVoltage = parser.Results.saturationVoltage;
    minClippedNSamples = parser.Results.minClippedNSamples;

    blankingNSamples = round(blankingPeriod * sampleRate);
    output = data;

    %% 1) Find peakIdx and adjust it if clipped
    peakIdx = findArtifactPeak(data, sampleRate, blankingPeriod, saturationVoltage, minClippedNSamples);
    blankingNSamples = max([blankingNSamples, peakIdx]);

    %% 2) Select interpolating points and extract the artifact shape
    startInterpXOffset = round(0.25 * 1e-3 * sampleRate);
    startInterpX = blankingNSamples - startInterpXOffset;
    interpX = exp(linspace(log(startInterpX), log(magicNPoints), round(magicNPoints * sFraction)));
    interpX = unique(round(interpX));

    interpX(interpX > length(output)) = [];
    
    largestIPI = interpX(end) - interpX(end-1);
    nExtraInterpX = floor((length(output) - interpX(end)) / largestIPI);
    if nExtraInterpX > 0
        interpX = [interpX, interpX(end) + largestIPI * (1:nExtraInterpX)];
    end

    IPI = [interpX(1), diff(interpX), length(output) - interpX(end)];

    minHalfInterval = 2;
    maxHalfInterval = 15;
    interpY = zeros(1, numel(interpX));
    for i = 1:numel(interpY)
        if floor(IPI(i) / 2) >= minHalfInterval && floor(IPI(i + 1) / 2) >= minHalfInterval
            intervalSamples = -min(maxHalfInterval, floor(IPI(i) / 2)):min(maxHalfInterval, floor(IPI(i + 1) / 2));
        else
            intervalSamples = 0;
        end
        interpY(i) = mean(output(intervalSamples + interpX(i)));
    end

    keyX = [blankingNSamples + 1, length(output)];
    keyY = output(keyX);

    [interpX, keptIdxs, ~] = unique([keyX, interpX]);   % Unique automatically sorts
    interpY = [keyY, interpY];
    interpY = interpY(keptIdxs);
    
    output = interp1(interpX, interpY, 1:length(output), 'linear');

    %% 3) Restore original data in the blanking period
    blankingSamples = 1:blankingNSamples;
    output(blankingSamples) = data(blankingSamples);

    %% 4) Return output values
    artifact = output;
    varargout{1} = blankingNSamples;
    varargout{2} = peakIdx;

    %% 5) Plot
    % t = 0:1/sampleRate:(length(data)/sampleRate - 1/sampleRate);
    % t = t*1e3;

    % fig = figure();
    % tiledlayout(2, 1);

    % ax = nexttile();
    % hold('on');
    % plot(t, data);
    % plot(t, artifact, 'Color', 'magenta')
    % scatter(1e3*(peakIdx/sampleRate - 1/sampleRate), artifact(peakIdx), 25, 'black', 'Marker', '*');
    % patch([0, blankingPeriod, blankingPeriod, 0] * 1e3, [min(ax.YLim), min(ax.YLim), max(ax.YLim), max(ax.YLim)], [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'LineStyle', 'none');
    % title('Raw Data');
    % xlabel('Time (ms)');
    % ylabel('Voltage (\mu{V})');

    % residuals = data-artifact;
    % ax = nexttile();
    % hold('on')
    % plot(t, residuals, 'Color', 'b')
    % scatter(1e3*(peakIdx/sampleRate - 1/sampleRate), residuals(peakIdx), 25, 'black', 'Marker', '*');
    % patch([0, blankingPeriod, blankingPeriod, 0] * 1e3, [min(ax.YLim), min(ax.YLim), max(ax.YLim), max(ax.YLim)], [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'LineStyle', 'none');
    % title('Residuals');
    % xlabel('Time (ms)');
    % ylabel('Voltage (\mu{V})');
    % set(gcf,'Visible','on');
    % uiwait(fig);

end