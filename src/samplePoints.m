% function [x, y] = samplePoints(data)
%     minInterval = 3;
%     maxInterval = 12;
% 
%     intervalValues = minInterval:maxInterval;
%     % m = (1-0) / (maxInterval-minInterval);
%     % q = (maxInterval*0 - minInterval*1) / (maxInterval - minInterval);
%     % cumulative = m * intervalValues + q;
%     probabilities = pdf('Gamma', intervalValues, 3, 6);
%     cumulative = cumsum(probabilities);
% 
%     % figure();
%     % tiledlayout(2, 1)
%     % nexttile()
%     % plot(intervalValues, probabilities)
%     % nexttile()
%     % plot(intervalValues, cumulative)
% 
%     intervals = 0;
% 
%     while sum(intervals) < length(data)
%         newInterval = interp1(cumulative, cumulative, rand(1, 1), 'nearest', 0);
%         newInterval = intervalValues(find(cumulative == newInterval, 1));
% 
%         intervals = [intervals, newInterval];
%     end
% 
%     intervals = intervals(intervals > 0);
%     intervals = sort(intervals);
%     intervals = intervals(1:(end - 1));
% 
%     x = unique([1, cumsum(intervals), length(data)]);
%     y = data(x);
% 
% end

function x = samplePoints(data, sampleRate)
    minInterval = 2;
    maxInterval = round(2e-3 * sampleRate);

    mu = [3; 20];
    sigma = cat(3, 30, 10);
    weights = [0.35, 0.65];

    gmd = gmdistribution(mu, sigma, weights);

    rng(33);
    intervals = round(random(gmd, round(length(data) / min(mu))))';
    intervals = intervals(intervals >= minInterval & intervals <= maxInterval);
    intervals = intervals(cumsum(intervals) < length(data));
    intervals = sort(intervals);

    x = unique(cumsum(intervals));

end