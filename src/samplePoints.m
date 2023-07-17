function [x, y] = samplePoints(data)
    minInterval = 3;
    maxInterval = 12;

    intervalValues = minInterval:maxInterval;
    % m = (1-0) / (maxInterval-minInterval);
    % q = (maxInterval*0 - minInterval*1) / (maxInterval - minInterval);
    % cumulative = m * intervalValues + q;
    probabilities = pdf('Gamma', intervalValues, 3, 6);
    cumulative = cumsum(probabilities);

    % figure();
    % tiledlayout(2, 1)
    % nexttile()
    % plot(intervalValues, probabilities)
    % nexttile()
    % plot(intervalValues, cumulative)

    intervals = 0;

    while sum(intervals) < length(data)
        newInterval = interp1(cumulative, cumulative, rand(1, 1), 'nearest', 0);
        newInterval = intervalValues(find(cumulative == newInterval, 1));
        
        intervals = [intervals, newInterval];
    end

    intervals = intervals(intervals > 0);
    intervals = sort(intervals);
    intervals = intervals(1:(end - 1));

    x = unique([1, cumsum(intervals), length(data)]);
    y = data(x);

end