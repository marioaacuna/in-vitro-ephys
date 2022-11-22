function outData = movingAverage(inData, windowSize)
% use a boxcar filter of length windowSize points on inData
% filteredData = movingAverage(rawData, windowSize);
% defaults:
%   windowSize = 10 points
if nargin < 2
    windowSize = 10;
end
if size(inData, 1) > size(inData, 2)
    longSide = 1;
    flatData = ones(windowSize, 1);
else
    longSide = 2;
    flatData = ones(1, windowSize);
end

cheatShift = int32(windowSize / 2);
outData = filter(flatData./(windowSize),1,cat(longSide, flatData.*inData(1), inData, flatData.*inData(end)));
outData = outData(windowSize + cheatShift:length(inData) + windowSize + cheatShift - 1);
end