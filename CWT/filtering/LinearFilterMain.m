function LinearFilterMain(dataPlotOrAxes)
%LINEARFILTERMAIN Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    dataPlotOrAxes = gca;
end

dataPlot = findLinesMenu(dataPlotOrAxes, 'Filter data selection menu');

if isempty(dataPlot)
    warning('no data selected');
    return
end

LinearFilterMenu(dataPlot);

end

