function fig = WaveletMenu(waveletPlot, varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
defaultParent = figure;
checkParent = @(f) isa(f, 'matlab.ui.Figure') || isa(f, 'matlab.ui.container.Panel')...
    || isa(f, 'matlab.ui.container.Tab') || isa(f, 'matlab.ui.container.ButtonGroup');
addRequired(p,'waveletPlot');
addParameter(p,'Parent', defaultParent, checkParent);
parse(p, waveletPlot, varargin{:})
fig = p.Results.Parent;



end

