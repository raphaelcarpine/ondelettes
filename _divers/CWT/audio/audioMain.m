function audioMain(dataPlotOrAxis)
%AUDIOMAIN Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0 % test
    Fs = 1000;
    T = (0:10*Fs) / Fs;
    X = 0.2 * randn(2, length(T));
    ax = axes(figure);
    dataPlotOrAxis = plot(ax, T, X);
end


%% data
plts = findLinesMenu(dataPlotOrAxis, 'Audio data selection menu');

if isempty(plts)
    warning('no line selected');
    return
end

%% audio
audioPlayer(plts);


end

