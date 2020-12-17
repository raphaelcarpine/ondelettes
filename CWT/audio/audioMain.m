function audioMain(ax)
%AUDIOMAIN Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0 % test
    Fs = 1000;
    T = (0:10*Fs) / Fs;
    X = 0.2 * randn(2, length(T));
    ax = axes(figure);
    plot(ax, T, X);
end


%% data
ax = findAxesAudio(ax);

linesAx = findobj(ax, 'Type', 'line');
if isempty(linesAx)
    warning('no line on axis');
    return
end
T = get(linesAx(1), 'XData');
X = get(linesAx(1), 'YData');
for k_l = 2:length(linesAx)
    if isequal(T, get(linesAx(k_l), 'XData'))
        X = [X; get(linesAx(k_l), 'YData')];
    end
end

%% axes
linkedAxes = [ax, findLinkedAxesAudio()];
[initFcn, updateFcn, closeFcn] = audioTimeOnAxes(linkedAxes);

%% audio
audioPlayer(T, X, initFcn, updateFcn, closeFcn);


end

