function Axes = findAxesAudio(Fs)
%FINDAXESAUDIO Summary of this function goes here
%   Detailed explanation goes here

allAxes = findall(0, 'type', 'axes');

Axes = [];
for kax = 1:length(allAxes)
    ax = allAxes(kax);
    allLines = findall(ax, 'type', 'line');
    for kline = 1:length(allLines)
        xdata = get(allLines(kline), 'XData');
        dx = mean(diff(xdata));
        if (Fs*dx - 1) < 1e-3 && all(abs(diff(xdata)/dx - 1) < 1e-3)
            Axes(end+1) = ax;
            continue
        end
    end
end

end

