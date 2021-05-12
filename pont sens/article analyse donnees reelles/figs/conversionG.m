%% lines Y

g = 9.81;

lines = findobj(gca, 'Type', 'line');

for k = 1:length(lines)
    lines(k).YData = g * lines(k).YData;
end

%% lines X

g = 9.81;

lines = findobj(gca, 'Type', 'line');

for k = 1:length(lines)
    lines(k).XData = g * lines(k).XData;
end

%% surf

g = 9.81;

surf = findobj(gca, 'Type', 'surface');

surf.CData = g * surf.CData;