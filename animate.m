function animate(object, t, x)
%animate Summary of this function goes here
%   Detailed explanation goes here
fig = figure;
ax = fig.CurrentAxes;
xlim([-0.5 1.5]);
ylim([-0.5 0.5]);

p1 = plot(ax, 0, 0, 'o');
p2 = plot(ax, 1, 0, 'o');
end

