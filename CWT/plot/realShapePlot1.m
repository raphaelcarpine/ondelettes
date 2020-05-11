function realShapePlot1(shape, figTitle)
%COMPLEXSHAPEPLOTDEFAULT Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    warning('shape mising');
    shape = randn(1, 5);
end

if nargin < 2
    figTitle = '';
end

n = length(shape);

fig = figure('Name', figTitle);
ax = axes(fig);
stem(ax, shape);
xlim(ax, [0, n+1]);
xlabel(ax, 'dof');
ylabel(ax, '\phi');
xticks(ax, 1:n);

end

