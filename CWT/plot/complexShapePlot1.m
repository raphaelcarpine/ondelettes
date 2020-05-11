function complexShapePlot1(shape, figTitle)
%COMPLEXSHAPEPLOTDEFAULT Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    warning('shape mising');
    shape = randn(1, 5) + 0.1i*randn(1, 5);
end

if nargin < 2
    figTitle = '';
end

fig = figure('Name', figTitle);
ax = polaraxes(fig);
hold(ax, 'on');
legendNames = {};
for k = 1:length(shape)
    polarplot(ax, [0, shape(k)], '-o');
    legendNames{end+1} = ['dof', num2str(k)];
end

legend(ax, legendNames);

end

