function fig = shapePlotPlate(plateDim, dofPos, shape, figTitle, displayChannels)
%COMPLEXSHAPEPLOTDEFAULT Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    warning('shape mising');
    plateDim = [3, 1];
    dofPos = rand(5, 2) .* plateDim;
    shape = randn(1, 5);
end

if nargin < 4
    figTitle = '';
end

if nargin < 5
    displayChannels = true;
end

%%
n = length(shape);

shape = 0.4 * max(plateDim)/max(abs(shape)) * shape;

fig = figure('Name', figTitle);
set(fig, 'Position', get(fig, 'Position') .* [1 1 1.2 0.9]);
ax = axes(fig);
hold(ax, 'on');
% legend(ax, 'AutoUpdate' ,'off');

% plate
Xplate = [0, 0; plateDim(1), plateDim(1)];
Yplate = [0, plateDim(2); 0, plateDim(2)];
Zplate = zeros(2);

surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', 'black', 'FaceColor', 0.5*[1 1 1], 'FaceAlpha', 0.3);

% dof
for k = 1:n
    X = dofPos(k, 1) * ones(1, 2);
    Y = dofPos(k, 2) * ones(1, 2);
    Z = [0, shape(k)];
    
    dofColor = 0.8 * [1 0 0];
    
    plot3(ax, X, Y, Z, 'LineWidth', 2, 'DisplayName', ['channel ', num2str(k)], 'Color', dofColor);
    
    %scatter3(ax, X(1), Y(1), Z(1), 2, 'X', 'LineWidth', 10, 'MarkerEdgeColor', dofColor);
    if shape(k) >= 0
        arrowPointer = '^';
    else
        arrowPointer = 'v';
    end
    scatter3(ax, X(2), Y(2), Z(2), 2, arrowPointer, 'LineWidth', 3, 'MarkerEdgeColor', dofColor);
    
    if displayChannels
        %text(ax, X(2), Y(2), Z(2) + 0.1*sign(shape(k))*max(abs(shape)), ['', num2str(k)]);
        text(ax, X(1), Y(1), Z(1), [' ch', num2str(k)]);
    end
end


% legendNames = {};
% for k = 1:length(shape)
%     legendNames{end+1} = ['dof', num2str(k)];
% end
% 
% legend(ax, legendNames);



view(ax, -30, 20);

pbaspect(ax, [plateDim, 0.3*max(plateDim)]);
set(ax,'visible','off');


end

