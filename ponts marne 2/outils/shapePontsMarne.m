function fig = shapePontsMarne(plateDim, sensorsPos, sensorsDir, shape, figTitle)
%COMPLEXSHAPEPLOTDEFAULT Summary of this function goes here
%   Detailed explanation goes here

displayChannelNb = false;
horizontalInterp = 'quad'; % 'none' 'lin', 'quad';
shapeCoeff = 0.1;


if nargin < 4
    warning('shape mising');
    
    L = 74 + 2*(0.51+0.125);
    l = 8 + 2*0;
    D = 12;
    plateDim = [L, l];
    sensorsPos = [L/2-2*D, 0;L/2-1*D, 0;L/2, 0;L/2+1*D, 0;L/2+2*D, 0;L/2+3*D, 0;L/2-3*D, l;L/2-2*D, l;L/2-1*D, l;L/2, l;L/2+1*D, l;L/2+2*D, l;L/2+3*D, l;L/2-3*D, 0;L/2, 0];
    sensorsDir = [0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 0, 1;0, 1, 0];
    shape = randn(1, 15);
    shape = shape / norm(shape);
end

if nargin < 5
    figTitle = '';
end

if any(~isreal(shape))
    error(' ');
end


if length(shape) ~= size(sensorsPos, 1)
    warning(sprintf('Mode shape dimensions: %u (%u expected)', [length(shape), size(sensorsPos, 1)]));
end

if size(shape, 1) == 1
    shape = shape.';
end


% shape = shapeCoeff*shape;
shape = shapeCoeff * max(plateDim)/(max(shape) - min(shape)) * shape;


%% mise en forme surface

indSensVert = sensorsDir(:, 3) == 1; % capteurs verticaux
indSensHor = find(sensorsDir(:, 2) == 1); % capteur horizontal

indSensLeft = find(indSensVert & sensorsPos(:, 2) <= mean(sensorsPos(indSensVert, 2))); % capteurs verticaux amont
indSensRight = find(indSensVert & sensorsPos(:, 2) > mean(sensorsPos(indSensVert, 2))); % capteurs verticaux aval

[~, I] = sort(sensorsPos(indSensLeft, 1)); % tri croissant
% sensorsPos(indSensLeft, 1) = sensorsPos(indSensLeft(I), 1);
indSensLeft = indSensLeft(I);
[~, I] = sort(sensorsPos(indSensRight, 1)); % tri croissant
% sensorsPos(indSensRight, 1) = sensorsPos(indSensRight(I), 1);
indSensRight = indSensRight(I);

posSensX = [sensorsPos(indSensLeft, 1), sensorsPos(indSensRight, 1)];
posSensY = [sensorsPos(indSensLeft, 2), sensorsPos(indSensRight, 2)];


%% interpolation deplacement horizontal

switch horizontalInterp
    case 'none'
        horInterpFunc = @(X) zeros(size(X));
    case 'lin'
        L = plateDim(1);
        x0 = sensorsPos(indSensHor, 1);
        y0 = shape(indSensHor);
        horInterpFunc = @(X) y0/x0 * X .* (X<=x0) - y0/(L-x0) * (X-L) .* (X>x0);
    case 'quad'
        L = plateDim(1);
        x0 = sensorsPos(indSensHor, 1);
        y0 = shape(indSensHor);
        coeffA = y0/(x0^2-L*x0);
        coeffB = -L*coeffA;
        horInterpFunc = @(X) coeffA*X.^2 + coeffB*X;
end

%% calcul surface

posSensZ = [shape(indSensLeft), shape(indSensRight)];

posSensY = posSensY + horInterpFunc(posSensX);

%% interpolation

n_interp = 100;
posSensXinterp = [linspace(posSensX(1, 1), posSensX(end, 1), n_interp).', linspace(posSensX(1, 2), posSensX(end, 2), n_interp).'];
posSensYinterp = nan(size(posSensXinterp));
posSensZinterp = nan(size(posSensXinterp));
for c = 1:2
    posSensYinterp(:, c) = interp1(posSensX(:, c), posSensY(:, c), posSensXinterp(:, c), 'pchip');
    posSensZinterp(:, c) = interp1(posSensX(:, c), posSensZ(:, c), posSensXinterp(:, c), 'pchip');
end


%% plot
n = length(shape);

fig = figure('Name', figTitle);
set(fig, 'Position', get(fig, 'Position') .* [1 1 1.2 0.9]);
ax = axes(fig);
hold(ax, 'on');
% legend(ax, 'AutoUpdate' ,'off');

% plate
Xplate = [0, 0; plateDim(1), plateDim(1)];
Yplate = [0, plateDim(2); 0, plateDim(2)];
Zplate = zeros(2);

surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', 'black', 'FaceColor', 0.3*[1 1 1], 'FaceAlpha', 0.3);


% dof
% surf(ax, posSensX, posSensY, posSensZ, 'EdgeColor', 'black', 'FaceColor', 0.5*[1 1 1],...
%     'FaceAlpha', 0., 'FaceLighting', 'gouraud', 'LineWidth', 2);
surf(ax, posSensXinterp, posSensYinterp, posSensZinterp, 'EdgeColor', 'none', 'FaceColor', 0.5*[1 1 1],...
    'FaceAlpha', 1, 'FaceLighting', 'gouraud', 'LineWidth', 2);
for c = 1:2
    plot3(ax, posSensXinterp(:, c), posSensYinterp(:, c), posSensZinterp(:, c), 'k', 'LineWidth', 2);
end
for l = 1:size(posSensX, 1)
    plot3(ax, posSensX(l, :), posSensY(l, :), posSensZ(l, :), 'k', 'LineWidth', 2);
end


% dof
for k = 1:n
    X = sensorsPos(k, 1) * [1 1];
    Y = sensorsPos(k, 2) * [1 1];
    Z = [0 0];
    X = X + shape(k) * [0, sensorsDir(k, 1)];
    Y = Y + shape(k) * [0, sensorsDir(k, 2)];
    Z = Z + shape(k) * [0, sensorsDir(k, 3)];
    
    dofColor = 0.8 * [1 0 0];
    
%     mArrow3(ax, [X(1), Y(1), Z(1)], [X(2), Y(2), Z(2)], 'color', dofColor, 'stemWidth', min(plateDim)/20);
    
    if displayChannelNb
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

daspect(ax, [1 1 1]);
set(ax, 'PositionConstraint', 'innerposition');
set(ax, 'XLim', [0, plateDim(1)]);
set(ax, 'YLim', [0, plateDim(2)]);
set(ax, 'ZLim', max(plateDim)*0.2*[-1 1]);
set(ax, 'Clipping', 'off');
% pbaspect(ax, [plateDim, 0.3*max(plateDim)]);
set(ax,'visible','off');
set(ax, 'InnerPosition', [0.02 0.02 0.96 0.96]);

% lightangle(ax, -70, 40);
% lighting gouraud


%% vecteurs capteurs (figures capteurs)

% Ksensors = [2:4, 9:11];
% Ksensors = [3, 10];
% Ksensors = [2, 4, 9, 11];
% Hsensors = 5;
% for ks = Ksensors
%     if ks == 2 || ks == 9
%         plot3(ax, sensorsPos(ks, 1)*[1 1], sensorsPos(ks, 2)*[1 1], [0 Hsensors], 'r', 'LineWidth', 3);
%         plot3(ax, sensorsPos(ks, 1), sensorsPos(ks, 2), Hsensors, 'r', 'LineWidth', 3, 'Marker', '^');
%     else
%         plot3(ax, sensorsPos(ks, 1)*[1 1], sensorsPos(ks, 2)*[1 1], [0 -Hsensors], 'r', 'LineWidth', 3);
%         plot3(ax, sensorsPos(ks, 1), sensorsPos(ks, 2), -Hsensors, 'r', 'LineWidth', 3, 'Marker', 'v');
%     end
% end


end

