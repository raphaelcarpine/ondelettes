function fig = plotModShape(shape, figTitle)
%PLOTMODSHAPE Plot the modale shape 
%   Detailed explanation goes here

interpExt = false;
printChanels = true;

%%

if nargin < 2
    figTitle = '';
end
if nargin < 1 % exemple
    shape = [1.0000, -0.0027, -0.9622, 0.6654, -0.0474, -0.7961, 0.4438, -0.0047, -0.0557];
    %shape = [1.0000, 0, 0, 0, 0, 0, 0, 0, 0];
    shape = shape / sqrt(shape * shape.');
end

%%

if size(shape, 1) > size(shape, 2)
    shape = transpose(shape);
end

if size(shape, 1) ~= 1 || size(shape, 2) ~= 9 || ~all(imag(shape) == 0)
    lastwarn('');
    fig = nan;
    return
end

maxAbsShape = max(abs(shape));

shape = reshape(shape, 3, 3).';

H = 150;
L = 131;

% sensors
Xch = [20, 65.5, 111; 20, 65.5, 111; 20, 65.5, 111];
Ych = [114.25, 114.25, 114.25; 75, 75, 75; 36.25, 36.25, 36.25];

nH = 100;
nL = 100;

% ground
X0 = [20, 65.5, 111];
Y0 = zeros(1, 3);
shape0 = zeros(1, 3);

% total shape points
X = [Xch; X0];
Y = [Ych; Y0];
shape = [shape; shape0];

if interpExt
    X = [X0; X];
    Y = [H*ones(1, 3); Y];
    shape = [shape(1,:); shape];
    X = [zeros(5, 1), X, L*ones(5, 1)];
    Y = [Y(:,1), Y, Y(:,end)];
    shape = [shape(:,1), shape, shape(:,end)];
end

% grid
x = linspace(20, 111, nL);
y = linspace(0, 114.25, nH);


if interpExt
    x = linspace(0, L, nL);
    y = linspace(0, H, nH);
end

[x, y] = meshgrid(x, y);

%% interpolation

shape2 = interp2(X, Y, shape, x, y);
% shape2 = interp2(X, Y, shape, x, y, 'makima'); %, 'spline'

%% affichage

fig = figure('Name', figTitle);
ax = axes(fig);
hold(ax, 'on');
plt = surf(ax, x, shape2, y, shape2, 'EdgeColor', 'none');

grey = 0.5;
xGrid = [20, 65.5, 111];
yGrid = [0, 36.25, 75, 114.25];
if interpExt
    xGrid = [0, 20, 65.5, 111, 131];
    yGrid = [0, 36.25, 75, 114.25, 150];
end
for x = xGrid
    plot3(ax, [x, x], [0, 0], [yGrid(1), yGrid(end)], 'Color', grey*[1 1 1]);
end
for y = yGrid
    plot3(ax, [xGrid(1), xGrid(end)], [0, 0], [y, y], 'Color', grey*[1 1 1]);
end

if printChanels
    xch = transpose(Xch);
    xch = xch(:);
    ych = transpose(Ych);
    ych = ych(:);
    for ch = 1:9
        text(ax, xch(ch), -1, ych(ch), ['ch', num2str(ch)]);
    end
end

load('mur silvia/customColorMap');

colormap(fig, customColorMap);
% colormap(fig, jet);
caxis(ax, [-maxAbsShape, maxAbsShape]);
colorbar(ax);

view(ax, 45, 20);

pbaspect(ax, [1 0.1 1]);
set(ax,'visible','off');

pos = get(fig, 'Position');
pos(3) = 380;
set(fig, 'Position', pos);

end

