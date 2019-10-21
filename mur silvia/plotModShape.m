function plt = plotModShape(shape)
%PLOTMODSHAPE Plot the modale shape 
%   Detailed explanation goes here


if nargin < 1 % exemple
    shape = [1.0000, -0.0027, -0.9622, 0.6654, -0.0474, -0.7961, 0.4438, -0.0047, -0.0557];
end

shape = reshape(shape, 3, 3).';

H = 150;
L = 131;

% sensors
X = [20, 65.5, 111; 20, 65.5, 111; 20, 65.5, 111];
Y = [114.25, 114.25, 114.25; 75, 75, 75; 36.25, 36.25, 36.25];

nH = 100;
nL = 100;

% ground
X0 = [20, 65.5, 111];
Y0 = zeros(1, 3);
shape0 = zeros(1, 3);

% total shape points
X = [X; X0];
Y = [Y; Y0];
shape = [shape; shape0];

% grid
x = linspace(0, L, nL);
y = linspace(0, H, nH);
[x, y] = meshgrid(x, y);

%% interpolation

shape2 = interp2(X, Y, shape, x, y, 'makima'); %, 'spline'

%% affichage

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plt = surf(ax, x, shape2, y, shape2, 'EdgeColor', 'none');
for x = [0, 20, 65.5, 111, 131]
    plot3(ax, [x, x], [0, 0], [0, 150], 'Color', 0.5*[1 1 1]);
end
for y = [0, 36.25, 75, 114.25, 150]
    plot3(ax, [0, 131], [0, 0], [y, y], 'Color', 0.5*[1 1 1]);
end

load('mur silvia/customColorMap');

colormap(fig, customColorMap);

view(ax, 45, 20);

pbaspect(ax, [1 0.1 1]);
set(ax,'visible','off');

end

