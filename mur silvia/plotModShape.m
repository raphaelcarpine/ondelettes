function plt = plotModShape(shape)
%PLOTMODSHAPE Plot the modale shape 
%   Detailed explanation goes here
shape = reshape(shape, 3, 3).';

H = 150;
L = 131;

% sensors
X = [20, 65.5, 111; 20, 65.5, 111; 20, 65.5, 111];
Y = [114.25, 114.25, 114.25; 75, 75, 75; 36.25, 36.25, 36.25];

nH = 40;
nL = 30;

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

shape2 = interp2(X, Y, shape, x, y); %, 'cubic'

%% affichage

fig = figure;
ax = axes(fig);
plt = surf(ax, x, y, shape2);

end

