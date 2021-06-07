function AltitudePlot(direction, XLabel, YLabel)
%ALTITUDEPLOT Summary of this function goes here
%   Detailed explanation goes here

switch direction
    case 'x'
        direction = 1;
    case 'y'
        direction = 2;
    otherwise
        error('Direction input error');
end

%% input

f = msgbox('Click on axis');
waitfor(f);
ax = gca;

if direction == 1
    prompt = 'Enter y:';
else
    prompt = 'Enter x:';
end
XY0 = inputdlg(prompt, 'Input');

%% surface extraction

S = findobj(ax, 'Type', 'Surface');

if isempty(S)
    error('No surface');
elseif length(S) > 1
    S = S(1);
    warning('Multiple surfaces');
end

X = S.XData;
Y = S.YData;
Z = S.ZData;
if all(Z == 0, 'all')
    Z = S.CData;
end

if all(X == X(1, :), 'all') && all(Y == Y(:, 1), 'all')
    x = X(1, :);
    y = Y(:, 1)';
elseif all(X == X(:, 1), 'all') && all(Y == Y(1, :), 'all')
    x = X(:, 1)';
    y = Y(1, :);
    Z = Z';
else
    error('Irregular pattern');
end


%% line

if direction == 1
    n = length(y);
    xy = y;
    xy2 = x;
else
    n = length(x);
    xy = x;
    xy2 = y;
end

k = 1;
while k <= n && XY0 < xy(k)
    k = k+1;
end
if k == 1 || k == n+1
    error('Outside limits');
end
    
r = (XY0 - xy(k-1)) / (xy(k) - xy(k-1));

    
end

