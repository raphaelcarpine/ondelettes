function fig = plotComplexModShape(shape, figTitle)
%PLOTMODSHAPE Plot the modale shape 
%   Detailed explanation goes here

printChanels = true;

%%

if nargin < 2
    figTitle = 'test';
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

if size(shape, 1) ~= 1 || size(shape, 2) ~= 9
    warning('');
    fig = nan;
    return
end


%% affichage

fig = figure('Position', [400 300 450 450], 'Name', figTitle);
ax = axes(fig);
hold(ax, 'on');

circle = exp(1i*linspace(0, 2*pi, 30));
circle2 = 1/2 * exp(1i*linspace(0, 2*pi, 30));
l0 = linspace(-1, 1, 2);
l30 = exp(1i*pi/6) * l0;
l60 = exp(1i*pi/3) * l0;
l90 = exp(1i*pi/2) * l0;
l120 = exp(2i*pi/3) * l0;
l150 = exp(5i*pi/6) * l0;

for k=1:9
    p0 = 2*mod(k-1, 3) + 2*1i*(3-fix((k-1)/3));
    p1 = p0 + shape(k);
    
    grey = 0.9;
    plot(ax, real(p0 + circle2), imag(p0 + circle2), 'color', grey+[0 0 0]);
    plot(ax, real(p0 + l0), imag(p0 + l0), 'color', grey+[0 0 0]);
    plot(ax, real(p0 + l30), imag(p0 + l30), 'color', grey+[0 0 0]);
    plot(ax, real(p0 + l60), imag(p0 + l60), 'color', grey+[0 0 0]);
    plot(ax, real(p0 + l90), imag(p0 + l90), 'color', grey+[0 0 0]);
    plot(ax, real(p0 + l120), imag(p0 + l120), 'color', grey+[0 0 0]);
    plot(ax, real(p0 + l150), imag(p0 + l150), 'color', grey+[0 0 0]);
    plot(ax, real(p0 + circle), imag(p0 + circle), 'black');
    
    plot(ax, real([p0 p1]), imag([p0, p1]), '-o');
    
    if printChanels
        text(ax, real(p0)-0.17, imag(p0)+0.5, ['ch', num2str(k)]);
    end
end

pbaspect(ax, [1 1 1]);
axis off

end

