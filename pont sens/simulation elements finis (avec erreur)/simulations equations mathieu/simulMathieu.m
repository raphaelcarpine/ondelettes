% params systeme
w0 = 3 * 2*pi;
zeta = 0.;

%params mathieu
eps = 0.01;
wm = 5 * 2*pi;


%% integration num

DY = @(t, Y) [Y(2); -2*zeta*w0*Y(2) - w0^2*(1 + eps*cos(wm*t))*Y(1)];

t = 0:0.01:1000;
Y0 = [1; 0];

[~, Y] = ode45(DY, t, Y0);

y = Y(:, 1).';

%% affichage

figure;
plt = plot(t, y);


fmin = 3;
fmax = 16;
Q = 10;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false,...
    'SignalUnit', 'm', 'SquaredSignalUnit', 'm²');