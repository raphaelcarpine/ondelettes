Fs = 100;
T = 10;
f1 = 5;
f2 = 3.16;


t = 0:(1/Fs):T;
x = sin(2*pi*f1*t);% .* exp(-0.03*2*pi*f1*t);
% x = [2; -1] * x;
% x = x + [0.5; -0.5] * sin(2*pi*f2*t+2) .* exp(-0.03*2*pi*f2*t);

figure;
plt = plot(t, x);

fmin = 1;
fmax = 10;
Q = 10;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q);