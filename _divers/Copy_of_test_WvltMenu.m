Fs = 100;
T = 1000*[-1 1];
f = 10;
d = pi/2;


t = T(1):(1/Fs):T(2);
A0 = 0.0;
A = ((A0 - t).*(t<0) + (A0 + t).*(t>=0));
A = 1 + exp(2*randn(size(t)));
% A = randn(size(t));
A = getSmoothSignal(t, A, 'gaussian', 0.2);
x = A .* sin(2*pi*f*t+d);

figure;
plt = plot(t, x);

fmin = 5;
fmax = 15;
Q = 3;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q);