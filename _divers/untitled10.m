f = 10;
z = 0.01;
l = 2*pi*f*(1i*sqrt(1-z^2) - z);

dt = 0.01;
T = 10.04;
t = 0:dt:T;

x = real(exp(l*t));

figure;
plt = plot(t, x);

WaveletMenu('WaveletPlot', plt, 'fmin', 5, 'fmax', 15);