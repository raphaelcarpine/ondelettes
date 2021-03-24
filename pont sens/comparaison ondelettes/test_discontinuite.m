%% tests discontinuité freq

f1 = 1;
f2 = 1.1;

t = linspace(-500, 500, 100000);
x1 = exp(2i*pi*f1*t) .* (t < 0);
x2 = exp(2i*pi*f2*t) .* (t >= 0);
% x = x1 + x2;
x = x1;


figure;
plt = plot(t, real(x));


fmin = 0.2;
fmax = 2;
NbFreq = 1000;
Q = 10;
MaxRidges = 1;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges);