t = linspace(0, 10, 10000);
x = [2; 1i] * exp(2i*pi*10*t);
x = real(x);

x = x .* (2 + cos(2*pi*t/10));

figure;
plt = plot(t, x);


WaveletMenu('WaveletPlot', plt, 'fmin', 5, 'fmax', 15);