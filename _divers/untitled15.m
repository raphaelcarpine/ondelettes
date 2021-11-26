dt = 0.001;
t = 0:dt:10;

x = sin(2*pi*10*t);
a = -(2*pi*10)^2*sin(2*pi*10*t);


figure('Name', 'Position');
pltx = plot(t, x);

figure('Name', 'Acceleration');
plta = plot(t, a);


fmin = 5;
fmax = 15;

WaveletMenu('WaveletPlot', pltx, 'fmin', fmin, 'fmax', fmax);
WaveletMenu('WaveletPlot', plta, 'fmin', fmin, 'fmax', fmax);