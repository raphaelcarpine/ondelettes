dt = 0.001;
t = -10:dt:10;

x1 = sin(2*pi*10*t);
x2 = sin(2*pi*20*t);
x = x1.*(t<=0) + x2.*(t>0);

figure;
plt = plot(t, x);

WaveletMenu('WaveletPlot', plt, 'fmin', 5, 'fmax', 25);