X = rand(10000, 9);
Fs = 600;

%%
X = transpose(X);
t = (0:size(X, 2)-1)/Fs;

figure;
plt = plot(t, X);

%%

fmin = 2;
fmax = 50;
Q = 10;
% XLim = [126, 132];
MotherWavelet =  'cauchy';

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', true, 'MotherWavelet', MotherWavelet);
