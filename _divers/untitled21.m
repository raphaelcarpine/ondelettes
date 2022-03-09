t = 0:0.005:10;
x = sin(2*pi*10*(t + 0*sin(t))) .* (-(t-t(end)/2+0.3)/300);

Q = 5;
freqs = linspace(7, 15, 100);
MotherWavelet = 'morlet';
ct = 3;

% figure;
% plot(t, x);

CWT = WvltComp(t, x, freqs, Q, 'MotherWavelet', MotherWavelet);

WvltPlot2(t, freqs, CWT, 'abs', Q, ct, MotherWavelet, 'log', '', '', 'lin');

ridge = SingleRidgeExtract(t, freqs, CWT, MotherWavelet, Q, ct, 'slope', 0.1);

figure;
plot(ridge.time, abs(ridge.val));
figure;
plot(ridge.time, ridge.freq);