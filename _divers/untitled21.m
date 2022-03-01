t = 0:0.01:100;
x = sin(2*pi*10*t);

Q = 5;
freqs = linspace(5, 15, 100);
MotherWavelet = 'morlet';
ct = 3;

CWT = WvltComp(t, x, freqs, Q, 'MotherWavelet', MotherWavelet);

ridge = SingleRidgeExtract(t, freqs, CWT, MotherWavelet, Q, ct, 'slope2', 10);

figure;
plot(ridge.time, abs(ridge.val));
figure;
plot(ridge.time, ridge.freq);