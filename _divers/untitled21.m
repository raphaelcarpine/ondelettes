t = 0:0.01:10;
x = sin(2*pi*10*t);

Q = 2;
freqs = linspace(5, 15, 300);
MotherWavelet = 'morlet';
ct = 3;
ridgeContinuity = 'slope';
ridgeTime = 3;

% figure;
% plot(t, x);

CWT = WvltComp(t, x, freqs, Q, 'MotherWavelet', MotherWavelet);

% WvltPlot2(t, freqs, CWT, 'abs', Q, ct, MotherWavelet, 'log', '', '', 'lin');

ridge = SingleRidgeExtract(t, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity, ridgeTime);

% figure;
% plot(ridge.time, abs(ridge.val));
figure('Name', [num2str(Q), ' ', ridgeContinuity, ' ', num2str(ridgeTime)]);
plot(ridge.time, ridge.freq);
ylim([1.8, 2.3]);