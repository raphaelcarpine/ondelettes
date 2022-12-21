S = systLin(1, (2*pi)^2, 0.02*2*pi);

dt = 0.05;
t = 0:dt:10000;

%%

f1 = zeros(size(t));
f1(1) = 1;
x1 = S.response(f1, dt, 1);
figure;
plt1 = plot(t, x1);
xlabel('Temps [s]');
ylabel('Déplacement');
yticks([]);
WaveletMenu('WaveletPlot', plt1, 'fmin', 0, 'fmax', 2);

%%

f2 = randn(size(t));
x2 = S.response(f2, dt, 1);
figure;
plt2 = plot(t, x2);
xlabel('Temps [s]');
ylabel('Déplacement');
yticks([]);
WaveletMenu('WaveletPlot', plt2, 'fmin', 0, 'fmax', 2);