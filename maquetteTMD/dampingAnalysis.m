load('maquetteTMD/donnees/mData.mat');

Fa = mesureAmortissement;

Fa = transpose(Fa);
Fa = Fa(1,:) + 1i*Fa(2,:);

a = real(ifft(Fa));
%a = [a, zeros(1, 9*length(a))];

t = 1/400 * (1:length(a));

%% t0

t0 = 0.62;

a = a(t>=t0);
t = t(t>=t0);

%% wavelet

fig = figure;
ax = axes(fig);
plt = plot(ax, t, a);


WaveletMenu('WaveletPlot', plt);












