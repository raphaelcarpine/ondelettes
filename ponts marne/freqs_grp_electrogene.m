[x, Fs] = audioread('ponts marne/data/groupe_electrogene.wav');

t = (1/Fs) * (0:(length(x)-1));

fig = figure;
ax = axes(fig);
plt = plot(ax, t, x);

WaveletMenu('WaveletPlot', plt);