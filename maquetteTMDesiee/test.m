load('maquetteTMDesiee/donnees/mData.mat');



%% affichage

mesure  = test2;

mesure = transpose(mesure);
t = mesure(1, :);
a = mesure(2, :);

fig  = figure;
ax = axes(fig);
plt = plot(ax, t, a);
xlabel(ax, 't');
ylabel(ax, 'a');

%% wvlt

fmin = 0.5;
fmax = 10;
NbFreq = 300;
Q = 10;


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq, 'Q', Q);



%%

getFreq(t, a, [0.8, 1.2], 1)





























