x = nan; % à completer
t = nan; % à compléter




%% plot

fig = figure;
ax = axes(fig);
plts = plot(ax, t, x);
xlabel(ax, 't');
ylabel(ax, 'a');




%% ondelette


Q = 5;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 1;
fmax = 100;
NbFreq = 300;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);