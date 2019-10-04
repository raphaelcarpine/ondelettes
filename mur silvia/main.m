%etape et transient
P = 0;
transient = 1;

[t, X] = getData(P, transient);



% fig = figure;
% ax = axes(fig);
% hold(ax, 'on');
% plts = nan(1, 9);
% for i = 1:9
%     plts(i) = plot(t, X(i,:), 'Parent', ax);
% end



fig = figure;
ax = axes(fig);
plts = plot(t, X, 'Parent', ax);
plts = transpose(plts);


%ondelette
Q = 15;
MaxParallelRidges = 1;
fmin = 10;
fmax = 20;
NbFreq = 300;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);