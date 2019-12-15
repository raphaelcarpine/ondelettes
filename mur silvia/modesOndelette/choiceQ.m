%etape et transient
P = 0;
transient = 3;

singleRidgeMode = false;

t0 = 0;
tf = inf;


[t, X] = getData(P, transient);

X = X(:, t>=t0 & t<tf);
t = t(t>=t0 & t<tf);




sensors = 1:9;



fig = figure;
ax = axes(fig);
plts = plot(t, X(sensors,:), 'Parent', ax);
plts = transpose(plts);


%ondelette
Q = 10;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 6;
fmax = 10;
NbFreq = 400;

ct = 3;
cf = 5;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges, 'CtEdgeEffects', ct);

