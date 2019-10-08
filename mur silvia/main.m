%etape et transient
P = 0;
transient = 3;

t0 = 0.25;
tf = inf;


[t, X] = getData(P, transient);

X = X(:, t>=t0 & t<tf);
t = t(t>=t0 & t<tf);





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
Q = 50;
MaxRidges = 10;
MaxParallelRidges = 2;
fmin = 50;
fmax = 70;
NbFreq = 200;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);


%RegressionMenu


%% test

ridges = {};
for k = 1:9
    ridges{end+1} = RidgeExtract(t, X(k,:), Q, fmin, fmax, NbFreq,...
        'NbMaxParallelRidges', 2, 'NbMaxRidges', 10);
end


[t, freqs, shapes] = getModes(ridges);













