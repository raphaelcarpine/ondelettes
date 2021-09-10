%% test threshold

load('piles mohamed\data\6_10cm_x_4096.mat');
[X, channelNames] = channelNamesConversion(X, channelNames);
X = X.';
T = T.';

figure;
plot(T(1:size(X, 2)), X);
legend(channelNames);
selectLine();

Kt0 = getShockIndex(T, X, channelNames);

for kt = Kt0
    xline(T(kt));
end

%% test autocorr

load('piles mohamed\data\7_10cm_x_4096.mat');
[X, channelNames] = channelNamesConversion(X, channelNames);
X = X.';
T = T.';

figure;
plt = plot(T(1:size(X, 2)), X);
legend(channelNames);
selectLine();

Kt0 = getShockIndex(T, X, channelNames);
T0 = T(Kt0)

fmin = 0;
fmax = 2000;
Q = 10;
XLim = T([Kt0(1), Kt0(2)-2]);
[fctDefModale, fctDefModaleAnimation] = defModales(6, channelNames, '', 5);

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'RemoveMean', true,...
    'XLim', XLim, 'RealShapePlot', fctDefModale, 'AutocorrelationMode', true,...
    'AutocorrelationNsvd', 3, 'AutocorrelationFourierSVDMode', true, 'FourierScale', 'log',...
    'FrequencyScale', 'log', 'AnimatedShapePlot', fctDefModaleAnimation);

%% test reponses libres

load('piles mohamed\data\2_5cm_x_4096.mat');
[X, channelNames] = channelNamesConversion(X, channelNames);
X = X.';
T = T.';

figure;
plt = plot(T(1:size(X, 2)), X);
legend(channelNames);
selectLine();

Kt0 = getShockIndex(T, X, channelNames);
T0 = T(Kt0)

fmin = 1;
fmax = 40;
Q = 7.14;
XLim = T([Kt0(1), Kt0(2)-2]);
[fctDefModale, fctDefModaleAnimation] = defModales(2, channelNames, '', 5);

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'RemoveMean', true,...
    'XLim', XLim, 'RealShapePlot', fctDefModale, 'FourierScale', 'log',...
    'AnimatedShapePlot', fctDefModaleAnimation, 'StopRidgeWhenIncreasing', true,...
    'MultiSignalMode', true, 'MultiSignalModeAbsValue', true);

%% test dephasage

load('piles mohamed\data\7_10cm_x_4096.mat');
[X, channelNames] = channelNamesConversion(X, channelNames);
X = X.';
T = T.';

X = X - mean(X, 2);

figure;
plt = plot(T(1:size(X, 2)), X);
legend(channelNames);
selectLine();

[synchro, resync, resyncVect] = testSynchro(T, X, channelNames)