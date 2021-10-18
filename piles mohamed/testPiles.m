dataFolder  = 'C:\Users\carpine\Documents\projets\experiences affouillement mohamed\data';

%% visualisation

pile = 4;
prof = 15;
direction = 'x';
g = 8;
load(fullfile(dataFolder, 'raw mat files', sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g])));
[X, channelNames] = channelNamesConversion(X, channelNames);

[X, T] = removeRedundantData(X, T);
X = X - mean(X, 2);

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();

% ylim(0.2*[-1 1]);

%% test threshold

pile = 1;
prof = 10;
direction = 'x';
g = 8;
load(fullfile(dataFolder, sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g])));

figure;
plot(T, X);
legend(channelNames);
selectLine();

Kt0 = getShockIndex(T, X, channelNames);

for kt = Kt0
    xline(T(kt));
end

%% test autocorr

pile = 1;
prof = 25;
direction = 'x';
g = 8;
load(fullfile(dataFolder, sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g])));

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();

Kt0 = getShockIndex(T, X, channelNames);
T0 = T(Kt0)

fmin = 0;
fmax = 2000;
Q = 10;
XLim = T([Kt0(1), Kt0(2)-2]);
[fctDefModale, fctDefModaleAnimation] = defModales(pile, channelNames,...
    sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g]), prof);

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'RemoveMean', true,...
    'XLim', XLim, 'RealShapePlot', fctDefModale, 'AutocorrelationMode', true,...
    'AutocorrelationNsvd', 3, 'AutocorrelationFourierSVDMode', true, 'FourierScale', 'log',...
    'FrequencyScale', 'log', 'AnimatedShapePlot', fctDefModaleAnimation);

%% test reponses libres

pile = 4; %3, 4, 4, 4
prof = 15; %15, 5, 10, 25
direction = 'x'; %y, x, y, x
g = 8;
load(fullfile(dataFolder, sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g])));

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();

Kt0 = getShockIndex(T, X, channelNames);
T0 = T(Kt0)

fmin = 1;
fmax = 40;
Q = 7.14;
XLim = T([Kt0(1), Kt0(2)-2]);
[fctDefModale, fctDefModaleAnimation] = defModales(pile, channelNames,...
    sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g]), prof);

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'RemoveMean', true,...
    'XLim', XLim, 'RealShapePlot', fctDefModale, 'FourierScale', 'log',...
    'AnimatedShapePlot', fctDefModaleAnimation, 'StopRidgeWhenIncreasing', true,...
    'MultiSignalMode', true, 'MultiSignalModeAbsValue', true);

%% test dephasage

pile = 6;
prof = 5;
direction = 'x';
g = 8;
load(fullfile(dataFolder, 'raw mat files', sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g])));

X = X - mean(X, 2);

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();

[synchro, maxDesync, resync, resyncVect] = testSynchro(T, X, channelNames)
if ~synchro
    [T, X] = resynchChannels(T, X, resyncVect);
end

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();