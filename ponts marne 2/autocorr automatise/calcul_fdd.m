clear all

removeNan = true;
plotTemp = 0;

plotFourier = 1;
Nsv = 6;

bridge = 0;

%% data

dataFilePath = choixData(bridge);

load(dataFilePath);

[~, dataFileName] = fileparts(dataFilePath);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp(dataFileName);
disp(startDate);
% T°
if plotTemp
    [TemperatureTime, TemperatureTemp] = getTemperature(startDate);
    figure('Name', dataFileName);
    plot(TemperatureTime, TemperatureTemp);
    ylabel('Temperature [°C]');
end

%% mise en forme

X = X.';
T = T.';

[X, T] = removeRedundantData(X, T);

if removeNan
    [X, T] = removeNanSignal(X, T);
end

X = -X; % capteurs vers le bas
X = X - mean(X, 2); % moyenne


%% deformee

dimensionsShapes2;

%% FDD

Taveraging = 50;

Naveraging = floor((T(end)-T(1))/Taveraging);
MaxLag = inf;
XcorrScale = 'unbiased';
HalfXcorr = false;

[f, SV, mod_shapes] = computeFDD(T, X, 'Naveraging', Naveraging, 'MaxLag',...
    MaxLag, 'XcorrScale', XcorrScale, 'HalfXcorr', HalfXcorr);

%%

if plotFourier
    figure;
    plot(f, SV(1:Nsv, :));
    set(gca, 'yscale', 'log');
    xlim([1 4]);
end













