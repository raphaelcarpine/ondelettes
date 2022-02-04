clear all

savePath = 'ponts marne 2\erreur amort\save';

MotherWavelet = 'cauchy';

% mode
manualMode = 1;
saveResults = 0;
saveResults = saveResults & ~manualMode;

% affichage
pauseModes = 0; % pause entre les modes
plotCWT = 0;
plotRidgeExtract = 0;
plotShapes = 0;
plotShapesTime = 0;
plotTemp = 0;

% filtrage passe haut
filtrage = 1;
fc_filtre = 3.5; % freq coupure
fmin_filtrage = 5; % min freq propre avec filtre

% séparation signal
Nsep = 1;

% pont et mode
bridge = 6;
Kf = 1;

% Meff
Meff0 = 1;

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

% mise en forme
X = X.';
T = T.';

[X, T] = removeRedundantData(X, T);

[X, T] = removeNanSignal(X, T);

X = -X; % capteurs vers le bas
X = X - mean(X, 2); % moyenne


% tri channels par ordre croissant
[channelNames, I] = sort(channelNames);
X = X(I, :);
I = 1:length(I);
for kc = 1:length(channelNames) % ch2 à la fin
    if strcmp(channelNames{kc}(end-2:end), 'ch2')
        channelNames = [channelNames(1:kc-1), channelNames(kc+1:end), channelNames(kc)];
        I = [I(1:kc-1), I(kc+1:end), I(kc)];
        break
    end
end
X = X(I, :);

Nt = length(T);
dt = (T(end) - T(1))/(Nt-1);
Ttot = Nt*dt;

%%

for kch = 1:size(X, 1)
    figure;
    Xch = X(kch, :);
    Xch = butterworthFilter(T, Xch, 0.1, 'high', 5);
    histogram(Xch, 'Normalization', 'pdf');
    hold on
    stdX = std(Xch);
    u = get(gca, 'XLim');
    u = linspace(u(1), u(2), 1000);
    loiNorm = 1/(stdX*sqrt(2*pi)) * exp(-u.^2/(2*stdX^2));
    plot(u, loiNorm);
end




















