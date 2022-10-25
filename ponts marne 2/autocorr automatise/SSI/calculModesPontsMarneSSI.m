clear all

savePath = 'ponts marne 2\autocorr automatise\SSI\save';

% mode
manualMode = 0;
saveResults = 0;
saveResults = saveResults & ~manualMode;

% affichage
pauseModes = 1; % pause entre les modes
plotRidgeExtract = 1;
plotShapes = 1;
plotTemp = 0;

% filtrage passe haut
filtrage = 0;
fc_filtre = 3.5; % freq coupure
fmin_filtrage = 5; % min freq propre avec filtre

bridge = 6; % annet
Kf = 0;

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

Ttot = T(end) - T(1);
dt = Ttot/(length(T)-1);

% deformees modales
dimensionsShapes2;

%% calcul SSI

freqsCWT = [2.0762    2.1761    2.8349    2.9441    5.4652    6.5129    7.8272   15.0428   16.7422   21.1347];
eps_freqs_cwt = 0.03;

Ts = 1.4326/(2*pi*2.08*0.009); % 1.4326/\mu_1 = Topt (k=1)

tic;
[fn0,zeta0,phi0,paraPlot] = SSICOV(X, dt, 'Ts', Ts, 'Nmin', 2, 'Nmax', 50, 'eps_cluster', 0.5);
toc

[h] = plotStabDiag(paraPlot.fn, X(2,:), 1/dt, paraPlot.status, paraPlot.Nmin, paraPlot.Nmax);





















