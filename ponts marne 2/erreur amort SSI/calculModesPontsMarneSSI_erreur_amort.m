clear all

savePath = 'ponts marne 2\erreur amort SSI\save';

% mode
saveResults = 1;

% affichage
dispSeps = 0;
plotShapes = 1;
plotTemp = 0;

% filtrage passe haut
filtrage = 0;
fc_filtre = 3.5; % freq coupure
fmin_filtrage = 5; % min freq propre avec filtre

% séparation signal
Nsep = 20;

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

% découpage sous signaux
Ntsep = floor(Nt/Nsep); % longueur sous-intervalles
Isep = cell(1, Nsep); % sous-intervalles
for ksep = 1:Nsep
    Isep{ksep} = (ksep-1)*Ntsep + (1:Ntsep);
end
Tsep = Ntsep*dt;

% deformees modales
dimensionsShapes2;

%% calcul SSI

freqsCWT = [2.0762    2.1761    2.8349    2.9441    5.4652    6.5129    7.8272   15.0428   16.7422   21.1347];
eps_freqs_cwt = 0.03;
Ts = 1.4326/(2*pi*2.08*0.009); % 1.4326/\mu_1 = Topt (k=1), temps pour calcul SSI

% save
Freqs = nan(length(freqsCWT), Nsep);
Damps = nan(length(freqsCWT), Nsep);
Shapes = nan(size(X, 1), length(freqsCWT), Nsep);

% waitbar
[initWaitBar, updateWaitBar, closeWaitBar] =...
    getWaitBar(Nsep, 'windowTitle', 'Extraction modes', 'displayTime', 0);
initWaitBar();

for ksep = 1:Nsep
    updateWaitBar(nan, sprintf('sub-signal %d/%d', [ksep, Nsep]));
    
    % calcul SSI
    tic;
    [fn0, zeta0, phi0, paraPlot] = SSICOV(X(:, Isep{ksep}), dt, 'Ts', Ts, 'Nmin', 2, 'Nmax', 50, 'eps_cluster', 1);
    toc
    
    % correspondance freq CWT
    fn = nan(length(freqsCWT), 1);
    zeta = nan(length(freqsCWT), 1);
    phi = nan(size(X, 1), length(freqsCWT));
    for kf = 1:length(freqsCWT)
        [minDiff, kf0] = min(abs(fn0-freqsCWT(kf)));
        if minDiff/freqsCWT(kf) <= eps_freqs_cwt
            fn(kf) = fn0(kf0);
            zeta(kf) = zeta0(kf0);
            phi(:, kf) = phi0(kf0, :).';
        end
    end
    
    % save
    Freqs(:, ksep) = fn;
    Damps(:, ksep) = zeta;
    Shapes(:, :, ksep) = phi;
    
    % waitbar
    updateWaitBar();
    
    % plot
    if dispSeps
        plotStabDiag(paraPlot.fn, X(2, Isep{ksep}), 1/dt, paraPlot.status, paraPlot.Nmin, paraPlot.Nmax);
        
        if plotShapes
            for kf = 1:length(fn)
                if ~isnan(fn(kf))
                    figName = sprintf('mode #%d: %.2fHz, %.2%%', [kf, fn(kf), 100*zeta(kf)]);
                    shapePlotBridge(phi(:, kf)*sign(phi(8, kf)), figName);
                end
            end
        end
        
        input('continue?');
    end
    
    close all
end

closeWaitBar();


%% results
disp(Freqs);
disp(Damps);

if saveResults
    saveFile = fullfile(savePath, ['modes_', dataFileName, '_Nsep', num2str(Nsep)]);
    save(saveFile, 'Freqs', 'Damps', 'Shapes');
    disp('saved');
end




















