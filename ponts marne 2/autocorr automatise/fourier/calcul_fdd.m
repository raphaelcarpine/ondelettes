clear all

removeNan = true;
plotTemp = 0;

Nsv = 6;

bridge = 0;

Taveraging = 10;

% filtrage passe haut
filtrage = 1;
fc_filtre = 3.5; % freq coupure

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


%% tri channels par ordre croissant

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

%% filtrage

if filtrage
    X = butterworthFilter(T, X, fc_filtre, 'high', 10);
end

%% FDD

Naveraging = floor((T(end)-T(1))/Taveraging);
MaxLag = inf;
XcorrScale = 'unbiased';
HalfXcorr = false;

[f, SV, mod_shapes] = computeFDD(T, X, 'Naveraging', Naveraging, 'MaxLag',...
    MaxLag, 'XcorrScale', XcorrScale, 'HalfXcorr', HalfXcorr);

%% plot fourier

figure('Name', dataFileName);
plot(f, SV(1:Nsv, :));
set(gca, 'yscale', 'log');
xlim([0 30]);

PeakPickingMenu(false, true);

%% deformees

dimensionsShapes2;

disp(' ');
disp('mode shapes');

while true
    f0 = input('f = ');
    if length(f0) ~= 1
        break
    end
    [~, kf0] = min(abs(f-f0));
    shape0 = mod_shapes(:, 1, kf0);
    shape0 = shape0 / sqrt(shape0.'*shape0);
    shape0 = shape0 * sign(max(real(shape0)) + min(real(shape0)));
    shapePlotBridgeAnim(shape0, [dataFileName, sprintf('; f = %.2fHz; I = %.1f%%', [f(kf0), 100*nonPropIndex(shape0)])]);
    fig = shapePlotBridge(real(shape0), [dataFileName, sprintf('; f = %.2fHz; I = %.1f%%', [f(kf0), 100*nonPropIndex(shape0)])]);
    
%     % saving
%     s = input('save? ', 's');
%     saveFolder = 'ponts marne 2\autocorr automatise\fourier\save';
%     switch s
%         case {'y', 'o', ''}
% %             savefig(fig, fullfile(saveFolder, ['mode_', dataFileName, sprintf('_%dcHz', 100*f(kf0))]));
%     end
end













