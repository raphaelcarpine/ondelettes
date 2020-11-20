dataFolder = 'pont sens\donnees reelles\data'; % dossier où les fichier .mat sont enregistrés
file = 'TGV3A';

load([dataFolder, '\', file, '.mat']);
eval(['X = ', file, ';']);
X = transpose(X);
t = X(1, :);
X = X(2:end, :);


%% information colonnes

colones = [3];

% temps_avant = dataArray{:, 1};
% c_02_avant = dataArray{:, 2};
% a_p12_avant = dataArray{:, 3};
% a_p15_avant = dataArray{:, 4};
% a_p13_avant = dataArray{:, 5};
% a_p10_avant = dataArray{:, 6};
% a_p07_avant = dataArray{:, 7};
% c_08_avant = dataArray{:, 8};
% a_p04_avant = dataArray{:, 9};
% a_p09_avant = dataArray{:, 10};
% c_11_avant = dataArray{:, 11};
% c_12_avant = dataArray{:, 12};
% c_13_avant = dataArray{:, 13};
% c_14_avant = dataArray{:, 14};
% a_p03_avant = dataArray{:, 15};
% c_16_avant = dataArray{:, 16};
% c_17_avant = dataArray{:, 17};

X = X(colones-1, :);


%% correction du temps

dt = diff(t);
dt = dt(abs(dt) < 0.25);
if any( abs(dt/mean(dt) - 1) > 1e-3)
    warning('pas de temps non constant');
end
dt = mean(dt);
t = dt * (0:length(t)-1);

%% affichage

figure;
plt = plot(t, X);


fmin = 1;
fmax = 10;
Q = 10;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
    'MultiSignalMode', false);


