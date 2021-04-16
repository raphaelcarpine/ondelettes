%% inputs : 2, 6, 8, 9, 13, 15
N = 1;
mode = 2; % 0: affichage, 1: t0, 2: ft, 3: passage train, 4 : apres train
saveResults = false;

channels = [3 4 5 6 7 9 10 15];
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

%% data
% longueurs
L_w = 18.7; % longueur wagon
L_t = 200.19 - 2*(5.02-1.5); % longueur train
L_loco = 1.5 + 14 + 1.5; % longueur locomotive
L_p = 17.5; % longueur pont

% geometrie pont et position capteurs
l_p = 4.84; % largeur pont
pos_capt = [0, 0; 1, 0; 2, 0; 3, 0; 4, 0; 5, 0; 6, 0; 7, 0];
realShapePlotPont = @(shape, figTitle) shapePlotPlate([L_p, l_p], pos_capt, shape, figTitle); % deformees pont

% fréquences
ft0 = 4.3; % estimation ft

%% paths
dataFolder = 'pont sens/donnees reelles/data';
dataName = ['TGV', num2str(N), 'A.mat'];
dataPath = fullfile(dataFolder, dataName);

T0Folder = 'pont sens/article analyse donnees reelles';
T0Name = 'T0.xlsx';
T0Path = fullfile(T0Folder, T0Name);

FtFolder = T0Folder;
FtName = 'Ft.xlsx';
FtPath = fullfile(FtFolder, FtName);

saveFolder = 'pont sens/article analyse donnees reelles/save';
saveName = ['TGV', num2str(N), 'save'];
savePath = fullfile(saveFolder, saveName);


%% data loading

% t0
T0 = readtable(T0Path);
t0 = T0.Var1(N);

% ft
Ft = readtable(FtPath);
ft = Ft.Var1(N);

% acceleration
load(dataPath);
X = transpose(X);
t = X(1, :);
X = X(channels, :);
X = X - mean(X, 2);

% correction du temps
dt = diff(t);
dt = dt(abs(dt) < 0.25);
if any( abs(dt/mean(dt) - 1) > 1e-3)
    warning('pas de temps non constant');
end
dt = mean(dt);
t = dt * (0:length(t)-1);


%% premier affichage

if mode == 0
    figure;
    plt = plot(t, X);
    WaveletMenu('WaveletPlot', plt, 'fmin', 2, 'fmax', 20, 'Q', 10,...
        'RealShapePlot', @realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true);
end


%% t0

if mode == 1
    figure;
    plt = plot(t, X);
end


%% ft

if mode == 2
    if ft == 0
        ft_estim = ft0;
    else
        ft_estim = ft;
    end
    disp(['ft_estim = ', num2str(ft_estim)]);
    v_t = L_w * ft_estim;
    Delta_t = (L_p + L_loco)/v_t;
    t1 = t0;
    t2 = t0 + (L_p + L_t)/v_t;
    
    figure;
    plt = plot(t, X);
    WaveletMenu('WaveletPlot', plt, 'fmin', 2, 'fmax', 7, 'Q', 7, 'XLim', [t1+Delta_t, t2-Delta_t],...
        'RealShapePlot', @realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true);
end






















