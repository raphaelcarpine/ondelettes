clear all

%% data

dataFolder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise operations 2021\donnees'; % dossier où les fichier .csv sont
dataFileName = 'testLabo.mat';

load(fullfile(dataFolder, dataFileName));

disp(startDate);

X = X.';
T = T.';

%% enlever nan

if true
    Nt = size(X, 2);
    for kt = size(X, 2):-1:1
        if any(isnan(X(:, kt)))
            Nt = kt-1;
        end
    end
    X = X(:, 1:Nt);
    T = T(:, 1:Nt);
end

%% deformee

dimensionsShapes2;


%% plots

fig = figure;
ax = axes(fig);
plts = plot(T, X);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Acceleration [m/s²]');


%% wavelet

Q = 10;
NbMaxRidges = 1;
NbMaxParallelRidges = 1;
fmin = 0.5;
fmax = 5;

WaveletMenu('WaveletPlot', plts, 'Q', Q, 'fmin', fmin, 'fmax', fmax,...
    'MaxRidges', NbMaxRidges, 'MaxParallelRidges', NbMaxParallelRidges, 'RealShapePlot', shapePlotBridge);





