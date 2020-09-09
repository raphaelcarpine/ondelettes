clear all

%%
kSimul = 1;
results_folder = 'testsMasseAjoutee';

folder_dir = 'pont sens/simulation elements finis/resultats';
if ~isempty(results_folder)
    folder_dir = [folder_dir, '/', results_folder];
end

% recherche du fichier
name0 = sprintf('simul%d', kSimul);
listing = dir(folder_dir);
listingNames = {listing.name};
fileName = [];
for kname = 1:length(listingNames)
    name1 = strsplit(listingNames{kname}, '_');
    name1 = name1{1};
    if strcmp(name1, name0)
        fileName = listingNames{kname};
        break
    end
end

if isempty(fileName)
    error('file not found');
end

load([folder_dir, '/', fileName]);

%% calcul freq propre, freq excitation

try
    for kfreq = 1:1
        fprintf('freq. propre %d : %.2fHz, th. %.2fHz (%.1f%% error)\n',...
            [kfreq, freqs(kfreq), freqsTh(kfreq), (freqs(kfreq)-freqsTh(kfreq))/freqsTh(kfreq)*100]);
    end
catch
end

fprintf('freq. excitation : %.2f\n', c/L_wagons);

%% temps

if ~exist('t1', 'var') || ~exist('t2', 'var')
    t1 = 0;
    t2 = ((N_bogies-1)*L_wagons + L)/c;
end
disp(sprintf('t1 = %.2f, t2 = %.2f', [t1, t2]));

%% resultats

% pos_capteurs = linspace(0, L, N);
% pos_capteurs = pos_capteurs(2:end-1);

% animation
moving_coeff = 1;
movingPlot(Ytot, t, L, essieux, c, pos_capteurs, moving_coeff);

% wavelet capteurs
Ycapt = getYcapt(Ytot, pos_capteurs, dx);
%Ycapt = Ycapt + 1e-6 * randn(size(Ycapt));


fig = figure;
ax = axes(fig);
plt = plot(ax, t, Ycapt);

fmin = 3;
fmax = 16;
Q = 10;
MaxRidges = 1;
RealShapePlot = deformeePont(L, pos_capteurs);
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', true,...
    'MaxRidges', MaxRidges, 'RealShapePlot', RealShapePlot);
