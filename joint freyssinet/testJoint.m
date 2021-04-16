clear all
% close all

results_folder = 'C:\Users\carpine\Documents\projets\donnees charles';
files_names = {'1_ZOOM0002_LR.mat', '2_ZOOM0003_LR.mat', '3_ZOOM0004_LR.mat', '4_ZOOM0001_LR.mat', '5_ZOOM0002_LR.mat'};

k = 1;
load(fullfile(results_folder, files_names{k}));


%%  rééchantillonage

t = (1/Fs) * (0:size(X, 2)-1);

if 0
    % fréquence de rééchantillonage
    Fs_new = 10000;
    
    n_resampling = floor(Fs/Fs_new);
    Fs_new = Fs / n_resampling;
    
    % filtrage
    Fc = 0.8 * Fs_new/2;
    Fc = 4000;
    X = butterworthFilter(t, X, Fc, 'low', 5);
    X = butterworthFilter(t, X, Fc, 'low', 5);
    
    % rééchantillonage
    X = X(:, 1:n_resampling:end);
    t = t(1:n_resampling:end);
end

%% selection

ti = 16;
tf = 972;

% ti = 20;
% tf = 200;

X = X(:, t >= ti & t <= tf);
t = t(t >= ti & t <= tf);


%% affichage

% X = X(2, :);

fig = figure;
ax = axes(fig);
plt = plot(ax, t, X);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Signal [V]');
legend(ax, {'channel 1', 'channel 2'});

XLim = [16, 972];
fmin = 10;
fmax = 24000;
Q = 30;
FrequencyScale = 'log';
WvltScale = 'log';
MotherWavelet =  'morlet';
SignalUnit = 'V';
SquaredSignalUnit = 'V²';

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false,...
    'WvltScale', WvltScale, 'FrequencyScale', FrequencyScale, 'MotherWavelet', MotherWavelet, 'XLim', XLim,...
    'SignalUnit', SignalUnit, 'SquaredSignalUnit', SquaredSignalUnit);



%% précision temporelle

% f = 2;
% [~, DeltaT] = getPrecision(f, Q);
% fprintf('DeltaT = %.2f s (for f = %.1f Hz & Q = %.1f)\n', [DeltaT, f, Q]);













