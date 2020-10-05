clear all
% close all

results_folder = 'joint freyssinet/donnees/';

load([results_folder, 'data11.mat']);
load([results_folder, 'data12.mat']);
load([results_folder, 'data21.mat']);
load([results_folder, 'data22.mat']);

X = [X11, X12; X21, X22];

%%  rééchantillonage

t = (1/Fs) * (0:size(X, 2)-1);

% fréquence de rééchantillonage
Fs_new = 200;

n_resampling = floor(Fs/Fs_new);
Fs_new = Fs / n_resampling;

% filtrage
Fc = 0.5 * Fs_new/2;
X = butterworthFilter(t, X, Fc, 'low', 5);
X = butterworthFilter(t, X, Fc, 'low', 5);

% rééchantillonage
X = X(:, 1:n_resampling:end);
t = t(1:n_resampling:end);


%% affichage

% X = X(2, :);

fig = figure;
ax = axes(fig);
plt = plot(ax, t, X(1, :));


fmin = 1;
fmax = 10;
Q = 3;
% MaxRidges = 1;
% XLimRidge = [t(kt0), t(end)];
% ctRidge = 1;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false);%,...
%     'MaxRidges', MaxRidges, 'XLim', XLimRidge, 'ctRidge', 1);



%% précision temporelle

% f = 2;
% [~, DeltaT] = getPrecision(f, Q);
% fprintf('DeltaT = %.2f s (for f = %.1f Hz & Q = %.1f)\n', [DeltaT, f, Q]);













