clear all
close all

kSimul = 12;

% w0 : 2, 8, 10, 15 ,16, 22, 33, 34, 43

%%  chargement fichiers
results_folder = 'TMD freyssinet/donnees';

% recherche du fichier
name0 = num2str(kSimul);
listing = dir(results_folder);
listingNames = {listing.name};
fileName = [];
for kname = 1:length(listingNames)
    name1 = strsplit(listingNames{kname}, ' ');
    name1 = name1{1};
    if strcmp(name1, name0)
        fileName = listingNames{kname};
        break
    end
end

if isempty(fileName)
    error('file not found');
end

load([results_folder, '/', fileName]);

disp(essai);


%% recherche t0

enveloppe = tableExcel(1:end, 6);

kt0 = 1;
while kt0 <= length(enveloppe) && (isnan(enveloppe(kt0)) || enveloppe(kt0) == 0)
    kt0 = kt0+1;
end

fprintf('t0 = %.2f\n', t(kt0));

%% affichage

X = X(2, :);

fig = figure('Name', fileName);
ax = axes(fig);
plt = plot(ax, t, X);


fmin = 1;
fmax = 10;
Q = 3;
MaxRidges = 1;
XLimRidge = [t(kt0), t(end)];
ctRidge = 1;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', true,...
    'MaxRidges', MaxRidges, 'XLim', XLimRidge, 'ctRidge', 1);



%% précision temporelle

f = 2;
[~, DeltaT] = getPrecision(f, Q);
fprintf('DeltaT = %.2f s (for f = %.1f Hz & Q = %.1f)\n', [DeltaT, f, Q]);













