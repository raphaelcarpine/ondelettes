kSimul = 1;

% recherche du fichier
name0 = sprintf('simul%d', kSimul);
listing = dir('pont sens/simulation elements finis/resultats');
listingNames = {listing.name};
fileName = [];
for kname = length(listingNames)
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

load(['pont sens/simulation elements finis/resultats/', fileName]);

%%
disp(['freq propre 1 : ', num2str(pi/(2*L^2) * sqrt(E*J/mu))]);

% animation
moving_coeff = 1;
movingPlot(Ytot_interp, t_interp, L, essieux, c, pos_capteurs, moving_coeff);

% wavelet capteurs
Ycapt_interp = getYcapt(Ytot_interp, pos_capteurs, dx);
Ycapt_interp = Ycapt_interp + 1e-6 * randn(size(Ycapt_interp));


fig = figure;
ax = axes(fig);
plt = plot(ax, t_interp, Ycapt_interp);

fmin = 3;
fmax = 16;
Q = 10;
MaxRidges = 1;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', true, 'MaxRidges', MaxRidges);
