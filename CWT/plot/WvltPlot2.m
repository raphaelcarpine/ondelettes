function plt = WvltPlot2(t, freqs, wavelet, plotQuantity, Q, ctEdgeEffects, MotherWavelet, scale, figureTitle, axesTitle)

if ~ismember(plotQuantity, {'abs', 'module', 'arg', 'phase'})
    error('');
end

if nargin < 7
    scale = 'log10';
end
if nargin < 8
    figureTitle = '';
end
if nargin < 9
    axesTitle = '';
end


if isequal(scale, 'lin')
    moduleScale = @(x) x;
elseif isequal(scale, 'log')
    moduleScale = @(x) log(x);
elseif isequal(scale, 'log10')
    moduleScale = @(x) log10(x);
else
    error("");
end


%% rééchantillonage
reechant = true;

Ntmax = 20000;

if reechant
    Nt = length(t);
    step = ceil(Nt/Ntmax);
    indices = 1:step:Nt;
    
    t = t(indices);
    wavelet = wavelet(:, indices);
end


%% figure

fig = figure('Name', figureTitle);
ax = axes(fig);
hold(ax, 'on');

%% plot

[T, Freqs] = meshgrid(t, freqs);

if isequal(plotQuantity, 'abs') || isequal(plotQuantity, 'module')
    wavelet = moduleScale(abs(wavelet));
%     plt = pcolor(ax, T, Freqs, wavelet);
    plt = surf(ax, T, Freqs, wavelet);
    colormap(ax, jet);
    
    ZedgeEffects = max(max(wavelet)) + 1;
    
elseif isequal(plotQuantity, 'arg') || isequal(plotQuantity, 'phase')
    wavelet = angle(wavelet);
    plt = pcolor(T, Freqs, wavelet);
    colormap(ax, hsv);
    
    ZedgeEffects = 1;
end

xlabel(ax, 'Time (s)');
ylabel(ax, 'Frequency (Hz)');
title(ax, axesTitle);
shading flat

%% zones d'effets de bord

[~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);

if ctEdgeEffects > 0    
    t1 = t(1) + ctEdgeEffects * DeltaT(freqs);
    t2 = t(end) - ctEdgeEffects * DeltaT(freqs);
    
    % zone grisée
    tEdge1 = [t1, t(1), t(1)];
    freqsEdge1 = [freqs, freqs(end), freqs(1)];
    
    tEdge2 = [t2, t(end), t(end)];
    freqsEdge2 = [freqs, freqs(end), freqs(1)];
    
    patch(tEdge1,freqsEdge1, ZedgeEffects * ones(size(tEdge1)),'white','EdgeColor','none','FaceAlpha',.35);
    patch(tEdge2,freqsEdge2, ZedgeEffects * ones(size(tEdge2)),'white','EdgeColor','none','FaceAlpha',.35);
    
    %ligne démarcation
    plot3(ax, t1, freqs, ZedgeEffects * ones(size(t1)), 'black--', 'LineWidth', 1.5);
    plot3(ax, t2, freqs, ZedgeEffects * ones(size(t2)), 'black--', 'LineWidth', 1.5);
    
    % limites axes
    axis tight
end

set(ax, 'XLim', [t(1), t(end)]);

end