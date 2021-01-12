function plt = WvltPlot2(t, freqs, wavelet, plotQuantity, Q, ctEdgeEffects, MotherWavelet, WvltScale, figureTitle, axesTitle, FreqScale)


plotScale = true;
autoScale = true;

if ~ismember(plotQuantity, {'abs', 'module', 'arg', 'phase'})
    error('');
end

if nargin < 7
    WvltScale = 'log';
end
if nargin < 8
    figureTitle = '';
end
if nargin < 9
    axesTitle = '';
end
if nargin < 10
    FreqScale = 'lin';
end


switch WvltScale
    case 'lin'
        moduleScale = @(x) x;
    case 'log'
        moduleScale = @(x) log10(x);
    otherwise
        error('wrong scale');
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
set(ax,'TickDir','out');
hold(ax, 'on');

%% plot

[T, Freqs] = meshgrid(t, freqs);

switch plotQuantity
    case {'abs', 'module'}
        wavelet = abs(wavelet);
%         plt = pcolor(ax, T, Freqs, wavelet);
        plt = surf(ax, T, Freqs, wavelet);
        shading(ax, 'flat');
        set(ax, 'ZScale', WvltScale);
        set(ax, 'ColorScale', WvltScale);
        colormap(ax, jet);
        if ~autoScale
            caxis(ax, [min(wavelet(wavelet ~= 0)), max(wavelet, [], 'all')]);
        end
        if plotScale
            c = colorbar(ax);
            set(c,'TickDir','out');
%             c.Label.String = '|CWT|';
        end
        
        ZedgeEffects = 1.1 * max(wavelet, [], 'all');
        
    case {'arg', 'phase'}
        wavelet = angle(wavelet);
        plt = pcolor(T, Freqs, wavelet);
        colormap(ax, hsv);
        
        ZedgeEffects = 1;
end

xlabel(ax, 'Time [s]');
ylabel(ax, 'Frequency [Hz]');
set(ax, 'YScale', FreqScale);
title(ax, axesTitle);


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
%     axis tight
end

set(ax, 'XLim', t([1, end]));
set(ax, 'YLim', freqs([1, end]));

end