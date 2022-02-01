function plt = SlicePlot(sliceAxis, ft0, t, freqs, wavelet, plotQuantity, Q, ctEdgeEffects, MotherWavelet, WvltScale, figureTitle, axesTitle, FreqScale, CWTunit)


if ~ismember(sliceAxis, {'freq', 'time'}) % freq = constant freq, time = constant time
    error('');
end
if ~ismember(plotQuantity, {'abs', 'module', 'arg', 'phase'})
    error('');
end

if nargin < 9
    WvltScale = 'log';
end
if nargin < 10
    figureTitle = '';
end
if nargin < 11
    axesTitle = '';
end
if nargin < 12
    FreqScale = 'lin';
end
if nargin < 13
    CWTunit = 'm/s²';
end


switch WvltScale
    case {'lin', 'log'}
    otherwise
        error('wrong scale');
end

%% slice

switch sliceAxis
    case 'freq'
        if length(freqs) == 1 && freqs == ft0
            % ok
        elseif length(freqs) > 1 && ft0 >= freqs(1) && ft0 <= freqs(end) % pas prévu normalement
            k = 2;
            while freqs(k) < ft0
                k = k+1;
            end
            r = (ft0 - freqs(k-1)) / (freqs(k) - freqs(k-1));
            wavelet = (1-r)*wavelet(k-1, :) + r*wavelet(k, :);
        else
            error('');
        end
    case 'time'
        if length(t) == 1 && t == ft0
            % ok mais pas prévu normalement
        elseif length(t) > 1 && ft0 >= t(1) && ft0 <= t(end)
            k = 2;
            while t(k) < ft0
                k = k+1;
            end
            r = (ft0 - t(k-1)) / (t(k) - t(k-1));
            wavelet = (1-r)*wavelet(:, k-1) + r*wavelet(:, k);
            wavelet = transpose(wavelet);
        else
            error('');
        end
end

switch plotQuantity
    case {'abs', 'module'}
        wavelet = abs(wavelet);
        YLab = ['Amplitude [', CWTunit, ']'];
        ordinateScale = WvltScale;
    case {'arg', 'phase'}
        wavelet = angle(wavelet);
        YLab = 'Phase [rad]';
        ordinateScale = 'lin';
end


%% zones d'effets de bord

waveletEdge = wavelet;
tEdge = t;
freqsEdge = freqs;

[~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);

switch sliceAxis
    case 'freq'
        dt = ctEdgeEffects * DeltaT(ft0);
        wavelet = wavelet(t >= t(1) + dt & t <= t(end) - dt);
        t = t(t >= t(1) + dt & t <= t(end) - dt);
    case 'time'
        for kf = 1:length(freqs)
            dt = ctEdgeEffects * DeltaT(freqs(kf));
            if ft0 < t(1) + dt || ft0 > t(end) - dt
                wavelet(kf) = nan;
            end
        end
        freqs = freqs(~isnan(wavelet));
        wavelet = wavelet(~isnan(wavelet));
end

%% figure

fig = figure('Name', figureTitle);
ax = axes(fig);
hold(ax, 'on');

%% plot

switch sliceAxis
    case 'freq'
        plot(ax, tEdge, waveletEdge, '--');
        ax.ColorOrderIndex = ax.ColorOrderIndex-1;
        plt = plot(ax, t, wavelet);
        xlabel(ax, 'Time [s]');
        ylabel(ax, YLab);
        set(ax, 'YScale', ordinateScale);
        set(ax, 'XLim', tEdge([1, end]));
    case 'time'
        plot(ax, freqsEdge, waveletEdge, '--');
        ax.ColorOrderIndex = ax.ColorOrderIndex-1;
        plt = plot(ax, freqs, wavelet);
        xlabel(ax, 'Frequency [Hz]');
        ylabel(ax, YLab);
        set(ax, 'YScale', ordinateScale);
        set(ax, 'XLim', freqsEdge([1, end]));
        set(ax, 'XScale', FreqScale);
end

title(ax, axesTitle);

end