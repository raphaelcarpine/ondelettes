function [fmax, zeta] = PeakPickingMax(varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here

defaultLine = [];
defaultFreqs = [];
defaultFFT = [];
defaultQuadraticMax = false;
defaultQuadraticFFT = false;
defaultBoundsFreq = [-inf, inf];
defaultPlotPeak = true;

p = inputParser();
addOptional(p, 'Line', defaultLine);
addOptional(p, 'Freqs', defaultFreqs);
addOptional(p, 'FFT', defaultFFT);
addOptional(p, 'QuadraticMax', defaultQuadraticMax);
addOptional(p, 'QuadraticFFT', defaultQuadraticFFT);
addOptional(p, 'BoundsFreq', defaultBoundsFreq);
addOptional(p, 'PlotPeak', defaultPlotPeak);
parse(p, varargin{:});

line = p.Results.Line;
freqs = p.Results.Freqs;
fft = p.Results.FFT;
quadraticMax = p.Results.QuadraticMax;
quadraticFFT = p.Results.QuadraticFFT;
boundsFreq = p.Results.BoundsFreq;
plotPeak = p.Results.PlotPeak;


%% data

if ~isempty(line)
    X = get(line, 'XData');
    Y = get(line, 'YData');
    Y = Y(boundsFreq(1) <= X & X <= boundsFreq(2));
    X = X(boundsFreq(1) <= X & X <= boundsFreq(2));
elseif ~isempty(freqs) && isequal(size(fft), size(freqs))
    X = freqs;
    Y = fft;
else
    error('no data');
end


%% find local maxs
localMax = [];
for k = 2:length(Y)-1
    if Y(k-1) < Y(k) && Y(k) >= Y(k+1)
        localMax(end+1) = k;
    end
end
[~, I] = sort(Y(localMax), 'descend');
localMax = localMax(I(1));


%% fmax
if quadraticMax
    [fmax, Ymax] = localMax3Points(X(localMax-1:localMax+1).', Y(localMax-1:localMax+1).');
else
    fmax = X(localMax);
    Ymax = Y(localMax);
end

%% zeta
if quadraticFFT
    H = Ymax/2;
else
    H = Ymax/sqrt(2);
end
k1 = localMax-1;
while k1 > 0 && Y(k1) > H
    k1 = k1 - 1;
end
k2 = localMax+1;
while k2 <= length(Y) && Y(k2) > H
    k2 = k2 + 1;
end
if k1 < 1 || k2 > length(Y)
    zeta = nan;
else
    fh1 = X(k1) + (H-Y(k1))*(X(k1+1)-X(k1))/(Y(k1+1)-Y(k1));
    fh2 = X(k2-1) + (H-Y(k2-1))*(X(k2)-X(k2-1))/(Y(k2)-Y(k2-1));
    Q = fmax/(fh2-fh1);
    zeta = 1/(2*Q);
end


%% plot

if plotPeak
    if isempty(line)
        figure;
        line = plot(freqs, fft);
        xlabel('Frequency [Hz]');
    end
    ax = get(line, 'Parent');
    hold(ax, 'on');
    
    scatter(ax, fmax, Ymax, '+', 'MarkerEdgeColor', 'red',...
        'LineWidth', 2, 'HandleVisibility', 'off');
    yline(ax, H, 'Color', 0.5*[1 1 1], 'HandleVisibility', 'off');
    text(ax, fmax, Ymax, [sprintf('  f_{max} = %.2f Hz ; ', fmax), '\zeta',...
        sprintf(' = %.2f %%', 100*zeta)]);
end


end

