%% levure

clear all
close all
load("article francoise/data/levure.mat");

x = Vert_defl;
t = dt * (0:(length(x)-1));

% filtrage
if true
    ordre = 3;
    f0 = 10;
    f1 = 200;
    n0 = floor(f0 * (t(end)-t(1)));
    n1 = ceil(f1 * (t(end)-t(1)));
    Fx = fft(x);
    f = (0:length(x)-1) / (t(end)-t(1));
    for k = 1:ordre
        Fx = Fx .* (1i*f/f0 ./ (1 + 1i*f/f0));
    end
    for k = 2:length(Fx)-1
        Fx(end+2-k) = conj(Fx(k));
    end
    x = ifft(Fx);
    x = real(x);
end

if false % filtrage
    f0 = 0.5;
    a = [1, 2*pi*f0];
    b = 1;
    x = filter(b, a, x);
end


% autocorrélation
Rx = xcorr(x, 'biased') / var(x);
Rx = Rx(ceil(length(Rx)/2):end);



% rééchantillonage
dt0 = dt;
fe0 = 1/dt0;
fe = 2000;
deltaN = ceil(fe0/fe);
dt = dt * deltaN;
fe = fe0 / deltaN;


x = x(1:deltaN:end);
Rx = Rx(1:deltaN:end);
t = dt * (0:(length(x)-1));


% limites
t0 = 0;
tf = inf;
% t0 = 35;
% tf = 37;
x = x(t>=t0 & t<tf);
t = t(t>=t0 & t<tf);

% plot
fig =  figure;
ax = axes(fig);
plt = plot(ax, t, x);
xlabel(ax, 't');
ylabel(ax, 'x');

fig =  figure;
ax = axes(fig);
plt2 = plot(ax, t, Rx);
xlabel(ax, 't');
ylabel(ax, 'Rx');


%% diagramme stabilisation

if false
    freqs = (0:(length(x)-1)) / (t(end)-t(1));
    freqs = freqs(1:floor(length(freqs)/2));
    fftX = fft(x);
    fftX = fftX(1:floor(length(fftX)/2));
    fftA = freqs .* fftX;
    
    % fftX = rand(1, length(freqs));
    
    figure;
    modalsd(fftA.', freqs, fe, 'MaxModes', 10);
end


%% ondelette
Q = 10;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 30;
fmax = 500;
NbFreq = 300;

ct = 3;
cf = 5;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct);

WaveletMenu('WaveletPlot', plt2, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct);

