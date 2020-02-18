%% levure

clear all
close all
load("article francoise/data/levure.mat");

% rééchantillonage
dt0 = dt;
fe0 = 1/dt0;
fe = 2000;
deltaN = ceil(fe0/fe);
dt = dt * deltaN;
fe = fe0 / deltaN;


x = Vert_defl(1:deltaN:end);
t = dt * (0:(length(x)-1));


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


%% diagramme stabilisation

freqs = (0:(length(x)-1)) / (t(end)-t(1));
freqs = freqs(1:floor(length(freqs)/2));
fftX = fft(x);
fftX = fftX(1:floor(length(fftX)/2));
fftA = freqs .* fftX;

% fftX = rand(1, length(freqs));

figure;
modalsd(fftA.', freqs, fe, 'MaxModes', 10);


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

