clear all
close all


%% système

w0 = 2*pi;
zeta = 0.01;


%% excitation

T = 100;
dt = 0.01;
fe = 1/dt;

t = 0:dt:T;
nt = length(t);

excitation = 'bruit'; % 'bruit' 'dirac'

if isequal(excitation, 'bruit')
    f = exp(2i*pi*rand(1, nt)); % bruit
    f(1) = 1;
    for k = 2:nt-1
        f(end+2-k) = conj(f(k));
    end
    f = ifft(f);
elseif isequal(excitation, 'bruit')
    f = ones(1, nt); % dirac
    f = ifft(f);
end

%% reponse

% ondelette
Q = 30;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 0.1;
fmax = 2;
NbFreq = 300;

ct = 3;
cf = 5;

% plots

if true % f plot
    fig = figure('Name', 'bruit blanc');
    ax = axes(fig);
    pltf = plot(ax, t, f);
    xlabel(ax, 't');
    ylabel(ax, 'x');
    
    
    WaveletMenu('WaveletPlot', pltf, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
        'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
        , 'CtEdgeEffects', ct);
end

x = reponseSyst1ddl(t, f, w0, zeta);

fig = figure('Name', ['zeta=', num2str(zeta)]);
ax = axes(fig);
plt = plot(ax, t, x);
xlabel(ax, 't');
ylabel(ax, 'x');

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct);


%% autocorrelation

Rx = autocorr(x);
fig = figure('Name', ['autocorrelation;zeta=', num2str(zeta)]);
ax = axes(fig);
plt = plot(ax, t, Rx);
xlabel(ax, 't');
ylabel(ax, 'Rx');



%% stabilisation

freqs = (0:(length(x)-1)) / (t(end)-t(1));
freqs = freqs(1:floor(length(freqs)/2));
fftX = fft(x);
fftX = fftX(1:floor(length(fftX)/2));

figure;
modalsd(fftX.', freqs, fe, 'MaxModes', 20);



